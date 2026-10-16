"""
Data Quality Control & Validation Engine for Population Pharmacogenomics.
Handles column normalization, missing genotype detection, duplicate sample detection,
gender standardization, and invalid allele detection.
"""

import pandas as pd
import numpy as np
from typing import Dict, Tuple, List, Set, Any
import logging

from config import SNP_CONFIG

logger = logging.getLogger(__name__)


class DataQCEngine:
    def __init__(self, df: pd.DataFrame, existing_sample_ids: Set[str] = None, custom_col_mapping: Dict[str, str] = None, custom_snp_rules: Dict[str, Dict] = None):
        """
        Initialize QC Engine with a raw DataFrame.
        """
        self.raw_df = df.copy()
        self.existing_sample_ids = set(existing_sample_ids) if existing_sample_ids else set()
        self.custom_col_mapping = custom_col_mapping or {}
        self.snp_rules = dict(SNP_CONFIG)
        if custom_snp_rules:
            for k, v in custom_snp_rules.items():
                if k in self.snp_rules:
                    self.snp_rules[k] = {**self.snp_rules[k], **v}
                else:
                    self.snp_rules[k] = v
        self.clean_df = pd.DataFrame()
        self.qc_report = {}

    def normalize_column_names(self) -> pd.DataFrame:
        """
        Maps raw input headers to standard column names.
        Deduplicates column names to ensure single Series access per variable.
        """
        df = self.raw_df.copy()
        col_mapping = {}
        
        # 1. Apply user custom mappings if provided
        for raw_col, target_col in self.custom_col_mapping.items():
            if raw_col in df.columns and target_col and target_col != "(None)":
                col_mapping[raw_col] = target_col
        
        for col in df.columns:
            if col in col_mapping:
                continue
            col_str = str(col).strip()
            col_lower = col_str.lower()
            
            # Specific SNP rules first
            if 'cyp2c19*2' in col_lower or 'rs4244285' in col_lower:
                col_mapping[col] = 'CYP2C19*2'
            elif 'cyp2c19*3' in col_lower or 'rs4986893' in col_lower:
                col_mapping[col] = 'CYP2C19*3'
            elif 'cyp2c19*17' in col_lower or 'rs12248560' in col_lower:
                col_mapping[col] = 'CYP2C19*17'
            elif '3 generation' in col_lower or '3_generation' in col_lower or 'past 3' in col_lower or 'generation' in col_lower:
                col_mapping[col] = 'Three_Generations_Residency'
            elif 'native' in col_lower:
                col_mapping[col] = 'Native_Place'
            elif 'sample' in col_lower and 'id' in col_lower:
                col_mapping[col] = 'Sample_ID'
            elif 'gender' in col_lower or 'sex' in col_lower:
                col_mapping[col] = 'Gender'
            elif 'birth' in col_lower or 'dob' in col_lower:
                col_mapping[col] = 'Date_of_Birth'
            elif 'state' in col_lower:
                col_mapping[col] = 'State'
            elif 'region' in col_lower:
                col_mapping[col] = 'Region'
            elif 'test' in col_lower:
                col_mapping[col] = 'Test_Requested'
        
        df = df.rename(columns=col_mapping)
        # Deduplicate columns if multiple mapped to the same name
        df = df.loc[:, ~df.columns.duplicated()].copy()
        return df

    def _get_series(self, df: pd.DataFrame, col_name: str) -> pd.Series:
        """Helper to safely retrieve a 1D Pandas Series."""
        if col_name not in df.columns:
            return pd.Series(dtype=str, index=df.index)
        val = df[col_name]
        if isinstance(val, pd.DataFrame):
            val = val.iloc[:, 0]
        return val

    def run_qc(self) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Executes full data quality control & cleaning pipeline.
        """
        df = self.normalize_column_names()
        total_input_samples = len(df)
        
        # Ensure Sample_ID column exists
        if 'Sample_ID' not in df.columns:
            df['Sample_ID'] = [f"SAMPLE_{i+1:05d}" for i in range(len(df))]
        
        sample_ids = self._get_series(df, 'Sample_ID').astype(str).str.strip()
        df['Sample_ID'] = sample_ids
        
        # Check duplicates within batch file
        duplicated_in_file = df[df['Sample_ID'].duplicated()]['Sample_ID'].tolist()
        all_duplicates = list(set(duplicated_in_file))
        
        # Clean Gender
        if 'Gender' in df.columns:
            gender_series = self._get_series(df, 'Gender').astype(str).str.strip().str.lower()
            gender_map = {
                'male': 'Male', 'm': 'Male',
                'female': 'Female', 'f': 'Female',
                'null': 'Unknown', 'not mentioned': 'Unknown', 'not specified': 'Unknown', 'nan': 'Unknown'
            }
            df['Gender_Clean'] = gender_series.map(lambda x: gender_map.get(x, 'Unknown'))
        else:
            df['Gender_Clean'] = 'Unknown'
            
        # Clean State & Native Place
        df['State_Clean'] = self._get_series(df, 'State').astype(str).str.strip() if 'State' in df.columns else ''
        df['Native_Place_Clean'] = self._get_series(df, 'Native_Place').astype(str).str.strip() if 'Native_Place' in df.columns else ''

        # SNP Validation
        snp_qc_stats = {}
        invalid_genotypes_log = {}
        
        for snp_name, snp_info in self.snp_rules.items():
            if snp_name in df.columns:
                raw_vals = self._get_series(df, snp_name)
                
                # Format genotype strings
                def clean_genotype(val):
                    if pd.isna(val):
                        return np.nan
                    v = str(val).strip().upper().replace('/', '')
                    if v in ['', 'NAN', 'NONE', 'NULL', 'NA', '?']:
                        return np.nan
                    
                    if v in snp_info['valid_genotypes']:
                        return v
                    return np.nan

                cleaned_vals = raw_vals.apply(clean_genotype)
                
                # Detect invalid non-empty values
                raw_str = raw_vals.astype(str).str.strip().str.upper()
                raw_not_na = raw_vals.notna() & (~raw_str.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?']))
                invalid_mask = raw_not_na & cleaned_vals.isna()
                invalid_samples = df.loc[invalid_mask, 'Sample_ID'].tolist()
                invalid_raw_vals = raw_vals.loc[invalid_mask].tolist()
                
                if invalid_samples:
                    invalid_genotypes_log[snp_name] = [
                        {'sample_id': sid, 'raw_value': rval} 
                        for sid, rval in zip(invalid_samples, invalid_raw_vals)
                    ]
                
                df[snp_name] = cleaned_vals
                
                # Stats per SNP
                valid_count = df[snp_name].notna().sum()
                missing_count = len(df) - valid_count - len(invalid_samples)
                missing_pct = (missing_count / total_input_samples) * 100.0 if total_input_samples > 0 else 0.0
                
                snp_qc_stats[snp_name] = {
                    'valid_samples': int(valid_count),
                    'missing_samples': int(missing_count),
                    'invalid_samples': int(len(invalid_samples)),
                    'missing_pct': round(missing_pct, 2)
                }

        total_invalid_calls = sum(len(items) for items in invalid_genotypes_log.values())
        fully_valid = max(0, total_input_samples - len(all_duplicates) - total_invalid_calls)

        self.clean_df = df
        self.qc_report = {
            'total_samples': total_input_samples,
            'fully_valid_samples': fully_valid,
            'duplicate_sample_count': len(all_duplicates),
            'duplicate_sample_ids': all_duplicates,
            'invalid_genotypes_count': total_invalid_calls,
            'invalid_genotypes_log': invalid_genotypes_log,
            'snp_qc_stats': snp_qc_stats
        }
        
        return self.clean_df, self.qc_report

