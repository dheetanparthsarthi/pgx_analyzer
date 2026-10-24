"""
Data Quality Control & Validation Engine for Population Pharmacogenomics.
Handles column normalization, missing genotype detection, duplicate sample detection,
gender standardization, and invalid allele detection.
"""

import pandas as pd
import numpy as np
from typing import Dict, Tuple, List, Set, Any
import logging

from config import SNP_CONFIG, CYP2C9_SNP_CONFIG

logger = logging.getLogger(__name__)


class DataQCEngine:
    def __init__(self, df: pd.DataFrame, existing_sample_ids: Set[str] = None, custom_col_mapping: Dict[str, str] = None, custom_snp_rules: Dict[str, Dict] = None, *args, **kwargs):
        """
        Initialize QC Engine with a raw DataFrame.
        """
        self.raw_df = df.copy()
        self.existing_sample_ids = set(existing_sample_ids) if existing_sample_ids else set()
        self.custom_col_mapping = custom_col_mapping or {}
        # Automatically include all standard pharmacogenetic SNP rules
        self.snp_rules = {**SNP_CONFIG, **CYP2C9_SNP_CONFIG}
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
            elif 'cyp2c9*2' in col_lower or 'rs1799853' in col_lower:
                col_mapping[col] = 'CYP2C9*2'
            elif 'cyp2c9*3' in col_lower or 'rs1057910' in col_lower:
                col_mapping[col] = 'CYP2C9*3'
            elif '3 generation' in col_lower or '3_generation' in col_lower or 'past 3' in col_lower or 'generation' in col_lower:
                col_mapping[col] = 'Three_Generations_Residency'
            elif 'native' in col_lower:
                col_mapping[col] = 'Native_Place'
            elif 'sample' in col_lower and 'id' in col_lower:
                col_mapping[col] = 'Sample_ID'
            elif 'lab' in col_lower and 'id' in col_lower:
                col_mapping[col] = 'Sample_ID'
            elif 'patient' in col_lower and 'id' in col_lower:
                col_mapping[col] = 'Sample_ID'
            elif col_lower == 'id':
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
        Calculates exact duplicate samples, missing genotypes, and invalid calls
        dynamically from the input dataset.
        """
        df = self.normalize_column_names()
        df = df.reset_index(drop=True)
        total_input_samples = len(df)
        
        # Ensure Sample_ID column exists
        if 'Sample_ID' not in df.columns:
            df['Sample_ID'] = [f"SAMPLE_{i+1:05d}" for i in range(len(df))]
        
        raw_sids = self._get_series(df, 'Sample_ID')
        sample_ids = [
            str(s).strip() if (pd.notna(s) and str(s).strip() not in ('', 'nan', 'None', 'NAN', 'NULL'))
            else f"SAMPLE_{i+1:05d}"
            for i, s in enumerate(raw_sids)
        ]
        df['Sample_ID'] = sample_ids
        
        # Check duplicates within batch file
        duplicated_in_file = df[df['Sample_ID'].duplicated()]['Sample_ID'].tolist()
        all_duplicates = sorted(list({str(s).strip() for s in duplicated_in_file}), key=lambda x: str(x))
        
        # Clean Gender
        if 'Gender' in df.columns:
            gender_series = self._get_series(df, 'Gender').astype(str).str.strip().str.lower()
            gender_map = {
                'male': 'Male', 'm': 'Male', 'man': 'Male',
                'female': 'Female', 'f': 'Female', 'woman': 'Female',
            }
            missing_gender = {'', 'null', 'none', 'not mentioned', 'not specified', 'nan', 'na', 'unknown'}
            df['Gender_Clean'] = gender_series.map(
                lambda value: 'Unknown' if value in missing_gender else gender_map.get(value, str(value).title())
            )
        else:
            df['Gender_Clean'] = 'Unknown'
            
        # Clean State & Native Place
        if 'State' in df.columns:
            state_raw_series = self._get_series(df, 'State')
            state_values = state_raw_series.astype(str).str.replace('\xa0', ' ', regex=False).str.strip()
            is_unclear_mask = (
                state_raw_series.isna()
                | state_values.str.lower().isin({'', 'nan', 'none', 'null', '<na>', 'unknown'})
                | state_values.str.lower().str.contains('needs review', regex=False)
                | (state_values.str.lower() == 'india')
            )
            df['State_Clean'] = state_values.mask(is_unclear_mask, '')
        else:
            df['State_Clean'] = ''
        df['Native_Place_Clean'] = self._get_series(df, 'Native_Place').astype(str).str.strip() if 'Native_Place' in df.columns else ''

        # Determine target genotype / marker columns dynamically
        non_marker_std = {
            'Sample_ID', 'Gender', 'Gender_Clean', 'State', 'State_Clean', 
            'Native_Place', 'Native_Place_Clean', 'Region', 'Date_of_Birth', 
            'Test_Requested', 'Three_Generations_Residency', 'Name'
        }
        target_cols = []
        for c in df.columns:
            if c in non_marker_std or str(c).startswith('Unnamed:'):
                continue
            if c in self.snp_rules or any(k in str(c).lower() for k in ['*', 'rs', 'cyp', 'vkorc', 'dpyd', 'tpmt', 'slco', 'snp', 'marker']):
                target_cols.append(c)
        if not target_cols:
            for raw_col, std in self.custom_col_mapping.items():
                if std not in non_marker_std and std in df.columns and std not in target_cols:
                    target_cols.append(std)

        # SNP Validation & Missing Call Tracking
        snp_qc_stats = {}
        invalid_genotypes_log = {}
        samples_with_missing_genotypes = set()
        samples_with_invalid_genotypes = set()
        total_missing_calls = 0

        for col_name in target_cols:
            raw_vals = self._get_series(df, col_name)
            raw_str = raw_vals.astype(str).str.strip().str.upper()
            is_missing_mask = raw_vals.isna() | raw_str.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN', 'NOT PERFORMED'])
            missing_count = int(is_missing_mask.sum())
            total_missing_calls += missing_count
            missing_sids = [str(s).strip() for s in df.loc[is_missing_mask, 'Sample_ID'] if pd.notna(s)]
            samples_with_missing_genotypes.update(missing_sids)

            # Clean and validate genotypes
            snp_info = self.snp_rules.get(col_name)
            valid_set = snp_info.get('valid_genotypes', ()) if snp_info else ()
            if valid_set:
                def clean_genotype(val):
                    if pd.isna(val): return np.nan
                    v = str(val).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
                    if v in ['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN']:
                        return np.nan
                    if v in valid_set:
                        return v
                    return np.nan
                cleaned_vals = raw_vals.apply(clean_genotype)
                invalid_mask = (~is_missing_mask) & cleaned_vals.isna()
            else:
                def clean_genotype(val):
                    if pd.isna(val): return np.nan
                    v = str(val).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
                    if v in ['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN']:
                        return np.nan
                    if v.isalpha() and len(v) in (1, 2):
                        return v
                    return np.nan
                cleaned_vals = raw_vals.apply(clean_genotype)
                invalid_mask = (~is_missing_mask) & cleaned_vals.isna()

            invalid_samples = [str(s).strip() for s in df.loc[invalid_mask, 'Sample_ID'] if pd.notna(s)]
            invalid_raw_vals = raw_vals.loc[invalid_mask].tolist()
            if invalid_samples:
                invalid_genotypes_log[col_name] = [
                    {'sample_id': sid, 'raw_value': rval}
                    for sid, rval in zip(invalid_samples, invalid_raw_vals)
                ]
                samples_with_invalid_genotypes.update(invalid_samples)

            df[col_name] = cleaned_vals
            valid_count = int(cleaned_vals.notna().sum())
            missing_pct = (missing_count / total_input_samples) * 100.0 if total_input_samples > 0 else 0.0
            snp_qc_stats[col_name] = {
                'valid_samples': valid_count,
                'missing_samples': missing_count,
                'invalid_samples': len(invalid_samples),
                'missing_pct': round(missing_pct, 2)
            }

        total_invalid_calls = sum(len(items) for items in invalid_genotypes_log.values())
        num_markers = len(target_cols)
        total_genotype_opportunities = total_input_samples * num_markers if num_markers > 0 else total_input_samples
        overall_missing_pct = round((total_missing_calls / total_genotype_opportunities) * 100.0, 2) if total_genotype_opportunities > 0 else 0.0

        # Fully valid samples: samples with unique ID, no invalid calls, and no missing calls
        problematic_samples = set(all_duplicates) | samples_with_missing_genotypes | samples_with_invalid_genotypes
        fully_valid_samples = max(0, total_input_samples - len(problematic_samples))

        # 4-Field Completeness Audit & Discarded Samples Tracking
        # 1. Sample ID Missing
        missing_sid_mask = [
            pd.isna(s) or str(s).strip().lower() in ('', 'nan', 'none', 'null')
            for s in raw_sids
        ]
        missing_sid_count = int(sum(missing_sid_mask))

        # 2. Patient Name Missing
        name_col = next((c for c in df.columns if c.lower() == 'name' or 'name' in c.lower()), None)
        if name_col:
            name_series = self._get_series(df, name_col)
            missing_name_mask = [
                pd.isna(n) or str(n).strip().lower() in ('', 'nan', 'none', 'null')
                for n in name_series
            ]
            missing_name_count = int(sum(missing_name_mask))
        else:
            missing_name_mask = [False] * total_input_samples
            missing_name_count = 0

        # 3. State Missing
        state_clean_series = df.get('State_Clean', pd.Series([''] * total_input_samples))
        missing_state_mask = [
            pd.isna(st) or str(st).strip().lower() in ('', 'nan', 'none', 'null')
            for st in state_clean_series
        ]
        missing_state_count = int(sum(missing_state_mask))

        # 4. Target Marker Missing
        missing_markers_count = len(samples_with_missing_genotypes)

        completeness_audit = [
            {
                'field': 'Sample ID',
                'description': 'Unique Patient Identifier',
                'total': total_input_samples,
                'missing': missing_sid_count,
                'valid': total_input_samples - missing_sid_count,
                'completion_pct': round(((total_input_samples - missing_sid_count) / total_input_samples) * 100, 1) if total_input_samples else 100.0,
                'status': 'PASS' if missing_sid_count == 0 else f'FLAGGED ({missing_sid_count} missing)'
            },
            {
                'field': 'Patient Name',
                'description': 'Subject Full Name' if name_col else 'Not Mapped in Dataset',
                'total': total_input_samples,
                'missing': missing_name_count,
                'valid': total_input_samples - missing_name_count,
                'completion_pct': round(((total_input_samples - missing_name_count) / total_input_samples) * 100, 1) if total_input_samples else 100.0,
                'status': 'PASS' if missing_name_count == 0 else f'FLAGGED ({missing_name_count} missing)'
            },
            {
                'field': 'State (Demographic)',
                'description': 'Indian State of Origin',
                'total': total_input_samples,
                'missing': missing_state_count,
                'valid': total_input_samples - missing_state_count,
                'completion_pct': round(((total_input_samples - missing_state_count) / total_input_samples) * 100, 1) if total_input_samples else 100.0,
                'status': 'PASS' if missing_state_count == 0 else f'FLAGGED ({missing_state_count} missing)'
            },
            {
                'field': 'Genetic Markers (Genotypes)',
                'description': f"Target Panel ({', '.join(target_cols)})" if target_cols else 'Genotype Calls',
                'total': total_input_samples,
                'missing': missing_markers_count,
                'valid': total_input_samples - missing_markers_count,
                'completion_pct': round(((total_input_samples - missing_markers_count) / total_input_samples) * 100, 1) if total_input_samples else 100.0,
                'status': 'PASS' if missing_markers_count == 0 else f'FLAGGED ({missing_markers_count} missing)'
            }
        ]

        # Row-by-row Discarded Samples Log
        discarded_samples_log = []
        for i in range(total_input_samples):
            sid = str(df.loc[i, 'Sample_ID']).strip()
            is_sid_bad = missing_sid_mask[i]
            is_name_bad = missing_name_mask[i]
            is_state_bad = missing_state_mask[i]

            marker_reasons = []
            row_marker_calls = {}
            for col_name in target_cols:
                c_val = df.loc[i, col_name] if col_name in df.columns else np.nan
                row_marker_calls[col_name] = str(c_val) if pd.notna(c_val) else 'Missing'
                if pd.isna(c_val) or str(c_val).strip().upper() in ('', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'):
                    marker_reasons.append(f"Missing {col_name}")

            reasons = []
            if is_sid_bad:
                reasons.append("Missing Sample ID")
            if is_name_bad:
                reasons.append("Missing Patient Name")
            if is_state_bad:
                raw_st_val = str(df.loc[i, 'State']).strip() if 'State' in df.columns and pd.notna(df.loc[i, 'State']) else ''
                raw_st_lower = raw_st_val.lower()
                if 'needs' in raw_st_lower and 'review' in raw_st_lower:
                    reasons.append("Ambiguous State ('Needs review')")
                elif raw_st_lower == 'india':
                    reasons.append("Unspecified State ('India')")
                else:
                    reasons.append("Missing State")
            if marker_reasons:
                reasons.extend(marker_reasons)

            if reasons:
                name_val = str(df.loc[i, name_col]).strip() if name_col and pd.notna(df.loc[i, name_col]) else 'N/A'
                raw_st = str(df.loc[i, 'State']).strip() if 'State' in df.columns and pd.notna(df.loc[i, 'State']) else 'Missing'
                discarded_samples_log.append({
                    'sample_id': sid,
                    'name': name_val,
                    'state': raw_st,
                    'reasons': '; '.join(reasons),
                    'marker_calls': row_marker_calls
                })

        self.clean_df = df
        self.qc_report = {
            'total_samples': total_input_samples,
            'fully_valid_samples': fully_valid_samples,
            'duplicate_sample_count': len(all_duplicates),
            'duplicate_sample_ids': sorted([str(s) for s in all_duplicates], key=lambda x: str(x)),
            'total_missing_genotypes_count': total_missing_calls,
            'overall_missing_pct': overall_missing_pct,
            'samples_with_missing_count': len(samples_with_missing_genotypes),
            'samples_with_missing_ids': sorted([str(s) for s in samples_with_missing_genotypes], key=lambda x: str(x)),
            'invalid_genotypes_count': total_invalid_calls,
            'invalid_genotypes_log': invalid_genotypes_log,
            'snp_qc_stats': snp_qc_stats,
            'completeness_audit': completeness_audit,
            'discarded_samples_count': len(discarded_samples_log),
            'discarded_samples_log': discarded_samples_log,
            'missing_sid_count': missing_sid_count,
            'missing_name_count': missing_name_count,
            'missing_state_count': missing_state_count,
            'missing_markers_count': missing_markers_count,
        }
        
        return self.clean_df, self.qc_report
