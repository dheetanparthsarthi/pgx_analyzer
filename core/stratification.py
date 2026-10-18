"""
Demographic & Geographic Stratification Engine for Population Pharmacogenomics.
Automatically groups datasets into Overall, Regional (North/South/East/West India),
State-Wise (Individual States), and Gender subgroups, computing full genetic and pharmacogenomic parameters for each.
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, List, Optional
import re

from config import STATE_TO_REGION, CITY_TO_STATE, SNP_CONFIG, CYP2C9_SNP_CONFIG
from core.stats import GeneticStatsEngine
from core.cyp2c19 import CYP2C19Translator
from core.cyp2c9 import CYP2C9Translator


class DemographicStratifier:
    @staticmethod
    def assign_region(row: pd.Series) -> str:
        """
        Assigns geographic region (North India, South India, East India, West India, Central India)
        to a sample based on State, Native Place, or pre-classified text.
        """
        # Detect state from all possible column aliases
        state_candidates = ['State_Clean', 'State', 'State / Country', 'State/Country', 'state', 'Region']
        state_raw = ''
        for sc in state_candidates:
            if sc in row and pd.notna(row[sc]) and str(row[sc]).strip() not in ('', 'nan', 'Null', 'None'):
                state_raw = str(row[sc]).strip()
                break
                
        # Detect native place from aliases
        native_candidates = ['Native_Place_Clean', 'Native_Place', 'Native place ', 'Native place', 'Native Place', 'native place']
        native_place_raw = ''
        for nc in native_candidates:
            if nc in row and pd.notna(row[nc]) and str(row[nc]).strip() not in ('', 'nan', 'Null', 'None'):
                native_place_raw = str(row[nc]).strip().lower()
                break

        # Clean non-breaking spaces and redundant whitespace
        state_raw = re.sub(r'[\s\u00a0]+', ' ', state_raw).strip()
        native_place_raw = re.sub(r'[\s\u00a0]+', ' ', native_place_raw).strip()
        region_raw = str(row.get('Region', '')).strip() if 'Region' in row else ''
        
        # 0. Check direct Region column match
        if region_raw:
            reg_clean = re.sub(r'[\s\u00a0]+', ' ', region_raw).title()
            if 'Central' in reg_clean: return 'Central India'
            if 'North' in reg_clean: return 'North India'
            if 'South' in reg_clean: return 'South India'
            if 'East' in reg_clean:  return 'East India'
            if 'West' in reg_clean:  return 'West India'

        # 1. Direct match on State
        if state_raw in STATE_TO_REGION:
            return STATE_TO_REGION[state_raw]
            
        # 2. Case & whitespace insensitive match
        state_title = state_raw.title()
        for k, v in STATE_TO_REGION.items():
            if k.lower() == state_raw.lower():
                return v

        # 3. Check for direct region string (e.g. "North India ", "South India ", "Central India ")
        state_clean = re.sub(r'[\s\u00a0]+', ' ', state_raw).title()
        if 'Central' in state_clean: return 'Central India'
        if 'North' in state_clean: return 'North India'
        if 'South' in state_clean: return 'South India'
        if 'East' in state_clean:  return 'East India'
        if 'West' in state_clean:  return 'West India'
        
        # 4. Predict State from Native Place (City/District gazetteer)
        for city, mapped_state in CITY_TO_STATE.items():
            if city in native_place_raw:
                return STATE_TO_REGION.get(mapped_state, 'Unclassified / Other')
                
        # 5. Check if native place directly matches a state name
        for state_name, reg_name in STATE_TO_REGION.items():
            if state_name.lower() in native_place_raw:
                return reg_name

        return 'Unclassified / Other'

    @classmethod
    def stratify_and_analyze(cls, df: pd.DataFrame, custom_snp_rules: Optional[Dict] = None, active_gene: Optional[str] = None) -> Dict[str, Any]:
        """
        Performs end-to-end demographic stratification and statistical calculations.
        Detects target gene (CYP2C19, CYP2C9, or custom) and computes full CPIC diplotypes,
        phenotypes, allele frequencies, HWE tests, and regional stratification.
        """
        data = df.copy()
        
        # Ensure Region column is populated
        data['Region'] = data.apply(cls.assign_region, axis=1)
        
        # Detect active gene from parameter or column names
        if not active_gene:
            cols_lower = [str(c).lower().replace(" ", "").replace("_", "").replace("-", "") for c in data.columns]
            is_cyp2c9 = any('cyp2c9' in c or 'rs1799853' in c or 'rs1057910' in c for c in cols_lower)
            is_cyp2c19 = any('cyp2c19' in c or 'rs4244285' in c or 'rs4986893' in c or 'rs12248560' in c for c in cols_lower)
            active_gene = 'CYP2C9' if is_cyp2c9 else ('CYP2C19' if is_cyp2c19 else 'Custom')
        else:
            active_gene = str(active_gene).upper().strip()
            is_cyp2c9 = (active_gene == 'CYP2C9')
            is_cyp2c19 = (active_gene == 'CYP2C19')

        if custom_snp_rules:
            active_snp_rules = dict(custom_snp_rules)
        elif is_cyp2c9:
            active_snp_rules = CYP2C9_SNP_CONFIG
        elif is_cyp2c19:
            active_snp_rules = SNP_CONFIG
        else:
            active_snp_rules = {}

        # Apply appropriate translator across full dataset
        if is_cyp2c9:
            data_processed, overall_pgx_summary = CYP2C9Translator.process_dataframe(data)
        elif is_cyp2c19:
            data_processed, overall_pgx_summary = CYP2C19Translator.process_dataframe(data)
        else:
            data_processed = data
            overall_pgx_summary = {'diplotypes': {}, 'phenotypes': {}, 'total_samples': len(data), 'valid_pharmacogenomic_samples': 0}
        
        # Ensure State_Clean and Gender_Clean exist for sub-group stratification
        if 'State_Clean' not in data_processed.columns:
            state_candidates = ['State', 'State / Country', 'State/Country', 'state']
            s_col = next((c for c in state_candidates if c in data_processed.columns), None)
            if s_col:
                data_processed['State_Clean'] = data_processed[s_col].apply(lambda x: re.sub(r'[\s\u00a0]+', ' ', str(x)).strip() if pd.notna(x) and str(x).strip() not in ('', 'nan', 'Null', 'None') else '')
            else:
                data_processed['State_Clean'] = ''

        if 'Gender_Clean' not in data_processed.columns:
            g_col = next((c for c in ['Gender', 'gender', 'sex', 'Sex'] if c in data_processed.columns), None)
            if g_col:
                gender_map = {'male': 'Male', 'm': 'Male', 'female': 'Female', 'f': 'Female'}
                data_processed['Gender_Clean'] = data_processed[g_col].apply(lambda x: gender_map.get(str(x).strip().lower(), 'Unknown') if pd.notna(x) else 'Unknown')
            else:
                data_processed['Gender_Clean'] = 'Unknown'

        # Candidate SNP list: prioritize mapped rules present in columns
        target_snps = [s for s in active_snp_rules.keys() if s in data_processed.columns]
        if not target_snps:
            target_snps = list(active_snp_rules.keys())
        if not target_snps or not any(s in data_processed.columns for s in target_snps):
            target_snps = [c for c in data_processed.columns if ('*' in c or 'rs' in c.lower() or any(g.lower() in c.lower() for g in ['cyp', 'tpmt', 'dpyd', 'vkorc1', 'nudt15', 'slco1b1'])) and c not in ('Sample_ID', 'Region', 'Gender', 'State', 'State_Clean', 'Gender_Clean', 'Diplotype', 'Phenotype', 'Activity_Score', 'Test_Requested', 'Date_of_Birth')]

        results = {
            'overall': {},
            'regional': {},
            'state_wise': {},
            'gender': {},
            'processed_dataframe': data_processed,
            'active_gene': active_gene,
            'target_gene': active_gene,
            'markers': target_snps
        }
        
        # Helper to compute full stats bundle for any dataframe slice
        def analyze_subgroup(sub_df: pd.DataFrame, group_name: str) -> Dict[str, Any]:
            n_sub = len(sub_df)
            
            missing_genotypes_sub = 0
            for snp_name in target_snps:
                matched_col = None
                if snp_name in sub_df.columns:
                    matched_col = snp_name
                else:
                    for c in sub_df.columns:
                        if snp_name.lower() in str(c).lower() or str(c).lower() in snp_name.lower():
                            matched_col = c
                            break
                if matched_col and matched_col in sub_df.columns:
                    missing_genotypes_sub += int(sub_df[matched_col].isna().sum())
                    
            denom_snps = len(target_snps) if target_snps else 1
            missing_pct_sub = (missing_genotypes_sub / (n_sub * denom_snps)) * 100.0 if n_sub > 0 else 0.0

            sub_results = {
                'sample_count': n_sub,
                'missing_count': missing_genotypes_sub,
                'missing_pct': round(missing_pct_sub, 2),
                'snps': {},
                'cyp2c19_summary': {},
                'cyp2c9_summary': {},
                'pgx_summary': {}
            }
            
            # SNP Stats
            for snp_name in target_snps:
                matched_col = None
                if snp_name in sub_df.columns:
                    matched_col = snp_name
                else:
                    for c in sub_df.columns:
                        if snp_name.lower() in str(c).lower() or str(c).lower() in snp_name.lower():
                            matched_col = c
                            break
                if matched_col and matched_col in sub_df.columns:
                    rule_info = active_snp_rules.get(snp_name)
                    sub_results['snps'][snp_name] = GeneticStatsEngine.analyze_snp(sub_df[matched_col], snp_name, custom_info=rule_info)
                    
            # Gene Translation Stats
            if is_cyp2c9:
                _, pgx_sum = CYP2C9Translator.process_dataframe(sub_df)
                sub_results['cyp2c9_summary'] = pgx_sum
                sub_results['cyp2c19_summary'] = pgx_sum  # For backward compatibility with existing tab views
                sub_results['pgx_summary'] = pgx_sum
            elif is_cyp2c19:
                _, pgx_sum = CYP2C19Translator.process_dataframe(sub_df)
                sub_results['cyp2c19_summary'] = pgx_sum
                sub_results['pgx_summary'] = pgx_sum
            else:
                sub_results['pgx_summary'] = {'phenotypes': {}, 'diplotypes': {}, 'total_samples': n_sub}
            
            return sub_results

        # 1. Overall Population Analysis
        results['overall'] = analyze_subgroup(data_processed, 'Overall Population')
        
        # 2. Regional Analysis (North, South, East, West, Central India)
        standard_regions = ['North India', 'South India', 'East India', 'West India', 'Central India']
        for region in standard_regions:
            reg_df = data_processed[data_processed['Region'] == region]
            results['regional'][region] = analyze_subgroup(reg_df, region)
            
        # Catch any unclassified regions
        unclass_df = data_processed[~data_processed['Region'].isin(standard_regions)]
        if len(unclass_df) > 0:
            results['regional']['Unclassified / Other'] = analyze_subgroup(unclass_df, 'Unclassified')

        # 3. State-Wise Analysis (Individual Indian States)
        unique_states = data_processed['State_Clean'].unique() if 'State_Clean' in data_processed.columns else []
        for state in unique_states:
            state_str = str(state).strip()
            if state_str and state_str.lower() != 'nan':
                st_df = data_processed[data_processed['State_Clean'] == state]
                if len(st_df) > 0:
                    results['state_wise'][state_str] = analyze_subgroup(st_df, state_str)

        # 4. Gender Analysis (Male, Female)
        genders = ['Male', 'Female']
        for gender in genders:
            gen_df = data_processed[data_processed['Gender_Clean'] == gender] if 'Gender_Clean' in data_processed.columns else pd.DataFrame()
            if not gen_df.empty:
                results['gender'][gender] = analyze_subgroup(gen_df, gender)
            
        return results
