"""
CYP2C9 Star Allele, Diplotype Call, and Phenotype Classifier.
Implements standard CPIC rules for CYP2C9 pharmacogenomics.
Biomarkers: CYP2C9*2 (rs1799853, 430C>T) and CYP2C9*3 (rs1057910, 1075A>C).
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, Tuple
from config import CYP2C9_DIPLOTYPE_PHENOTYPE_MAP, CYP2C9_PHENOTYPE_ORDER, CYP2C9_CLINICAL_GUIDELINES


class CYP2C9Translator:
    @staticmethod
    def assign_sample_diplotype(row: pd.Series) -> Tuple[str, str]:
        """
        Translates a single sample's SNP genotypes (CYP2C9*2, CYP2C9*3) 
        into a Diplotype (*1/*1, *1/*2, *1/*3, *2/*2, *2/*3, *3/*3) and CPIC Metabolizer Phenotype.
        """
        # Retrieve values searching for standard names or flexible patterns
        g_cyp2 = None
        g_cyp3 = None
        for col in row.index:
            c_norm = str(col).lower().replace(" ", "").replace("_", "").replace("-", "")
            if "*2" in c_norm or "rs1799853" in c_norm or "cyp2c9_2" in c_norm or "cyp2c9*2" in c_norm:
                v = row.get(col)
                if pd.notna(v) and str(v).strip().upper() not in ['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN']:
                    g_cyp2 = str(v).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
            elif "*3" in c_norm or "rs1057910" in c_norm or "cyp2c9_3" in c_norm or "cyp2c9*3" in c_norm:
                v = row.get(col)
                if pd.notna(v) and str(v).strip().upper() not in ['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN']:
                    g_cyp3 = str(v).strip().upper().replace('/', '').replace('|', '').replace(' ', '')

        # Check missingness
        if not g_cyp2 and not g_cyp3:
            return 'Indeterminate', 'Indeterminate'

        # *2 (rs1799853): Ref = C, Var = T. CC = 0 copies, CT = 1 copy, TT = 2 copies
        star2_count = 0
        if g_cyp2 == 'TT':
            star2_count = 2
        elif g_cyp2 in ('CT', 'TC'):
            star2_count = 1
        elif g_cyp2 == 'CC':
            star2_count = 0

        # *3 (rs1057910): Ref = A, Var = C. AA = 0 copies, AC = 1 copy, CC = 2 copies
        star3_count = 0
        if g_cyp3 == 'CC':
            star3_count = 2
        elif g_cyp3 in ('AC', 'CA'):
            star3_count = 1
        elif g_cyp3 == 'AA':
            star3_count = 0

        # CPIC Diplotype Calling
        if star2_count == 0 and star3_count == 0:
            diplotype = '*1/*1'
        elif star2_count == 1 and star3_count == 0:
            diplotype = '*1/*2'
        elif star2_count == 0 and star3_count == 1:
            diplotype = '*1/*3'
        elif star2_count == 2 and star3_count == 0:
            diplotype = '*2/*2'
        elif star2_count == 1 and star3_count == 1:
            diplotype = '*2/*3'
        elif star2_count == 0 and star3_count == 2:
            diplotype = '*3/*3'
        elif star2_count == 2 and star3_count == 1:
            diplotype = '*2/*3'
        elif star2_count == 1 and star3_count == 2:
            diplotype = '*2/*3'
        elif star2_count == 2 and star3_count == 2:
            diplotype = '*2/*3'
        else:
            diplotype = '*1/*1'

        phenotype = CYP2C9_DIPLOTYPE_PHENOTYPE_MAP.get(diplotype, 'Indeterminate')
        return diplotype, phenotype

    @classmethod
    def process_dataframe(cls, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Applies CYP2C9 diplotype and phenotype classification across the entire dataset.
        Computes diplotype & phenotype counts, frequencies, and percentages.
        """
        res_df = df.copy()

        diplotypes = []
        phenotypes = []

        for _, row in res_df.iterrows():
            d, p = cls.assign_sample_diplotype(row)
            diplotypes.append(d)
            phenotypes.append(p)

        res_df['Diplotype'] = diplotypes
        res_df['Phenotype'] = phenotypes

        # Compute population summary
        total_valid = len(res_df[res_df['Phenotype'] != 'Indeterminate'])
        total_samples = len(res_df)

        # Diplotype Summary
        dip_counts = res_df[res_df['Diplotype'] != 'Indeterminate']['Diplotype'].value_counts().to_dict()
        dip_summary = {}
        denom = total_valid if total_valid > 0 else total_samples
        for dip, count in sorted(dip_counts.items(), key=lambda x: x[1], reverse=True):
            freq = count / denom if denom > 0 else 0.0
            dip_summary[dip] = {
                'count': int(count),
                'frequency': round(freq, 4),
                'percentage': round(freq * 100, 2)
            }

        # Phenotype Summary
        pheno_counts = res_df['Phenotype'].value_counts().to_dict()
        pheno_summary = {}
        for pheno in CYP2C9_PHENOTYPE_ORDER:
            count = pheno_counts.get(pheno, 0)
            freq = count / denom if (pheno != 'Indeterminate' and denom > 0) else (count / total_samples if total_samples > 0 else 0.0)
            pheno_summary[pheno] = {
                'count': int(count),
                'frequency': round(freq, 4),
                'percentage': round(freq * 100, 2)
            }

        summary = {
            'total_samples': total_samples,
            'valid_pharmacogenomic_samples': total_valid,
            'diplotypes': dip_summary,
            'phenotypes': pheno_summary,
            'gene': 'CYP2C9'
        }

        return res_df, summary
