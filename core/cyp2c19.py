"""
CYP2C19 Star Allele, Diplotype Call, and Phenotype Classifier.
Implements standard CPIC rules for CYP2C19 pharmacogenomics.
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, Tuple
from config import DIPLOTYPE_PHENOTYPE_MAP, PHENOTYPE_ORDER


class CYP2C19Translator:
    @staticmethod
    def assign_sample_diplotype(row: pd.Series) -> Tuple[str, str]:
        """
        Translates a single sample's SNP genotypes (CYP2C19*2, *3, *17) 
        into a Diplotype (*1/*1, *1/*2, etc.) and CPIC Metabolizer Phenotype.
        """
        g_cyp2 = str(row.get('CYP2C19*2', '')).strip().upper() if pd.notna(row.get('CYP2C19*2')) else None
        g_cyp3 = str(row.get('CYP2C19*3', '')).strip().upper() if pd.notna(row.get('CYP2C19*3')) else None
        g_cyp17 = str(row.get('CYP2C19*17', '')).strip().upper() if pd.notna(row.get('CYP2C19*17')) else None

        # Check missingness
        if not g_cyp2 and not g_cyp3 and not g_cyp17:
            return 'Indeterminate', 'Indeterminate'
            
        # Count variant alleles for *2 (GA=1, AA=2), *3 (GA=1, AA=2), *17 (CT=1, TT=2)
        star2_count = 2 if g_cyp2 == 'AA' else (1 if g_cyp2 == 'GA' else 0)
        star3_count = 2 if g_cyp3 == 'AA' else (1 if g_cyp3 == 'GA' else 0)
        star17_count = 2 if g_cyp17 == 'TT' else (1 if g_cyp17 == 'CT' else 0)
        
        # Diplotype calling rules
        if star2_count == 2:
            diplotype = '*2/*2'
        elif star3_count == 2:
            diplotype = '*3/*3'
        elif star17_count == 2:
            diplotype = '*17/*17'
        elif star2_count == 1 and star3_count == 1:
            diplotype = '*2/*3'
        elif star2_count == 1 and star17_count == 1:
            diplotype = '*2/*17'
        elif star3_count == 1 and star17_count == 1:
            diplotype = '*3/*17'
        elif star2_count == 1:
            diplotype = '*1/*2'
        elif star3_count == 1:
            diplotype = '*1/*3'
        elif star17_count == 1:
            diplotype = '*1/*17'
        else:
            diplotype = '*1/*1'

        phenotype = DIPLOTYPE_PHENOTYPE_MAP.get(diplotype, 'Indeterminate')
        return diplotype, phenotype

    @classmethod
    def process_dataframe(cls, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Applies diplotype and phenotype classification across the entire dataset.
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
        dip_counts = res_df['Diplotype'].value_counts().to_dict()
        dip_summary = {}
        for dip, count in dip_counts.items():
            freq = count / total_samples if total_samples > 0 else 0.0
            dip_summary[dip] = {
                'count': int(count),
                'frequency': round(freq, 4),
                'percentage': round(freq * 100, 2)
            }
            
        # Phenotype Summary
        pheno_counts = res_df['Phenotype'].value_counts().to_dict()
        pheno_summary = {}
        for pheno in PHENOTYPE_ORDER:
            count = pheno_counts.get(pheno, 0)
            freq = count / total_samples if total_samples > 0 else 0.0
            pheno_summary[pheno] = {
                'count': int(count),
                'frequency': round(freq, 4),
                'percentage': round(freq * 100, 2)
            }
            
        summary = {
            'total_samples': total_samples,
            'valid_pharmacogenomic_samples': total_valid,
            'diplotypes': dip_summary,
            'phenotypes': pheno_summary
        }
        
        return res_df, summary
