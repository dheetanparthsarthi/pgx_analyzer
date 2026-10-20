"""
Population Genetics Engine for Genotype, Allele, and Hardy-Weinberg Equilibrium Analysis.
Implements Chi-square goodness-of-fit test and exact HWE calculation.
"""

import pandas as pd
import numpy as np
from scipy import stats
from typing import Dict, Any, Optional
import math

from config import SNP_CONFIG, CYP2C9_SNP_CONFIG, ALPHA_SIGNIFICANCE


class GeneticStatsEngine:
    @staticmethod
    def _normalize_genotype_call(value: Any) -> str:
        """Normalize a diploid call without assuming alleles are one character."""
        if pd.isna(value):
            return ""
        return str(value).strip().upper().replace("/", "").replace("|", "").replace(" ", "")

    @staticmethod
    def exact_hwe_pvalue(n_AA: int, n_Aa: int, n_aa: int) -> float:
        """
        Calculates Haldane / Wigginton Exact HWE P-value for biallelic SNPs.
        Useful when sample size is small or expected counts are < 5.
        """
        n = n_AA + n_Aa + n_aa
        if n == 0:
            return 1.0

        n_A = 2 * n_AA + n_Aa
        n_a = 2 * n_aa + n_Aa
        
        # Possible count of heterozygotes with same allele counts
        # n_Aa must have same parity as n_A (or n_a)
        start_het = n_A % 2
        
        # Calculate probability distribution of het count
        def het_prob(h):
            # h is number of heterozygotes
            n_AA_h = (n_A - h) // 2
            n_aa_h = (n_a - h) // 2
            if n_AA_h < 0 or n_aa_h < 0:
                return 0.0
            # log probability via log gamma/factorial: (N! * 2^h / (n_AA! * h! * n_aa!)) * (n_A! * n_a! / (2N)!)
            log_p = (math.lgamma(n + 1) + math.lgamma(n_A + 1) + math.lgamma(n_a + 1) + (h * math.log(2.0))
                     - (math.lgamma(2*n + 1) + math.lgamma(n_AA_h + 1) + math.lgamma(h + 1) + math.lgamma(n_aa_h + 1)))
            return math.exp(log_p)

        observed_prob = het_prob(n_Aa)
        
        p_val = 0.0
        max_het = min(n_A, n_a)
        for h in range(start_het, max_het + 1, 2):
            prob = het_prob(h)
            if prob <= observed_prob + 1e-12:
                p_val += prob
                
        return min(1.0, p_val)

    @classmethod
    def analyze_snp(cls, series: pd.Series, snp_name: str, custom_info: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Performs full Genotype, Allele, and HWE analysis for a single SNP column.
        Supports configured and generic/custom biallelic markers.
        """
        clean_series = series.dropna().map(cls._normalize_genotype_call)
        clean_series = clean_series[~clean_series.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN', 'NOTPERFORMED'])]
        
        all_configs = {**SNP_CONFIG, **CYP2C9_SNP_CONFIG}
        snp_info = custom_info or all_configs.get(snp_name, None)
        
        if not snp_info:
            # Dynamically infer alleles if not pre-configured
            counts = clean_series.value_counts()
            alleles = sorted({a for call in counts.index.tolist() if len(call) == 2 for a in call if a.isalpha()})
            ref = alleles[0] if len(alleles) > 0 else 'A'
            var = alleles[1] if len(alleles) > 1 else 'G'
            snp_info = {
                'ref_allele': ref,
                'var_allele': var,
                'wildtype_genotype': f"{ref}{ref}",
                'het_genotype': f"{ref}{var}",
                'hom_var_genotype': f"{var}{var}"
            }
        
        ref_allele = snp_info.get('ref_allele')
        var_allele = snp_info.get('var_allele')
        if not ref_allele or not var_allele:
            counts = clean_series.value_counts()
            alleles = sorted({a for call in counts.index.tolist() if len(call) == 2 for a in call if a.isalpha()})
            ref_allele = ref_allele or (alleles[0] if len(alleles) > 0 else 'A')
            var_allele = var_allele or (alleles[1] if len(alleles) > 1 else 'G')

        ref_allele = str(ref_allele).strip().upper()
        var_allele = str(var_allele).strip().upper()
        wt_gen = cls._normalize_genotype_call(snp_info.get('wildtype_genotype') or f"{ref_allele}{ref_allele}")
        het_gen = cls._normalize_genotype_call(snp_info.get('het_genotype') or f"{ref_allele}{var_allele}")
        var_gen = cls._normalize_genotype_call(snp_info.get('hom_var_genotype') or f"{var_allele}{var_allele}")

        ref_calls = {cls._normalize_genotype_call(f"{ref_allele}{ref_allele}"), wt_gen}
        het_calls = {cls._normalize_genotype_call(f"{ref_allele}{var_allele}"), cls._normalize_genotype_call(f"{var_allele}{ref_allele}"), het_gen}
        var_calls = {cls._normalize_genotype_call(f"{var_allele}{var_allele}"), var_gen}

        # Count only calls that match the configured biallelic genotype pairs.
        # This supports multi-character allele labels such as DEL and INS.
        call_counts = clean_series.value_counts()
        c_wt = int(sum(call_counts.get(call, 0) for call in ref_calls))
        c_het = int(sum(call_counts.get(call, 0) for call in (het_calls - ref_calls - var_calls)))
        c_var = int(sum(call_counts.get(call, 0) for call in (var_calls - ref_calls)))
        n_valid = c_wt + c_het + c_var
        if n_valid == 0:
            return {
                'snp_name': snp_name,
                'valid_samples': 0,
                'genotype_counts': {wt_gen: 0, het_gen: 0, var_gen: 0},
                'genotype_freqs': {wt_gen: 0.0, het_gen: 0.0, var_gen: 0.0},
                'genotype_pcts': {wt_gen: 0.0, het_gen: 0.0, var_gen: 0.0},
                'allele_counts': {ref_allele: 0, var_allele: 0},
                'allele_freqs': {ref_allele: 0.0, var_allele: 0.0},
                'allele_pcts': {ref_allele: 0.0, var_allele: 0.0},
                'hwe': {
                    'expected_counts': {wt_gen: 0.0, het_gen: 0.0, var_gen: 0.0},
                    'expected_freqs': {wt_gen: 0.0, het_gen: 0.0, var_gen: 0.0},
                    'chi2_stat': 0.0,
                    'df': 1,
                    'p_value': 1.0,
                    'exact_p_value': 1.0,
                    'interpretation': 'No Data'
                }
            }
        
        f_wt = c_wt / n_valid
        f_het = c_het / n_valid
        f_var = c_var / n_valid
        
        # 2. Allele Counts
        c_ref_allele = (2 * c_wt) + c_het
        c_var_allele = (2 * c_var) + c_het
        total_alleles = 2 * n_valid
        
        p = c_ref_allele / total_alleles if total_alleles > 0 else 0.0
        q = c_var_allele / total_alleles if total_alleles > 0 else 0.0
        
        # 3. Hardy-Weinberg Equilibrium
        exp_wt = n_valid * (p ** 2)
        exp_het = n_valid * (2 * p * q)
        exp_var = n_valid * (q ** 2)
        
        # Chi-square calculation with safety for division by zero
        chi2_stat = 0.0
        if exp_wt > 0:  chi2_stat += ((c_wt - exp_wt) ** 2) / exp_wt
        if exp_het > 0: chi2_stat += ((c_het - exp_het) ** 2) / exp_het
        if exp_var > 0: chi2_stat += ((c_var - exp_var) ** 2) / exp_var
        
        p_val = stats.chi2.sf(chi2_stat, df=1) if chi2_stat > 0 else 1.0
        exact_p_val = cls.exact_hwe_pvalue(c_wt, c_het, c_var)
        
        # Interpretation
        final_p = exact_p_val if (exp_wt < 5 or exp_het < 5 or exp_var < 5) else p_val
        interpretation = "In HWE" if final_p >= ALPHA_SIGNIFICANCE else "Departure from HWE"
        
        return {
            'snp_name': snp_name,
            'valid_samples': n_valid,
            'genotype_counts': {wt_gen: c_wt, het_gen: c_het, var_gen: c_var},
            'genotype_freqs': {wt_gen: round(f_wt, 4), het_gen: round(f_het, 4), var_gen: round(f_var, 4)},
            'genotype_pcts': {wt_gen: round(f_wt * 100, 2), het_gen: round(f_het * 100, 2), var_gen: round(f_var * 100, 2)},
            'allele_counts': {ref_allele: c_ref_allele, var_allele: c_var_allele},
            'allele_freqs': {ref_allele: round(p, 4), var_allele: round(q, 4)},
            'allele_pcts': {ref_allele: round(p * 100, 2), var_allele: round(q * 100, 2)},
            'hwe': {
                'expected_counts': {wt_gen: round(exp_wt, 2), het_gen: round(exp_het, 2), var_gen: round(exp_var, 2)},
                'expected_freqs': {wt_gen: round(p ** 2, 4), het_gen: round(2 * p * q, 4), var_gen: round(q ** 2, 4)},
                'chi2_stat': round(chi2_stat, 4),
                'df': 1,
                'p_value': round(p_val, 4),
                'exact_p_value': round(exact_p_val, 4),
                'interpretation': interpretation
            }
        }
