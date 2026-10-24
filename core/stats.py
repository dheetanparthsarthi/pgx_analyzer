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

        # An explicitly identical reference/variant label means this cohort's
        # marker is being analyzed as monomorphic. Report its observed genotype
        # and allele once (dict keys cannot represent duplicate allele labels),
        # and do not present a meaningless biallelic HWE test.
        if ref_allele == var_allele:
            mono_gen = cls._normalize_genotype_call(f"{ref_allele}{ref_allele}")
            mono_count = int(clean_series.eq(mono_gen).sum())
            if mono_count == 0:
                return {
                    'snp_name': snp_name,
                    'valid_samples': 0,
                    'genotype_counts': {mono_gen: 0},
                    'genotype_freqs': {mono_gen: 0.0},
                    'genotype_pcts': {mono_gen: 0.0},
                    'allele_counts': {ref_allele: 0},
                    'allele_freqs': {ref_allele: 0.0},
                    'allele_pcts': {ref_allele: 0.0},
                    'hwe': {
                        'expected_counts': {mono_gen: 0.0},
                        'expected_freqs': {mono_gen: 0.0},
                        'chi2_stat': 0.0,
                        'df': 0,
                        'p_value': 1.0,
                        'decision_p_value': 1.0,
                        'exact_p_value': 1.0,
                        'interpretation': 'Monomorphic; HWE not informative'
                    }
                }
            return {
                'snp_name': snp_name,
                'valid_samples': mono_count,
                'genotype_counts': {mono_gen: mono_count},
                'genotype_freqs': {mono_gen: 1.0},
                'genotype_pcts': {mono_gen: 100.0},
                'allele_counts': {ref_allele: 2 * mono_count},
                'allele_freqs': {ref_allele: 1.0},
                'allele_pcts': {ref_allele: 100.0},
                'hwe': {
                    'expected_counts': {mono_gen: float(mono_count)},
                    'expected_freqs': {mono_gen: 1.0},
                    'chi2_stat': 0.0,
                    'df': 0,
                    'p_value': 1.0,
                    'decision_p_value': 1.0,
                    'exact_p_value': 1.0,
                    'interpretation': 'Monomorphic; HWE not informative'
                }
            }

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
                    'decision_p_value': 1.0,
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
                'decision_p_value': round(final_p, 4),
                'exact_p_value': round(exact_p_val, 4),
                'interpretation': interpretation
            }
        }

    @staticmethod
    def calculate_odds_ratio_haldane(a: int, b: int, c: int, d: int) -> Dict[str, Any]:
        """
        Calculate Haldane-Anscombe corrected Odds Ratio, 95% Confidence Interval (Woolf logit),
        and Fisher's Exact Test two-sided P-value for a 2x2 contingency table:
            [[a, b],
             [c, d]]
        Haldane correction adds 0.5 to all cells if any cell is 0 to avoid division by zero
        and reduce small-sample bias.
        """
        table_raw = [[a, b], [c, d]]
        try:
            _, fisher_p = stats.fisher_exact(table_raw)
        except Exception:
            fisher_p = 1.0

        haldane_applied = (a == 0 or b == 0 or c == 0 or d == 0)
        adj_a = a + (0.5 if haldane_applied else 0.0)
        adj_b = b + (0.5 if haldane_applied else 0.0)
        adj_c = c + (0.5 if haldane_applied else 0.0)
        adj_d = d + (0.5 if haldane_applied else 0.0)

        if adj_b * adj_c == 0:
            or_val = 1.0
            se_ln_or = 0.0
            ci_lower = 0.0
            ci_upper = float('inf')
        else:
            or_val = (adj_a * adj_d) / (adj_b * adj_c)
            se_ln_or = math.sqrt((1.0 / adj_a) + (1.0 / adj_b) + (1.0 / adj_c) + (1.0 / adj_d))
            ci_lower = math.exp(math.log(or_val) - 1.96 * se_ln_or)
            ci_upper = math.exp(math.log(or_val) + 1.96 * se_ln_or)

        return {
            'odds_ratio': round(or_val, 3),
            'ci_lower': round(ci_lower, 3),
            'ci_upper': round(ci_upper, 3),
            'p_value': round(fisher_p, 4),
            'se': round(se_ln_or, 4),
            'haldane_applied': haldane_applied,
            'table': table_raw
        }

    @staticmethod
    def calculate_wilson_ci(k: int, n: int, confidence: float = 0.95) -> Tuple[float, float]:
        """
        Calculate Wilson Score confidence interval for a binomial proportion k/n.
        Essential for rare allele frequencies where Wald normal approximation fails.
        """
        if n <= 0:
            return 0.0, 0.0
        p = k / n
        z = stats.norm.ppf(1.0 - (1.0 - confidence) / 2.0)
        z2 = z ** 2
        denominator = 1.0 + (z2 / n)
        center = (p + (z2 / (2.0 * n))) / denominator
        half_width = (z / denominator) * math.sqrt((p * (1.0 - p) / n) + (z2 / (4.0 * (n ** 2))))
        lower = max(0.0, center - half_width)
        upper = min(1.0, center + half_width)
        return round(lower, 4), round(upper, 4)

    @staticmethod
    def calculate_fst(subgroup_counts: Dict[str, Dict[str, int]]) -> Dict[str, Any]:
        """
        Calculate Wright's Fixation Index (F_ST) across subpopulations for a biallelic marker.
        subgroup_counts: {subgroup_name: {ref_allele: count, var_allele: count}}
        F_ST = (H_T - H_S) / H_T, measuring population genetic substructure.
        F_ST < 0.05 indicates little genetic differentiation (pan-Indian homogeneity).
        """
        valid_groups = {
            k: v for k, v in subgroup_counts.items()
            if sum(v.values()) > 0
        }
        if len(valid_groups) < 2:
            return {'fst': 0.0, 'interpretation': 'Insufficient subgroups (>=2 required)', 'h_s': 0.0, 'h_t': 0.0}

        total_var = 0
        total_alleles = 0
        sub_props = []
        sub_weights = []

        for grp, counts in valid_groups.items():
            ref_c = list(counts.values())[0] if len(counts) > 0 else 0
            var_c = list(counts.values())[1] if len(counts) > 1 else 0
            tot_sub = ref_c + var_c
            if tot_sub > 0:
                p_sub = var_c / tot_sub
                sub_props.append(p_sub)
                sub_weights.append(tot_sub)
                total_var += var_c
                total_alleles += tot_sub

        if total_alleles == 0:
            return {'fst': 0.0, 'interpretation': 'No allele data', 'h_s': 0.0, 'h_t': 0.0}

        # Subpopulation expected heterozygosity H_S = sum(w_i * 2 * p_i * (1 - p_i))
        h_s = sum(
            (w / total_alleles) * (2.0 * p_i * (1.0 - p_i))
            for p_i, w in zip(sub_props, sub_weights)
        )

        # Total expected heterozygosity H_T = 2 * p_total * (1 - p_total)
        p_total = total_var / total_alleles
        h_t = 2.0 * p_total * (1.0 - p_total)

        if h_t <= 0.0:
            fst = 0.0
            interp = "Monomorphic marker across all populations"
        else:
            fst = max(0.0, (h_t - h_s) / h_t)
            if fst < 0.005:
                interp = "Very low differentiation (<0.005) - Uniform pan-population genetic homogeneity"
            elif fst < 0.05:
                interp = "Low genetic differentiation (<0.05) - Consistent with single gene pool"
            elif fst < 0.15:
                interp = "Moderate genetic differentiation (0.05 - 0.15)"
            else:
                interp = "High genetic differentiation (>=0.15)"

        return {
            'fst': round(fst, 5),
            'h_s': round(h_s, 5),
            'h_t': round(h_t, 5),
            'total_alleles': total_alleles,
            'subgroups_analyzed': len(valid_groups),
            'interpretation': interp
        }

    @staticmethod
    def calculate_statistical_power(n1: int, n2: int, p1: float, p2: float, alpha: float = 0.05) -> float:
        """
        Calculate post-hoc statistical power for comparing two proportions (two-sided z-test).
        Used to assess whether small cohorts (e.g. Central India N=7) had adequate power.
        """
        if n1 <= 0 or n2 <= 0:
            return 0.0
        if p1 == p2:
            return alpha

        p_pool = (n1 * p1 + n2 * p2) / (n1 + n2)
        se0 = math.sqrt(max(1e-9, p_pool * (1.0 - p_pool) * ((1.0 / n1) + (1.0 / n2))))
        se1 = math.sqrt(max(1e-9, (p1 * (1.0 - p1) / n1) + (p2 * (1.0 - p2) / n2)))

        z_alpha = stats.norm.ppf(1.0 - alpha / 2.0)
        z_beta = (abs(p1 - p2) - z_alpha * se0) / se1
        power = float(stats.norm.cdf(z_beta))
        return round(max(0.0, min(1.0, power)), 4)

    @classmethod
    def calculate_linkage_disequilibrium(
        cls,
        df: pd.DataFrame,
        marker1: str,
        marker2: str,
        custom_info1: Optional[Dict] = None,
        custom_info2: Optional[Dict] = None
    ) -> Dict[str, Any]:
        """
        Calculate Linkage Disequilibrium (D, D', r, r^2) between two biallelic markers.
        Demonstrates whether markers (e.g. CYP2C9*2 and *3) segregate independently.
        """
        if marker1 not in df.columns or marker2 not in df.columns:
            return {'error': f"Markers {marker1} or {marker2} not found in dataset."}

        s1 = df[marker1].dropna().map(cls._normalize_genotype_call)
        s2 = df[marker2].dropna().map(cls._normalize_genotype_call)

        valid_idx = s1[~s1.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING'])].index.intersection(
            s2[~s2.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING'])].index
        )
        if len(valid_idx) == 0:
            return {'error': "No overlapping valid samples between markers."}

        # Encode dosage (0 = homozygous ref, 1 = het, 2 = hom var)
        all_configs = {**SNP_CONFIG, **CYP2C9_SNP_CONFIG}
        info1 = custom_info1 or all_configs.get(marker1, {})
        info2 = custom_info2 or all_configs.get(marker2, {})

        wt1 = cls._normalize_genotype_call(info1.get('wildtype_genotype', ''))
        het1 = cls._normalize_genotype_call(info1.get('het_genotype', ''))
        var1 = cls._normalize_genotype_call(info1.get('hom_var_genotype', ''))

        wt2 = cls._normalize_genotype_call(info2.get('wildtype_genotype', ''))
        het2 = cls._normalize_genotype_call(info2.get('het_genotype', ''))
        var2 = cls._normalize_genotype_call(info2.get('hom_var_genotype', ''))

        def get_dosage(call, wt, het, var):
            if call == wt: return 0
            if call == het: return 1
            if call == var: return 2
            return None

        dosages1 = []
        dosages2 = []
        for idx in valid_idx:
            d1 = get_dosage(s1.loc[idx], wt1, het1, var1)
            d2 = get_dosage(s2.loc[idx], wt2, het2, var2)
            if d1 is not None and d2 is not None:
                dosages1.append(d1)
                dosages2.append(d2)

        n_samples = len(dosages1)
        if n_samples < 5:
            return {'error': "Insufficient valid sample pairs (<5)."}

        # Pearson correlation r between dosage vectors
        arr1 = np.array(dosages1, dtype=float)
        arr2 = np.array(dosages2, dtype=float)

        std1 = np.std(arr1)
        std2 = np.std(arr2)

        if std1 == 0 or std2 == 0:
            r = 0.0
            r2 = 0.0
            chi2_stat = 0.0
            p_val = 1.0
            d_val = 0.0
            d_prime = 0.0
        else:
            r, _ = stats.pearsonr(arr1, arr2)
            r2 = r ** 2
            chi2_stat = n_samples * r2
            p_val = stats.chi2.sf(chi2_stat, df=1)

            # Allele frequencies
            p_a = np.mean(arr1) / 2.0
            p_b = np.mean(arr2) / 2.0
            # Covariance approximation for D: Cov(X1, X2) / 2
            cov = np.cov(arr1, arr2)[0, 1]
            d_val = cov / 2.0

            # D_prime = D / D_max
            if d_val >= 0:
                d_max = min(p_a * (1.0 - p_b), (1.0 - p_a) * p_b)
            else:
                d_max = min(p_a * p_b, (1.0 - p_a) * (1.0 - p_b))

            d_prime = (d_val / d_max) if d_max > 0 else 0.0
            d_prime = max(-1.0, min(1.0, d_prime))

        return {
            'marker1': marker1,
            'marker2': marker2,
            'samples_analyzed': n_samples,
            'r': round(r, 4),
            'r2': round(r2, 4),
            'd': round(d_val, 5),
            'd_prime': round(d_prime, 4),
            'chi2_stat': round(chi2_stat, 4),
            'p_value': round(p_val, 4),
            'interpretation': (
                "Independent segregation (No Linkage Disequilibrium, r^2 ~ 0)"
                if r2 < 0.1 else (
                    "Moderate Linkage Disequilibrium"
                    if r2 < 0.8 else "Strong Linkage Disequilibrium (Co-inherited)"
                )
            )
        }

    @staticmethod
    def calculate_monte_carlo_chi2(contingency_table: np.ndarray, n_sim: int = 2000, seed: int = 42) -> Dict[str, Any]:
        """
        Monte Carlo simulated Chi-square test for R x C contingency tables with sparse cells (<5).
        Permutes row and column marginals to estimate the empirical null distribution.
        """
        table = np.array(contingency_table, dtype=float)
        row_sums = table.sum(axis=1)
        col_sums = table.sum(axis=0)
        total = table.sum()

        if total == 0 or len(row_sums) < 2 or len(col_sums) < 2:
            return {'chi2': 0.0, 'p_asymptotic': 1.0, 'p_monte_carlo': 1.0, 'cochran_warning': False}

        # Expected counts under independence
        expected = np.outer(row_sums, col_sums) / total
        low_cells = (expected < 5).sum()
        cochran_violated = (low_cells / expected.size) > 0.20

        # Observed Chi2 statistic
        with np.errstate(divide='ignore', invalid='ignore'):
            chi2_cells = np.where(expected > 0, ((table - expected) ** 2) / expected, 0.0)
            obs_chi2 = np.sum(chi2_cells)

        dof = (table.shape[0] - 1) * (table.shape[1] - 1)
        p_asymp = stats.chi2.sf(obs_chi2, df=dof) if dof > 0 else 1.0

        # Monte Carlo permutation
        rng = np.random.default_rng(seed)
        sim_chi2 = []
        p_col = col_sums / total
        for _ in range(n_sim):
            sim_table = np.zeros_like(table)
            for i, r_sum in enumerate(row_sums):
                if r_sum > 0:
                    sim_table[i] = rng.multinomial(int(r_sum), p_col)
            with np.errstate(divide='ignore', invalid='ignore'):
                sim_stat = np.sum(np.where(expected > 0, ((sim_table - expected) ** 2) / expected, 0.0))
                sim_chi2.append(sim_stat)

        p_mc = float(np.mean(np.array(sim_chi2) >= obs_chi2))

        return {
            'chi2_stat': round(float(obs_chi2), 4),
            'df': int(dof),
            'p_asymptotic': round(float(p_asymp), 4),
            'p_monte_carlo': round(float(p_mc), 4),
            'expected_counts': np.round(expected, 2).tolist(),
            'cells_below_5': int(low_cells),
            'pct_cells_below_5': round(float((low_cells / expected.size) * 100), 1),
            'cochran_violated': cochran_violated
        }

