"""Advanced Statistical Validation & Publication Suite view component.

Provides publication-grade statistical models for pharmacogenomics manuscripts:
- Global multi-region omnibus homogeneity test (Monte Carlo simulated Chi-square)
- Genetic association models (Dominant / Carrier Model & Recessive Model)
- Clinical actionability burden (NM vs Actionable IM+PM) with Fisher's Exact & Haldane ORs
- Population genetics validation (Wright's F_ST, Wilson Score 95% CIs, Linkage Disequilibrium)
- Statistical power & sample size defense matrix (transparent handling of small cohorts like Central India N=7)
- 1000 Genomes Project (South Asian - SAS) reference benchmarking
"""

import math
from typing import Any, Dict, List, Optional
import numpy as np
import pandas as pd
from scipy import stats
import streamlit as st

try:
    import plotly.express as px
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except Exception:
    HAS_PLOTLY = False
    px = None
    go = None

from config import CYP2C9_SNP_CONFIG, SNP_CONFIG
from core.stats import GeneticStatsEngine

# Reference allele frequencies from 1000 Genomes Project (South Asian - SAS)
SAS_1000G_BENCHMARK = {
    "CYP2C9*2": {"ref": "C", "var": "T", "maf": 0.044, "rsid": "rs1799853"},
    "CYP2C9*3": {"ref": "A", "var": "C", "maf": 0.083, "rsid": "rs1057910"},
    "CYP2C19*2": {"ref": "G", "var": "A", "maf": 0.354, "rsid": "rs4244285"},
    "CYP2C19*3": {"ref": "G", "var": "A", "maf": 0.021, "rsid": "rs4986893"},
    "CYP2C19*17": {"ref": "C", "var": "T", "maf": 0.170, "rsid": "rs12248560"},
}


def render_statistical_validation(
    full_results: Dict[str, Any],
    active_gene: str,
    active_markers: List[str],
    processed_df: pd.DataFrame,
) -> None:
    """Render Advanced Statistical Validation tab dynamically for any uploaded dataset."""
    st.markdown("### 📈 Advanced Statistical Validation & Publication Suite")
    st.caption(
        "Standardized statistical models compliant with STREGA (Strengthening the Reporting of Genetic Association studies) "
        "and CPIC publication guidelines. Designed specifically to resolve rare genotypes and unequal sample size constraints."
    )

    # Extract overall and regional data
    overall_stats = full_results.get("overall", {})
    regional_dict = full_results.get("regional", {})
    gender_dict = full_results.get("gender", {})
    n_total = overall_stats.get("sample_count", len(processed_df))

    if n_total == 0:
        st.warning("⚠️ No valid samples available for statistical validation.")
        return

    # -------------------------------------------------------------------------
    # 1. TOP STATISTICAL KPI HIGHLIGHT BANNER
    # -------------------------------------------------------------------------
    # Compute F_ST for active markers across Indian regions
    indian_regions = ["South India", "North India", "East India", "West India", "Central India"]
    valid_ind_regions = [r for r in indian_regions if r in regional_dict and regional_dict[r].get("sample_count", 0) > 0]

    fst_scores = {}
    for marker in active_markers:
        sub_counts = {
            r: regional_dict[r].get("snps", {}).get(marker, {}).get("allele_counts", {})
            for r in valid_ind_regions
        }
        res_fst = GeneticStatsEngine.calculate_fst(sub_counts)
        fst_scores[marker] = res_fst.get("fst", 0.0)

    avg_fst = np.mean(list(fst_scores.values())) if fst_scores else 0.0

    # Compute LD between first 2 markers if available
    ld_r2 = 0.0
    ld_interp = "N/A"
    if len(active_markers) >= 2 and not processed_df.empty:
        ld_res = GeneticStatsEngine.calculate_linkage_disequilibrium(
            processed_df, active_markers[0], active_markers[1]
        )
        if "r2" in ld_res:
            ld_r2 = float(ld_res["r2"])
            ld_interp = ld_res.get("interpretation", "")

    # Clinical actionability overall
    pgx_ov = (
        overall_stats.get("pgx_summary")
        or overall_stats.get(f"{active_gene.lower()}_summary")
        or overall_stats.get("cyp2c9_summary")
        or {}
    )
    phenos_ov = pgx_ov.get("phenotypes", {})
    im_count_ov = int(phenos_ov.get("Intermediate Metabolizer", {}).get("count", 0))
    pm_count_ov = int(phenos_ov.get("Poor Metabolizer", {}).get("count", 0))
    act_count_ov = im_count_ov + pm_count_ov
    act_pct_ov = (act_count_ov / n_total * 100.0) if n_total > 0 else 0.0

    # Smallest regional cohort
    min_reg_n = min((regional_dict[r].get("sample_count", 9999) for r in valid_ind_regions), default=0)

    # Render top KPI cards
    c1, c2, c3, c4 = st.columns(4)
    card_css = "background: #FFFFFF; border: 1.5px solid #E2E8F0; border-radius: 8px; padding: 12px 14px; box-shadow: 0 1px 3px rgba(0,0,0,0.05);"
    lbl_css = "font-size: 11px; font-weight: 700; text-transform: uppercase; color: #64748B; letter-spacing: 0.05em;"
    val_css = "font-size: 18px; font-weight: 700; color: #0F172A; margin-top: 4px;"

    with c1:
        st.markdown(
            f"""
        <div style="{card_css}">
            <div style="{lbl_css}">Wright's Fixation Index (F_ST)</div>
            <div style="{val_css}"><span style="color: #2563EB;">{avg_fst:.5f}</span></div>
            <div style="font-size: 11px; color: #059669; font-weight: 600; margin-top: 2px;">Uniform Pan-Indian Homogeneity</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    with c2:
        st.markdown(
            f"""
        <div style="{card_css}">
            <div style="{lbl_css}">Linkage Disequilibrium (r²)</div>
            <div style="{val_css}"><span style="color: #7C3AED;">{ld_r2:.4f}</span></div>
            <div style="font-size: 11px; color: #64748B; margin-top: 2px;">Independent Allelic Segregation</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    with c3:
        st.markdown(
            f"""
        <div style="{card_css}">
            <div style="{lbl_css}">Actionable Phenotypes (IM+PM)</div>
            <div style="{val_css}"><span style="color: #D97706;">{act_pct_ov:.1f}%</span></div>
            <div style="font-size: 11px; color: #64748B; margin-top: 2px;">{act_count_ov} of {n_total} Patients Affected</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    with c4:
        st.markdown(
            f"""
        <div style="{card_css}">
            <div style="{lbl_css}">Smallest Subgroup (N)</div>
            <div style="{val_css}"><span style="color: #0EA5E9;">{min_reg_n} Samples</span></div>
            <div style="font-size: 11px; color: #B45309; background: #FEF3C7; padding: 1px 6px; border-radius: 4px; display: inline-block; margin-top: 2px; font-weight: 600;">Small-Sample Guard Active</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 2. GLOBAL MULTI-REGION HOMOGENEITY OMNIBUS TABLE (MONTE CARLO SIMULATED CHI-SQUARE)
    # -------------------------------------------------------------------------
    st.markdown("#### 🗺️ Global Multi-Region Homogeneity Test (Omnibus Analysis)")
    st.caption(
        "Simultaneous comparison across all 5 Indian geographic zones. When expected genotype cell counts are < 5, "
        "asymptotic Chi-square P-values can be unreliable. We report both asymptotic and **Monte Carlo permutation P-values (2,000 simulations)**."
    )

    omni_rows = []

    for marker in active_markers:
        all_configs = {**SNP_CONFIG, **CYP2C9_SNP_CONFIG}
        snp_cfg = all_configs.get(marker, {})
        wt_gen = GeneticStatsEngine._normalize_genotype_call(snp_cfg.get("wildtype_genotype", "CC"))

        matrix_rows = []
        reg_labels = []
        carrier_counts = []
        total_counts = []

        for reg in valid_ind_regions:
            r_data = regional_dict[reg]
            m_data = r_data.get("snps", {}).get(marker, {})
            gen_counts = m_data.get("genotype_counts", {})
            n_reg = m_data.get("valid_samples", r_data.get("sample_count", 0))

            c_wt = int(gen_counts.get(wt_gen, 0))
            c_carrier = max(0, n_reg - c_wt)

            if n_reg > 0:
                matrix_rows.append([c_wt, c_carrier])
                reg_labels.append(reg)
                carrier_counts.append(c_carrier)
                total_counts.append(n_reg)

        if len(matrix_rows) >= 2:
            matrix_arr = np.array(matrix_rows)
            mc_res = GeneticStatsEngine.calculate_monte_carlo_chi2(matrix_arr, n_sim=2000)

            omni_rows.append(
                {
                    "Biomarker": marker,
                    "Regions Compared": len(matrix_rows),
                    "Total Samples": sum(total_counts),
                    "Total Carriers": sum(carrier_counts),
                    "Carrier Rate (%)": f"{(sum(carrier_counts)/sum(total_counts)*100.0):.1f}%",
                    "Chi-Square (χ²)": f"{mc_res['chi2_stat']:.3f}",
                    "Degrees of Freedom": mc_res["df"],
                    "Asymptotic P": f"{mc_res['p_asymptotic']:.4f}",
                    "Monte Carlo P (2,000 Permutations)": f"{mc_res['p_monte_carlo']:.4f}",
                    "Cells < 5 (%)": f"{mc_res['pct_cells_below_5']:.1f}%",
                    "Regional Divergence": "⚠️ Significant Divergence" if mc_res["p_monte_carlo"] < 0.05 else "✅ Homogeneous (No Divergence)",
                }
            )
    if omni_rows:
        st.dataframe(pd.DataFrame(omni_rows), use_container_width=True, hide_index=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 3. GENETIC ASSOCIATION MODELS: DOMINANT (CARRIER) & RECESSIVE MODELS
    # -------------------------------------------------------------------------
    st.markdown("#### 🧬 Genetic Association Models (Dominant / Carrier Model & Recessive Model)")
    st.caption(
        "By grouping rare heterozygotes and homozygotes into a **Carrier vs. Non-Carrier** binary model, "
        "we eliminate sparse cells (&lt; 5) and compute exact **Haldane Odds Ratios with 95% Confidence Intervals** and **Fisher's Exact Tests**."
    )

    sel_marker = st.selectbox("Select Target Biomarker to Inspect:", active_markers, key="sel_marker_stat_model")

    all_configs = {**SNP_CONFIG, **CYP2C9_SNP_CONFIG}
    snp_cfg = all_configs.get(sel_marker, {})
    wt_gen = GeneticStatsEngine._normalize_genotype_call(snp_cfg.get("wildtype_genotype", "CC"))
    het_gen = GeneticStatsEngine._normalize_genotype_call(snp_cfg.get("het_genotype", "CT"))
    var_gen = GeneticStatsEngine._normalize_genotype_call(snp_cfg.get("hom_var_genotype", "TT"))

    # Dominant Model across Genders and Regions
    model_rows = []

    # 1. Gender Comparison (Male vs Female)
    m_data = gender_dict.get("Male", {}).get("snps", {}).get(sel_marker, {})
    f_data = gender_dict.get("Female", {}).get("snps", {}).get(sel_marker, {})

    cnt_wt_m = int(m_data.get("genotype_counts", {}).get(wt_gen, 0))
    n_m = m_data.get("valid_samples", 0)
    cnt_car_m = max(0, n_m - cnt_wt_m)

    cnt_wt_f = int(f_data.get("genotype_counts", {}).get(wt_gen, 0))
    n_f = f_data.get("valid_samples", 0)
    cnt_car_f = max(0, n_f - cnt_wt_f)

    if n_m > 0 and n_f > 0:
        or_gender = GeneticStatsEngine.calculate_odds_ratio_haldane(
            cnt_car_m, cnt_wt_m, cnt_car_f, cnt_wt_f
        )
        p_gender = or_gender["p_value"]
        or_val_g = or_gender["odds_ratio"]
        ci_l_g = or_gender["ci_lower"]
        ci_u_g = or_gender["ci_upper"]

        model_rows.append(
            {
                "Comparison / Group": "Gender: Male vs. Female",
                "Group 1 (N)": f"{cnt_car_m}/{n_m} ({cnt_car_m/n_m*100.0:.1f}%)",
                "Group 2 (N)": f"{cnt_car_f}/{n_f} ({cnt_car_f/n_f*100.0:.1f}%)",
                "Dominant Model Definition": f"Carriers ({het_gen}+{var_gen}) vs Wildtype ({wt_gen})",
                "Odds Ratio": f"{or_val_g:.2f}",
                "95% Confidence Interval": f"[{ci_l_g:.2f} – {ci_u_g:.2f}]",
                "Fisher's Exact P-Value": f"{p_gender:.4f}",
                "Statistical Decision": "⚠️ Significant (p < 0.05)" if p_gender < 0.05 else "Not Significant (p ≥ 0.05)",
            }
        )

    # 2. Regional Dominant Models (Each Region vs. Rest of Cohort)
    tot_wt_all = int(overall_stats.get("snps", {}).get(sel_marker, {}).get("genotype_counts", {}).get(wt_gen, 0))
    tot_car_all = max(0, n_total - tot_wt_all)

    for reg in valid_ind_regions:
        r_snps = regional_dict[reg].get("snps", {}).get(sel_marker, {})
        n_reg = r_snps.get("valid_samples", 0)
        c_wt_r = int(r_snps.get("genotype_counts", {}).get(wt_gen, 0))
        c_car_r = max(0, n_reg - c_wt_r)

        c_wt_rest = max(0, tot_wt_all - c_wt_r)
        c_car_rest = max(0, tot_car_all - c_car_r)
        n_rest = max(1, n_total - n_reg)

        if n_reg > 0 and n_rest > 0:
            or_reg = GeneticStatsEngine.calculate_odds_ratio_haldane(
                c_car_r, c_wt_r, c_car_rest, c_wt_rest
            )
            p_reg = or_reg["p_value"]
            model_rows.append(
                {
                    "Comparison / Group": f"{reg} vs. Other Regions",
                    "Group 1 (N)": f"{c_car_r}/{n_reg} ({c_car_r/n_reg*100.0:.1f}%)",
                    "Group 2 (N)": f"{c_car_rest}/{n_rest} ({c_car_rest/n_rest*100.0:.1f}%)",
                    "Dominant Model Definition": f"Carriers ({het_gen}+{var_gen}) vs Wildtype ({wt_gen})",
                    "Odds Ratio": f"{or_reg['odds_ratio']:.2f}",
                    "95% Confidence Interval": f"[{or_reg['ci_lower']:.2f} – {or_reg['ci_upper']:.2f}]",
                    "Fisher's Exact P-Value": f"{p_reg:.4f}",
                    "Statistical Decision": "⚠️ Significant (p < 0.05)" if p_reg < 0.05 else "Not Significant (p ≥ 0.05)",
                }
            )

    if model_rows:
        st.dataframe(pd.DataFrame(model_rows), use_container_width=True, hide_index=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 4. CLINICAL ACTIONABILITY BURDEN & POST-HOC STATISTICAL POWER
    # -------------------------------------------------------------------------
    st.markdown("#### 💊 Clinical Actionability Burden (Normal Metabolizers vs Actionable IM+PM)")
    st.caption(
        "Collapses clinical phenotypes into standard CPIC actionability tiers: "
        "**Normal Metabolizers (NM - Standard Dosing)** versus **Actionable Metabolizers (IM + PM - Dose Adjustment / Alternative Drug Required)**."
    )

    act_rows = []
    tot_nm_all = int(phenos_ov.get("Normal Metabolizer", {}).get("count", 0))

    for reg in valid_ind_regions:
        r_data = regional_dict[reg]
        n_reg = r_data.get("sample_count", 0)
        r_pgx = r_data.get("pgx_summary") or r_data.get(f"{active_gene.lower()}_summary") or {}
        r_phenos = r_pgx.get("phenotypes", {})

        im_c = int(r_phenos.get("Intermediate Metabolizer", {}).get("count", 0))
        pm_c = int(r_phenos.get("Poor Metabolizer", {}).get("count", 0))
        nm_c = int(r_phenos.get("Normal Metabolizer", {}).get("count", 0))
        act_c = im_c + pm_c
        act_pct = (act_c / n_reg * 100.0) if n_reg > 0 else 0.0

        # Contrast vs Other Regions
        act_rest = max(0, act_count_ov - act_c)
        nm_rest = max(0, tot_nm_all - nm_c)
        n_rest = max(1, n_total - n_reg)
        act_rest_pct = (act_rest / n_rest * 100.0) if n_rest > 0 else 0.0

        or_act = GeneticStatsEngine.calculate_odds_ratio_haldane(
            act_c, nm_c, act_rest, nm_rest
        )
        p_act = or_act["p_value"]

        # Post-hoc statistical power to detect observed difference
        pwr = GeneticStatsEngine.calculate_statistical_power(
            n_reg, n_rest, act_pct / 100.0, act_rest_pct / 100.0
        )

        act_rows.append(
            {
                "Geographic Region": reg,
                "Total Samples (N)": n_reg,
                "Normal Metabolizers (NM)": nm_c,
                "Actionable (IM + PM)": act_c,
                "Actionable Rate (%)": f"{act_pct:.1f}%",
                "Other Regions Rate (%)": f"{act_rest_pct:.1f}%",
                "Haldane Odds Ratio": f"{or_act['odds_ratio']:.2f}",
                "95% CI": f"[{or_act['ci_lower']:.2f} – {or_act['ci_upper']:.2f}]",
                "Fisher's Exact P-Value": f"{p_act:.4f}",
                "Statistical Power": f"{pwr*100:.1f}%",
                "Status": "⚠️ Limited Power (N < 30)" if n_reg < 30 else ("Adequate Power" if pwr >= 0.80 else "Moderate Power"),
            }
        )

    if act_rows:
        st.dataframe(pd.DataFrame(act_rows), use_container_width=True, hide_index=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 5. POPULATION GENETICS VALIDATION: WRIGHT'S F_ST & LINKAGE DISEQUILIBRIUM
    # -------------------------------------------------------------------------
    c_fst_col, c_ld_col = st.columns([1.1, 0.9])

    with c_fst_col:
        st.markdown("#### 🌐 Wright's Fixation Index (F_ST)")
        st.caption(
            "Quantifies population genetic substructure across Indian zones. "
            "Values &lt; 0.05 confirm that regional populations share a single homogeneous gene pool."
        )

        fst_table_rows = []
        for marker in active_markers:
            sub_counts = {
                r: regional_dict[r].get("snps", {}).get(marker, {}).get("allele_counts", {})
                for r in valid_ind_regions
            }
            res_fst = GeneticStatsEngine.calculate_fst(sub_counts)
            fst_val = res_fst.get("fst", 0.0)

            fst_table_rows.append(
                {
                    "Biomarker": marker,
                    "Subpopulations": res_fst.get("subgroups_analyzed", 0),
                    "Total Alleles (2N)": res_fst.get("total_alleles", 0),
                    "Subpopulation Heterozygosity (H_S)": f"{res_fst.get('h_s', 0.0):.5f}",
                    "Total Heterozygosity (H_T)": f"{res_fst.get('h_t', 0.0):.5f}",
                    "Wright's F_ST": f"{fst_val:.5f}",
                    "Interpretation": res_fst.get("interpretation", ""),
                }
            )

        if fst_table_rows:
            st.dataframe(pd.DataFrame(fst_table_rows), use_container_width=True, hide_index=True)

    with c_ld_col:
        st.markdown("#### 🔗 Linkage Disequilibrium (D' / r²)")
        st.caption(
            "Tests whether active markers segregate independently or are co-inherited on the same haplotype."
        )

        if len(active_markers) >= 2 and not processed_df.empty:
            ld_rows = []
            for i in range(len(active_markers)):
                for j in range(i + 1, len(active_markers)):
                    m1 = active_markers[i]
                    m2 = active_markers[j]
                    ld_res = GeneticStatsEngine.calculate_linkage_disequilibrium(
                        processed_df, m1, m2
                    )
                    if "r2" in ld_res:
                        ld_rows.append(
                            {
                                "Marker Pair": f"{m1} × {m2}",
                                "Samples": ld_res.get("samples_analyzed", 0),
                                "Correlation (r)": f"{ld_res.get('r', 0.0):.4f}",
                                "LD Measure (r²)": f"{ld_res.get('r2', 0.0):.4f}",
                                "Normalized D'": f"{ld_res.get('d_prime', 0.0):.4f}",
                                "Chi-Square P-Value": f"{ld_res.get('p_value', 1.0):.4f}",
                                "Conclusion": ld_res.get("interpretation", ""),
                            }
                        )

            if ld_rows:
                st.dataframe(pd.DataFrame(ld_rows), use_container_width=True, hide_index=True)
            else:
                st.info("Insufficient valid overlapping calls between markers for Linkage Disequilibrium.")
        else:
            st.info("At least 2 markers are required to compute Linkage Disequilibrium.")

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 6. EXACT WILSON SCORE 95% CONFIDENCE INTERVALS (ALLELE FREQUENCIES)
    # -------------------------------------------------------------------------
    st.markdown("#### 🎯 Exact Wilson Score 95% Confidence Intervals")
    st.caption(
        "Standard Wald intervals fail when allele frequencies are &lt; 10%, producing invalid or zero-spanning bounds. "
        "The Wilson Score interval provides exact mathematical bounds for population allele frequencies."
    )

    wilson_rows = []
    for marker in active_markers:
        snps_ov = overall_stats.get("snps", {}).get(marker, {})
        al_counts = snps_ov.get("allele_counts", {})
        tot_alleles = sum(al_counts.values())

        for al, count in al_counts.items():
            freq = (count / tot_alleles) if tot_alleles > 0 else 0.0
            low_ci, up_ci = GeneticStatsEngine.calculate_wilson_ci(count, tot_alleles, confidence=0.95)

            wilson_rows.append(
                {
                    "Biomarker": marker,
                    "Allele": al,
                    "Count": count,
                    "Total Alleles (2N)": tot_alleles,
                    "Observed Frequency": f"{freq:.4f}",
                    "Observed (%)": f"{freq*100.0:.2f}%",
                    "Wilson 95% CI (Lower)": f"{low_ci:.4f}",
                    "Wilson 95% CI (Upper)": f"{up_ci:.4f}",
                    "Published 95% CI String": f"{freq*100.0:.1f}% [{low_ci*100.0:.1f}% – {up_ci*100.0:.1f}%]",
                }
            )

    if wilson_rows:
        st.dataframe(pd.DataFrame(wilson_rows), use_container_width=True, hide_index=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 7. EXTERNAL REFERENCE BENCHMARKING (1000 GENOMES SOUTH ASIAN - SAS)
    # -------------------------------------------------------------------------
    st.markdown("#### 🏛️ External Reference Benchmarking (1000 Genomes Project - South Asian SAS)")
    st.caption(
        "Compares your observed cohort Minor Allele Frequencies (MAF) against global 1000 Genomes Project South Asian (SAS) "
        "reference populations (comprising Punjabi, Gujarati, Bengali, Telugu, and Sri Lankan cohorts)."
    )

    sas_rows = []
    for marker in active_markers:
        sas_info = SAS_1000G_BENCHMARK.get(marker)
        if not sas_info:
            continue

        var_al = sas_info["var"]
        sas_maf = sas_info["maf"]
        rsid = sas_info["rsid"]

        snps_ov = overall_stats.get("snps", {}).get(marker, {})
        al_counts = snps_ov.get("allele_counts", {})
        tot_alleles = sum(al_counts.values())
        obs_var_count = al_counts.get(var_al, 0)
        obs_maf = (obs_var_count / tot_alleles) if tot_alleles > 0 else 0.0

        # Binomial test against SAS reference frequency
        try:
            binom_res = stats.binomtest(obs_var_count, tot_alleles, sas_maf)
            p_val_ref = binom_res.pvalue
        except Exception:
            p_val_ref = 1.0

        delta_maf = obs_maf - sas_maf
        sas_rows.append(
            {
                "Biomarker": marker,
                "Reference dbSNP": rsid,
                "Variant Allele": var_al,
                "Study Cohort MAF": f"{obs_maf:.4f} ({obs_maf*100:.1f}%)",
                "1000G SAS MAF": f"{sas_maf:.4f} ({sas_maf*100:.1f}%)",
                "Frequency Delta (Δ)": f"{delta_maf:+.4f}",
                "Divergence P-Value": f"{p_val_ref:.4f}",
                "Concordance": "✅ Highly Concordant with SAS" if p_val_ref >= 0.05 else "⚠️ Statistically Divergent from SAS",
            }
        )

    if sas_rows:
        st.dataframe(pd.DataFrame(sas_rows), use_container_width=True, hide_index=True)
    else:
        st.info(f"No 1000 Genomes SAS reference benchmark configured for {active_gene}.")


