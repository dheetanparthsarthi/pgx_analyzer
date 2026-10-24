"""Hardy-Weinberg Equilibrium analytical view component.

Provides comprehensive HWE testing:
1. Haldane / Wigginton Exact Test (gold standard for rare alleles < 5)
2. Inbreeding Coefficient (F_IS / Wright's Fixation Index)
3. Cross-Cohort Consistency Matrix (testing across 5 Indian regions and sexes for the Wahlund Effect)
"""

from typing import Any, Callable, Dict, Optional, Tuple
import pandas as pd
import plotly.express as px
import streamlit as st


def render_hardy_weinberg(
    render_subcohort_pills: Callable[[str], Tuple[str, Dict[str, Any]]],
    full_results: Optional[Dict[str, Any]] = None,
) -> None:
    """Render the Hardy-Weinberg Equilibrium Audit view with Inbreeding Coefficients and Cross-Cohort testing."""
    sel_key, sel_cohort = render_subcohort_pills("tab_hwe")
    snps_data = sel_cohort.get("snps", {})

    cohort_total_n = sel_cohort.get("sample_count", 0)

    st.markdown(f"### ⚖️ Hardy-Weinberg Equilibrium Audit ({sel_key} Cohort)")
    st.caption(
        f"Cohort Total: **{cohort_total_n:,} samples**. HWE uses valid genotyped calls. "
        "The Decision P-Value uses the **Haldane/Wigginton Exact Test** when expected genotype counts are < 5 "
        "and standard Chi-square (χ²) otherwise. Inbreeding coefficients (F_IS) quantify heterozygote excess or deficiency."
    )

    hwe_table_rows = []
    for snp_name, s_res in snps_data.items():
        hw = s_res.get("hwe", {})
        gc = s_res.get("genotype_counts", {})
        exp_c = hw.get("expected_counts", {})
        chi2_val = hw.get("chi2_stat", 0.0)
        v_n = s_res.get("valid_samples", 0)
        m_n = max(0, cohort_total_n - v_n)
        m_pct = (m_n / cohort_total_n * 100.0) if cohort_total_n > 0 else 0.0
        f_is = hw.get("inbreeding_coeff", 0.0)

        hwe_table_rows.append(
            {
                "SNP Variant": snp_name,
                "Cohort Total (N)": cohort_total_n,
                "Valid Analyzed (N)": v_n,
                "Missing Calls": f"{m_n} ({m_pct:.1f}%)" if m_n > 0 else "0 (0.0%)",
                "Observed Genotypes": " / ".join([f"{k}:{v}" for k, v in gc.items()]),
                "Expected Genotypes": " / ".join([f"{k}:{exp_c.get(k, 0)}" for k in gc.keys()]),
                "_raw_chi2": chi2_val,
                "Degrees of Freedom": hw.get("df", 1),
                "P-Value": f"{hw.get('p_value', 1.0):.4f}",
                "Decision P-Value": f"{hw.get('decision_p_value', hw.get('p_value', 1.0)):.4f}",
                "Exact Test P-Value": f"{hw.get('exact_p_value', 1.0):.4f}",
                "Inbreeding Coeff (F_IS)": f"{f_is:+.4f}",
                "HWE Interpretation": hw.get("interpretation", "Not available"),
            }
        )

    c_hwe_tbl, c_hwe_fig = st.columns([1.15, 0.85])
    with c_hwe_tbl:
        HWE_HEADERS = [
            "SNP Variant",
            "Cohort Total (N)",
            "Valid Analyzed (N)",
            "Missing Calls",
            "Observed Genotypes",
            "Expected Genotypes",
            "Chi-Square (χ²)",
            "df",
            "P-Value",
            "Decision P-Value",
            "Exact P-Val",
            "Inbreeding (F_IS)",
            "HWE Interpretation",
        ]
        if not hwe_table_rows:
            hwe_df = pd.DataFrame(columns=HWE_HEADERS)
        else:
            hwe_df = pd.DataFrame(
                [
                    {
                        "SNP Variant": r["SNP Variant"],
                        "Cohort Total (N)": r["Cohort Total (N)"],
                        "Valid Analyzed (N)": r["Valid Analyzed (N)"],
                        "Missing Calls": r["Missing Calls"],
                        "Observed Genotypes": r["Observed Genotypes"],
                        "Expected Genotypes": r["Expected Genotypes"],
                        "Chi-Square (χ²)": f"⚠️ {r['_raw_chi2']:.4f}"
                        if r["HWE Interpretation"] == "Departure from HWE"
                        else f"{r['_raw_chi2']:.4f}",
                        "df": r["Degrees of Freedom"],
                        "P-Value": r["P-Value"],
                        "Decision P-Value": r["Decision P-Value"],
                        "Exact P-Val": r["Exact Test P-Value"],
                        "Inbreeding (F_IS)": r["Inbreeding Coeff (F_IS)"],
                        "HWE Interpretation": r["HWE Interpretation"],
                    }
                    for r in hwe_table_rows
                ]
            )

        def style_hwe_chi2(val: Any) -> str:
            if "⚠️" in str(val):
                return "background-color: #FEE2E2; color: #DC2626; font-weight: bold; border-radius: 4px;"
            return "color: #166534; font-weight: 600;"

        if not hwe_df.empty and "Chi-Square (χ²)" in hwe_df.columns:
            styled_hwe = hwe_df.style.map(style_hwe_chi2, subset=["Chi-Square (χ²)"])
            st.dataframe(styled_hwe, use_container_width=True, hide_index=True)
        else:
            st.dataframe(hwe_df, use_container_width=True, hide_index=True)

    with c_hwe_fig:
        chi2_rows = []
        for snp_name, s_res in snps_data.items():
            raw_c = s_res.get("hwe", {}).get("chi2_stat", 0.0)
            try:
                c_num = float(raw_c) if (raw_c is not None and str(raw_c).strip() not in ["", "nan", "none", "-"]) else 0.0
            except (ValueError, TypeError):
                c_num = 0.0
            chi2_rows.append({"SNP": str(snp_name), "Chi2 Stat": c_num})

        chi2_chart_df = pd.DataFrame(chi2_rows, columns=["SNP", "Chi2 Stat"])
        if not chi2_chart_df.empty and len(chi2_rows) > 0:
            fig_hwe = px.bar(
                chi2_chart_df,
                x="SNP",
                y="Chi2 Stat",
                title=f"HWE Chi-Square (χ²) by Marker ({sel_key})",
                color="SNP",
                color_discrete_sequence=["#2563EB", "#7C3AED", "#DB2777"],
            )
            if {s_res.get("hwe", {}).get("df", 1) for s_res in snps_data.values()} == {1}:
                fig_hwe.add_hline(
                    y=3.841,
                    line_dash="dash",
                    line_color="red",
                    annotation_text="α = 0.05 (df = 1)",
                )
            fig_hwe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_hwe, use_container_width=True)
        else:
            st.info("No Chi-Square data available to plot for this cohort.")

    # -------------------------------------------------------------------------
    # 2. CROSS-COHORT HWE CONSISTENCY & WAHLUND EFFECT MATRIX
    # -------------------------------------------------------------------------
    if full_results:
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("#### 🌐 Cross-Cohort HWE Consistency & Wahlund Effect Matrix")
        st.caption(
            "Tests whether HWE holds consistently within each individual geographic region and sex. "
            "When individual sub-cohorts are in HWE (p ≥ 0.05) and inbreeding coefficients (F_IS) are ~ 0, "
            "this rules out the **Wahlund Effect** (hidden population stratification), confirming that the cohort is free from genetic artifacts."
        )

        subgroup_list = []
        # Overall
        ov_data = full_results.get("overall", {})
        if ov_data.get("sample_count", 0) > 0:
            subgroup_list.append(("Overall Population", "National", ov_data))

        # Regions
        reg_dict = full_results.get("regional", {})
        for r_name in ["South India", "North India", "East India", "West India", "Central India"]:
            if r_name in reg_dict and reg_dict[r_name].get("sample_count", 0) > 0:
                subgroup_list.append((r_name, "Geographic Region", reg_dict[r_name]))

        # Genders
        gen_dict = full_results.get("gender", {})
        for g_name in ["Male", "Female"]:
            if g_name in gen_dict and gen_dict[g_name].get("sample_count", 0) > 0:
                subgroup_list.append((f"Gender: {g_name}", "Sex-Stratified", gen_dict[g_name]))

        wahlund_rows = []
        for sub_name, sub_type, sub_data in subgroup_list:
            sub_n = sub_data.get("sample_count", 0)
            sub_snps = sub_data.get("snps", {})
            for m_name, m_res in sub_snps.items():
                m_hw = m_res.get("hwe", {})
                m_gc = m_res.get("genotype_counts", {})
                m_exp = m_hw.get("expected_counts", {})
                f_is = m_hw.get("inbreeding_coeff", 0.0)
                exact_p = m_hw.get("exact_p_value", 1.0)
                chi2_v = m_hw.get("chi2_stat", 0.0)

                wahlund_rows.append(
                    {
                        "Subpopulation Cohort": sub_name,
                        "Cohort Level": sub_type,
                        "Biomarker": m_name,
                        "Sample N": sub_n,
                        "Observed Calls": " / ".join([f"{k}:{v}" for k, v in m_gc.items()]),
                        "Expected Calls": " / ".join([f"{k}:{m_exp.get(k, 0)}" for k in m_gc.keys()]),
                        "Chi-Square (χ²)": f"{chi2_v:.3f}",
                        "Exact P-Value": f"{exact_p:.4f}",
                        "Inbreeding (F_IS)": f"{f_is:+.4f}",
                        "HWE Status": "✅ In HWE" if exact_p >= 0.05 else "⚠️ Departure (p < 0.05)",
                    }
                )

        if wahlund_rows:
            st.dataframe(pd.DataFrame(wahlund_rows), use_container_width=True, hide_index=True)
            st.markdown(
                """
            <div style="background: #F8FAFC; border-left: 4px solid #10B981; padding: 10px 14px; border-radius: 4px; font-size: 12px; color: #065F46; margin-top: 6px;">
                <b>💡 STREGA Population Genetics Insight:</b><br>
                When inbreeding coefficients (F_IS) are close to zero across all regional zones and biological sexes, 
                it mathematically confirms that <b>random mating (panmixia)</b> applies across India for this gene, with no evidence of endogamous inbreeding distortions or sex-biased genotyping error.
            </div>
            """,
                unsafe_allow_html=True,
            )

    # Reference cutoffs expander
    with st.expander(
        "📖 Standard Chi-Square Distribution Cutoff Reference Table (Theoretical Distribution, df = 1)",
        expanded=False,
    ):
        st.caption(
            "ℹ️ **Theoretical Reference:** These critical values apply to a chi-square test with df = 1. The decision p-value uses the exact HWE test when expected genotype counts are small."
        )
        chi2_ref_df = pd.DataFrame(
            [
                {
                    "Significance Level (α)": "0.50 (50%)",
                    "Critical Cutoff Value (χ², df=1)": "0.455",
                    "Interpretation & Decision Rule": "Consistent with HWE (Random variation)",
                },
                {
                    "Significance Level (α)": "0.10 (10%)",
                    "Critical Cutoff Value (χ², df=1)": "2.706",
                    "Interpretation & Decision Rule": "Non-significant variation (In HWE)",
                },
                {
                    "Significance Level (α)": "0.05 (5%) [Standard HWE Cutoff]",
                    "Critical Cutoff Value (χ², df=1)": "3.841",
                    "Interpretation & Decision Rule": "⚠️ Statistically Significant HWE Departure (p < 0.05)",
                },
                {
                    "Significance Level (α)": "0.01 (1%)",
                    "Critical Cutoff Value (χ², df=1)": "6.635",
                    "Interpretation & Decision Rule": "⚠️ Very Strong HWE Departure (p < 0.01)",
                },
                {
                    "Significance Level (α)": "0.001 (0.1%)",
                    "Critical Cutoff Value (χ², df=1)": "10.828",
                    "Interpretation & Decision Rule": "⚠️ Extremely Significant HWE Departure (p < 0.001)",
                },
            ]
        )
        st.dataframe(chi2_ref_df, use_container_width=True, hide_index=True)
