"""Hardy-Weinberg Equilibrium analytical view component."""

from typing import Any, Callable, Dict, Tuple
import pandas as pd
import plotly.express as px
import streamlit as st


def render_hardy_weinberg(
    render_subcohort_pills: Callable[[str], Tuple[str, Dict[str, Any]]]
) -> None:
    """Render the Hardy-Weinberg Equilibrium Audit view."""
    sel_key, sel_cohort = render_subcohort_pills("tab_hwe")
    snps_data = sel_cohort.get("snps", {})

    cohort_total_n = sel_cohort.get("sample_count", 0)

    st.markdown(f"### Hardy-Weinberg Equilibrium Audit ({sel_key} Cohort)")
    st.caption(
        f"Cohort Total: **{cohort_total_n:,} samples**. HWE uses valid, called genotypes only; the decision p-value uses the exact test when expected genotype counts are small and the chi-square test otherwise."
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

        hwe_table_rows.append(
            {
                "SNP Variant": snp_name,
                "Cohort Total (N)": cohort_total_n,
                "Valid Analyzed (N)": v_n,
                "Missing Calls": f"{m_n} ({m_pct:.1f}%)"
                if m_n > 0
                else "0 (0.0%)",
                "Observed Genotypes": " / ".join(
                    [f"{k}:{v}" for k, v in gc.items()]
                ),
                "Expected Genotypes": " / ".join(
                    [f"{k}:{exp_c.get(k, 0)}" for k in gc.keys()]
                ),
                "_raw_chi2": chi2_val,
                "Degrees of Freedom": hw.get("df", 1),
                "P-Value": f"{hw.get('p_value', 1.0):.4f}",
                "Decision P-Value": f"{hw.get('decision_p_value', hw.get('p_value', 1.0)):.4f}",
                "Exact Test P-Value": f"{hw.get('exact_p_value', 1.0):.4f}",
                "HWE Interpretation": hw.get(
                    "interpretation", "Not available"
                ),
            }
        )

    c_hwe_tbl, c_hwe_fig = st.columns([1.1, 0.9])
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
            styled_hwe = hwe_df.style.map(
                style_hwe_chi2, subset=["Chi-Square (χ²)"]
            )
            st.dataframe(styled_hwe, use_container_width=True, hide_index=True)
        else:
            st.dataframe(hwe_df, use_container_width=True, hide_index=True)

    with c_hwe_fig:
        chi2_rows = []
        for snp_name, s_res in snps_data.items():
            raw_c = s_res.get("hwe", {}).get("chi2_stat", 0.0)
            try:
                c_num = (
                    float(raw_c)
                    if (
                        raw_c is not None
                        and str(raw_c).strip() not in ["", "nan", "none", "-"]
                    )
                    else 0.0
                )
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
            if {
                s_res.get("hwe", {}).get("df", 1)
                for s_res in snps_data.values()
            } == {1}:
                fig_hwe.add_hline(
                    y=3.841,
                    line_dash="dash",
                    line_color="red",
                    annotation_text="α = 0.05 (df = 1)",
                )
            else:
                st.caption(
                    "Critical values vary by degrees of freedom. Monomorphic markers have df = 0 and are reported as not informative for HWE."
                )
            fig_hwe.update_layout(
                height=260, margin=dict(l=20, r=20, t=35, b=20)
            )
            st.plotly_chart(fig_hwe, use_container_width=True)
        else:
            st.info("No Chi-Square data available to plot for this cohort.")

    with st.expander(
        "📖 Standard Chi-Square Distribution Cutoff Reference Table (Theoretical Distribution, df = 1)",
        expanded=False,
    ):
        st.caption(
            "ℹ️ **Theoretical Reference:** These critical values apply to a chi-square test with df = 1. The decision p-value shown above uses the exact HWE test when expected genotype counts are small; monomorphic markers have df = 0 and HWE is not informative."
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
                    "Significance Level (α)": "0.02 (2%)",
                    "Critical Cutoff Value (χ², df=1)": "5.412",
                    "Interpretation & Decision Rule": "⚠️ Strong HWE Departure (p < 0.02)",
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
        st.info(
            "📌 Use the displayed HWE interpretation and decision p-value for each marker. The 3.841 cutoff applies only to the asymptotic df = 1 chi-square reference."
        )
