"""Gender-Wise Population Pharmacogenomics analytical view component."""

from typing import Any, Dict, List
import pandas as pd
import plotly.express as px
import streamlit as st


def render_gender_analysis(
    active_gene: str,
    active_markers: List[str],
    full_results: Dict[str, Any],
) -> None:
    """Render gender-wise population pharmacogenomics analysis."""
    st.markdown(
        f"### 👫 Gender-Wise Population Pharmacogenomics Analysis ({active_gene})"
    )
    gender_data = full_results.get("gender", {})
    ov = full_results.get("overall", {})
    tot_count = ov.get("sample_count", 0)
    gender_summary = pd.DataFrame(
        [
            {
                "Population Group": group,
                "Sample Count": group_data.get("sample_count", 0),
                "Cohort Percentage": f"{group_data.get('sample_count', 0) / tot_count * 100:.1f}%"
                if tot_count
                else "0.0%",
            }
            for group, group_data in gender_data.items()
        ]
    )
    st.dataframe(gender_summary, use_container_width=True, hide_index=True)

    pheno_order = list(
        dict.fromkeys(
            phenotype
            for group_data in gender_data.values()
            for phenotype in (
                (
                    group_data.get("pgx_summary")
                    or group_data.get(f"{active_gene.lower()}_summary")
                    or {}
                ).get("phenotypes", {})
            )
            if phenotype != "Not configured"
        )
    )
    pheno_chart_df = []
    for gender_name, group_data in gender_data.items():
        group_phenotypes = (
            group_data.get("pgx_summary")
            or group_data.get(f"{active_gene.lower()}_summary")
            or {}
        ).get("phenotypes", {})
        for p_name in pheno_order:
            metric = group_phenotypes.get(p_name, {})
            if metric.get("count", 0) > 0:
                pheno_chart_df.append(
                    {
                        "Population Group": gender_name,
                        "Phenotype": p_name.split(" (")[0],
                        "Count": metric["count"],
                    }
                )

    st.markdown("<br>", unsafe_allow_html=True)
    cg1, cg2 = st.columns(2)
    with cg1:
        if pheno_chart_df:
            fig_g2 = px.bar(
                pd.DataFrame(pheno_chart_df),
                x="Phenotype",
                y="Count",
                color="Population Group",
                barmode="group",
                title=f"Phenotype Comparison by Population Group ({active_gene})",
                color_discrete_sequence=px.colors.qualitative.Safe,
            )
            fig_g2.update_layout(
                height=280, margin=dict(l=20, r=20, t=35, b=20)
            )
            st.plotly_chart(fig_g2, use_container_width=True)
        else:
            configured_by_group = [
                (
                    group_data.get("pgx_summary")
                    or group_data.get(f"{active_gene.lower()}_summary")
                    or {}
                ).get("translation_configured", False)
                for group_data in gender_data.values()
            ]
            if not any(configured_by_group):
                st.info(
                    "No gender-stratified phenotype calls are available because this gene and marker panel have no configured translation rules. Configure validated diplotype and phenotype rules in the Gene Calling tab; marker genotype and allele comparisons remain available here."
                )
            else:
                st.info(
                    "Translation rules are configured, but no gender-group samples matched a determinate diplotype pattern. Check marker coverage and missing or invalid genotype calls."
                )

    with cg2:
        sel_gender_marker = (
            st.selectbox(
                "Select Biomarker to Compare:",
                active_markers,
                key="gender_biomarker_sel",
            )
            if len(active_markers) > 1
            else (active_markers[0] if active_markers else "Biomarker")
        )
        gen_gt_df = []
        for gender_name, group_data in gender_data.items():
            group_genotypes = (
                group_data.get("snps", {})
                .get(sel_gender_marker, {})
                .get("genotype_counts", {})
            )
            for genotype, count in group_genotypes.items():
                if count:
                    gen_gt_df.append(
                        {
                            "Population Group": gender_name,
                            "Genotype": genotype,
                            "Count": count,
                        }
                    )
        if gen_gt_df:
            fig_gen_gt = px.bar(
                pd.DataFrame(gen_gt_df),
                x="Genotype",
                y="Count",
                color="Population Group",
                barmode="group",
                title=f"{sel_gender_marker} Genotypes by Population Group",
                color_discrete_sequence=px.colors.qualitative.Safe,
            )
            fig_gen_gt.update_layout(
                height=280, margin=dict(l=20, r=20, t=35, b=20)
            )
            st.plotly_chart(fig_gen_gt, use_container_width=True)
        else:
            st.info(f"No genotype calls to compare for {sel_gender_marker}.")
