"""Dedicated State Profiles analytical view component."""

from typing import Any, Dict, List
import pandas as pd
import plotly.express as px
import streamlit as st


def render_state_profiles(
    active_gene: str,
    active_markers: List[str],
    active_snp_conf: Dict[str, Any],
    full_results: Dict[str, Any],
    processed_df: pd.DataFrame,
) -> None:
    """Render dedicated individual state profiles, marker statistics, and state-wise matrices."""
    st.markdown(f"### 📍 Dedicated Individual State Profiles ({active_gene})")
    state_dict = full_results.get("state_wise", {})
    state_names = sorted(list(state_dict.keys()))

    if state_names:
        selected_st = st.selectbox(
            "Select Dedicated State Profile:",
            state_names,
            index=0,
            key="st_tab_selectbox",
        )
        st_data = state_dict[selected_st]

        if "Unclassified" in selected_st:
            st.info(
                "ℹ️ **Unclassified / Missing State Cohort**: These samples have missing or unmapped state information in your dataset. "
                "Their genotypes, allele frequencies, and phenotypes are preserved and displayed here for full visibility."
            )

        primary_marker = (
            st.selectbox(
                "Select marker for this state profile:",
                active_markers,
                key="state_profile_marker",
            )
            if active_markers
            else "Biomarker"
        )
        m_conf = active_snp_conf.get(primary_marker, {})
        var_al = m_conf.get("var_allele")
        al_freqs = (
            st_data.get("snps", {})
            .get(primary_marker, {})
            .get("allele_freqs", {})
        )
        if not var_al:
            ref_al = m_conf.get("ref_allele")
            non_refs = [a for a in al_freqs.keys() if a != ref_al]
            var_al = (
                non_refs[0]
                if non_refs
                else (list(al_freqs.keys())[0] if al_freqs else "Var")
            )
        marker_freq = al_freqs.get(
            var_al, st_data.get("snps", {}).get(primary_marker, {}).get("maf", 0)
        )

        pgx_st = (
            st_data.get("pgx_summary")
            or st_data.get(f"{active_gene.lower()}_summary")
            or st_data.get("cyp2c19_summary", {})
        )
        state_phenotypes = {
            k: v
            for k, v in pgx_st.get("phenotypes", {}).items()
            if k != "Not configured"
        }
        top_state_phenotype = (
            max(
                state_phenotypes.items(),
                key=lambda item: item[1].get("count", 0),
            )
            if state_phenotypes
            else (None, {})
        )

        sk1, sk2, sk3, sk4 = st.columns(4)
        sk1.metric(
            "State Cohort Size (N)",
            f"{st_data.get('sample_count', 0):,} samples",
            help=f"Total number of patients/samples from {selected_st} present in this cohort",
        )
        sk2.metric(
            "Missing Data Rate",
            f"{st_data.get('missing_pct', 0):.1f}%",
            help=f"Percentage of uncalled / missing genotypes in {selected_st}",
        )
        allele_metric_label = (
            f"{primary_marker} Single Allele Frequency ({var_al})"
            if var_al == m_conf.get("ref_allele")
            else f"{primary_marker} Variant Allele Frequency ({var_al})"
        )
        sk3.metric(
            allele_metric_label,
            f"{marker_freq:.4f}"
            if isinstance(marker_freq, (int, float))
            else str(marker_freq),
            help=f"Frequency of configured allele {var_al} in {selected_st}",
        )
        sk4.metric(
            "Most Frequent Phenotype",
            f"{top_state_phenotype[0]} ({top_state_phenotype[1].get('percentage', 0):.1f}%)"
            if top_state_phenotype[0]
            else "Rules needed",
        )

        st_phenos = {
            k: v
            for k, v in pgx_st.get("phenotypes", {}).items()
            if k != "Not configured"
        }
        st_snps = (
            st_data.get("snps", {})
            .get(primary_marker, {})
            .get("genotype_counts", {})
        )

        cs1, cs2 = st.columns(2)
        with cs1:
            st_gt_df = pd.DataFrame(
                [{"Genotype": k, "Count": v} for k, v in st_snps.items()]
            )
            if not st_gt_df.empty:
                fig_st_gt = px.bar(
                    st_gt_df,
                    x="Genotype",
                    y="Count",
                    title=f"{selected_st} — {primary_marker} Genotypes",
                    color="Genotype",
                    color_discrete_sequence=[
                        "#16A34A",
                        "#2563EB",
                        "#DC2626",
                        "#D97706",
                    ],
                )
                fig_st_gt.update_layout(
                    height=260, margin=dict(l=20, r=20, t=35, b=20)
                )
                st.plotly_chart(fig_st_gt, use_container_width=True)
            else:
                st.info(f"No {primary_marker} genotype data for {selected_st}.")

        with cs2:
            st_phe_df = pd.DataFrame(
                [
                    {
                        "Phenotype": k.split(" (")[0],
                        "Count": v["count"],
                    }
                    for k, v in st_phenos.items()
                    if v["count"] > 0
                ]
            )
            if not st_phe_df.empty:
                fig_st_phe = px.pie(
                    st_phe_df,
                    names="Phenotype",
                    values="Count",
                    hole=0.45,
                    title=f"{selected_st} — {active_gene} Metabolizer Phenotypes",
                    color_discrete_sequence=[
                        "#166534",
                        "#1E40AF",
                        "#4B5563",
                        "#991B1B",
                        "#D97706",
                        "#0D9488",
                    ],
                )
                fig_st_phe.update_layout(
                    height=260, margin=dict(l=20, r=20, t=35, b=20)
                )
                st.plotly_chart(fig_st_phe, use_container_width=True)
            else:
                st.info(f"No {active_gene} phenotype calls for {selected_st}.")

        state_marker_rows = []
        for marker in active_markers:
            marker_stats = st_data.get("snps", {}).get(marker, {})
            hwe_stats = marker_stats.get("hwe", {})
            state_marker_rows.append(
                {
                    "Marker": marker,
                    "Valid N": marker_stats.get("valid_samples", 0),
                    "Genotype Counts": ", ".join(
                        f"{gt}: {count}"
                        for gt, count in marker_stats.get(
                            "genotype_counts", {}
                        ).items()
                    ),
                    "Allele Counts": ", ".join(
                        f"{allele}: {count}"
                        for allele, count in marker_stats.get(
                            "allele_counts", {}
                        ).items()
                    ),
                    "HWE χ²": hwe_stats.get("chi2_stat", 0.0),
                    "HWE df": hwe_stats.get("df", 1),
                    "Asymptotic p-value": hwe_stats.get("p_value", 1.0),
                    "Decision p-value": hwe_stats.get(
                        "decision_p_value", hwe_stats.get("p_value", 1.0)
                    ),
                    "HWE status": hwe_stats.get(
                        "interpretation", "Not available"
                    ),
                }
            )
        st.markdown(
            f"#### Marker Statistics for {selected_st} (All Selected Markers)"
        )
        st.dataframe(
            pd.DataFrame(state_marker_rows),
            use_container_width=True,
            hide_index=True,
        )

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(
            f"#### 💊 State-Wise Metabolizer Phenotype Distribution Matrix ({active_gene})"
        )
        st.caption(
            f"Observed {active_gene} phenotype categories across the states present in this dataset."
        )

        pheno_order = list(
            dict.fromkeys(
                name
                for group in state_dict.values()
                for name in (
                    (
                        group.get("pgx_summary")
                        or group.get(f"{active_gene.lower()}_summary")
                        or {}
                    ).get("phenotypes", {})
                )
                if name != "Not configured"
            )
        )
        st_dict = full_results.get("state_wise", {})
        state_pheno_matrix = []
        for s_name, s_data in sorted(
            st_dict.items(),
            key=lambda x: x[1].get("sample_count", 0),
            reverse=True,
        ):
            if s_data.get("sample_count", 0) > 0:
                p_dict = (
                    s_data.get("pgx_summary")
                    or s_data.get(f"{active_gene.lower()}_summary")
                    or s_data.get("cyp2c19_summary", {})
                ).get("phenotypes", {})
                row = {
                    "State": s_name,
                    "Total Sample N": s_data.get("sample_count", 0),
                }
                for p_cat in pheno_order:
                    cnt = p_dict.get(p_cat, {}).get("count", 0)
                    pct = p_dict.get(p_cat, {}).get("percentage", 0.0)
                    p_label = p_cat.split(" (")[0]
                    row[p_label] = f"{cnt} ({pct:.1f}%)"
                state_pheno_matrix.append(row)

        ov = full_results.get("overall", {})
        tot_p = (
            ov.get("pgx_summary")
            or ov.get(f"{active_gene.lower()}_summary")
            or ov.get("cyp2c19_summary", {})
        ).get("phenotypes", {})
        tot_row = {
            "State": "Total Cohort (All States Combined)",
            "Total Sample N": ov.get("sample_count", 0),
        }
        for p_cat in pheno_order:
            cnt = tot_p.get(p_cat, {}).get("count", 0)
            pct = tot_p.get(p_cat, {}).get("percentage", 0.0)
            p_label = p_cat.split(" (")[0]
            tot_row[p_label] = f"{cnt} ({pct:.1f}%)"
        state_pheno_matrix.append(tot_row)

        st.dataframe(
            pd.DataFrame(state_pheno_matrix),
            use_container_width=True,
            hide_index=True,
        )

        st.markdown(f"#### Diplotype Distribution by State ({active_gene})")
        state_diplotype_categories = list(
            dict.fromkeys(
                name
                for group in state_dict.values()
                for name in (
                    (
                        group.get("pgx_summary")
                        or group.get(f"{active_gene.lower()}_summary")
                        or {}
                    ).get("diplotypes", {})
                )
                if name != "Not configured"
            )
        )
        state_diplotype_matrix = []
        for s_name, s_data in sorted(
            state_dict.items(),
            key=lambda item: item[1].get("sample_count", 0),
            reverse=True,
        ):
            dip_data = (
                s_data.get("pgx_summary")
                or s_data.get(f"{active_gene.lower()}_summary")
                or {}
            ).get("diplotypes", {})
            row = {
                "State": s_name,
                "Total Sample N": s_data.get("sample_count", 0),
            }
            for category in state_diplotype_categories:
                metric = dip_data.get(category, {})
                row[category] = (
                    f"{metric.get('count', 0)} ({metric.get('percentage', 0.0):.1f}%)"
                )
            state_diplotype_matrix.append(row)
        st.dataframe(
            pd.DataFrame(state_diplotype_matrix),
            use_container_width=True,
            hide_index=True,
        )
    else:
        mapped_state = (
            "State_Clean" in processed_df.columns
            and processed_df["State_Clean"]
            .astype(str)
            .str.strip()
            .ne("")
            .any()
        )
        if mapped_state:
            st.info(
                "No state profiles could be summarized from the mapped state values. Check that the selected state/locality column contains non-empty values."
            )
        else:
            st.info(
                "No state or locality column is currently mapped. Map a state/locality column during Upload or Map Columns, then recalculate to populate individual profiles. Geographic-group summaries are available in the Geographic Groups tab when a region column is mapped."
            )
