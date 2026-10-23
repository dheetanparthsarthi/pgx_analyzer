"""Population groups and Genotype PCA analytical view component."""

from typing import Any, Dict, List
import pandas as pd
import plotly.express as px
import streamlit as st
from analysis_core import genotype_pca


def render_population_groups(
    full_results: Dict[str, Any],
    active_markers: List[str],
    active_snp_conf: Dict[str, Any],
) -> None:
    """Render population subgroup sizes and genotype PCA analysis."""
    st.markdown("### Group sizes")
    ov = full_results.get("overall", {})
    reg_dict = full_results.get("regional", {})
    gen_dict = full_results.get("gender", {})
    tot_n = ov.get("sample_count", 0)

    group_rows = [{"GROUP": "Overall", "N": tot_n, "% OF TOTAL": "100.0%"}]
    for group_name, group_data in reg_dict.items():
        group_n = group_data.get("sample_count", 0)
        group_rows.append(
            {
                "GROUP": f"Region — {group_name}",
                "N": group_n,
                "% OF TOTAL": f"{group_n / tot_n * 100:.1f}%"
                if tot_n
                else "0.0%",
            }
        )
    for group_name, group_data in gen_dict.items():
        group_n = group_data.get("sample_count", 0)
        group_rows.append(
            {
                "GROUP": f"Gender — {group_name}",
                "N": group_n,
                "% OF TOTAL": f"{group_n / tot_n * 100:.1f}%"
                if tot_n
                else "0.0%",
            }
        )

    c_pop1, c_pop2 = st.columns([1.0, 1.0])
    with c_pop1:
        st.table(pd.DataFrame(group_rows))

    with c_pop2:
        pop_chart_df = pd.DataFrame(
            [
                {"Group": r["GROUP"], "N": r["N"]}
                for r in group_rows
                if r["GROUP"] != "Overall"
            ]
        )
        fig_pop = px.bar(
            pop_chart_df,
            x="N",
            y="Group",
            orientation="h",
            title="Population Subgroup Sizes (N)",
            color="Group",
            color_discrete_sequence=px.colors.qualitative.Bold,
        )
        fig_pop.update_layout(
            height=320,
            margin=dict(l=20, r=20, t=35, b=20),
            yaxis=dict(autorange="reversed"),
        )
        st.plotly_chart(fig_pop, use_container_width=True)

    st.markdown("### 🧬 Genotype Principal Component Analysis")
    st.caption(
        "PCA summarizes genetic variation across the selected markers. Missing calls are mean-imputed for this visualization only."
    )
    pca_source = full_results.get("processed_dataframe", pd.DataFrame())
    pca_markers = [
        marker for marker in active_markers if marker in pca_source.columns
    ]
    pca_scores, pca_info = genotype_pca(
        pca_source, pca_markers, active_snp_conf
    )
    if pca_info.get("reason"):
        st.info(pca_info["reason"])
    elif not pca_scores.empty:
        pca_plot_df = pca_scores.copy()
        for metadata_col in ("Sample_ID", "Region", "Gender_Clean"):
            if metadata_col in pca_source.columns:
                pca_plot_df[metadata_col] = pca_source[metadata_col]
        color_col = next(
            (
                col
                for col in ("Region", "Gender_Clean")
                if col in pca_plot_df.columns
                and pca_plot_df[col].nunique(dropna=True) > 1
            ),
            None,
        )
        variance = pca_info["explained_variance_percent"]
        pca_plot_df["PC1"] = pca_plot_df["PC1"].astype(float)
        pca_plot_df["PC2"] = pca_plot_df["PC2"].astype(float)
        fig_pca_2d = px.scatter(
            pca_plot_df,
            x="PC1",
            y="PC2",
            color=color_col,
            hover_name="Sample_ID"
            if "Sample_ID" in pca_plot_df.columns
            else None,
            title="Genotype PCA — PC1 vs PC2",
            labels={
                "PC1": f"PC1 ({variance[0]:.1f}% variance)",
                "PC2": f"PC2 ({variance[1]:.1f}% variance)",
                "Region": "Population Region",
                "Gender_Clean": "Gender",
            },
            color_discrete_sequence=px.colors.qualitative.Safe,
        )
        fig_pca_2d.update_layout(
            height=430, margin=dict(l=20, r=20, t=55, b=20)
        )
        st.plotly_chart(fig_pca_2d, use_container_width=True)
        if len(variance) >= 3:
            pca_plot_df["PC3"] = pca_plot_df["PC3"].astype(float)
            fig_pca_3d = px.scatter_3d(
                pca_plot_df,
                x="PC1",
                y="PC2",
                z="PC3",
                color=color_col,
                hover_name="Sample_ID"
                if "Sample_ID" in pca_plot_df.columns
                else None,
                title="Genotype PCA — PC1, PC2, and PC3",
                labels={
                    "PC1": f"PC1 ({variance[0]:.1f}% variance)",
                    "PC2": f"PC2 ({variance[1]:.1f}% variance)",
                    "PC3": f"PC3 ({variance[2]:.1f}% variance)",
                    "Region": "Population Region",
                    "Gender_Clean": "Gender",
                },
                color_discrete_sequence=px.colors.qualitative.Safe,
            )
            fig_pca_3d.update_layout(
                height=500, margin=dict(l=20, r=20, t=55, b=20)
            )
            st.plotly_chart(fig_pca_3d, use_container_width=True)
        else:
            st.caption(
                "A 3D PCA view needs at least three informative components; the available data supports 2D only."
            )
        st.caption(f"Markers used: {', '.join(pca_info['markers_used'])}")
        if pca_info.get("skipped_markers"):
            st.caption(
                f"Skipped markers without usable biallelic variation: {', '.join(pca_info['skipped_markers'])}"
            )
