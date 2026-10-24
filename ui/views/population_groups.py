"""Population Stratification & Genotype PCA interactive view component.

Provides interactive demographic and genetic stratification options:
1. Subgroup Sample Sizes & Proportions (Region, State, Gender)
2. Interactive Cross-Tabulation Matrix (Region x Phenotype, Region x Gender, Gender x Phenotype)
3. Genotype PCA Clustering with interactive color-by dimensions (Region, Gender, Diplotype, Phenotype)
"""

from typing import Any, Dict, List
import numpy as np
import pandas as pd
import plotly.express as px
from scipy import stats
import streamlit as st
from core.pca import genotype_pca


def render_population_groups(
    full_results: Dict[str, Any],
    active_markers: List[str],
    active_snp_conf: Dict[str, Any],
) -> None:
    """Render population stratification with interactive view options and PCA grouping controls."""
    st.markdown("### 🌐 Population Stratification Analysis")
    st.caption(
        "Analyze demographic stratification, multi-variable cross-tabulations, "
        "and genetic clustering across geographic regions, states, and genders."
    )

    ov = full_results.get("overall", {})
    reg_dict = full_results.get("regional", {})
    gen_dict = full_results.get("gender", {})
    state_dict = full_results.get("state_wise", {})
    proc_df = full_results.get("processed_dataframe", pd.DataFrame())
    tot_n = ov.get("sample_count", len(proc_df))

    if tot_n == 0:
        st.warning("⚠️ No valid samples available for population stratification.")
        return

    # -------------------------------------------------------------------------
    # 0. INTERACTIVE STRATIFICATION VIEW SELECTOR
    # -------------------------------------------------------------------------
    strat_view_mode = st.radio(
        "Select Population Stratification View:",
        [
            "📊 Subgroup Sizes & Proportions",
            "🔀 Multi-Variable Cross-Tabulation",
            "🔬 Genotype PCA Clustering",
        ],
        horizontal=True,
        key="pop_strat_view_mode_radio",
    )
    st.markdown("<div style='margin-bottom: 12px;'></div>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # MODE 1: SUBGROUP SIZES & PROPORTIONS
    # -------------------------------------------------------------------------
    if strat_view_mode == "📊 Subgroup Sizes & Proportions":
        st.markdown("#### 👥 Subgroup Sample Sizes & Cohort Representation")
        st.caption("Distribution of total analyzed cohort across demographic dimensions.")

        group_rows = [{"GROUP": "Overall Population", "N": tot_n, "% OF TOTAL": "100.0%"}]
        for group_name, group_data in reg_dict.items():
            group_n = group_data.get("sample_count", 0)
            if group_n > 0:
                group_rows.append(
                    {
                        "GROUP": f"Region — {group_name}",
                        "N": group_n,
                        "% OF TOTAL": f"{group_n / tot_n * 100:.1f}%",
                    }
                )
        for group_name, group_data in gen_dict.items():
            group_n = group_data.get("sample_count", 0)
            if group_n > 0:
                group_rows.append(
                    {
                        "GROUP": f"Gender — {group_name}",
                        "N": group_n,
                        "% OF TOTAL": f"{group_n / tot_n * 100:.1f}%",
                    }
                )

        c_pop1, c_pop2 = st.columns([1.0, 1.0])
        with c_pop1:
            st.dataframe(pd.DataFrame(group_rows), use_container_width=True, hide_index=True)

        with c_pop2:
            pop_chart_df = pd.DataFrame(
                [
                    {"Group": r["GROUP"], "N": r["N"]}
                    for r in group_rows
                    if r["GROUP"] != "Overall Population"
                ]
            )
            if not pop_chart_df.empty:
                fig_pop = px.bar(
                    pop_chart_df,
                    x="N",
                    y="Group",
                    orientation="h",
                    title="Population Subgroup Representation (N)",
                    color="Group",
                    color_discrete_sequence=px.colors.qualitative.Bold,
                )
                fig_pop.update_layout(
                    height=320,
                    margin=dict(l=20, r=20, t=35, b=20),
                    yaxis=dict(autorange="reversed"),
                    showlegend=False,
                )
                st.plotly_chart(fig_pop, use_container_width=True)

    # -------------------------------------------------------------------------
    # MODE 2: MULTI-VARIABLE CROSS-TABULATION MATRIX
    # -------------------------------------------------------------------------
    elif strat_view_mode == "🔀 Multi-Variable Cross-Tabulation":
        st.markdown("#### 🔀 Demographic × Clinical Phenotype Cross-Tabulation")
        st.caption(
            "Select demographic and clinical dimensions to explore multi-variable cross-stratified distributions with automated contingency Chi-square testing."
        )

        c_opt1, c_opt2 = st.columns([1.2, 0.8])
        with c_opt1:
            cross_mode = st.selectbox(
                "Select Stratification Cross-Tabulation:",
                [
                    "🗺️ Region × 💊 Clinical Metabolizer Phenotype",
                    "🗺️ Region × 👫 Gender",
                    "👫 Gender × 💊 Clinical Metabolizer Phenotype",
                    "📍 Top Indian States × 💊 Clinical Metabolizer Phenotype",
                ],
                key="sel_cross_strat_mode",
            )

        with c_opt2:
            display_pct = st.radio(
                "Display Values As:",
                ["Count + Row Percentage (%)", "Raw Count Only", "Percentage (%) Only"],
                horizontal=True,
                key="cross_display_fmt_radio",
            )

        if not proc_df.empty:
            if "Region × 💊 Clinical Metabolizer" in cross_mode:
                row_col, col_col = "Region", "Phenotype"
                title_lbl = "Geographic Region vs Clinical Phenotype"
            elif "Region × 👫 Gender" in cross_mode:
                row_col, col_col = "Region", "Gender_Clean"
                title_lbl = "Geographic Region vs Gender"
            elif "Gender × 💊 Clinical Metabolizer" in cross_mode:
                row_col, col_col = "Gender_Clean", "Phenotype"
                title_lbl = "Gender vs Clinical Phenotype"
            else:  # Top States
                row_col, col_col = "State_Clean", "Phenotype"
                title_lbl = "Top Indian States vs Clinical Phenotype"

            # Filter valid rows
            valid_cross_df = proc_df[
                proc_df[row_col].notna()
                & proc_df[col_col].notna()
                & ~proc_df[row_col].astype(str).str.lower().isin(["", "nan", "none", "null", "unknown", "unclassified / other state"])
                & ~proc_df[col_col].astype(str).str.lower().isin(["", "nan", "none", "null", "unknown", "not configured"])
            ].copy()

            # For states, keep top 10 states by sample size
            if row_col == "State_Clean":
                top_states = valid_cross_df["State_Clean"].value_counts().head(10).index.tolist()
                valid_cross_df = valid_cross_df[valid_cross_df["State_Clean"].isin(top_states)]

            if not valid_cross_df.empty and valid_cross_df[row_col].nunique() >= 2 and valid_cross_df[col_col].nunique() >= 2:
                # Raw contingency table
                contingency_raw = pd.crosstab(valid_cross_df[row_col], valid_cross_df[col_col])

                # Chi-Square test of independence
                try:
                    chi2_val, p_val, dof, _ = stats.chi2_contingency(contingency_raw.values)
                    chi2_str = f"χ² = {chi2_val:.3f}, df = {dof}, p = {p_val:.4f}"
                    sig_badge = "⚠️ Statistically Significant (p < 0.05)" if p_val < 0.05 else "✅ Homogeneous / Not Significant (p ≥ 0.05)"
                except Exception:
                    chi2_str = "N/A"
                    sig_badge = "N/A"

                st.markdown(
                    f"""
                <div style="background: #F8FAFC; border: 1.5px solid #CBD5E1; border-radius: 6px; padding: 10px 14px; margin-bottom: 12px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;">
                    <div style="font-weight: 700; color: #1E293B; font-size: 13px;">Contingency Independence Test ({title_lbl}):</div>
                    <div style="font-size: 12px; color: #475569;"><b>{chi2_str}</b> · <span style="font-weight: 600;">{sig_badge}</span></div>
                </div>
                """,
                    unsafe_allow_html=True,
                )

                # Format matrix based on display selection
                row_totals = contingency_raw.sum(axis=1)
                formatted_matrix = pd.DataFrame(index=contingency_raw.index)

                for col_name in contingency_raw.columns:
                    col_vals = contingency_raw[col_name]
                    if display_pct == "Count + Row Percentage (%)":
                        formatted_matrix[col_name] = [
                            f"{cnt} ({(cnt / tot * 100.0):.1f}%)" if tot > 0 else "0 (0.0%)"
                            for cnt, tot in zip(col_vals, row_totals)
                        ]
                    elif display_pct == "Percentage (%) Only":
                        formatted_matrix[col_name] = [
                            f"{(cnt / tot * 100.0):.1f}%" if tot > 0 else "0.0%"
                            for cnt, tot in zip(col_vals, row_totals)
                        ]
                    else:
                        formatted_matrix[col_name] = col_vals

                formatted_matrix["Total Samples (N)"] = row_totals
                st.dataframe(formatted_matrix.reset_index(), use_container_width=True, hide_index=True)

                # Grouped visual breakdown chart
                plot_cross = valid_cross_df.groupby([row_col, col_col]).size().reset_index(name="Sample Count")
                fig_cross = px.bar(
                    plot_cross,
                    x=row_col,
                    y="Sample Count",
                    color=col_col,
                    barmode="group",
                    title=f"Cross-Stratified Distribution ({title_lbl})",
                    color_discrete_sequence=px.colors.qualitative.Safe,
                )
                fig_cross.update_layout(height=320, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_cross, use_container_width=True)
            else:
                st.info("Insufficient multi-category data for this cross-tabulation.")
        else:
            st.info("Processed dataset not available for cross-tabulation.")

    # -------------------------------------------------------------------------
    # MODE 3: GENOTYPE PCA CLUSTERING
    # -------------------------------------------------------------------------
    else:
        st.markdown("#### 🔬 Genotype Principal Component Analysis (PCA)")
        st.caption(
            "PCA reduces multi-marker genotypes into continuous spatial coordinates. "
            "Use the interactive controls below to color and inspect genetic clusters by region, gender, phenotype, or diplotype."
        )

        pca_markers = [m for m in active_markers if m in proc_df.columns]
        pca_scores, pca_info = genotype_pca(proc_df, pca_markers, active_snp_conf)

        if pca_info.get("reason"):
            st.info(pca_info["reason"])
        elif not pca_scores.empty:
            pca_plot_df = pca_scores.copy()
            for metadata_col in ("Sample_ID", "Region", "State_Clean", "Gender_Clean", "Diplotype", "Phenotype"):
                if metadata_col in proc_df.columns:
                    pca_plot_df[metadata_col] = proc_df[metadata_col]

            # Add Carrier status column dynamically
            all_configs = {**active_snp_conf}
            wt_calls = [
                all_configs.get(m, {}).get("wildtype_genotype", "")
                for m in pca_markers
            ]
            if len(pca_markers) > 0:
                is_all_wt = pd.Series(True, index=proc_df.index)
                for m, wt_g in zip(pca_markers, wt_calls):
                    if m in proc_df.columns and wt_g:
                        norm_m = proc_df[m].astype(str).str.strip().str.upper().str.replace("/", "")
                        norm_wt = str(wt_g).strip().upper().replace("/", "")
                        is_all_wt = is_all_wt & (norm_m == norm_wt)
                pca_plot_df["Carrier_Status"] = np.where(is_all_wt, "Wildtype Non-Carrier", "Variant Carrier")

            c_ctrl1, c_ctrl2 = st.columns([1.2, 0.8])
            with c_ctrl1:
                color_option = st.selectbox(
                    "🎨 Color PCA Clusters By:",
                    [
                        "🗺️ Geographic Region (Region)",
                        "💊 Clinical Metabolizer (Phenotype)",
                        "🧬 Diplotype Assignment",
                        "👫 Gender",
                        "🎯 Variant Carrier Status (Dominant Model)",
                    ],
                    key="pca_color_dimension_sel",
                )

            with c_ctrl2:
                plot_dim = st.radio(
                    "Display Dimensions:",
                    ["2D Scatter (PC1 vs PC2)", "3D Scatter (PC1, PC2, PC3)"],
                    horizontal=True,
                    key="pca_plot_dim_radio",
                )

            # Map selection to column
            if "Geographic Region" in color_option:
                color_col = "Region"
            elif "Clinical Metabolizer" in color_option:
                color_col = "Phenotype"
            elif "Diplotype" in color_option:
                color_col = "Diplotype"
            elif "Gender" in color_option:
                color_col = "Gender_Clean"
            else:
                color_col = "Carrier_Status"

            variance = pca_info["explained_variance_percent"]
            pca_plot_df["PC1"] = pca_plot_df["PC1"].astype(float)
            pca_plot_df["PC2"] = pca_plot_df["PC2"].astype(float)

            hover_cols = [c for c in ["Sample_ID", "Region", "State_Clean", "Gender_Clean", "Diplotype", "Phenotype"] if c in pca_plot_df.columns]

            if plot_dim == "2D Scatter (PC1 vs PC2)":
                fig_pca_2d = px.scatter(
                    pca_plot_df,
                    x="PC1",
                    y="PC2",
                    color=color_col if color_col in pca_plot_df.columns else None,
                    hover_name="Sample_ID" if "Sample_ID" in pca_plot_df.columns else None,
                    hover_data=hover_cols,
                    title=f"Genotype PCA (PC1 vs PC2) — Grouped by {color_col}",
                    labels={
                        "PC1": f"PC1 ({variance[0]:.1f}% variance)",
                        "PC2": f"PC2 ({variance[1]:.1f}% variance)",
                    },
                    color_discrete_sequence=px.colors.qualitative.Safe,
                )
                fig_pca_2d.update_layout(height=450, margin=dict(l=20, r=20, t=50, b=20))
                st.plotly_chart(fig_pca_2d, use_container_width=True)
            else:
                if len(variance) >= 3:
                    pca_plot_df["PC3"] = pca_plot_df["PC3"].astype(float)
                    fig_pca_3d = px.scatter_3d(
                        pca_plot_df,
                        x="PC1",
                        y="PC2",
                        z="PC3",
                        color=color_col if color_col in pca_plot_df.columns else None,
                        hover_name="Sample_ID" if "Sample_ID" in pca_plot_df.columns else None,
                        hover_data=hover_cols,
                        title=f"Genotype PCA (3D) — Grouped by {color_col}",
                        labels={
                            "PC1": f"PC1 ({variance[0]:.1f}%)",
                            "PC2": f"PC2 ({variance[1]:.1f}%)",
                            "PC3": f"PC3 ({variance[2]:.1f}%)",
                        },
                        color_discrete_sequence=px.colors.qualitative.Safe,
                    )
                    fig_pca_3d.update_layout(height=520, margin=dict(l=20, r=20, t=50, b=20))
                    st.plotly_chart(fig_pca_3d, use_container_width=True)
                else:
                    st.caption("A 3D PCA view requires at least three informative components; data currently supports 2D.")

            st.caption(f"Biomarkers in PCA model: **{', '.join(pca_info['markers_used'])}**")
            if pca_info.get("skipped_markers"):
                st.caption(f"Skipped monomorphic markers: **{', '.join(pca_info['skipped_markers'])}**")
