"""Unified results dashboard view component coordinating all 9 analysis tabs."""

import io
import logging
from typing import Any, Dict
import pandas as pd
import streamlit as st
from config import CYP2C9_SNP_CONFIG, SNP_CONFIG
from marker_catalog import (
    detect_genes_in_columns,
    find_best_column_match,
    get_default_markers_for_gene,
    load_marker_catalog,
)
from reports.report_generator import ReportGenerator
from ui.session_manager import (
    has_active_data,
    save_active_dataset,
    save_active_session,
)
from ui.views.cohort_comparison import render_cohort_comparison
from ui.views.diplotype_phenotype import render_diplotype_phenotype
from ui.views.executive_overview import render_executive_overview
from ui.views.gender_analysis import render_gender_analysis
from ui.views.genotype_allele import render_genotype_allele
from ui.views.geographic_regions import render_geographic_regions
from ui.views.hardy_weinberg import render_hardy_weinberg
from ui.views.population_groups import render_population_groups
from ui.views.quality_audit import render_quality_audit
from ui.views.state_profiles import render_state_profiles
from ui.views.statistical_validation import render_statistical_validation


logger = logging.getLogger(__name__)


def get_df_download_excel(df: pd.DataFrame) -> bytes:
    """Serialize DataFrame to Excel bytes."""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df.to_excel(writer, index=False)
    return output.getvalue()


def get_df_download_csv(df: pd.DataFrame) -> bytes:
    """Serialize DataFrame to CSV bytes."""
    return df.to_csv(index=False).encode("utf-8")


def render_results_dashboard(
    df_raw: pd.DataFrame,
    processed_df: pd.DataFrame,
    full_results: Dict[str, Any],
    qc_report: Dict[str, Any],
) -> None:
    """Render Step 4 Results Dashboard coordinating the 9 analysis tabs and export section."""
    if not has_active_data() or qc_report is None:
        st.info(
            "ℹ️ Workspace is empty (0 samples). All existing data was deleted. Upload a new Excel/CSV dataset file or add sample records on Page 1."
        )
        return

    if len(st.session_state.get("datasets_dict", {})) > 1:
        st.markdown(
            f"""
        <div style="background: #F8FAFC; border: 1.5px solid #3B82F6; border-radius: 8px; padding: 12px 16px; margin-bottom: 12px;">
            <div style="display: flex; align-items: center; justify-content: space-between;">
                <div style="font-weight: 700; color: #1E40AF; font-size: 14px;">📂 Multi-Dataset Profile Active ({len(st.session_state['datasets_dict'])} Datasets Uploaded)</div>
                <div style="font-size: 12px; color: #64748B;">Select dataset from dropdown to calculate and view profile</div>
            </div>
        </div>
        """,
            unsafe_allow_html=True,
        )
        ds_list = list(st.session_state["datasets_dict"].keys())
        active_current = st.session_state.get("active_ds_name", ds_list[0])
        idx_cur = (
            ds_list.index(active_current) if active_current in ds_list else 0
        )
        c_ds1, c_ds2 = st.columns([1.5, 1.5])
        with c_ds1:
            sel_ds_name = st.selectbox(
                "📊 Select Active Dataset Profile to View:",
                ds_list,
                index=idx_cur,
                key="sel_multi_ds_switcher",
            )
            if (
                sel_ds_name
                and st.session_state.get("active_ds_name") != sel_ds_name
            ):
                st.session_state["active_ds_name"] = sel_ds_name
                new_active_df = st.session_state["datasets_dict"][sel_ds_name]
                st.session_state["uploaded_df"] = new_active_df
                # Dynamically synchronize target gene to match the switched dataset
                ds_detected_genes = detect_genes_in_columns(
                    list(new_active_df.columns)
                )
                if ds_detected_genes:
                    target_g_switched = ds_detected_genes[0]
                    st.session_state["selected_gene"] = target_g_switched
                    all_g_m = get_default_markers_for_gene(
                        target_g_switched, catalog_df=load_marker_catalog()
                    )
                    matched_m = [
                        m
                        for m in all_g_m
                        if find_best_column_match(
                            m,
                            m.split(" ")[0],
                            target_g_switched,
                            list(new_active_df.columns),
                        )
                    ]
                    st.session_state["selected_markers"] = (
                        matched_m if matched_m else all_g_m
                    )
                save_active_dataset(new_active_df)
                save_active_session()
                st.rerun()
        with c_ds2:
            st.markdown(
                "<div style='margin-top: 28px;'></div>", unsafe_allow_html=True
            )
            active_df_preview = st.session_state["datasets_dict"].get(
                sel_ds_name, pd.DataFrame()
            )
            st.info(
                f"Viewing: **{sel_ds_name}** ({len(active_df_preview):,} samples, {len(active_df_preview.columns)} columns)"
            )
        st.divider()

    cur_uploaded = st.session_state.get("uploaded_df")
    detected_cols_gene = (
        detect_genes_in_columns(list(cur_uploaded.columns))
        if cur_uploaded is not None and not cur_uploaded.empty
        else []
    )
    active_gene = (
        st.session_state.get("selected_gene")
        or (detected_cols_gene[0] if detected_cols_gene else None)
        or full_results.get("active_gene")
        or full_results.get("target_gene")
        or "Custom"
    )
    active_snp_conf = {
        **SNP_CONFIG,
        **CYP2C9_SNP_CONFIG,
        **st.session_state.get("custom_snp_rules", {}),
    }
    ov = full_results["overall"]
    active_markers = list(ov.get("snps", {}).keys())
    if not active_markers and "markers" in full_results:
        active_markers = [m for m in full_results["markers"] if m]
    if (
        not active_markers
        and "uploaded_df" in st.session_state
        and st.session_state["uploaded_df"] is not None
    ):
        active_markers = [
            c
            for c in st.session_state["uploaded_df"].columns
            if ("*" in c or "rs" in c.lower() or active_gene.lower() in c.lower())
            and c
            not in (
                "Sample_ID",
                "Region",
                "Gender",
                "State",
                "State_Clean",
                "Gender_Clean",
                "Diplotype",
                "Phenotype",
                "Activity_Score",
                "Test_Requested",
                "Date_of_Birth",
            )
        ]
    if not active_markers:
        active_markers = get_default_markers_for_gene(
            active_gene, catalog_df=load_marker_catalog()
        )

    # 11 Main Tabs Consolidating All Analysis Views, Statistical Validation & Comparison Studio
    (
        tab_overview,
        tab_gen,
        tab_stat,
        tab_comparison,
        tab_hwe,
        tab_pop,
        tab_pgx,
        tab_qc,
        tab_geo,
        tab_gender,
        tab_state,
    ) = st.tabs(
        [
            "📊 Executive Overview",
            "🧬 Genotype & Allele",
            "📈 Statistical Validation",
            "⚖️ Geographic & State Comparison",
            "⚖️ Hardy-Weinberg",
            "🌐 Population Stratification",
            f"💊 {active_gene} Calling",
            "🛡️ Data Quality Audit",
            "🗺️ Geographic Groups",
            "👫 Gender-Wise",
            "📍 Dedicated State Profiles",
        ]
    )

    # Helper for Sub-cohort & State-Wise Filter UI
    def render_subcohort_pills(key_prefix: str):
        st_dict = full_results.get("state_wise", {})
        reg_dict = full_results.get("regional", {})

        c_mode, c_sub = st.columns([1.1, 0.9])
        with c_mode:
            filter_mode = st.radio(
                "Filter Population By:",
                [
                    "🌐 Overall Cohort",
                    "🗺️ Geographic Region",
                    "📍 Individual State",
                ],
                horizontal=True,
                key=f"{key_prefix}_filter_mode",
            )

        with c_sub:
            if filter_mode == "🗺️ Geographic Region":
                reg_names = sorted(
                    r
                    for r, item in reg_dict.items()
                    if item.get("sample_count", 0) > 0
                )
                reg_options = [
                    f"{r} (n={reg_dict.get(r, {}).get('sample_count', 0):,})"
                    for r in reg_names
                ]
                if reg_options:
                    selected_reg_str = st.selectbox(
                        "Select Geographic Region:",
                        reg_options,
                        key=f"{key_prefix}_sel_region",
                    )
                    r_name = selected_reg_str.split(" (")[0]
                    sel_key = r_name
                    sel_cohort = reg_dict.get(
                        r_name, full_results.get("overall", {})
                    )
                else:
                    st.info("No regional data available.")
                    sel_key = "Overall Population"
                    sel_cohort = full_results.get("overall", {})
            elif filter_mode == "📍 Individual State":
                state_items = sorted(
                    [
                        (s, data.get("sample_count", 0))
                        for s, data in st_dict.items()
                        if data.get("sample_count", 0) > 0
                    ],
                    key=lambda x: x[1],
                    reverse=True,
                )
                if state_items:
                    state_options = [
                        f"{s} (n={cnt:,})" for s, cnt in state_items
                    ]
                    selected_state_str = st.selectbox(
                        "Select State Profile:",
                        state_options,
                        key=f"{key_prefix}_sel_state",
                    )
                    s_name = selected_state_str.split(" (")[0]
                    sel_key = f"State: {s_name}"
                    sel_cohort = st_dict.get(
                        s_name, full_results.get("overall", {})
                    )
                else:
                    st.info("No state data available.")
                    sel_key = "Overall Population"
                    sel_cohort = full_results.get("overall", {})
            else:
                sel_key = "Overall Population"
                sel_cohort = full_results.get("overall", {})
                st.caption(
                    f"Displaying full analyzed cohort (**{sel_cohort.get('sample_count', 0):,}** samples)."
                )

        return sel_key, sel_cohort

    # 1. Executive Overview
    with tab_overview:
        render_executive_overview(
            full_results, qc_report, active_gene, active_markers
        )

    # 2. Genotype & Allele
    with tab_gen:
        render_genotype_allele(
            active_gene,
            active_markers,
            active_snp_conf,
            full_results,
            render_subcohort_pills,
        )

    # 3. Advanced Statistical Validation
    with tab_stat:
        render_statistical_validation(
            full_results,
            active_gene,
            active_markers,
            processed_df,
        )

    # 4. Hardy-Weinberg
    with tab_hwe:
        render_hardy_weinberg(render_subcohort_pills)

    # 4. Population Stratification
    with tab_pop:
        render_population_groups(
            full_results, active_markers, active_snp_conf
        )

    # 5. Gene Calling & Translation
    with tab_pgx:
        render_diplotype_phenotype(
            active_gene,
            active_markers,
            active_snp_conf,
            processed_df,
            full_results,
            render_subcohort_pills,
        )

    # 6. Quality Audit
    with tab_qc:
        render_quality_audit(qc_report)

    # 7. Geographic Groups
    with tab_geo:
        render_geographic_regions(
            active_gene, active_markers, active_snp_conf, full_results
        )

    # 8. Gender Analysis
    with tab_gender:
        render_gender_analysis(active_gene, active_markers, full_results)

    # 9. State Profiles
    with tab_state:
        render_state_profiles(
            active_gene,
            active_markers,
            active_snp_conf,
            full_results,
            processed_df,
        )

    # 10. Geographic & State Comparison & Predictions
    with tab_comparison:
        render_cohort_comparison(
            active_gene,
            active_markers,
            full_results,
            processed_df,
        )

    # Export Section
    st.divider()
    st.markdown("### Export")
    st.caption("Download the full result set for the current group selection.")

    st.markdown('<div class="export-nav-anchor"></div>', unsafe_allow_html=True)
    e_prev, e1, e2, e3, e_remap, e_retarget = st.columns(
        [0.18, 0.18, 0.18, 0.18, 0.14, 0.14]
    )

    rep_gen = (
        ReportGenerator(full_results, qc_report) if has_active_data() else None
    )

    pdf_bytes = None
    if rep_gen:
        try:
            pdf_bytes = rep_gen.generate_pdf()
        except Exception as ex:
            logger.warning("Error generating PDF preview bytes: %s", ex)

    with e_prev:
        if pdf_bytes:
            preview_label = (
                "👁️ Hide Preview"
                if st.session_state.get("show_pdf_preview", False)
                else "👁️ Preview PDF"
            )
            if st.button(
                preview_label,
                use_container_width=True,
                key="btn_toggle_pdf_prev",
            ):
                st.session_state["show_pdf_preview"] = not st.session_state.get(
                    "show_pdf_preview", False
                )
                st.rerun()
        else:
            st.button("👁️ Preview PDF", disabled=True, use_container_width=True)

    with e1:
        if rep_gen:
            excel_bytes = rep_gen.generate_excel()
            st.download_button(
                "⬇ Excel (.xlsx)",
                data=excel_bytes,
                file_name=f"{active_gene}_Population_Report.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
                key="btn_exp_excel",
            )
        else:
            st.button(
                "⬇ Excel (.xlsx)", disabled=True, use_container_width=True
            )

    with e2:
        if rep_gen:
            csv_bytes = rep_gen.generate_csv()
            st.download_button(
                "⬇ CSV (summary)",
                data=csv_bytes,
                file_name=f"{active_gene}_Population_Summary.csv",
                mime="text/csv",
                use_container_width=True,
                key="btn_exp_csv",
            )
        else:
            st.button(
                "⬇ CSV (summary)", disabled=True, use_container_width=True
            )

    with e3:
        if pdf_bytes:
            st.download_button(
                "⬇ PDF Report",
                data=pdf_bytes,
                file_name=f"{active_gene}_Population_Report.pdf",
                mime="application/pdf",
                use_container_width=True,
                key="btn_exp_pdf",
            )
        else:
            st.button("⬇ PDF Report", disabled=True, use_container_width=True)

    with e_remap:
        if st.button(
            "← Re-map", use_container_width=True, key="btn_res_remap"
        ):
            st.session_state["workflow_step"] = "3 · Map columns"
            save_active_session()
            st.rerun()

    with e_retarget:
        if st.button(
            "🎯 Target",
            use_container_width=True,
            key="btn_res_retarget",
            help="Switch target gene or select different biomarkers",
        ):
            st.session_state["workflow_step"] = "1 · Select Target"
            save_active_session()
            st.rerun()

    # Render Inline Interactive PDF Preview if toggled
    if st.session_state.get("show_pdf_preview", False) and pdf_bytes:
        st.markdown(
            "<div style='margin-top: 14px;'></div>", unsafe_allow_html=True
        )
        with st.expander(
            f"📑 Live Document Preview: {active_gene} Population Report",
            expanded=True,
        ):
            p_col1, p_col2 = st.columns([0.8, 0.2])
            with p_col1:
                st.caption(
                    "Inspect formatted tables, HWE statistics, allele distributions, and CPIC phenotype summaries before downloading."
                )
            with p_col2:
                st.download_button(
                    "📥 Download PDF Directly",
                    data=pdf_bytes,
                    file_name=f"{active_gene}_Population_Report.pdf",
                    mime="application/pdf",
                    use_container_width=True,
                    key="btn_prev_dl_pdf",
                )

            overall_preview = full_results.get("overall", {})
            preview_snps = overall_preview.get("snps", {})
            preview_pgx = (
                overall_preview.get("pgx_summary")
                or overall_preview.get(f"{active_gene.lower()}_summary")
                or overall_preview.get("cyp2c19_summary")
                or overall_preview.get("cyp2c9_summary")
                or {}
            )

            st.markdown(
                f"#### Cohort summary — {overall_preview.get('sample_count', 0):,} samples"
            )
            hwe_preview_rows = []
            for snp_name, snp_data in preview_snps.items():
                hwe = snp_data.get("hwe", {})
                allele_freqs = snp_data.get("allele_freqs", {})
                hwe_preview_rows.append(
                    {
                        "Biomarker": snp_name,
                        "Valid N": snp_data.get("valid_samples", 0),
                        "Allele frequencies": ", ".join(
                            f"{a}: {f}" for a, f in allele_freqs.items()
                        ),
                        "Chi-square": hwe.get("chi2_stat", 0),
                        "P-value": hwe.get("p_value", 1),
                        "Exact P-value": hwe.get("exact_p_value", 1),
                        "HWE status": hwe.get("interpretation", "No Data"),
                    }
                )
            st.markdown("##### Hardy–Weinberg and allele summary")
            st.dataframe(
                pd.DataFrame(hwe_preview_rows),
                use_container_width=True,
                hide_index=True,
            )

            phenotype_rows = [
                {
                    "Phenotype": name,
                    "Count": info.get("count", 0),
                    "Percentage": f"{info.get('percentage', 0):.1f}%",
                }
                for name, info in preview_pgx.get("phenotypes", {}).items()
            ]
            st.markdown("##### Metabolizer phenotype summary")
            st.dataframe(
                pd.DataFrame(phenotype_rows),
                use_container_width=True,
                hide_index=True,
            )

            region_preview_rows = []
            for region_name, region_data in full_results.get(
                "regional", {}
            ).items():
                region_pgx = (
                    region_data.get("pgx_summary")
                    or region_data.get(f"{active_gene.lower()}_summary")
                    or region_data.get("cyp2c19_summary")
                    or region_data.get("cyp2c9_summary")
                    or {}
                )
                poor = region_pgx.get("phenotypes", {}).get(
                    "Poor Metabolizer", {}
                )
                region_preview_rows.append(
                    {
                        "Region": region_name,
                        "Samples": region_data.get("sample_count", 0),
                        "Poor metabolizers": poor.get("count", 0),
                        "Poor metabolizer %": f"{poor.get('percentage', 0):.1f}%",
                    }
                )
            st.markdown("##### Regional summary")
            st.dataframe(
                pd.DataFrame(region_preview_rows),
                use_container_width=True,
                hide_index=True,
            )

    if df_raw is not None and not df_raw.empty:
        st.markdown(
            "<div style='margin-top: 10px;'></div>", unsafe_allow_html=True
        )
        target_full_df = (
            processed_df
            if (processed_df is not None and not processed_df.empty)
            else df_raw
        )
        c_full_dl1, c_full_dl2 = st.columns(2)
        with c_full_dl1:
            label_desc = (
                "⬇️ Download Complete Classified Dataset (.xlsx)"
                if (processed_df is not None and not processed_df.empty)
                else "⬇️ Download Complete Dataset (.xlsx)"
            )
            st.download_button(
                label=label_desc,
                data=get_df_download_excel(target_full_df),
                file_name=f"complete_pgx_dataset_{len(target_full_df)}_samples.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
                key="btn_res_full_dl_xlsx",
            )
        with c_full_dl2:
            st.download_button(
                label="⬇️ Download Complete Dataset (.csv)",
                data=get_df_download_csv(target_full_df),
                file_name=f"complete_pgx_dataset_{len(target_full_df)}_samples.csv",
                mime="text/csv",
                use_container_width=True,
                key="btn_res_full_dl_csv",
            )

    target_gene_lbl = (
        st.session_state.get("selected_gene")
        or full_results.get("active_gene")
        or full_results.get("target_gene")
        or "Custom"
    )
    st.markdown(
        f"""
    <div class="footer-disclaimer-text">
        PGx Pop MVP · {target_gene_lbl} calling uses configured translation rules. Validate clinical interpretation independently before clinical use.
    </div>
    """,
        unsafe_allow_html=True,
    )
