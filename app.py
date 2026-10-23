"""Automated Population Pharmacogenomics Analysis Platform

Modularized Streamlit application orchestrating target selection, dataset upload,
column mapping, and the comprehensive 9-tab results dashboard.
"""

import json
import logging
import pandas as pd
import streamlit as st
from core.qc import DataQCEngine
from core.stratification import DemographicStratifier
from marker_catalog import detect_genes_in_columns
from ui.column_mapping import render_column_mapping
from ui.dataset_upload import render_dataset_upload
from ui.session_manager import (
    clear_deletion_flag,
    has_active_data,
    init_session_state,
    reset_all_data,
    save_active_session,
)
from ui.styles import apply_custom_styles, render_sidebar_controls
from ui.target_selection import render_target_selection
from ui.views.results_dashboard import render_results_dashboard

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# 1. PAGE CONFIG & STYLING
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Population Pharmacogenomics Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)
apply_custom_styles()

# -----------------------------------------------------------------------------
# 2. SESSION STATE & PERSISTENCE INITIALIZATION
# -----------------------------------------------------------------------------
init_session_state()


# -----------------------------------------------------------------------------
# 3. ANALYSIS EXECUTION PIPELINE
# -----------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def run_fast_pipeline(
    df: pd.DataFrame,
    mapped_cols_tuple: tuple,
    snp_rules_tuple: tuple = (),
    active_gene: str = None,
    translation_rules_json: str = "{}",
):
    """Cached execution of data quality control and stratification engine."""
    mapped_dict = dict(mapped_cols_tuple)
    snp_rules = (
        {k: dict(v) for k, v in snp_rules_tuple} if snp_rules_tuple else None
    )
    translation_rules = (
        json.loads(translation_rules_json) if translation_rules_json else {}
    )

    if df is None or len(df) == 0:
        clean_df = pd.DataFrame()
        qc_report = {
            "total_samples": 0,
            "fully_valid_samples": 0,
            "duplicate_sample_count": 0,
            "duplicate_sample_ids": [],
            "invalid_genotypes_count": 0,
            "invalid_genotypes_log": {},
            "snp_qc_stats": {},
        }
        full_results = {
            "overall": {
                "sample_count": 0,
                "snps": {},
                "pgx_summary": {"phenotypes": {}, "diplotypes": {}},
                "cyp2c19_summary": {"phenotypes": {}, "diplotypes": {}},
                "cyp2c9_summary": {"phenotypes": {}, "diplotypes": {}},
            },
            "regional": {},
            "state_wise": {},
            "gender": {},
            "processed_dataframe": pd.DataFrame(),
            "active_gene": active_gene,
            "target_gene": active_gene,
        }
        return clean_df, qc_report, full_results

    try:
        qc_engine = DataQCEngine(
            df, custom_col_mapping=mapped_dict, custom_snp_rules=snp_rules
        )
    except TypeError:
        qc_engine = DataQCEngine(df, custom_col_mapping=mapped_dict)
        if snp_rules:
            try:
                if hasattr(qc_engine, "snp_rules"):
                    for k, v in snp_rules.items():
                        if k in qc_engine.snp_rules:
                            qc_engine.snp_rules[k] = {
                                **qc_engine.snp_rules[k],
                                **v,
                            }
                        else:
                            qc_engine.snp_rules[k] = v
            except Exception:
                pass

    clean_df, qc_report = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(
        clean_df,
        custom_snp_rules=snp_rules,
        active_gene=active_gene,
        custom_translation_rules=translation_rules,
    )
    return clean_df, qc_report, full_results


def serialize_rules(rules_dict):
    """Serialize custom rules dictionary to hashable tuple for Streamlit caching."""
    if not rules_dict:
        return ()
    serialized = []
    for k, v in sorted(rules_dict.items()):
        inner_items = []
        for sk, sv in sorted(v.items()):
            if isinstance(sv, (set, list, tuple)):
                inner_items.append((sk, tuple(sorted(sv))))
            else:
                inner_items.append((sk, sv))
        serialized.append((k, tuple(inner_items)))
    return tuple(serialized)


# Execute Pipeline
df_raw = st.session_state["uploaded_df"]
mapped_tuple = tuple(sorted(st.session_state["mapped_cols"].items()))
snp_rules_items = serialize_rules(st.session_state.get("custom_snp_rules", {}))
df_cols_detected = (
    detect_genes_in_columns(list(df_raw.columns))
    if df_raw is not None and not df_raw.empty
    else []
)
active_gene_to_run = st.session_state.get("selected_gene") or (
    df_cols_detected[0] if df_cols_detected else "Custom"
)
translation_rules_json = json.dumps(
    st.session_state.get("custom_translation_rules", {}), sort_keys=True
)

clean_df, qc_report, full_results = run_fast_pipeline(
    df_raw,
    mapped_tuple,
    snp_rules_items,
    active_gene=active_gene_to_run,
    translation_rules_json=translation_rules_json,
)
processed_df = full_results.get("processed_dataframe", pd.DataFrame())

# -----------------------------------------------------------------------------
# 4. SIDEBAR CONTROLS
# -----------------------------------------------------------------------------
render_sidebar_controls(df_raw, reset_all_data)

# -----------------------------------------------------------------------------
# 5. MAIN WORKFLOW STEP NAVIGATION
# -----------------------------------------------------------------------------
VALID_STEPS = ["1 · Select Target", "2 · Upload", "3 · Map columns", "4 · Results"]
curr_step = st.session_state.get("workflow_step", "")
if curr_step not in VALID_STEPS:
    if "Target" in curr_step:
        st.session_state["workflow_step"] = "1 · Select Target"
    elif "Upload" in curr_step:
        st.session_state["workflow_step"] = "2 · Upload"
    elif "Map" in curr_step:
        st.session_state["workflow_step"] = "3 · Map columns"
    elif "Result" in curr_step:
        st.session_state["workflow_step"] = "4 · Results"
    else:
        st.session_state["workflow_step"] = "1 · Select Target"

main_tab_home, main_tab_about = st.tabs(["🏠 Home Workspace", "ℹ️ About Platform"])

with main_tab_about:
    st.markdown(
        """
    <div style="background: linear-gradient(135deg, #0F172A 0%, #1E1B4B 45%, #2E1065 100%); border-radius: 16px; padding: 2.5rem 2.8rem; color: #FFFFFF; margin-top: 0.5rem; margin-bottom: 2rem; box-shadow: 0 12px 30px rgba(15, 23, 42, 0.4); border: 1px solid #4338CA;">
        <div style="background: linear-gradient(90deg, #4F46E5 0%, #7C3AED 100%); border: 1px solid #A78BFA; border-radius: 12px; padding: 0.8rem 1.5rem; display: inline-block; margin-bottom: 1.2rem; box-shadow: 0 4px 14px rgba(124, 58, 237, 0.35);">
            <h1 style="color: #FFFFFF; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; font-weight: 800; font-size: 1.85rem; margin: 0; letter-spacing: -0.01em;">
                🧬 Automated Population Pharmacogenomics Analysis Platform
            </h1>
        </div>
        <p style="color: #F1F5F9; font-size: 1.1rem; line-height: 1.7; font-weight: 400; max-width: 900px; margin-bottom: 1.8rem;">
            <strong>PGx Analytics Pro</strong> is a standardized, high-performance pharmacogenomics engine engineered for instant quality control, Hardy–Weinberg Equilibrium (HWE) evaluation, regional demographic stratification, and CPIC metabolizer phenotype translation.
        </p>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 1rem;">
            <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #34D399; backdrop-filter: blur(8px);">
                <div style="font-size: 0.78rem; text-transform: uppercase; color: #A7F3D0; font-weight: 800; letter-spacing: 0.06em;">Target Gene & SNPs</div>
                <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">Selected gene and markers</div>
                <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">Configured from each uploaded dataset</div>
            </div>
            <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #60A5FA; backdrop-filter: blur(8px);">
                <div style="font-size: 0.78rem; text-transform: uppercase; color: #93C5FD; font-weight: 800; letter-spacing: 0.06em;">Statistical & HWE Engine</div>
                <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">Chi-Square (χ²) & Haldane</div>
                <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">Test and interpretation shown per marker</div>
            </div>
            <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #FBBF24; backdrop-filter: blur(8px);">
                <div style="font-size: 0.78rem; text-transform: uppercase; color: #FDE68A; font-weight: 800; letter-spacing: 0.06em;">Demographic Matrix</div>
                <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">Uploaded population groups</div>
                <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">Region and locality values from your data</div>
            </div>
            <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #F472B6; backdrop-filter: blur(8px);">
                <div style="font-size: 0.78rem; text-transform: uppercase; color: #F472B6; font-weight: 800; letter-spacing: 0.06em;">CPIC Clinical Calling</div>
                <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">Configured phenotype classes</div>
                <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">Based on supplied interpretation rules</div>
            </div>
        </div>
    </div>
    """,
        unsafe_allow_html=True,
    )
    st.info(
        "🔒 **100% Local Browser Privacy**: All dataset parsing, statistical math, and SQLite indexing (`pgx_analysis.db`) run locally on your device with zero cloud data transmission."
    )
    st.success(
        "📊 **Multi-Format Publication Export**: One-click download of ready-to-present multi-tab Excel (`.xlsx`), raw summary CSV (`.csv`), and formatted PDF (`.pdf`) reports."
    )
    st.markdown("<br>", unsafe_allow_html=True)

with main_tab_home:
    step_cols = st.columns(4)
    t1 = (
        "primary"
        if st.session_state["workflow_step"] == "1 · Select Target"
        else "secondary"
    )
    t2 = (
        "primary"
        if st.session_state["workflow_step"] == "2 · Upload"
        else "secondary"
    )
    t3 = (
        "primary"
        if st.session_state["workflow_step"] == "3 · Map columns"
        else "secondary"
    )
    t4 = (
        "primary"
        if st.session_state["workflow_step"] == "4 · Results"
        else "secondary"
    )

    if step_cols[0].button(
        "1 · Select Target",
        use_container_width=True,
        type=t1,
        key="nav_btn_step1",
    ):
        st.session_state["workflow_step"] = "1 · Select Target"
        save_active_session()
        st.rerun()

    if step_cols[1].button(
        "2 · Upload", use_container_width=True, type=t2, key="nav_btn_step2"
    ):
        st.session_state["workflow_step"] = "2 · Upload"
        save_active_session()
        st.rerun()

    if step_cols[2].button(
        "3 · Map columns",
        use_container_width=True,
        type=t3,
        key="nav_btn_step3",
    ):
        st.session_state["workflow_step"] = "3 · Map columns"
        save_active_session()
        st.rerun()

    if step_cols[3].button(
        "4 · Results", use_container_width=True, type=t4, key="nav_btn_step4"
    ):
        st.session_state["workflow_step"] = "4 · Results"
        save_active_session()
        st.rerun()

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # WORKFLOW PIPELINE ROUTER
    # -------------------------------------------------------------------------
    if st.session_state["workflow_step"] == "1 · Select Target":
        render_target_selection(df_raw)

    elif st.session_state["workflow_step"] == "2 · Upload":
        render_dataset_upload(df_raw)

    elif st.session_state["workflow_step"] == "3 · Map columns":
        render_column_mapping(df_raw)

    elif st.session_state["workflow_step"] == "4 · Results":
        render_results_dashboard(df_raw, processed_df, full_results, qc_report)

save_active_session()

# Footer
st.markdown(
    """
<br><hr>
<div style="font-size: 0.8rem; color: #94A3B8; text-align: center;">
    Automated Population Pharmacogenomics Analysis Platform — Built for Reproducible Genotype, Allele, HWE, and CPIC Phenotype Analysis.
</div>
""",
    unsafe_allow_html=True,
)
