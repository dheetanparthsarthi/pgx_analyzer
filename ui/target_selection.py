"""Step 1 Target Selection component (Gene and Biomarkers)."""

import pandas as pd
import streamlit as st
from marker_catalog import get_default_markers_for_gene, load_marker_catalog
from ui.session_manager import has_active_data, save_active_session, set_active_target_gene


def render_target_selection(df_raw: pd.DataFrame) -> None:
    """Render Step 1: Select Target (Gene & Biomarkers)."""
    st.markdown("## 1️⃣ Select Pharmacogenetic Target (Gene & Biomarkers)")
    st.caption(
        "Choose which gene and specific biomarker(s) to analyze. Only your selected markers will be uploaded, mapped, and visualized."
    )

    catalog_df = load_marker_catalog()

    if has_active_data():
        c_wb1, c_wb2 = st.columns([1.5, 0.5])
        with c_wb1:
            st.info(
                f"📂 **Active Dataset Loaded**: `{len(df_raw):,} samples` in workspace. You can change your target gene and re-map columns anytime, or view your current results."
            )
        with c_wb2:
            if st.button(
                "📊 View Results (Step 4) →",
                use_container_width=True,
                key="btn_step1_jump_res",
            ):
                st.session_state["workflow_step"] = "4 · Results"
                save_active_session()
                st.rerun()

    catalog_genes = (
        sorted(catalog_df["Gene"].unique()) if not catalog_df.empty else []
    )
    custom_genes = st.session_state.get("custom_genes", [])
    gene_list = sorted(list(set(catalog_genes + custom_genes)))
    gene_options = ["(Select Target Gene)"] + gene_list
    curr_gene = st.session_state.get("selected_gene")
    g_idx = (
        gene_options.index(curr_gene)
        if (curr_gene and curr_gene in gene_options)
        else 0
    )

    if (
        "step1_sel_gene" in st.session_state
        and st.session_state["step1_sel_gene"] not in gene_options
    ):
        del st.session_state["step1_sel_gene"]

    c_g1, c_g2 = st.columns([1, 1.8])
    with c_g1:
        sel_gene_choice = st.selectbox(
            "1️⃣ Select Target Gene:",
            gene_options,
            index=g_idx,
            key="step1_sel_gene",
        )
        if sel_gene_choice != "(Select Target Gene)":
            if sel_gene_choice != st.session_state.get("selected_gene"):
                set_active_target_gene(
                    sel_gene_choice, markers=None, catalog_df=catalog_df
                )
                st.rerun()
        else:
            if st.session_state.get("selected_gene") is not None:
                set_active_target_gene(
                    None, markers=[], catalog_df=catalog_df
                )
                st.rerun()

    with c_g2:
        sel_gene = st.session_state.get("selected_gene")
        if sel_gene:
            base_markers = get_default_markers_for_gene(sel_gene, catalog_df)
            custom_m_for_gene = st.session_state.get("custom_markers", {}).get(
                sel_gene, []
            )
            available_markers = []
            for m in base_markers + custom_m_for_gene:
                if m not in available_markers:
                    available_markers.append(m)

            curr_selected = st.session_state.get("selected_markers", [])
            valid_defaults = [m for m in curr_selected if m in available_markers]
            if not valid_defaults:
                valid_defaults = list(available_markers)
                st.session_state["selected_markers"] = valid_defaults

            c_btn_m1, c_btn_m2, _ = st.columns([0.3, 0.3, 0.4])
            with c_btn_m1:
                if st.button(
                    "⚡ Select All",
                    key="btn_sel_all_markers_step1",
                    help="Select all available biomarkers",
                ):
                    st.session_state["selected_markers"] = list(
                        available_markers
                    )
                    if "step1_sel_markers" in st.session_state:
                        del st.session_state["step1_sel_markers"]
                    st.rerun()
            with c_btn_m2:
                if st.button(
                    "🧹 Clear All",
                    key="btn_desel_all_markers_step1",
                    help="Deselect all biomarkers",
                ):
                    st.session_state["selected_markers"] = []
                    if "step1_sel_markers" in st.session_state:
                        del st.session_state["step1_sel_markers"]
                    st.rerun()

            if "step1_sel_markers" in st.session_state:
                if any(
                    m not in available_markers
                    for m in st.session_state["step1_sel_markers"]
                ):
                    del st.session_state["step1_sel_markers"]

            sel_markers = st.multiselect(
                f"2️⃣ Select Biomarker(s) for {sel_gene}:",
                available_markers,
                default=valid_defaults,
                key="step1_sel_markers",
                placeholder="Choose biomarkers to analyze or click 'Select All'",
            )
            if sel_markers != st.session_state.get("selected_markers"):
                st.session_state["selected_markers"] = sel_markers
                for k in list(st.session_state.keys()):
                    if k.startswith(
                        (
                            "w_map_",
                            "sel_ref_",
                            "sel_var_",
                            "custom_ref_",
                            "custom_var_",
                            "_mapped_",
                            "_col_tracker_",
                        )
                    ):
                        del st.session_state[k]
        else:
            st.info(
                "👆 Please select a target gene to configure biomarkers, or proceed to **2 · Upload** to upload your dataset and auto-detect your target gene."
            )

    # Expandable Section: Add Custom Gene / Biomarker Manually
    with st.expander(
        "➕ Can't find your Gene or Biomarker? Add Custom Targets Manually",
        expanded=False,
    ):
        st.markdown(
            "Add custom target genes or biomarkers not present in the default catalog:"
        )
        t_mode = st.radio(
            "Target Addition Type:",
            [
                f"Add Biomarker(s) to Selected Gene ({sel_gene})",
                "Add a Completely New Gene & Biomarker(s)",
            ],
            horizontal=True,
            key="step1_custom_type_radio",
        )
        if "Completely New" in t_mode:
            c_cg1, c_cg2 = st.columns([1, 1.8])
            with c_cg1:
                in_new_gene = (
                    st.text_input(
                        "New Gene Symbol / Name:",
                        placeholder="e.g. CYP2D6, DPYD, TPMT",
                        key="in_new_custom_gene_step1",
                    )
                    .strip()
                    .upper()
                )
            with c_cg2:
                in_new_markers = st.text_input(
                    "Biomarker Name(s) (comma-separated):",
                    placeholder="e.g. *4 (rs3892097), *10 (rs1065852)",
                    key="in_new_custom_markers_step1",
                ).strip()
        else:
            in_new_gene = sel_gene
            in_new_markers = st.text_input(
                f"New Biomarker Name(s) for {sel_gene} (comma-separated):",
                placeholder="e.g. *4 (rs28399504), *6",
                key="in_existing_custom_markers_step1",
            ).strip()

        c_btn_add, c_btn_reset = st.columns([0.35, 0.65])
        with c_btn_add:
            if st.button(
                "➕ Add to Target List",
                type="secondary",
                key="btn_add_custom_target_step1",
            ):
                if not in_new_gene:
                    st.error("Please provide a valid Gene name.")
                else:
                    c_genes = st.session_state.get("custom_genes", [])
                    if in_new_gene not in c_genes and (
                        catalog_df.empty
                        or in_new_gene not in catalog_df["Gene"].unique()
                    ):
                        c_genes.append(in_new_gene)
                        st.session_state["custom_genes"] = c_genes

                    c_markers_dict = st.session_state.get(
                        "custom_markers", {}
                    )
                    current_m = c_markers_dict.get(in_new_gene, [])
                    if in_new_markers:
                        parsed_m = [
                            m.strip()
                            for m in in_new_markers.split(",")
                            if m.strip()
                        ]
                        for m in parsed_m:
                            if m not in current_m:
                                current_m.append(m)
                    else:
                        if not current_m:
                            current_m.append(f"{in_new_gene}_marker1")
                    c_markers_dict[in_new_gene] = current_m
                    st.session_state["custom_markers"] = c_markers_dict

                    set_active_target_gene(
                        in_new_gene, markers=current_m, catalog_df=catalog_df
                    )
                    st.success(
                        f"Successfully added **{in_new_gene}** ({', '.join(current_m)}) to your active targets!"
                    )
                    st.rerun()

        with c_btn_reset:
            custom_genes_list = st.session_state.get("custom_genes", [])
            custom_markers_dict = st.session_state.get("custom_markers", {})
            if custom_genes_list or custom_markers_dict:
                if st.button(
                    "🗑️ Clear Custom Targets",
                    key="btn_clear_custom_targets_step1",
                ):
                    st.session_state["custom_genes"] = []
                    st.session_state["custom_markers"] = {}
                    set_active_target_gene(None, catalog_df=catalog_df)
                    st.rerun()

    # Active Target Summary Card
    st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
    active_gene = st.session_state.get("selected_gene") or "None Selected"
    active_markers = st.session_state.get("selected_markers", [])
    num_m = len(active_markers)

    st.markdown(
        f"""
    <div style="background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 12px 16px;">
        <div style="font-size: 11px; font-weight: 700; text-transform: uppercase; color: #64748B;">🎯 Active Configuration</div>
        <div style="font-size: 14px; font-weight: 600; color: #0F172A; margin-top: 4px;">
            Gene: <span style="color: #2563EB;">{active_gene}</span> &nbsp;|&nbsp; 
            Selected Markers ({num_m}): <span style="color: #059669;">{', '.join(active_markers) if active_markers else 'None'}</span>
        </div>
    </div>
    """,
        unsafe_allow_html=True,
    )

    st.markdown("<br>", unsafe_allow_html=True)
    c_sp, c_next1 = st.columns([0.7, 0.3])
    with c_next1:
        if st.button(
            "Continue to Upload Dataset →",
            type="primary",
            use_container_width=True,
            key="btn_step0_cont",
            disabled=num_m == 0,
        ):
            st.session_state["workflow_step"] = "2 · Upload"
            st.rerun()
