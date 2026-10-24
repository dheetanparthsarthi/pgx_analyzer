"""Step 3 Column Mapping component."""

from collections import Counter
import re
import numpy as np
import pandas as pd
import streamlit as st
from config import SNP_CONFIG, CYP2C9_SNP_CONFIG
from marker_catalog import (
    detect_genes_in_columns,
    find_best_column_match,
    get_column_gene,
    get_default_markers_for_gene,
    load_marker_catalog,
)
from ui.session_manager import save_active_dataset, save_active_session, set_active_target_gene


def render_column_mapping(df_raw: pd.DataFrame) -> None:
    """Render Step 3: Map Columns, allele configuration, and permitted genotypes."""
    st.markdown("## 3️⃣ Map Columns: Gene & Marker Selection")
    st.caption("First select or add your target gene and biomarker(s), then map each marker to your spreadsheet columns. You can manually customize Reference and Variant alleles, define permitted genotypes, and enter multiple rows.")

    # -------------------------------------------------------------------------
    # 1. CASCADE DROPDOWNS: SELECT GENE FIRST, THEN SELECT MARKER(S)
    # -------------------------------------------------------------------------
    catalog_df = load_marker_catalog()
    catalog_genes = sorted(catalog_df['Gene'].unique()) if not catalog_df.empty else []
    custom_genes = st.session_state.get('custom_genes', [])
    all_catalog_genes = sorted(list(set(catalog_genes + custom_genes)))
    gene_options = ["(Select Target Gene)"] + all_catalog_genes

    curr_gene = st.session_state.get('selected_gene')
    g_idx = gene_options.index(curr_gene) if (curr_gene and curr_gene in gene_options) else 0

    if 'map_step_gene_sel' in st.session_state and st.session_state['map_step_gene_sel'] not in gene_options:
        del st.session_state['map_step_gene_sel']

    c_g1, c_g2 = st.columns([1, 1.8])
    with c_g1:
        sel_map_gene_choice = st.selectbox("1️⃣ Select Target Gene:", gene_options, index=g_idx, key="map_step_gene_sel")
        if sel_map_gene_choice != "(Select Target Gene)":
            if sel_map_gene_choice != st.session_state.get('selected_gene'):
                set_active_target_gene(sel_map_gene_choice, markers=None, catalog_df=catalog_df)
                st.rerun()
        else:
            if st.session_state.get('selected_gene') is not None:
                set_active_target_gene(None, markers=[], catalog_df=catalog_df)
                st.rerun()

    with c_g2:
        sel_map_gene = st.session_state.get('selected_gene')
        if sel_map_gene:
            base_markers = get_default_markers_for_gene(sel_map_gene, catalog_df)
            custom_m_for_gene = st.session_state.get('custom_markers', {}).get(sel_map_gene, [])
            available_markers = []
            for m in base_markers + custom_m_for_gene:
                if m not in available_markers:
                    available_markers.append(m)

            curr_selected = st.session_state.get('selected_markers', [])
            valid_defaults = [m for m in curr_selected if m in available_markers]
            if not valid_defaults:
                valid_defaults = list(available_markers)
                st.session_state['selected_markers'] = valid_defaults

            c_btn_m1, c_btn_m2, _ = st.columns([0.3, 0.3, 0.4])
            with c_btn_m1:
                if st.button("⚡ Select All", key="btn_sel_all_markers_step3", help="Select all available biomarkers"):
                    st.session_state['selected_markers'] = list(available_markers)
                    if 'map_step_marker_sel' in st.session_state:
                        del st.session_state['map_step_marker_sel']
                    st.rerun()
            with c_btn_m2:
                if st.button("🧹 Clear All", key="btn_desel_all_markers_step3", help="Deselect all biomarkers"):
                    st.session_state['selected_markers'] = []
                    if 'map_step_marker_sel' in st.session_state:
                        del st.session_state['map_step_marker_sel']
                    st.rerun()

            if 'map_step_marker_sel' in st.session_state:
                if any(m not in available_markers for m in st.session_state['map_step_marker_sel']):
                    del st.session_state['map_step_marker_sel']

            sel_map_markers = st.multiselect(
                f"2️⃣ Select Biomarker(s) for {sel_map_gene}:",
                available_markers,
                default=valid_defaults,
                key="map_step_marker_sel",
                placeholder="Choose biomarkers to analyze or click 'Select All'"
            )
            if sel_map_markers != st.session_state.get('selected_markers'):
                st.session_state['selected_markers'] = sel_map_markers
                for k in list(st.session_state.keys()):
                    if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', '_mapped_', '_col_tracker_')):
                        del st.session_state[k]
        else:
            sel_map_markers = []
            st.info("👆 Please select a target gene to configure biomarkers.")

    # Expandable Section: Add Custom Gene / Biomarker Manually directly in Step 3
    with st.expander("➕ Can't find your Gene or Biomarker? Add Custom Targets Manually", expanded=False):
        st.markdown("Add custom target genes or biomarkers not present in the default catalog:")
        t_mode_step3 = st.radio(
            "Target Addition Type:",
            [f"Add Biomarker(s) to Selected Gene ({sel_map_gene})", "Add a Completely New Gene & Biomarker(s)"],
            horizontal=True,
            key="step3_custom_type_radio"
        )
        if "Completely New" in t_mode_step3:
            c_cg1_3, c_cg2_3 = st.columns([1, 1.8])
            with c_cg1_3:
                in_new_gene_3 = st.text_input("New Gene Symbol / Name:", placeholder="e.g. CYP2D6, DPYD, TPMT", key="in_new_custom_gene_step3").strip().upper()
            with c_cg2_3:
                in_new_markers_3 = st.text_input("Biomarker Name(s) (comma-separated):", placeholder="e.g. *4 (rs3892097), *10 (rs1065852)", key="in_new_custom_markers_step3").strip()
        else:
            in_new_gene_3 = sel_map_gene
            in_new_markers_3 = st.text_input(f"New Biomarker Name(s) for {sel_map_gene} (comma-separated):", placeholder="e.g. *4 (rs28399504), *6", key="in_existing_custom_markers_step3").strip()

        c_btn_add_3, c_btn_reset_3 = st.columns([0.35, 0.65])
        with c_btn_add_3:
            if st.button("➕ Add to Target List", type="secondary", key="btn_add_custom_target_step3"):
                if not in_new_gene_3:
                    st.error("Please provide a valid Gene name.")
                else:
                    c_genes = st.session_state.get('custom_genes', [])
                    if in_new_gene_3 not in c_genes and (catalog_df.empty or in_new_gene_3 not in catalog_df['Gene'].unique()):
                        c_genes.append(in_new_gene_3)
                        st.session_state['custom_genes'] = c_genes

                    c_markers_dict = st.session_state.get('custom_markers', {})
                    current_m = c_markers_dict.get(in_new_gene_3, [])
                    if in_new_markers_3:
                        parsed_m = [m.strip() for m in in_new_markers_3.split(',') if m.strip()]
                        for m in parsed_m:
                            if m not in current_m:
                                current_m.append(m)
                    else:
                        if not current_m:
                            current_m.append(f"{in_new_gene_3}_marker1")
                    c_markers_dict[in_new_gene_3] = current_m
                    st.session_state['custom_markers'] = c_markers_dict

                    set_active_target_gene(in_new_gene_3, markers=current_m, catalog_df=catalog_df)
                    st.success(f"Successfully added **{in_new_gene_3}** ({', '.join(current_m)}) to your active targets!")
                    st.rerun()

        with c_btn_reset_3:
            custom_genes_list_3 = st.session_state.get('custom_genes', [])
            custom_markers_dict_3 = st.session_state.get('custom_markers', {})
            if custom_genes_list_3 or custom_markers_dict_3:
                if st.button("🗑️ Clear Custom Targets", key="btn_clear_custom_targets_step3"):
                    st.session_state['custom_genes'] = []
                    st.session_state['custom_markers'] = {}
                    set_active_target_gene(None, catalog_df=catalog_df)
                    st.rerun()

    if not sel_map_markers:
        st.warning(f"⚠️ Please select at least one biomarker for **{sel_map_gene}** from the dropdown above to map columns.")
    elif df_raw is not None and len(df_raw.columns) > 0:
        all_cols = list(df_raw.columns)
        demo_mapping = st.session_state.get('demographic_mapping', {})
        confirmed_demo_cols = {str(v).strip().lower() for k, v in demo_mapping.items() if v}
        cur_insp = st.session_state.get('upload_inspection', {})
        if cur_insp:
            for k in ['sample_id_col', 'gender_col', 'state_col']:
                if cur_insp.get(k):
                    confirmed_demo_cols.add(str(cur_insp[k]).strip().lower())

        NON_GENOTYPE_KEYWORDS = [
            'sample', 'id', 'gender', 'sex', 'state', 'region', 'place', 
            'date', 'birth', 'dob', 'test', 'family', 'generation', 'native', 
            'unnamed', 'patient', 'name', 'age', 'center', 'hospital', 'clinic', 'india'
        ]

        # Cached extraction of column genotype calls to eliminate UI lag during reruns
        _col_calls_cache = {}
        def get_col_calls_cached(c_name: str):
            if not c_name or c_name == "(Select Column)" or c_name not in df_raw.columns:
                return pd.Series(dtype=str), {}, []
            if c_name not in _col_calls_cache:
                raw_s = df_raw[c_name].dropna().astype(str).str.strip().str.upper()
                valid_s = raw_s[~raw_s.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])]
                norm_calls = valid_s.map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                obs_counts = dict(norm_calls.value_counts())
                det_calls = list(obs_counts.keys())
                _col_calls_cache[c_name] = (norm_calls, obs_counts, det_calls)
            return _col_calls_cache[c_name]

        def is_demographic_or_metadata(col_name: str) -> bool:
            clean = str(col_name).strip().lower()
            if clean.startswith('unnamed:'):
                return True
            if clean in confirmed_demo_cols:
                return True
            has_marker_symbol = any(sym in clean for sym in ['*', 'rs', 'cyp', 'dpyd', 'tpmt', 'vkorc', 'nudt', 'ugt', 'slco', 'hla', 'g6pd'])
            if not has_marker_symbol:
                for kw in NON_GENOTYPE_KEYWORDS:
                    if kw in clean:
                        return True
            return False

        # Clean candidate genotype columns
        genotype_cols = [c for c in all_cols if not is_demographic_or_metadata(c)]
        if not genotype_cols:
            genotype_cols = [c for c in all_cols if not str(c).strip().startswith('Unnamed:')]
        if not genotype_cols:
            genotype_cols = all_cols

        # Check if uploaded columns suggest a specific target gene (e.g. CYP2C9, CYP2D6)
        detected_gene_matches = []
        for c in all_cols:
            c_clean = str(c).upper().replace(' ', '').replace('_', '').replace('-', '')
            for cat_gene in (catalog_df['Gene'].unique() if not catalog_df.empty else []):
                g_clean = str(cat_gene).upper().replace(' ', '').replace('_', '').replace('-', '')
                if g_clean != sel_map_gene.upper().replace(' ', '') and len(g_clean) >= 4 and g_clean in c_clean:
                    if cat_gene not in detected_gene_matches:
                        detected_gene_matches.append(cat_gene)

        # Define helper functions and callbacks for flexible marker mapping
        def extract_alleles_from_calls(calls):
            # Only infer alleles automatically for ordinary two-character,
            # single-base calls. Multi-character labels (e.g. DEL/INS) and
            # numeric calls require the user to enter the actual alleles.
            normalized_calls = [
                str(call).replace('/', '').replace('|', '').replace(' ', '').upper()
                for call in calls
            ]
            if not normalized_calls or any(len(call) != 2 or not call.isalpha() for call in normalized_calls):
                return '', ''
            al_counts = Counter()
            for call in normalized_calls:
                for allele in call:
                    al_counts[allele] += 1
            most = [al for al, _ in al_counts.most_common()]
            r = most[0] if len(most) > 0 else 'A'
            v = most[1] if len(most) > 1 else ''
            return r, v


        def get_marker_defaults(marker_label: str, sel_gene: str, detected_calls: list = None, col_name: str = ""):
            """
            Generates clean, user-driven biomarker configuration without hardcoding.
            Permitted genotypes are NOT pre-filled with catalog assumptions,
            leaving full flexibility to the user to type or click observed genotypes.
            """
            clean_label = str(marker_label).strip()
            if ' ' in clean_label and '(' in clean_label:
                base_name = clean_label.split(' ')[0]
            else:
                base_name = clean_label

            auto_r, auto_v = '', ''
            if detected_calls and len(detected_calls) > 0:
                cand_r, cand_v = extract_alleles_from_calls(detected_calls)
                if cand_r:
                    auto_r = cand_r
                if cand_v:
                    auto_v = cand_v

            badge_short = f"Column: {col_name}" if col_name else "Spreadsheet Mapping"
            if detected_calls:
                badge_detail = f"Observed calls in data: <b>{', '.join(detected_calls)}</b>"
            else:
                badge_detail = f"Configure permitted genotypes and alleles for <b>{base_name}</b>"

            # DO NOT pre-fill with hardcoded catalog strings - let the user enter or click chips!
            def_permitted = ""
            return base_name, auto_r, auto_v, def_permitted, badge_short, badge_detail

        # Define callbacks to update widget session state safely before widgets instantiate
        def cb_switch_target_gene(new_g, new_markers):
            set_active_target_gene(new_g, markers=new_markers, catalog_df=catalog_df)

        def cb_remove_marker_card(idx_to_remove):
            cur_markers = list(st.session_state.get('selected_markers', []))
            if 0 <= idx_to_remove < len(cur_markers):
                cur_markers.pop(idx_to_remove)
                st.session_state['selected_markers'] = cur_markers
                prefixes = ['w_map_name_', 'w_map_col_', 'w_map_ref_', 'w_map_var_', 'w_map_perm_', 'sel_ref_choice_', 'sel_var_choice_', 'custom_ref_in_', 'custom_var_in_']
                for i in range(idx_to_remove, len(cur_markers) + 1):
                    for prefix in prefixes:
                        cur_k = f"{prefix}{i}"
                        next_k = f"{prefix}{i+1}"
                        if next_k in st.session_state:
                            st.session_state[cur_k] = st.session_state[next_k]
                            del st.session_state[next_k]
                        elif cur_k in st.session_state:
                            del st.session_state[cur_k]

        def cb_add_new_marker_card(gene_prefix):
            cur_markers = list(st.session_state.get('selected_markers', []))
            new_num = len(cur_markers) + 1
            new_marker_lbl = f"{gene_prefix}_Marker_{new_num}"
            while new_marker_lbl in cur_markers:
                new_num += 1
                new_marker_lbl = f"{gene_prefix}_Marker_{new_num}"
            cur_markers.append(new_marker_lbl)
            st.session_state['selected_markers'] = cur_markers

        def cb_toggle_permitted_call(p_k, call, r_k, v_k, sel_r_k, sel_v_k):
            cur_p = str(st.session_state.get(p_k, "")).strip()
            calls = [g.strip().upper() for g in cur_p.split(',') if g.strip()]
            call_u = str(call).strip().upper()
            if call_u in calls:
                calls.remove(call_u)
            else:
                calls.append(call_u)
            new_p = ", ".join(calls)
            st.session_state[p_k] = new_p
            if calls:
                r_c, v_c = extract_alleles_from_calls(calls)
                if r_c:
                    st.session_state[r_k] = r_c
                    if sel_r_k:
                        st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                if v_c:
                    st.session_state[v_k] = v_c
                    if sel_v_k:
                        st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'

        def cb_allow_all_detected(p_k, calls_list, r_k, v_k, sel_r_k, sel_v_k):
            st.session_state[p_k] = ", ".join(calls_list)
            r_c, v_c = extract_alleles_from_calls(calls_list)
            if r_c:
                st.session_state[r_k] = r_c
                if sel_r_k:
                    st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
            if v_c:
                st.session_state[v_k] = v_c
                if sel_v_k:
                    st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'

        def cb_clear_permitted(p_k):
            st.session_state[p_k] = ""

        def cb_normalize_allele_input(state_key):
            st.session_state[state_key] = str(st.session_state.get(state_key, "")).strip().upper()

        def cb_automap_detected(p_k, calls_list, r_k, v_k, sel_r_k, sel_v_k):
            st.session_state[p_k] = ", ".join(calls_list)
            r_c, v_c = extract_alleles_from_calls(calls_list)
            if r_c:
                st.session_state[r_k] = r_c
                if sel_r_k:
                    st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                    if r_c not in ['A', 'C', 'G', 'T']:
                        idx_str = sel_r_k.split('_')[-1]
                        st.session_state[f"custom_ref_in_{idx_str}"] = r_c
            if v_c:
                st.session_state[v_k] = v_c
                if sel_v_k:
                    st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'
                    if v_c not in ['A', 'C', 'G', 'T', 'None']:
                        idx_str = sel_v_k.split('_')[-1]
                        st.session_state[f"custom_var_in_{idx_str}"] = v_c

        def cb_add_unmapped(p_k, unexp_list, r_k, v_k, sel_r_k, sel_v_k):
            cur_p = str(st.session_state.get(p_k, "")).strip()
            add_str = ", ".join(unexp_list)
            if cur_p:
                if not cur_p.endswith(','):
                    cur_p += ", "
                new_p = cur_p + add_str
            else:
                new_p = add_str
            st.session_state[p_k] = new_p
            all_c = [g.strip().upper() for g in new_p.split(',') if g.strip()]
            r_c, v_c = extract_alleles_from_calls(all_c)
            if r_c:
                st.session_state[r_k] = r_c
                if sel_r_k:
                    st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                    if r_c not in ['A', 'C', 'G', 'T']:
                        idx_str = sel_r_k.split('_')[-1]
                        st.session_state[f"custom_ref_in_{idx_str}"] = r_c
            if v_c:
                st.session_state[v_k] = v_c
                if sel_v_k:
                    st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'
                    if v_c not in ['A', 'C', 'G', 'T', 'None']:
                        idx_str = sel_v_k.split('_')[-1]
                        st.session_state[f"custom_var_in_{idx_str}"] = v_c

        def cb_align_alleles(r_k, v_k, r_c, v_c, sel_r_k, sel_v_k):
            st.session_state[r_k] = r_c
            st.session_state[v_k] = v_c
            if sel_r_k:
                st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                if r_c not in ['A', 'C', 'G', 'T']:
                    idx_str = sel_r_k.split('_')[-1]
                    st.session_state[f"custom_ref_in_{idx_str}"] = r_c
            if sel_v_k:
                st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'
                if v_c not in ['A', 'C', 'G', 'T', 'None']:
                    idx_str = sel_v_k.split('_')[-1]
                    st.session_state[f"custom_var_in_{idx_str}"] = v_c

        def cb_extra_automap(p_k, r_k, v_k, calls_list):
            st.session_state[p_k] = ", ".join(calls_list)
            r_c, v_c = extract_alleles_from_calls(calls_list)
            if r_c:
                st.session_state[r_k] = r_c
            if v_c:
                st.session_state[v_k] = v_c

        def cb_extra_add_unmapped(p_k, unexp_list, r_k, v_k):
            cur_p = str(st.session_state.get(p_k, "")).strip()
            add_str = ", ".join(unexp_list)
            if cur_p:
                if not cur_p.endswith(','):
                    cur_p += ", "
                new_p = cur_p + add_str
            else:
                new_p = add_str
            st.session_state[p_k] = new_p
            all_c = [g.strip().upper() for g in new_p.split(',') if g.strip()]
            r_c, v_c = extract_alleles_from_calls(all_c)
            if r_c:
                st.session_state[r_k] = r_c
            if v_c:
                st.session_state[v_k] = v_c


        c_opt_filter, c_mode_sel = st.columns([0.45, 0.55])
        with c_opt_filter:
            show_all_cols = st.checkbox("Show all spreadsheet columns in options", value=False, key="cb_show_all_cols_step3", help="Enable if your genotype column was auto-classified as metadata")
        candidate_cols = [c for c in (all_cols if show_all_cols else genotype_cols) if c not in confirmed_demo_cols and not str(c).startswith('Unnamed:')]
        if not candidate_cols and all_cols:
            candidate_cols = [c for c in all_cols if not str(c).startswith('Unnamed:')]
        col_options = ["(Select Column)"] + candidate_cols

        with c_mode_sel:
            mapping_view_mode = st.radio(
                "Mapping Interface:",
                ["🗂️ Interactive Marker Cards", "📋 Multi-Row Table Editor"],
                horizontal=True,
                key="step3_mapping_view_mode"
            )

        mapped_marker_dict = {}
        mapping_has_errors = False
        extra_rows = st.session_state.get('extra_mapping_rows', [])

        # -------------------------------------------------------------------------
        # -------------------------------------------------------------------------
        # MODE A: INTERACTIVE MARKER CARDS (WITH LIVE DATA QC & MULTI-VALUE COMBINATIONS)
        # -------------------------------------------------------------------------
        if mapping_view_mode == "🗂️ Interactive Marker Cards":
            st.markdown("#### 🧬 Biomarker Mapping & Genotype Combination Rules")
            st.caption("Map spreadsheet columns and configure any combination of permitted genotypes (e.g. GA, GG, AG or CC, CT, TT). Check live data breakdown and verify cohort calls instantly.")

            detected_file_genes = detect_genes_in_columns(all_cols)
            if detected_file_genes and sel_map_gene.upper() not in [g.upper() for g in detected_file_genes]:
                found_genes_str = ", ".join(detected_file_genes)
                st.markdown(f"""
                <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 12px 16px; margin-bottom: 14px;">
                    <div style="display: flex; align-items: center; gap: 8px;">
                        <span style="font-size: 18px;">🚨</span>
                        <div style="font-weight: 800; color: #991B1B; font-size: 14px;">TARGET GENE & COLUMN HEADER DISCREPANCY DETECTED</div>
                    </div>
                    <div style="color: #7F1D1D; font-size: 13px; margin-top: 4px;">
                        You are mapping markers for <b>{sel_map_gene}</b>, but spreadsheet column headers contain <b>{found_genes_str}</b>.<br>
                        If you uploaded data for <b>{detected_file_genes[0]}</b>, click the button below to switch your target gene:
                    </div>
                </div>
                """, unsafe_allow_html=True)
                c_sw_map, _ = st.columns([0.4, 0.6])
                with c_sw_map:
                    if st.button(f"👉 Switch Target Gene to {detected_file_genes[0]}", key="btn_switch_gene_step3_banner", type="primary"):
                        set_active_target_gene(detected_file_genes[0], catalog_df=catalog_df)
                        st.rerun()

            for idx, marker_label in enumerate(sel_map_markers):
                col_key = f"w_map_col_{idx}"
                name_key = f"w_map_name_{idx}"
                ref_key = f"w_map_ref_{idx}"
                var_key = f"w_map_var_{idx}"
                perm_key = f"w_map_perm_{idx}"
                tracker_key = f"_mapped_tracker_{idx}"
                col_tracker_key = f"_col_tracker_{idx}"

                # If marker at this index changed (e.g. user switched gene or marker selection), reset stale keys
                if st.session_state.get(tracker_key) != marker_label:
                    for k in [col_key, name_key, ref_key, var_key, perm_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}", f"custom_ref_in_{idx}", f"custom_var_in_{idx}", col_tracker_key]:
                        if k in st.session_state:
                            del st.session_state[k]
                    st.session_state[tracker_key] = marker_label

                # Determine best matching column
                best_match = find_best_column_match(marker_label, marker_label.split(' ')[0], sel_map_gene, col_options)
                if col_key not in st.session_state or st.session_state[col_key] not in col_options:
                    st.session_state[col_key] = best_match if (best_match and best_match in col_options) else "(Select Column)"

                sel_col = st.session_state[col_key]

                # Check if user changed the column selection in the dropdown
                if col_tracker_key in st.session_state and st.session_state[col_tracker_key] != sel_col:
                    for k in [ref_key, var_key, perm_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}", f"custom_ref_in_{idx}", f"custom_var_in_{idx}"]:
                        if k in st.session_state:
                            del st.session_state[k]
                    st.session_state[col_tracker_key] = sel_col
                else:
                    st.session_state[col_tracker_key] = sel_col

                # Live Data Inspection for this column (fast cached lookup)
                norm_calls, observed_counts, detected_calls_list = get_col_calls_cached(sel_col)
                tot_samples_col = len(norm_calls)

                # Determine defaults prioritizing detected calls in the data
                std_name, def_ref, def_var, def_permitted, badge_short, badge_detail = get_marker_defaults(marker_label, sel_map_gene, detected_calls_list, col_name=(sel_col if sel_col != "(Select Column)" else ""))

                # Initialize session state for this card
                if name_key not in st.session_state:
                    st.session_state[name_key] = std_name
                if ref_key not in st.session_state:
                    st.session_state[ref_key] = def_ref
                if var_key not in st.session_state:
                    st.session_state[var_key] = def_var
                if perm_key not in st.session_state:
                    st.session_state[perm_key] = def_permitted

                # Migrate values saved by the old dropdown + "Other" controls.
                if st.session_state.get(ref_key) == 'Other (Manual)':
                    st.session_state[ref_key] = st.session_state.get(f"custom_ref_in_{idx}", "") or def_ref
                if st.session_state.get(var_key) in {'Other (Manual)', 'None'}:
                    st.session_state[var_key] = st.session_state.get(f"custom_var_in_{idx}", "") or def_var

                cur_marker_name = st.session_state.get(name_key, std_name)
                if sel_col == "(Select Column)":
                    mapping_has_errors = True
                    col_badge_html = '<span style="background-color: #FEF3C7; color: #92400E; font-size: 0.8rem; font-weight: 600; padding: 3px 10px; border-radius: 9999px; border: 1px solid #FDE68A;">Column: <i>Not Selected</i></span>'
                    obs_str = '<span style="color: #92400E; font-style: italic;">Select spreadsheet column below to inspect observed genotypes</span>'
                    card_border = "#F59E0B"
                else:
                    col_badge_html = f'<span style="background-color: #EFF6FF; color: #1D4ED8; font-size: 0.8rem; font-weight: 600; padding: 3px 10px; border-radius: 9999px; border: 1px solid #BFDBFE;">Column: <code>{sel_col}</code></span>'
                    obs_str = ', '.join([f"<b>{c}</b> ({cnt:,})" for c, cnt in observed_counts.items()]) if observed_counts else "<i>No calls detected in this column</i>"
                    card_border = "#3B82F6"

                st.markdown(f"""
                <div class="clean-card-box" style="margin-bottom: 8px; border-left: 4px solid {card_border}; padding: 12px 16px; background: #F8FAFC; border-radius: 8px; border: 1px solid #E2E8F0;">
                    <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;">
                        <h5 style="margin:0; color:#0F172A; font-weight:700;">🧬 Biomarker #{idx+1}: {cur_marker_name}</h5>
                        {col_badge_html}
                    </div>
                    <div style="margin-top: 6px; font-size: 0.85rem; color: #475569;">
                        Observed Genotypes in column: {obs_str}
                    </div>
                </div>
                """, unsafe_allow_html=True)

                if sel_col != "(Select Column)":
                    # Check for header discrepancies (different gene, star allele, or rsID)
                    col_gene = get_column_gene(sel_col)
                    is_diff_gene = bool(col_gene and sel_map_gene and col_gene.upper() != sel_map_gene.upper())

                    m_star = re.search(r'\*\d+[a-zA-Z]?', cur_marker_name)
                    c_star = re.search(r'\*\d+[a-zA-Z]?', str(sel_col))
                    is_diff_star = bool(m_star and c_star and m_star.group(0).lower() != c_star.group(0).lower())

                    m_rs = re.search(r'rs\d+', cur_marker_name.lower())
                    c_rs = re.search(r'rs\d+', str(sel_col).lower())
                    is_diff_rs = bool(m_rs and c_rs and m_rs.group(0) != c_rs.group(0))

                    if is_diff_gene:
                        st.markdown(f"""
                        <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 10px 14px; margin-top: 6px; margin-bottom: 8px;">
                            <div style="font-weight: 800; color: #991B1B; font-size: 13px;">🚨 GENE HEADER MISMATCH: Column <code>{sel_col}</code> belongs to {col_gene}, but active gene is {sel_map_gene}!</div>
                            <div style="color: #7F1D1D; font-size: 12px; margin-top: 3px;">
                                Target marker is <b>{cur_marker_name}</b> ({sel_map_gene}), but selected column <code>{sel_col}</code> represents gene <b>{col_gene}</b>.<br>
                                Please verify if you selected the wrong column from the dropdown below.
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
                    elif is_diff_star:
                        st.markdown(f"""
                        <div style="background: #FFFBEB; border: 2px solid #F59E0B; border-radius: 8px; padding: 10px 14px; margin-top: 6px; margin-bottom: 8px;">
                            <div style="font-weight: 800; color: #B45309; font-size: 13px;">⚠️ STAR ALLELE MISMATCH: Target is {m_star.group(0)}, but column header is {c_star.group(0)}!</div>
                            <div style="color: #92400E; font-size: 12px; margin-top: 3px;">
                                Target marker <b>{cur_marker_name}</b> is for <b>{m_star.group(0)}</b>, but column header <code>{sel_col}</code> indicates <b>{c_star.group(0)}</b>.
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
                    elif is_diff_rs:
                        st.markdown(f"""
                        <div style="background: #FFFBEB; border: 2px solid #F59E0B; border-radius: 8px; padding: 10px 14px; margin-top: 6px; margin-bottom: 8px;">
                            <div style="font-weight: 800; color: #B45309; font-size: 13px;">⚠️ RSID MISMATCH: Target is {m_rs.group(0)}, but column header has {c_rs.group(0)}!</div>
                            <div style="color: #92400E; font-size: 12px; margin-top: 3px;">
                                Please verify if column <code>{sel_col}</code> corresponds to marker <b>{cur_marker_name}</b>.
                            </div>
                        </div>
                        """, unsafe_allow_html=True)

                c_top1, c_top2, c_top3 = st.columns([1.1, 1.5, 0.4])
                with c_top1:
                    user_marker_name = st.text_input(
                        "Biomarker Name:",
                        value=st.session_state.get(name_key, std_name),
                        key=name_key,
                        help="Edit marker name (e.g. CYP2C9*2, CYP2C19*2, rs1799853)"
                    ).strip()
                    if not user_marker_name:
                        user_marker_name = std_name
                with c_top2:
                    def_col_idx = col_options.index(st.session_state[col_key]) if st.session_state[col_key] in col_options else 0
                    new_sel_col = st.selectbox(
                        "Spreadsheet Column:",
                        col_options,
                        index=def_col_idx,
                        key=col_key,
                        help=f"Spreadsheet column for {user_marker_name}"
                    )
                    if new_sel_col != sel_col:
                        sel_col = new_sel_col
                        st.session_state[col_key] = new_sel_col
                        for k in [ref_key, var_key, perm_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}", f"custom_ref_in_{idx}", f"custom_var_in_{idx}"]:
                            if k in st.session_state:
                                del st.session_state[k]
                        st.session_state[col_tracker_key] = new_sel_col
                        st.rerun()
                with c_top3:
                    st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
                    st.button("🗑️ Remove", key=f"btn_remove_marker_{idx}", on_click=cb_remove_marker_card, args=(idx,), help=f"Remove marker from mapping")

                # Reference and variant alleles are directly editable.
                c2_m, c3_m = st.columns(2)
                with c2_m:
                    ref_val = st.text_input(
                        "Reference allele:",
                        key=ref_key,
                        max_chars=16,
                        on_change=cb_normalize_allele_input,
                        args=(ref_key,),
                        help="Enter the reference allele for this marker, for example A, C, G, T, DEL, or INS.",
                    ).strip().upper()

                with c3_m:
                    var_val = st.text_input(
                        "Variant allele:",
                        key=var_key,
                        max_chars=16,
                        on_change=cb_normalize_allele_input,
                        args=(var_key,),
                        help="Enter the variant allele for this marker, for example A, C, G, T, DEL, or INS.",
                    ).strip().upper()

                # Permitted Genotypes Input Row
                c_perm_in, c_perm_info = st.columns([1.6, 1.4])
                with c_perm_in:
                    perm_val = st.text_input(
                        "Permitted Genotypes (comma-separated):",
                        value=st.session_state.get(perm_key, def_permitted),
                        placeholder="e.g. CC, CT, TT or GA, GG, AA (or click chips below)",
                        key=perm_key,
                        help="Type any combinations e.g. CC, CT or GG, GA or AA or TT, or click the observed genotype chips below"
                    )
                with c_perm_info:
                    hom_ref_str = f"{ref_val}{ref_val}" if ref_val else "N/A"
                    het_str = f"{ref_val}{var_val}" if (ref_val and var_val and var_val != 'None') else "N/A"
                    hom_var_str = f"{var_val}{var_val}" if (var_val and var_val != 'None') else "N/A"
                    st.markdown(f"""
                    <div style="background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px; margin-top: 5px;">
                        <div style="font-size: 11px; font-weight: 700; color: #64748B; text-transform: uppercase;">🧬 Genotype-to-Allele Mapping:</div>
                        <div style="font-size: 13px; color: #0F172A; margin-top: 4px;">
                            Hom Ref = <b>{hom_ref_str}</b> &nbsp;|&nbsp; 
                            Het = <b>{het_str}</b> &nbsp;|&nbsp; 
                            Hom Var = <b>{hom_var_str}</b>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                # Interactive Call Toggle Chips (Manual Selection directly from Observed Data)
                if detected_calls_list:
                    st.markdown("<div style='font-size: 11px; font-weight: 700; color: #475569; margin-top: 6px; margin-bottom: 4px; text-transform: uppercase;'>Observed Calls in this Column — Click chip to toggle in/out of permitted rules:</div>", unsafe_allow_html=True)
                    cur_perm_tokens_set = {x.strip().upper().replace('/', '').replace('|', '').replace(' ', '') for x in str(st.session_state.get(perm_key, '')).split(',') if x.strip()}

                    n_items = len(detected_calls_list)
                    chunk_size = 5
                    for c_start in range(0, n_items, chunk_size):
                        c_chunk = detected_calls_list[c_start:c_start + chunk_size]
                        is_first_row = (c_start == 0)

                        col_widths = [1] * len(c_chunk)
                        if is_first_row:
                            col_widths.extend([1.3, 0.9])
                            chip_cols = st.columns(col_widths)
                            for c_i, call_item in enumerate(c_chunk):
                                cnt = observed_counts.get(call_item, 0)
                                is_active = call_item in cur_perm_tokens_set
                                btn_lbl = f"✅ {call_item} ({cnt:,})" if is_active else f"➕ {call_item} ({cnt:,})"
                                chip_cols[c_i].button(
                                    btn_lbl,
                                    key=f"btn_chip_{idx}_{call_item}",
                                    type="primary" if is_active else "secondary",
                                    use_container_width=True,
                                    on_click=cb_toggle_permitted_call,
                                    args=(perm_key, call_item, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}"),
                                    help=f"Click to {'remove' if is_active else 'permit'} {call_item}"
                                )
                            chip_cols[-2].button("⚡ Allow All", key=f"btn_chip_all_{idx}", use_container_width=True,
                                                on_click=cb_allow_all_detected, args=(perm_key, detected_calls_list, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}"), help="Add all observed calls to permitted list")
                            chip_cols[-1].button("🧹 Clear", key=f"btn_chip_clr_{idx}", use_container_width=True,
                                                on_click=cb_clear_permitted, args=(perm_key,), help="Clear permitted list")
                        else:
                            chip_cols = st.columns(col_widths)
                            for c_i, call_item in enumerate(c_chunk):
                                cnt = observed_counts.get(call_item, 0)
                                is_active = call_item in cur_perm_tokens_set
                                btn_lbl = f"✅ {call_item} ({cnt:,})" if is_active else f"➕ {call_item} ({cnt:,})"
                                chip_cols[c_i].button(
                                    btn_lbl,
                                    key=f"btn_chip_{idx}_{call_item}",
                                    type="primary" if is_active else "secondary",
                                    use_container_width=True,
                                    on_click=cb_toggle_permitted_call,
                                    args=(perm_key, call_item, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}"),
                                    help=f"Click to {'remove' if is_active else 'permit'} {call_item}"
                                )

                # Parse permitted genotypes & expanded set
                raw_perm_list = [g.strip().upper() for g in perm_val.split(',') if g.strip()]
                valid_clean_set = {g.replace('/', '').replace('|', '').replace(' ', '') for g in raw_perm_list}
                expanded_set = set(valid_clean_set)
                for g in valid_clean_set:
                    if len(g) == 2:
                        expanded_set.add(g[1] + g[0])

                # LIVE DATA AUDIT & STRICT DETECTION / FLAGGING ENGINE
                # CHECK 1: Detect Non-Genomic / Demographic / Text Column
                METADATA_WORDS = {
                    'MALE', 'FEMALE', 'MAN', 'WOMAN', 'UNKNOWN', 'NULL', 'NONE', 'NAN',
                    'CASE', 'CONTROL', 'YES', 'NO', 'TRUE', 'FALSE', 'POS', 'NEG',
                    'POSITIVE', 'NEGATIVE', 'NORMAL', 'MUTANT', 'HETEROZYGOUS', 'HOMOZYGOUS',
                    'INDIA', 'DELHI', 'PUNJAB', 'MAHARASHTRA', 'TAMIL', 'GUJARAT', 'KERALA',
                    'PATIENT', 'SAMPLE', 'CONTROL_GROUP', 'CASE_GROUP'
                }
                suspicious_calls = []
                for c in detected_calls_list:
                    c_u = str(c).strip().upper()
                    if c_u in METADATA_WORDS or len(c_u) > 6 or any(ch.isdigit() for ch in c_u) or any(ch not in 'ACGTID/*|- ' for ch in c_u):
                        suspicious_calls.append(c)

                is_metadata_col = bool(suspicious_calls and (len(suspicious_calls) >= len(detected_calls_list) * 0.5 or any(c.strip().upper() in METADATA_WORDS for c in detected_calls_list)))
                if is_metadata_col:
                    mapping_has_errors = True
                    st.markdown(f"""
                    <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 12px; margin-top: 10px; margin-bottom: 8px;">
                        <div style="font-weight: 800; color: #991B1B; font-size: 14px;">🚨 WRONG COLUMN DETECTED: Non-Genotype Data in `{sel_col}`</div>
                        <div style="color: #7F1D1D; font-size: 13px; margin-top: 4px;">
                            Column <code>{sel_col}</code> appears to contain demographic or text values (e.g. <b>{suspicious_calls[:4]}</b>) rather than genetic marker calls. Please select a valid genotype column above.
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                # CHECK 2: Severe Allele Disconnect / Wrong Marker Column
                if raw_perm_list:
                    col_nucleotides = {ch for call in detected_calls_list for ch in str(call).upper() if ch in 'ACGT'}
                    configured_alleles = {ch for ch in [ref_val, var_val] if ch in 'ACGT'}
                    is_severe_disconnect = bool(col_nucleotides and configured_alleles and not col_nucleotides.intersection(configured_alleles))
                    if is_severe_disconnect and not is_metadata_col:
                        st.markdown(f"""
                        <div style="background: #FFF7ED; border: 2px solid #F97316; border-radius: 8px; padding: 12px; margin-top: 10px; margin-bottom: 8px;">
                            <div style="font-weight: 800; color: #9A3412; font-size: 14px;">🚨 SEVERE ALLELE MISMATCH: Possible Wrong Column Selected!</div>
                            <div style="color: #7C2D12; font-size: 13px; margin-top: 4px;">
                                Observed calls in <code>{sel_col}</code> are <b>{detected_calls_list}</b> (nucleotides: <b>{sorted(col_nucleotides)}</b>).<br>
                                These have <b>ZERO overlap</b> with configured alleles for <b>{user_marker_name}</b> (Ref: <b>{ref_val}</b>, Var: <b>{var_val}</b>).<br>
                                Please verify if you selected the wrong column, or click below to align your alleles to this column if intended.
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
                        cand_col_r, cand_col_v = extract_alleles_from_calls(detected_calls_list)
                        st.button(
                            f"⚡ Align Alleles to this Column (Ref: {cand_col_r}, Var: {cand_col_v})",
                            key=f"btn_align_severe_{idx}",
                            on_click=cb_align_alleles,
                            args=(ref_key, var_key, cand_col_r, cand_col_v, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}")
                        )

                # CHECK 3: Extra / Unmapped Calls in Dataset
                unexpected_calls = [k for k in observed_counts.keys() if k not in expanded_set] if (observed_counts and raw_perm_list) else []
                tot_unmapped_cnt = sum(observed_counts[u] for u in unexpected_calls) if unexpected_calls else 0
                unmapped_pct = (tot_unmapped_cnt / tot_samples_col * 100) if tot_samples_col > 0 else 0

                cand_r, cand_v = extract_alleles_from_calls(raw_perm_list)
                perm_chars = {ch for g in raw_perm_list for ch in g if ch.isalpha()}
                ref_mismatch = bool(ref_val and ref_val not in perm_chars) if raw_perm_list else False
                var_mismatch = bool(var_val and var_val != 'None' and var_val not in perm_chars) if raw_perm_list else False

                missing_allele_definition = not ref_val or not var_val
                if missing_allele_definition:
                    st.error("⛔ Enter both a reference allele and a variant allele for this marker.")
                    mapping_has_errors = True

                identical_alleles = bool(ref_val == var_val and ref_val != '')

                st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)
                c_flg1, c_flg2 = st.columns([1.5, 1.5])
                if sel_col == "(Select Column)":
                    with c_flg1:
                        st.info(f"👆 **Action Required:** Please select the spreadsheet column for **{user_marker_name}** from the dropdown above.")
                    with c_flg2:
                        st.warning(f"⚪ **Status:** Biomarker **{user_marker_name}** is waiting for a column assignment. Cohort verification and QC checks will activate once selected.")
                else:
                    with c_flg1:
                        st.markdown(f"**Data Audit for `{sel_col}` ({tot_samples_col:,} valid samples):**")
                        if observed_counts:
                            check_rows = []
                            for call, cnt in sorted(observed_counts.items(), key=lambda x: x[1], reverse=True):
                                pct = (cnt / tot_samples_col * 100) if tot_samples_col > 0 else 0
                                is_covered = call in expanded_set
                                if is_covered:
                                    status_str = "✅ Permitted"
                                elif not raw_perm_list:
                                    status_str = "⚪ Pending User Entry"
                                else:
                                    status_str = "🚨 Flagged: Unmapped (Extra)"
                                check_rows.append({
                                    'Observed Genotype': call,
                                    'Sample Count': cnt,
                                    'Cohort %': f"{pct:.1f}%",
                                    'Mapping Status': status_str
                                })
                            st.dataframe(pd.DataFrame(check_rows), use_container_width=True, hide_index=True)
                        else:
                            st.info(f"ℹ️ No valid genotype calls found in column `{sel_col}`.")

                    with c_flg2:
                        if not raw_perm_list:
                            st.info(f"ℹ️ **Action Needed:** Please enter permitted genotypes for **{user_marker_name}** above (type calls or click **[⚡ Allow All]** / chip buttons).")
                            mapping_has_errors = True
                        elif unexpected_calls:
                            st.markdown(f"""
                            <div style="background: #FEF2F2; border: 1.5px solid #F87171; border-radius: 8px; padding: 12px; margin-bottom: 8px;">
                                <div style="font-weight: 700; color: #991B1B; font-size: 13px;">🚨 FLAGGED: Extra / Unmapped Entry in Dataset</div>
                                <div style="color: #7F1D1D; font-size: 12px; margin-top: 4px;">
                                    Found <b>{len(unexpected_calls)} unmapped call(s)</b>: <code>{unexpected_calls}</code> ({tot_unmapped_cnt:,} sample(s), {unmapped_pct:.1f}% of cohort).<br>
                                    Permitted genotypes are: <code>{perm_val}</code>.<br>
                                    <i>These samples will be treated as invalid/missing during QC unless you explicitly allow them.</i>
                                </div>
                            </div>
                            """, unsafe_allow_html=True)
                            unexp_str = ", ".join(unexpected_calls)
                            st.button(
                                f"➕ Allow Flagged Call(s) ({unexp_str}) into Permitted Genotypes",
                                key=f"btn_add_unexp_{idx}",
                                on_click=cb_add_unmapped,
                                args=(perm_key, unexpected_calls, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}")
                            )
                        elif tot_samples_col > 0:
                            st.markdown(f"""
                            <div style="background: #F0FDF4; border: 1.5px solid #86EFAC; border-radius: 8px; padding: 12px; margin-bottom: 8px;">
                                <div style="font-weight: 700; color: #166534; font-size: 13px;">✅ 100% Cohort Verified: All {tot_samples_col:,} Samples Permitted</div>
                                <div style="color: #14532D; font-size: 12px; margin-top: 3px;">
                                    All observed calls in <code>{sel_col}</code> match your configured genotype rules (<code>{', '.join(raw_perm_list)}</code>).
                                </div>
                            </div>
                            """, unsafe_allow_html=True)

                        # Flag: Allele Inconsistency
                        if (ref_mismatch or var_mismatch) and cand_r:
                            st.markdown(f"""
                            <div style="background: #FFFBEB; border: 1.5px solid #FCD34D; border-radius: 8px; padding: 10px; margin-bottom: 8px;">
                                <div style="font-weight: 700; color: #92400E; font-size: 13px;">⚠️ FLAG: Allele Rule Inconsistency</div>
                                <div style="color: #78350F; font-size: 12px; margin-top: 3px;">
                                    Configured Ref (<b>{ref_val}</b>) or Var (<b>{var_val}</b>) allele is absent from permitted genotypes. Expected alleles: Ref = <b>{cand_r}</b>, Var = <b>{cand_v}</b>.
                                </div>
                            </div>
                            """, unsafe_allow_html=True)
                            st.button(
                                f"👉 Align Alleles to Ref: {cand_r}, Var: {cand_v}",
                                key=f"btn_align_alleles_{idx}",
                                on_click=cb_align_alleles,
                                args=(ref_key, var_key, cand_r, cand_v, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}")
                            )

                        # A marker can be monomorphic in this cohort. Keep
                        # the mapping and report its single observed allele;
                        # statistical analysis handles it as monomorphic.
                        if identical_alleles:
                            st.info(f"ℹ️ Monomorphic marker: both allele fields are {ref_val}. Genotype/allele counts will be reported for this single allele; HWE and clinical translation are not informative for this marker.")

                if sel_col != "(Select Column)" and sel_col in df_raw.columns:
                    v_al_stored = 'A' if var_val == 'None' else var_val.strip().upper()
                    mapped_marker_dict[sel_col] = {
                        'std_name': user_marker_name,
                        'ref_allele': ref_val.strip().upper(),
                        'var_allele': v_al_stored,
                        'marker_label': user_marker_name,
                        'permitted_genotypes': raw_perm_list
                    }
                else:
                    mapping_has_errors = True

            # Button to add new marker card row
            st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
            c_add_card, _ = st.columns([0.4, 0.6])
            with c_add_card:
                st.button("➕ Add Marker Mapping Row", key="btn_add_marker_card_row", on_click=cb_add_new_marker_card, args=(sel_map_gene,))

            # Render Extra Custom Mapping Rows
            if extra_rows:
                st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                st.markdown("#### ➕ Custom Marker Mapping Rows")
                rows_to_remove = []
                for r_idx, row_item in enumerate(extra_rows):
                    r_id = row_item.get('id', f"r_{r_idx}")
                    st.markdown(f"""
                    <div class="clean-card-box" style="margin-bottom: 8px; border-left: 4px solid #8B5CF6;">
                        <div style="display: flex; justify-content: space-between; align-items: center;">
                            <h5 style="margin:0; color:#0F172A; font-weight:700;">🧬 Custom Row #{r_idx+1}: {row_item.get('marker_name', 'Custom Marker')}</h5>
                            <span style="background-color: #F3E8FF; color: #7E22CE; font-size: 0.8rem; font-weight: 600; padding: 2px 8px; border-radius: 9999px;">Manual Entry</span>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                    c_r1, c_r2, c_r3, c_r4, c_del = st.columns([1.4, 1.4, 0.6, 0.6, 0.4])
                    with c_r1:
                        c_name = st.text_input("Marker Name", value=row_item.get('marker_name', f"Marker_{r_idx+1}"), key=f"ex_name_{r_id}")
                        row_item['marker_name'] = c_name
                    with c_r2:
                        def_col_ex = row_item.get('col', "(Select Column)")
                        col_idx_ex = col_options.index(def_col_ex) if def_col_ex in col_options else 0
                        c_col = st.selectbox("Genotype Column", col_options, index=col_idx_ex, key=f"ex_col_{r_id}")
                        row_item['col'] = c_col
                    with c_r3:
                        c_ref = st.text_input("Ref Allele", value=row_item.get('ref_allele', 'A'), max_chars=4, key=f"ex_ref_{r_id}").strip().upper()
                        row_item['ref_allele'] = c_ref
                    with c_r4:
                        c_var = st.text_input("Var Allele", value=row_item.get('var_allele', 'G'), max_chars=4, key=f"ex_var_{r_id}").strip().upper()
                        row_item['var_allele'] = c_var
                    with c_del:
                        st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
                        if st.button("🗑️", key=f"ex_del_{r_id}", help="Delete this row"):
                            rows_to_remove.append(r_idx)

                    # Live column detection for custom row
                    ex_observed_counts = {}
                    ex_detected_list = []
                    ex_norm = pd.Series(dtype=str)
                    if c_col in df_raw.columns:
                        raw_ex = df_raw[c_col].dropna().astype(str).str.strip().str.upper()
                        valid_ex = raw_ex[~raw_ex.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])]
                        ex_norm = valid_ex.map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                        ex_observed_counts = dict(ex_norm.value_counts())
                        ex_detected_list = list(ex_observed_counts.keys())

                    c_p_in, c_p_desc = st.columns([1.6, 1.4])
                    with c_p_in:
                        def_p_ex = row_item.get('permitted', f"{c_ref}{c_ref}, {c_ref}{c_var}, {c_var}{c_var}")
                        c_perm = st.text_input("Permitted Combinations (comma-separated):", value=def_p_ex, key=f"ex_perm_{r_id}")
                        row_item['permitted'] = c_perm

                        if ex_detected_list:
                            st.button(
                                f"⚡ Map Detected Calls ({', '.join(ex_detected_list)})",
                                key=f"btn_automap_ex_{r_id}",
                                on_click=cb_extra_automap,
                                args=(f"ex_perm_{r_id}", f"ex_ref_{r_id}", f"ex_var_{r_id}", ex_detected_list)
                            )

                    with c_p_desc:
                        st.markdown(f"""
                        <div style="font-size: 12px; color: #475569; padding-top: 28px;">
                            🧬 <b>Genotype Mapping:</b> Hom Ref: <code>{c_ref}{c_ref}</code> | Het: <code>{c_ref}{c_var}</code> | Hom Var: <code>{c_var}{c_var}</code>
                        </div>
                        """, unsafe_allow_html=True)

                    # Custom row live check table
                    if c_col in df_raw.columns and ex_observed_counts:
                        allowed_c = {g.strip().upper().replace('/', '').replace('|', '').replace(' ', '') for g in c_perm.split(',') if g.strip()}
                        exp_c = set(allowed_c)
                        for g in allowed_c:
                            if len(g) == 2:
                                exp_c.add(g[1] + g[0])
                        c_unexp = [k for k in ex_observed_counts.keys() if k not in exp_c]

                        ex_chk_rows = []
                        for call, cnt in sorted(ex_observed_counts.items(), key=lambda x: x[1], reverse=True):
                            pct = (cnt / len(ex_norm) * 100) if len(ex_norm) > 0 else 0
                            ex_chk_rows.append({
                                'Observed Genotype': call,
                                'Sample Count': cnt,
                                'Cohort %': f"{pct:.1f}%",
                                'Mapping Status': "✅ Permitted" if call in exp_c else "🚨 Flagged: Unmapped"
                            })
                        st.markdown(f"**Data Check for `{c_col}` ({len(ex_norm):,} valid samples):**")
                        st.dataframe(pd.DataFrame(ex_chk_rows), use_container_width=True, hide_index=True)

                        if c_unexp:
                            st.warning(f"🚨 Detected unmapped calls: `{c_unexp}`.")
                            st.button(
                                f"➕ Add `{', '.join(c_unexp)}` to Permitted List",
                                key=f"btn_add_unexp_ex_{r_id}",
                                on_click=cb_extra_add_unmapped,
                                args=(f"ex_perm_{r_id}", c_unexp, f"ex_ref_{r_id}", f"ex_var_{r_id}")
                            )
                        else:
                            st.success(f"✅ All {len(ex_norm):,} samples verified for `{c_name}`.")

                    if c_col and c_col != "(Select Column)" and c_col in df_raw.columns:
                        mapped_marker_dict[c_col] = {
                            'std_name': c_name,
                            'ref_allele': c_ref,
                            'var_allele': c_var,
                            'marker_label': c_name,
                            'permitted_genotypes': [g.strip().upper() for g in c_perm.split(',') if g.strip()]
                        }
                    else:
                        mapping_has_errors = True

                if rows_to_remove:
                    for idx_to_del in sorted(rows_to_remove, reverse=True):
                        extra_rows.pop(idx_to_del)
                    st.session_state['extra_mapping_rows'] = extra_rows
                    st.rerun()

            # Button to add extra mapping row
            st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)
            if st.button("➕ Add Custom Marker Mapping Row", key="btn_add_extra_mapping_row"):
                new_id = len(extra_rows) + 1
                extra_rows.append({
                    'id': f"extra_{len(extra_rows)}_{pd.Timestamp.now().strftime('%H%M%S')}",
                    'marker_name': f"Custom_Marker_{new_id}",
                    'col': "(Select Column)",
                    'ref_allele': 'A',
                    'var_allele': 'G',
                    'permitted': ''
                })
                st.session_state['extra_mapping_rows'] = extra_rows
                st.rerun()

        # -------------------------------------------------------------------------
        # MODE B: MULTI-ROW TABLE EDITOR (FAST BATCH ENTRY WITH LIVE VERIFICATION)
        # -------------------------------------------------------------------------
        else:
            st.markdown("#### 📋 Batch Multi-Row Mapping Table")
            st.caption("Enter and edit all mapping rows in an interactive table. Directly edit cells, select columns, configure alleles, and add rows with '➕ Add row'.")

            editor_data = []
            for m_idx, marker_label in enumerate(sel_map_markers):
                s_name, d_ref, d_var, d_perm, _, _ = get_marker_defaults(marker_label, sel_map_gene)
                is_stale = (st.session_state.get(f"_mapped_tracker_{m_idx}") != marker_label)
                cur_c = find_best_column_match(marker_label, s_name, sel_map_gene, col_options) if is_stale else st.session_state.get(f"w_map_col_{m_idx}", find_best_column_match(marker_label, s_name, sel_map_gene, col_options))
                cur_r = d_ref if is_stale else st.session_state.get(f"w_map_ref_{m_idx}", d_ref)
                cur_v = d_var if is_stale else st.session_state.get(f"w_map_var_{m_idx}", d_var)
                cur_p = d_perm if is_stale else st.session_state.get(f"w_map_perm_{m_idx}", d_perm)
                cur_n = s_name if is_stale else st.session_state.get(f"w_map_name_{m_idx}", s_name)
                editor_data.append({
                    'Marker Name': cur_n,
                    'Spreadsheet Column': cur_c if (cur_c and cur_c in col_options) else "(Select Column)",
                    'Reference Allele': cur_r,
                    'Variant Allele': cur_v,
                    'Permitted Genotypes': cur_p
                })

            for r_idx, row_item in enumerate(extra_rows):
                cur_c_ex = row_item.get('col', "(Select Column)")
                editor_data.append({
                    'Marker Name': row_item.get('marker_name', f"Custom_{r_idx+1}"),
                    'Spreadsheet Column': cur_c_ex if (cur_c_ex and cur_c_ex in col_options) else "(Select Column)",
                    'Reference Allele': row_item.get('ref_allele', 'A'),
                    'Variant Allele': row_item.get('var_allele', 'G'),
                    'Permitted Genotypes': row_item.get('permitted', '')
                })

            init_table_df = pd.DataFrame(editor_data)

            edited_table = st.data_editor(
                init_table_df,
                num_rows="dynamic",
                use_container_width=True,
                column_config={
                    'Marker Name': st.column_config.TextColumn("Marker Name", required=True),
                    'Spreadsheet Column': st.column_config.SelectboxColumn("Spreadsheet Column", options=col_options, required=True),
                    'Reference Allele': st.column_config.TextColumn("Reference Allele", required=True),
                    'Variant Allele': st.column_config.TextColumn("Variant Allele", required=True),
                    'Permitted Genotypes': st.column_config.TextColumn("Permitted Genotypes (comma-separated)", required=True)
                },
                key="batch_data_editor_step3"
            )

            c_app1, c_app2 = st.columns([0.3, 0.7])
            with c_app1:
                if st.button("💾 Apply Table Mappings", type="secondary", key="btn_apply_table_step3"):
                    new_extra = []
                    for t_idx, row_val in edited_table.iterrows():
                        m_name = str(row_val.get('Marker Name', '')).strip()
                        m_col = str(row_val.get('Spreadsheet Column', '')).strip()
                        m_ref = str(row_val.get('Reference Allele', 'A')).strip().upper()
                        m_var = str(row_val.get('Variant Allele', 'G')).strip().upper()
                        m_perm = str(row_val.get('Permitted Genotypes', '')).strip().upper()

                        if t_idx < len(sel_map_markers):
                            st.session_state[f"w_map_name_{t_idx}"] = m_name
                            st.session_state[f"w_map_col_{t_idx}"] = m_col
                            st.session_state[f"w_map_ref_{t_idx}"] = m_ref
                            st.session_state[f"w_map_var_{t_idx}"] = m_var
                            st.session_state[f"w_map_perm_{t_idx}"] = m_perm
                        else:
                            new_extra.append({
                                'id': f"extra_{t_idx}_{pd.Timestamp.now().strftime('%H%M%S')}",
                                'marker_name': m_name,
                                'col': m_col,
                                'ref_allele': m_ref,
                                'var_allele': m_var,
                                'permitted': m_perm
                            })
                    st.session_state['extra_mapping_rows'] = new_extra
                    st.success("All table mappings saved and synchronized!")
                    st.rerun()

            for _, row_val in edited_table.iterrows():
                m_name = str(row_val.get('Marker Name', '')).strip()
                m_col = str(row_val.get('Spreadsheet Column', '')).strip()
                m_ref = str(row_val.get('Reference Allele', 'A')).strip().upper()
                m_var = str(row_val.get('Variant Allele', 'G')).strip().upper()
                m_perm = str(row_val.get('Permitted Genotypes', '')).strip().upper()
                if m_col and m_col != "(Select Column)" and m_col in df_raw.columns:
                    mapped_marker_dict[m_col] = {
                        'std_name': m_name,
                        'ref_allele': m_ref,
                        'var_allele': 'A' if m_var == 'None' else m_var,
                        'marker_label': m_name,
                        'permitted_genotypes': [g.strip().upper() for g in m_perm.split(',') if g.strip()]
                    }
                else:
                    mapping_has_errors = True

            # Live Batch Data Verification Table
            st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
            st.markdown("##### 📊 Batch Live Data Verification Across All Mapped Columns")
            batch_check_rows = []
            for _, r_val in edited_table.iterrows():
                m_nm = str(r_val.get('Marker Name', '')).strip()
                m_cl = str(r_val.get('Spreadsheet Column', '')).strip()
                m_rf = str(r_val.get('Reference Allele', 'A')).strip().upper()
                m_vr = str(r_val.get('Variant Allele', 'G')).strip().upper()
                m_pm = str(r_val.get('Permitted Genotypes', '')).strip().upper()

                if m_cl and m_cl != "(Select Column)" and m_cl in df_raw.columns:
                    norm_c, c_cnts, _ = get_col_calls_cached(m_cl)
                    tot_c = len(norm_c)

                    pm_list = [g.strip().upper() for g in m_pm.split(',') if g.strip()]
                    pm_set = {g.replace('/', '').replace('|', '').replace(' ', '') for g in pm_list}
                    for g in list(pm_set):
                        if len(g) == 2:
                            pm_set.add(g[1] + g[0])

                    unexp = [k for k in c_cnts.keys() if k not in pm_set]
                    unexp_cnt = sum(c_cnts[k] for k in unexp)
                    dist_str = ", ".join([f"{k} ({v})" for k, v in sorted(c_cnts.items())])

                    status_str = f"✅ Verified ({tot_c:,} samples)" if not unexp else f"🚨 Flagged: {len(unexp)} Unmapped ({unexp_cnt} samples)"
                    batch_check_rows.append({
                        'Marker': m_nm,
                        'Column': m_cl,
                        'Ref / Var': f"{m_rf} / {m_vr}",
                        'Permitted Combinations': m_pm,
                        'Valid Samples': tot_c,
                        'Observed Call Breakdown': dist_str,
                        'QC Status': status_str
                    })
                else:
                    batch_check_rows.append({
                        'Marker': m_nm,
                        'Column': '(Not Selected)',
                        'Ref / Var': f"{m_rf} / {m_vr}",
                        'Permitted Combinations': m_pm,
                        'Valid Samples': 0,
                        'Observed Call Breakdown': 'N/A',
                        'QC Status': '⚪ Pending Column Selection'
                    })

            if batch_check_rows:
                st.dataframe(pd.DataFrame(batch_check_rows), use_container_width=True, hide_index=True)

        # Check for Duplicate Column Mapping Conflicts
        mapped_cols_list = [c for c in mapped_marker_dict.keys() if c and c != "(Select Column)"]
        dup_cols = {c for c in mapped_cols_list if mapped_cols_list.count(c) > 1}
        if dup_cols:
            st.error(f"⛔ **Duplicate Column Conflict:** Spreadsheet column(s) `{dup_cols}` are mapped to multiple markers! Each biomarker must be mapped to a distinct column.")
            mapping_has_errors = True

        total_configured_targets = len(sel_map_markers) + len(extra_rows)
        if len(mapped_marker_dict) < total_configured_targets:
            mapping_has_errors = True

        # Check for unmapped calls across all mapped columns to display summary advisory (using cached calls)
        total_unmapped_samples = 0
        for col_name, m_info in mapped_marker_dict.items():
            if col_name in df_raw.columns:
                p_list = m_info.get('permitted_genotypes', [])
                p_set = {g.replace('/', '').replace('|', '').replace(' ', '') for g in p_list if g.strip()}
                for g in list(p_set):
                    if len(g) == 2: p_set.add(g[1] + g[0])
                _, c_counts, _ = get_col_calls_cached(col_name)
                unm_cnt = sum(cnt for call, cnt in c_counts.items() if call not in p_set)
                total_unmapped_samples += unm_cnt

        if total_unmapped_samples > 0:
            st.warning(f"⚠️ **Advisory:** Detected {total_unmapped_samples:,} sample call(s) across your dataset that are not permitted by your genotype rules. Any unpermitted calls will be flagged during QC audit and treated as missing.")

        st.markdown("<br>", unsafe_allow_html=True)
        c_b1, c_b2, c_b3 = st.columns([0.25, 0.5, 0.25])
        with c_b1:
            if st.button("← Back to Upload", use_container_width=True, key="btn_map_back"):
                st.session_state['workflow_step'] = '2 · Upload'
                st.rerun()
        with c_b3:
            if mapping_has_errors:
                st.button("Continue to Results →", type="primary", disabled=True, use_container_width=True, key="btn_map_cont_dis", help="Resolve mapping errors above to proceed")
            elif st.button("Continue to Results →", type="primary", use_container_width=True, key="btn_map_cont"):
                # BUILD FINAL MAPPING WITH CONFIRMED METADATA & ALL MARKERS
                final_mapping = {}
                demo = st.session_state.get('demographic_mapping', {})
                if demo.get('Sample_ID') and demo['Sample_ID'] in all_cols:
                    final_mapping[demo['Sample_ID']] = 'Sample_ID'
                if demo.get('Gender') and demo['Gender'] in all_cols:
                    final_mapping[demo['Gender']] = 'Gender'
                if demo.get('State') and demo['State'] in all_cols:
                    final_mapping[demo['State']] = 'State'
                if demo.get('Region') and demo['Region'] in all_cols and demo['Region'] not in {demo.get('Sample_ID'), demo.get('Gender'), demo.get('State')}:
                    final_mapping[demo['Region']] = 'Region'
                custom_rules_for_qc = {}
                for m_col, m_info in mapped_marker_dict.items():
                    final_mapping[m_col] = m_info['std_name']
                    perm_list = m_info.get('permitted_genotypes', [])
                    r_a = m_info['ref_allele']
                    v_a = m_info['var_allele']
                    if perm_list:
                        clean_perms = set(g.replace('/', '').replace('|', '').replace(' ', '') for g in perm_list if g.strip())
                        augmented_p = set(clean_perms)
                        for g in clean_perms:
                            if len(g) == 2:
                                augmented_p.add(g[1] + g[0])
                        custom_rules_for_qc[m_info['std_name']] = {
                            'ref_allele': r_a,
                            'var_allele': v_a,
                            'valid_genotypes': tuple(sorted(augmented_p)),
                            'wildtype_genotype': f"{r_a}{r_a}",
                            'het_genotype': f"{r_a}{v_a}",
                            'hom_var_genotype': f"{v_a}{v_a}"
                        }
                    else:
                        dyn_perms = {r_a + r_a, r_a + v_a, v_a + r_a, v_a + v_a}
                        custom_rules_for_qc[m_info['std_name']] = {
                            'ref_allele': r_a,
                            'var_allele': v_a,
                            'valid_genotypes': tuple(sorted(dyn_perms)),
                            'wildtype_genotype': f"{r_a}{r_a}",
                            'het_genotype': f"{r_a}{v_a}",
                            'hom_var_genotype': f"{v_a}{v_a}"
                        }

                st.session_state['mapped_cols'] = final_mapping
                st.session_state['mapped_marker_details'] = mapped_marker_dict
                st.session_state['custom_snp_rules'] = custom_rules_for_qc

                # RESTRICT DATASET: RETAIN ONLY MAPPED COLUMNS & ENFORCE PERMITTED GENOTYPES
                if df_raw is not None and not df_raw.empty:
                    cols_to_keep = [c for c in df_raw.columns if c in final_mapping]
                    if cols_to_keep:
                        filtered_df = df_raw[cols_to_keep].copy()

                        for m_col, m_info in mapped_marker_dict.items():
                            if m_col in filtered_df.columns:
                                user_permitted = m_info.get('permitted_genotypes', [])
                                if user_permitted:
                                    valid_set = {str(g).strip().upper().replace('/', '').replace('|', '').replace(' ', '') for g in user_permitted if str(g).strip()}
                                    augmented = set(valid_set)
                                    for g in valid_set:
                                        if len(g) == 2:
                                            augmented.add(g[1] + g[0])
                                    def clean_user_calls(v):
                                        if pd.isna(v): return np.nan
                                        clean = str(v).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
                                        if clean in valid_set: return clean
                                        if clean in augmented:
                                            rev = clean[1] + clean[0] if len(clean) == 2 else clean
                                            if rev in valid_set: return rev
                                            return clean
                                        return np.nan
                                    filtered_df[m_col] = filtered_df[m_col].apply(clean_user_calls)
                                else:
                                    r_a = m_info['ref_allele']
                                    v_a = m_info['var_allele']
                                    valid_dyn = {r_a + r_a, r_a + v_a, v_a + r_a, v_a + v_a}
                                    def clean_dyn(v):
                                        if pd.isna(v): return np.nan
                                        clean = str(v).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
                                        if clean in valid_dyn: return clean
                                        return np.nan
                                    filtered_df[m_col] = filtered_df[m_col].apply(clean_dyn)

                        st.session_state['uploaded_df'] = filtered_df
                        cur_ds_name = st.session_state.get('active_ds_name', 'Active Dataset')
                        st.session_state.setdefault('datasets_dict', {})[cur_ds_name] = filtered_df
                        save_active_dataset(filtered_df)

                st.session_state['workflow_step'] = '4 · Results'
                save_active_session()
                st.rerun()

