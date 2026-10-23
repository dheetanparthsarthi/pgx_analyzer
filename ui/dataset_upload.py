"""Step 2 Dataset Upload component (Validation, multi-sheet handling, and dataset management)."""

import io
from typing import Any, Callable, Dict
import pandas as pd
import streamlit as st
from data_loader import inspect_uploaded_file, smart_merge_datasets
from marker_catalog import (
    detect_genes_in_columns,
    find_best_column_match,
    get_default_markers_for_gene,
    load_marker_catalog,
)
from ui.session_manager import (
    clear_deletion_flag,
    has_active_data,
    reset_all_data,
    save_active_dataset,
    set_active_target_gene,
)
from ui.views.results_dashboard import (
    get_df_download_csv,
    get_df_download_excel,
)


def render_dataset_upload(df_raw: pd.DataFrame) -> None:
    """Render Step 2: Upload Genotype Dataset."""
    st.markdown("## 2️⃣ Upload Genotype Dataset")
    st.caption(
        "Upload your spreadsheet (.xlsx, .xls, .csv). Compulsory fields (Sample ID, Gender, State) will be inspected immediately."
    )

    catalog_df = load_marker_catalog()

    st.markdown(
        """
    <div class="upload-container-box">
        <div class="upload-icon-style">⇧</div>
        <div class="upload-title-style">Click to choose a file, or drag it here</div>
        <div class="upload-sub-style">.xlsx, .xls, or .csv — first row must be headers</div>
    </div>
    """,
        unsafe_allow_html=True,
    )

    # Display active workspace status if data is loaded
    if has_active_data():
        c_status, c_act1, c_down_top, c_act2 = st.columns([1.1, 0.7, 0.7, 0.5])
        with c_status:
            st.info(
                f"📂 **Active Workspace Data:** `{len(df_raw):,} samples` loaded."
            )
        with c_act1:
            if st.button(
                "📊 View Results →",
                use_container_width=True,
                key="btn_step1_quick_res",
            ):
                st.session_state["workflow_step"] = "4 · Results"
                st.rerun()
        with c_down_top:
            st.download_button(
                label="⬇️ Download (.xlsx)",
                data=get_df_download_excel(df_raw),
                file_name=f"complete_dataset_{len(df_raw)}_samples.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
                key="btn_step2_top_download",
            )
        with c_act2:
            if st.button(
                "🗑️ Delete", use_container_width=True, key="btn_step1_quick_del"
            ):
                reset_all_data()
                st.rerun()

    file_upload = st.file_uploader(
        "Select dataset file",
        type=["xlsx", "xls", "csv"],
        key="w_step1_uploader",
        label_visibility="collapsed",
    )
    if file_upload is not None:
        uploaded_file_id = f"{file_upload.name}_{file_upload.size}"
        if st.session_state.get("last_uploaded_file_id") != uploaded_file_id:
            try:
                file_bytes = io.BytesIO(file_upload.getvalue())
                if file_upload.name.lower().endswith(".csv"):
                    df_load = pd.read_csv(file_bytes)
                else:
                    excel_file = pd.ExcelFile(file_bytes)
                    sheet_names = excel_file.sheet_names
                    if len(sheet_names) > 1:
                        best_sheet = sheet_names[0]
                        max_density = -1
                        for s in sheet_names:
                            s_df_prev = pd.read_excel(
                                excel_file, sheet_name=s, nrows=50
                            )
                            gt_cols = [
                                c
                                for c in s_df_prev.columns
                                if any(
                                    k in str(c).lower()
                                    for k in [
                                        "*2",
                                        "*3",
                                        "*17",
                                        "cyp",
                                        "rs",
                                        "snp",
                                        "marker",
                                        "gene",
                                    ]
                                )
                            ]
                            score = (
                                int(s_df_prev[gt_cols].notna().sum().sum())
                                if gt_cols
                                else int(s_df_prev.notna().sum().sum())
                            )
                            if score > max_density:
                                max_density = score
                                best_sheet = s
                        df_load = pd.read_excel(
                            excel_file, sheet_name=best_sheet
                        )
                        st.session_state["uploaded_sheet_names"] = sheet_names
                        st.session_state["active_sheet_name"] = best_sheet
                        st.session_state["uploaded_excel_bytes"] = (
                            file_upload.getvalue()
                        )
                    else:
                        df_load = pd.read_excel(file_bytes)

                chk_res = inspect_uploaded_file(df_load)

                if not chk_res["has_sample_id"]:
                    st.session_state["upload_error"] = (
                        "🚫 **Upload Blocked: Mandatory 'Sample ID' column is missing!**\n\n"
                        "A unique Sample ID (e.g. `Sample ID`, `sample_id`, `Patient ID`, `ID`) is strictly compulsory to track samples and calculate frequencies.\n\n"
                        "**You cannot proceed without a Sample ID column.** Please add a Sample ID column to your spreadsheet and re-upload."
                    )
                    st.session_state["pending_new_df"] = None
                else:
                    st.session_state["upload_error"] = None
                    st.session_state["upload_inspection"] = chk_res
                    st.session_state["demographic_mapping"] = {
                        "Sample_ID": chk_res["sample_id_col"],
                        "Gender": chk_res["gender_col"],
                        "State": chk_res["state_col"],
                    }

                    if has_active_data():
                        st.session_state["pending_new_df"] = df_load
                        st.session_state["pending_file_id"] = uploaded_file_id
                        st.session_state["pending_file_name"] = file_upload.name
                        st.warning(
                            f"⚠️ Dataset **{file_upload.name}** ({len(df_load):,} samples) uploaded! Please choose below: **Merge with Active Dataset**, **Keep as Separate Dataset**, or **Replace Active Dataset**."
                        )
                    else:
                        clear_deletion_flag()
                        st.session_state["last_uploaded_file_id"] = (
                            uploaded_file_id
                        )
                        st.session_state["uploaded_df"] = df_load
                        active_label = f"Dataset 1 ({file_upload.name})"
                        st.session_state["datasets_dict"] = {
                            active_label: df_load
                        }
                        st.session_state["active_ds_name"] = active_label
                        st.session_state["mapped_cols"] = {}
                        st.session_state["pending_new_df"] = None
                        st.session_state["pending_file_id"] = None
                        st.session_state["pending_file_name"] = None
                        save_active_dataset(df_load)

                        detected_file_genes = detect_genes_in_columns(
                            df_load.columns
                        )
                        curr_g = st.session_state.get("selected_gene")
                        if not curr_g and detected_file_genes:
                            primary_g = detected_file_genes[0]
                            all_g_markers = get_default_markers_for_gene(
                                primary_g, catalog_df
                            )
                            matched_f_markers = [
                                m
                                for m in all_g_markers
                                if find_best_column_match(
                                    m,
                                    m.split(" ")[0],
                                    primary_g,
                                    list(df_load.columns),
                                )
                            ]
                            chosen_markers = (
                                matched_f_markers
                                if matched_f_markers
                                else all_g_markers
                            )
                            set_active_target_gene(
                                primary_g,
                                markers=chosen_markers,
                                catalog_df=catalog_df,
                            )
                        elif curr_g:
                            all_g_markers = get_default_markers_for_gene(
                                curr_g, catalog_df
                            )
                            matched_f_markers = [
                                m
                                for m in all_g_markers
                                if find_best_column_match(
                                    m,
                                    m.split(" ")[0],
                                    curr_g,
                                    list(df_load.columns),
                                )
                            ]
                            if matched_f_markers:
                                st.session_state["selected_markers"] = (
                                    matched_f_markers
                                )

                        st.session_state["workflow_step"] = "3 · Map columns"
                        st.success(
                            f"✅ Loaded **{file_upload.name}** ({len(df_load):,} samples, {len(df_load.columns)} columns). Advancing to Map Columns..."
                        )
                        st.rerun()
            except Exception as e:
                st.error(f"Error loading file: {e}")

    # Multi-Sheet Switcher if uploaded workbook has multiple sheets
    if (
        len(st.session_state.get("uploaded_sheet_names", [])) > 1
        and st.session_state.get("uploaded_excel_bytes")
    ):
        all_sheets = st.session_state["uploaded_sheet_names"]
        cur_s = st.session_state.get("active_sheet_name", all_sheets[0])
        s_idx = all_sheets.index(cur_s) if cur_s in all_sheets else 0
        c_sh1, c_sh2 = st.columns([1.2, 1.8])
        with c_sh1:
            sel_s = st.selectbox(
                "📑 Multi-Sheet Detected — Switch Active Sheet:",
                all_sheets,
                index=s_idx,
                key="sel_active_excel_sheet",
            )
            if sel_s != cur_s:
                st.session_state["active_sheet_name"] = sel_s
                new_df = pd.read_excel(
                    io.BytesIO(st.session_state["uploaded_excel_bytes"]),
                    sheet_name=sel_s,
                )
                chk_res = inspect_uploaded_file(new_df)
                st.session_state["uploaded_df"] = new_df
                st.session_state["upload_inspection"] = chk_res
                st.session_state["demographic_mapping"] = {
                    "Sample_ID": chk_res["sample_id_col"],
                    "Gender": chk_res["gender_col"],
                    "State": chk_res["state_col"],
                }
                detected_file_genes = detect_genes_in_columns(new_df.columns)
                curr_g = st.session_state.get("selected_gene")
                if not curr_g and detected_file_genes:
                    primary_g = detected_file_genes[0]
                    all_g_markers = get_default_markers_for_gene(
                        primary_g, catalog_df
                    )
                    matched_f_markers = [
                        m
                        for m in all_g_markers
                        if find_best_column_match(
                            m,
                            m.split(" ")[0],
                            primary_g,
                            list(new_df.columns),
                        )
                    ]
                    chosen_markers = (
                        matched_f_markers
                        if matched_f_markers
                        else all_g_markers
                    )
                    set_active_target_gene(
                        primary_g,
                        markers=chosen_markers,
                        catalog_df=catalog_df,
                    )
                elif curr_g:
                    all_g_markers = get_default_markers_for_gene(
                        curr_g, catalog_df
                    )
                    matched_f_markers = [
                        m
                        for m in all_g_markers
                        if find_best_column_match(
                            m,
                            m.split(" ")[0],
                            curr_g,
                            list(new_df.columns),
                        )
                    ]
                    if matched_f_markers:
                        st.session_state["selected_markers"] = matched_f_markers
                save_active_dataset(new_df)
                st.rerun()

    # Show Error Banner if Sample ID Missing
    if st.session_state.get("upload_error"):
        st.error(st.session_state["upload_error"])

    st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)
    st.markdown("#### 📋 Sample and Population Group Mapping")
    st.caption(
        "Sample ID is required. Gender, state, and an optional geographic/population-group column enable their respective stratified views."
    )

    cur_insp = st.session_state.get("upload_inspection")
    if not cur_insp and df_raw is not None and not df_raw.empty:
        cur_insp = inspect_uploaded_file(df_raw)
        st.session_state["upload_inspection"] = cur_insp

    raw_cols = (
        list(df_raw.columns) if df_raw is not None and not df_raw.empty else []
    )
    demo_map = st.session_state.get("demographic_mapping", {})
    region_tokens = (
        "region",
        "geographic group",
        "population group",
        "ancestry",
        "ethnicity",
        "continent",
        "country",
    )
    reserved_demo_cols = {
        demo_map.get("Sample_ID")
        or (cur_insp.get("sample_id_col") if cur_insp else None),
        demo_map.get("Gender")
        or (cur_insp.get("gender_col") if cur_insp else None),
        demo_map.get("State")
        or (cur_insp.get("state_col") if cur_insp else None),
    }
    region_guess = next(
        (
            col
            for col in raw_cols
            if col not in reserved_demo_cols
            and any(token in str(col).lower() for token in region_tokens)
        ),
        None,
    )
    if not demo_map and cur_insp:
        demo_map = {
            "Sample_ID": cur_insp.get("sample_id_col"),
            "Gender": cur_insp.get("gender_col"),
            "State": cur_insp.get("state_col"),
            "Region": region_guess,
        }
        st.session_state["demographic_mapping"] = demo_map

    sid_mapped = (
        demo_map.get("Sample_ID")
        if "Sample_ID" in demo_map
        else (cur_insp["sample_id_col"] if cur_insp else None)
    )
    gen_mapped = (
        demo_map.get("Gender")
        if "Gender" in demo_map
        else (cur_insp["gender_col"] if cur_insp else None)
    )
    state_mapped = (
        demo_map.get("State")
        if "State" in demo_map
        else (cur_insp["state_col"] if cur_insp else None)
    )
    region_mapped = (
        demo_map.get("Region") if "Region" in demo_map else region_guess
    )
    if region_mapped and "Region" not in demo_map:
        demo_map["Region"] = region_mapped
        st.session_state["demographic_mapping"] = demo_map

    c_aud1, c_aud2, c_aud3, c_aud4 = st.columns(4)
    with c_aud1:
        if sid_mapped and sid_mapped in raw_cols:
            st.success(f"✅ **Sample ID (Required):** `{sid_mapped}`")
        else:
            st.error("❌ **Sample ID (Required):** Missing")
    with c_aud2:
        if gen_mapped and gen_mapped in raw_cols:
            st.success(f"✅ **Gender (Required):** `{gen_mapped}`")
        else:
            st.warning("⚠️ **Gender (Required):** Not detected in columns")
    with c_aud3:
        if state_mapped and state_mapped in raw_cols:
            st.success(f"✅ **State:** `{state_mapped}`")
        else:
            st.info("State column not mapped")
    with c_aud4:
        if region_mapped and region_mapped in raw_cols:
            st.success(f"✅ **Geographic Group:** `{region_mapped}`")
        else:
            st.info("Optional geographic group not mapped")

    if raw_cols:
        with st.expander(
            "🔧 Confirm / Adjust Required Metadata Columns", expanded=False
        ):
            c_adj1, c_adj2, c_adj3, c_adj4 = st.columns(4)
            sid_opt = raw_cols
            gen_opt = ["(None)"] + raw_cols
            state_opt = ["(None)"] + raw_cols

            idx_sid = sid_opt.index(sid_mapped) if sid_mapped in sid_opt else 0
            idx_gen = gen_opt.index(gen_mapped) if gen_mapped in gen_opt else 0
            idx_state = (
                state_opt.index(state_mapped)
                if state_mapped in state_opt
                else 0
            )

            chosen_sid = c_adj1.selectbox(
                "Sample ID Column *", sid_opt, index=idx_sid, key="u_map_sid"
            )
            chosen_gen = c_adj2.selectbox(
                "Gender Column *", gen_opt, index=idx_gen, key="u_map_gen"
            )
            chosen_state = c_adj3.selectbox(
                "State / Locality Column",
                state_opt,
                index=idx_state,
                key="u_map_state",
            )
            region_opt = ["(None)"] + [
                col
                for col in raw_cols
                if col not in {chosen_sid, chosen_gen, chosen_state}
            ]
            idx_region = (
                region_opt.index(region_mapped)
                if region_mapped in region_opt
                else 0
            )
            chosen_region = c_adj4.selectbox(
                "Geographic Group Column",
                region_opt,
                index=idx_region,
                key="u_map_region",
            )

            st.session_state["demographic_mapping"] = {
                "Sample_ID": chosen_sid,
                "Gender": chosen_gen if chosen_gen != "(None)" else None,
                "State": chosen_state if chosen_state != "(None)" else None,
                "Region": chosen_region if chosen_region != "(None)" else None,
            }

    if raw_cols:
        detected_file_genes = detect_genes_in_columns(raw_cols)
        target_gene_active = st.session_state.get("selected_gene")
        if not target_gene_active and detected_file_genes:
            target_gene_active = detected_file_genes[0]
            st.session_state["selected_gene"] = target_gene_active

        target_markers_active = st.session_state.get("selected_markers", [])

        if target_gene_active:
            all_g_markers = get_default_markers_for_gene(
                target_gene_active, catalog_df
            )
            custom_m = st.session_state.get("custom_markers", {}).get(
                target_gene_active, []
            )
            valid_m = all_g_markers + custom_m
            filtered_markers = [
                m
                for m in target_markers_active
                if any(
                    m == vm or m.split(" ")[0] == vm.split(" ")[0]
                    for vm in valid_m
                )
            ]
            if not filtered_markers:
                matched_f = [
                    m
                    for m in all_g_markers
                    if find_best_column_match(
                        m, m.split(" ")[0], target_gene_active, raw_cols
                    )
                ]
                target_markers_active = (
                    matched_f if matched_f else all_g_markers
                )
                st.session_state["selected_markers"] = target_markers_active
            else:
                target_markers_active = filtered_markers
                st.session_state["selected_markers"] = target_markers_active

        st.markdown(
            "<div style='margin-top: 14px;'></div>", unsafe_allow_html=True
        )
        st.markdown("#### 🎯 Target Biomarker & Header Audit")
        if target_gene_active:
            st.caption(
                f"Verifying uploaded column headers against your selection: Target Gene **{target_gene_active}** ({len(target_markers_active)} biomarker(s))."
            )
        else:
            st.caption(
                "Verifying uploaded column headers against pharmacogenetic targets."
            )

        gene_mismatch = bool(
            detected_file_genes
            and target_gene_active
            and target_gene_active.upper()
            not in [g.upper() for g in detected_file_genes]
        )

        if gene_mismatch:
            found_genes_str = ", ".join(detected_file_genes)
            st.markdown(
                f"""
            <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 14px; margin-bottom: 12px;">
                <div style="display: flex; align-items: center; gap: 8px;">
                    <span style="font-size: 20px;">🚨</span>
                    <div style="font-weight: 800; color: #991B1B; font-size: 15px;">TARGET GENE & COLUMN HEADER MISMATCH DETECTED</div>
                </div>
                <div style="color: #7F1D1D; font-size: 13px; margin-top: 6px; line-height: 1.5;">
                    In Step 1 you selected target gene: <b>{target_gene_active}</b>.<br>
                    However, the column headers in your uploaded spreadsheet belong to: <b>{found_genes_str}</b>.<br>
                    If you intended to analyze <b>{found_genes_str}</b>, click the button below to switch your target gene to match this dataset.
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )

            c_sw_btn, _ = st.columns([0.4, 0.6])
            with c_sw_btn:
                primary_detected = detected_file_genes[0]
                if st.button(
                    f"👉 Switch Target Gene to {primary_detected}",
                    type="primary",
                    use_container_width=True,
                    key="btn_switch_gene_step2",
                ):
                    set_active_target_gene(
                        primary_detected, catalog_df=catalog_df
                    )
                    st.rerun()

        audit_cols_display = []
        for m in target_markers_active:
            m_base = m.split(" ")[0] if " " in m else m
            matched_col = find_best_column_match(
                m, m_base, target_gene_active or "", raw_cols
            )

            audit_cols_display.append(
                {
                    "Biomarker (Target)": m,
                    "Header in File": matched_col
                    if matched_col
                    else "❌ Not found in file",
                    "Match Status": "✅ Header Match"
                    if matched_col
                    else "⚠️ Different / Missing Header",
                }
            )

        c_m_aud1, c_m_aud2 = st.columns([1.6, 1.4])
        with c_m_aud1:
            st.dataframe(
                pd.DataFrame(audit_cols_display),
                use_container_width=True,
                hide_index=True,
            )
        with c_m_aud2:
            n_matched = sum(
                1 for a in audit_cols_display if "✅" in a["Match Status"]
            )
            if n_matched == len(target_markers_active) and not gene_mismatch:
                st.success(
                    f"✅ All {n_matched} selected biomarker header(s) matched in uploaded file."
                )
            elif gene_mismatch:
                st.error(
                    f"🚨 **Target Gene Mismatch:** Uploaded file headers contain `{found_genes_str}`, but target gene is `{target_gene_active}`."
                )
            else:
                st.warning(
                    f"⚠️ **Header Notice:** {n_matched} of {len(target_markers_active)} biomarker header(s) match file columns. You can manually map different headers in Step 3."
                )

    if st.session_state.get("pending_new_df") is not None:
        pending_df = st.session_state["pending_new_df"]
        pending_fname = st.session_state.get(
            "pending_file_name", "Uploaded File"
        )
        st.markdown(
            f"""
        <div style="background: #EFF6FF; border: 2px solid #3B82F6; border-radius: 8px; padding: 14px; margin-top: 12px; margin-bottom: 12px;">
            <div style="display: flex; align-items: center; gap: 8px;">
                <span style="font-size: 20px;">📂</span>
                <div style="font-weight: 800; color: #1E40AF; font-size: 15px;">Multiple Datasets Detected in Workspace</div>
            </div>
            <div style="color: #1E3A8A; font-size: 13px; margin-top: 6px; line-height: 1.5;">
                You have an active dataset with <b>{len(st.session_state['uploaded_df']):,} samples</b>, and you just uploaded <b>{pending_fname}</b> ({len(pending_df):,} samples).<br>
                Choose how you would like to handle this data:
            </div>
        </div>
        """,
            unsafe_allow_html=True,
        )

        c_mrg1, c_mrg2, c_mrg3 = st.columns(3)
        with c_mrg1:
            if st.button(
                "🔄 Merge with Active Dataset",
                type="primary",
                use_container_width=True,
                key="btn_merge_dataset",
                help="Allowed when both datasets share the same target gene",
            ):
                base_df = st.session_state["uploaded_df"]
                merged_df, align_log = smart_merge_datasets(base_df, pending_df)
                is_blocked = any(
                    "⛔ Cannot merge" in str(msg) for msg in align_log
                )
                if is_blocked:
                    st.session_state["merge_alignment_log"] = align_log
                    st.error(align_log[0])
                else:
                    clear_deletion_flag()
                    st.session_state["last_uploaded_file_id"] = st.session_state[
                        "pending_file_id"
                    ]
                    st.session_state["uploaded_df"] = merged_df
                    merged_name = (
                        f"Master Merged Dataset ({len(merged_df):,} samples)"
                    )
                    st.session_state.setdefault("datasets_dict", {})[
                        merged_name
                    ] = merged_df
                    st.session_state["active_ds_name"] = merged_name
                    st.session_state["merge_alignment_log"] = align_log
                    st.session_state["pending_new_df"] = None
                    st.session_state["pending_file_id"] = None
                    st.session_state["pending_file_name"] = None
                    save_active_dataset(merged_df)
                    st.session_state["workflow_step"] = "3 · Map columns"
                    st.success(
                        f"✅ Successfully merged {len(pending_df):,} samples! Total dataset now contains {len(merged_df):,} samples."
                    )
                    st.rerun()
        with c_mrg2:
            if st.button(
                "📁 Keep as Separate Dataset",
                type="primary",
                use_container_width=True,
                key="btn_separate_dataset",
                help="Preserves both datasets so you can switch between them in Results via dropdown",
            ):
                clear_deletion_flag()
                d_dict = st.session_state.setdefault("datasets_dict", {})
                cur_name = st.session_state.get(
                    "active_ds_name",
                    f"Dataset 1 ({len(st.session_state['uploaded_df']):,} samples)",
                )
                d_dict[cur_name] = st.session_state["uploaded_df"]
                new_ds_lbl = f"Dataset {len(d_dict) + 1} ({pending_fname} - {len(pending_df):,} samples)"
                d_dict[new_ds_lbl] = pending_df
                st.session_state["last_uploaded_file_id"] = st.session_state[
                    "pending_file_id"
                ]
                st.session_state["uploaded_df"] = pending_df
                st.session_state["active_ds_name"] = new_ds_lbl
                st.session_state["mapped_cols"] = {}
                st.session_state["pending_new_df"] = None
                st.session_state["pending_file_id"] = None
                st.session_state["pending_file_name"] = None
                save_active_dataset(pending_df)
                st.session_state["workflow_step"] = "3 · Map columns"
                st.success(
                    f"📁 Added '{pending_fname}' as a separate dataset for independent analysis!"
                )
                st.rerun()
        with c_mrg3:
            if st.button(
                "🆕 Replace Active Dataset",
                type="secondary",
                use_container_width=True,
                key="btn_replace_dataset",
                help="Replaces the current active workspace",
            ):
                clear_deletion_flag()
                st.session_state["last_uploaded_file_id"] = st.session_state[
                    "pending_file_id"
                ]
                st.session_state["uploaded_df"] = pending_df
                single_lbl = f"Active Dataset ({pending_fname})"
                st.session_state["datasets_dict"] = {single_lbl: pending_df}
                st.session_state["active_ds_name"] = single_lbl
                st.session_state["mapped_cols"] = {}
                st.session_state["pending_new_df"] = None
                st.session_state["pending_file_id"] = None
                st.session_state["pending_file_name"] = None
                save_active_dataset(pending_df)
                st.session_state["workflow_step"] = "3 · Map columns"
                st.success(
                    f"🆕 Replaced active workspace with '{pending_fname}' ({len(pending_df):,} samples)."
                )
                st.rerun()

    # Show Gene Alignment Log if merge occurred
    if st.session_state.get("merge_alignment_log"):
        with st.expander(
            "🧬 Gene & Column Alignment Verification Log", expanded=True
        ):
            for log_entry in st.session_state["merge_alignment_log"]:
                st.markdown(f"• {log_entry}")

    # Show active dataset preview and download options if loaded
    if has_active_data():
        st.markdown(
            "<div style='margin-top: 14px;'></div>", unsafe_allow_html=True
        )
        c_prev_head, c_dl_xl, c_dl_csv = st.columns([1.2, 0.9, 0.9])
        with c_prev_head:
            st.markdown(
                f"##### 📋 Active Dataset Preview ({len(df_raw):,} samples, {len(df_raw.columns)} columns):"
            )
        with c_dl_xl:
            st.download_button(
                label="⬇️ Download Complete Data (.xlsx)",
                data=get_df_download_excel(df_raw),
                file_name=f"complete_pgx_dataset_{len(df_raw)}_samples.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
                key="btn_step2_prev_download_xlsx",
            )
        with c_dl_csv:
            st.download_button(
                label="⬇️ Download Complete Data (.csv)",
                data=get_df_download_csv(df_raw),
                file_name=f"complete_pgx_dataset_{len(df_raw)}_samples.csv",
                mime="text/csv",
                use_container_width=True,
                key="btn_step2_prev_download_csv",
            )
        st.dataframe(df_raw.head(5), use_container_width=True)

    st.markdown(
        """
    <div class="privacy-info-banner">
        No data leaves your browser. Parsing and all statistics run locally.
    </div>
    """,
        unsafe_allow_html=True,
    )

    c_back0, c_space, c_next = st.columns([0.25, 0.45, 0.3])
    with c_back0:
        if st.button(
            "← Back to Select Target",
            use_container_width=True,
            key="btn_step1_back",
        ):
            st.session_state["workflow_step"] = "1 · Select Target"
            st.rerun()
    with c_next:
        can_proceed = has_active_data() and not st.session_state.get(
            "upload_error"
        )
        if st.button(
            "Continue to Map Columns →",
            type="primary",
            use_container_width=True,
            key="btn_step1_cont",
            disabled=not can_proceed,
        ):
            st.session_state["workflow_step"] = "3 · Map columns"
            st.rerun()
