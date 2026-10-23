"""Pharmacogenomic Calling, Diplotype and Phenotype analytical view component."""

from collections import Counter
import re
from typing import Any, Callable, Dict, List, Tuple
import pandas as pd
import plotly.express as px
import streamlit as st
from config import CYP2C9_SNP_CONFIG, SNP_CONFIG
from marker_catalog import load_marker_catalog
from ui.session_manager import save_active_session


def render_cyp2c19_combination_studio(
    processed_df: pd.DataFrame,
    full_results: dict,
    view_mode: str = "All Views",
) -> None:
    """Renders genotype combination and geographic/state summaries for the markers

    selected from the active dataset.
    """
    if processed_df is None or processed_df.empty:
        st.info("No pharmacogenomic data available for combination analysis.")
        return

    # -------------------------------------------------------------------------
    # STUDIO CONTAINER & HEADER (CLEAN, COMPACT TYPOGRAPHY)
    # -------------------------------------------------------------------------
    target_gene = (
        st.session_state.get("selected_gene")
        or full_results.get("target_gene")
        or full_results.get("active_gene")
        or "Custom"
    )
    st.markdown(
        f"""
    <div style="border-top: 1px solid #E2E8F0; margin-top: 1.5rem; padding-top: 1rem; margin-bottom: 0.75rem;">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div>
                <span style="font-size: 16px; font-weight: 700; color: #0F172A;">🧩 {target_gene} Multi-Variant & Geographic Studio</span>
                <span style="margin-left: 8px; font-size: 11px; font-weight: 600; background: #EEF2F6; color: #334155; padding: 2px 7px; border-radius: 4px; border: 1px solid #CBD5E1;">Dataset Rules</span>
            </div>
        </div>
        <p style="margin: 3px 0 0 0; font-size: 12px; color: #64748B;">
            Genotype combination and sub-population distribution analysis for the selected markers.
        </p>
    </div>
    """,
        unsafe_allow_html=True,
    )

    # -------------------------------------------------------------------------
    # 1. VARIANT COMBINATION FOCUS SELECTOR (DYNAMIC ACROSS GENES)
    # -------------------------------------------------------------------------
    mapped_snps_in_df = [
        marker
        for marker in full_results.get("overall", {}).get("snps", {})
        if marker in processed_df.columns
    ]
    if not mapped_snps_in_df:
        mapped_snps_in_df = [
            marker
            for marker in full_results.get("markers", [])
            if marker in processed_df.columns
        ]
    if not mapped_snps_in_df:
        st.info(
            "ℹ️ Multi-Locus Combination Studio requires at least one mapped marker column with valid genotype calls."
        )
        return

    active_snps = st.multiselect(
        "Markers to include in the combination analysis:",
        options=mapped_snps_in_df,
        default=mapped_snps_in_df,
        key=f"studio_active_markers_{re.sub(r'[^A-Za-z0-9]+', '_', str(target_gene))}",
    )
    if not active_snps:
        st.info(
            "Select at least one mapped marker to calculate genotype combinations."
        )
        return
    combo_title = " + ".join(active_snps)
    std_snps_label = (
        f"{len(active_snps)} selected marker{'s' if len(active_snps) != 1 else ''}"
    )

    active_snps = [c for c in active_snps if c in processed_df.columns]
    if not active_snps:
        st.info(
            "ℹ️ No valid marker columns found in the active dataset to build a combination studio matrix."
        )
        return

    state_col = (
        "State_Clean"
        if "State_Clean" in processed_df.columns
        else ("State" if "State" in processed_df.columns else None)
    )

    # -------------------------------------------------------------------------
    # 2. GEOGRAPHIC LEVEL, SPECIFIC DRILL-DOWN, AND NOMENCLATURE FORMAT
    # -------------------------------------------------------------------------
    c_geo, c_specific = st.columns([1.1, 1.4])

    with c_geo:
        default_scope_idx = 1 if "Geographic" in view_mode else 0
        geo_scope = st.selectbox(
            "🌍 Geographic Level:",
            [
                "🌐 Overall Cohort",
                "🗺️ Geographic Groups",
                "📍 State / Locality Groups",
            ],
            index=default_scope_idx,
            key="studio_geo_scope",
        )

    selected_region = "All Regions (Comparative Overview)"
    selected_state = "All States (Comparative Overview)"

    with c_specific:
        if "Geographic Groups" in geo_scope:
            region_list = ["All Regions (Comparative Overview)"]
            if "Region" in processed_df.columns:
                unique_regs = sorted(
                    [
                        str(r)
                        for r in processed_df["Region"].dropna().unique()
                        if str(r).strip() != ""
                    ]
                )
                region_list.extend(unique_regs)
            selected_region = st.selectbox(
                "🗺️ Specific Geographical Region:",
                region_list,
                key="studio_sel_region",
            )
        elif "State / Locality Groups" in geo_scope:
            state_list = ["All States (Comparative Overview)"]
            if state_col and state_col in processed_df.columns:
                unique_states = sorted(
                    [
                        str(s)
                        for s in processed_df[state_col].dropna().unique()
                        if str(s).strip() != ""
                    ]
                )
                state_list.extend(unique_states)
            selected_state = st.selectbox(
                "📍 Specific State Profile:",
                state_list,
                key="studio_sel_state",
            )
        else:
            st.selectbox(
                "📍 Stratification Status:",
                ["All uploaded samples included"],
                disabled=True,
                key="studio_disabled_geo",
            )

    # -------------------------------------------------------------------------
    # 3. DYNAMIC DATA FILTERING BASED ON GEOGRAPHY / STATE
    # -------------------------------------------------------------------------
    active_df = processed_df.copy()
    scope_label = "Overall Cohort"

    if "Geographic Groups" in geo_scope:
        if selected_region != "All Regions (Comparative Overview)":
            active_df = active_df[
                active_df["Region"].astype(str) == str(selected_region)
            ]
            scope_label = f"Region: {selected_region}"
        else:
            scope_label = "All Regions (Comparative)"
    elif "State / Locality Groups" in geo_scope:
        if selected_state != "All States (Comparative Overview)" and state_col:
            active_df = active_df[
                active_df[state_col].astype(str) == str(selected_state)
            ]
            scope_label = f"State: {selected_state}"
        else:
            scope_label = "All States (Comparative)"

    tot_cohort = len(active_df)
    if tot_cohort == 0:
        st.warning(
            f"⚠️ No patient records found for the selected geographic filter: **{scope_label}**."
        )
        return

    # Compute combination breakdown across the active filtered dataset
    combo_grp = (
        active_df.groupby(active_snps, dropna=False)
        .size()
        .reset_index(name="Sample Count (N)")
    )
    combo_grp = combo_grp.sort_values(by="Sample Count (N)", ascending=False)
    combo_grp["Frequency (%)"] = combo_grp["Sample Count (N)"].apply(
        lambda x: f"{(x / tot_cohort * 100):.1f}%" if tot_cohort > 0 else "0.0%"
    )
    combo_grp["pct_val"] = combo_grp["Sample Count (N)"].apply(
        lambda x: round(x / tot_cohort * 100, 1) if tot_cohort > 0 else 0.0
    )

    # Add Combination Pattern string
    if len(active_snps) == 1:
        combo_grp["Combination Pattern"] = combo_grp[active_snps[0]].astype(str)
    elif len(active_snps) == 2:
        combo_grp["Combination Pattern"] = combo_grp.apply(
            lambda r: f"{r[active_snps[0]]} + {r[active_snps[1]]}", axis=1
        )
    else:
        combo_grp["Combination Pattern"] = combo_grp.apply(
            lambda r: " + ".join([str(r[c]) for c in active_snps]), axis=1
        )

    # Associated Diplotypes and Phenotypes for each combination in active_df
    def get_associated_dips(row):
        if "Diplotype" not in active_df.columns:
            return "N/A"
        cond = True
        for col in active_snps:
            cond = cond & (active_df[col] == row[col])
        sub = active_df[cond]
        dips = ", ".join(
            sorted([str(d) for d in sub["Diplotype"].dropna().unique()])
        )
        return dips if dips else "N/A"

    def get_associated_phes(row):
        if "Phenotype" not in active_df.columns:
            return "N/A"
        cond = True
        for col in active_snps:
            cond = cond & (active_df[col] == row[col])
        sub = active_df[cond]
        phes = ", ".join(
            sorted([str(p) for p in sub["Phenotype"].dropna().unique()])
        )
        return phes if phes else "N/A"

    combo_grp["Associated Diplotypes"] = combo_grp.apply(
        get_associated_dips, axis=1
    )
    combo_grp["CPIC Phenotype(s)"] = combo_grp.apply(
        get_associated_phes, axis=1
    )

    # Top KPI Metrics reflecting selected geography
    top_pattern = (
        combo_grp.iloc[0]["Combination Pattern"]
        if not combo_grp.empty
        else "N/A"
    )
    top_pct = (
        combo_grp.iloc[0]["Frequency (%)"] if not combo_grp.empty else "0%"
    )
    num_patterns = len(combo_grp)

    # -------------------------------------------------------------------------
    # COMPACT, ELEGANT KPI SUMMARY CARDS (NO OVERSIZED FONTS)
    # -------------------------------------------------------------------------
    st.markdown("<div style='margin-top: 8px;'></div>", unsafe_allow_html=True)
    k1, k2, k3, k4 = st.columns(4)

    card_style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 12px; height: 100%; box-shadow: 0 1px 2px rgba(0,0,0,0.03);"
    label_style = "font-size: 11px; font-weight: 600; text-transform: uppercase; color: #64748B; letter-spacing: 0.04em;"
    val_style = "font-size: 15px; font-weight: 700; color: #0F172A; margin-top: 3px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
    sub_style = "font-size: 11px; color: #475569; margin-top: 2px;"

    with k1:
        st.markdown(
            f"""
        <div style="{card_style}">
            <div style="{label_style}">Variant Target</div>
            <div style="{val_style}">{combo_title}</div>
            <div style="{sub_style}"><span style="color: #2563EB;">{std_snps_label}</span></div>
        </div>
        """,
            unsafe_allow_html=True,
        )
    with k2:
        st.markdown(
            f"""
        <div style="{card_style}">
            <div style="{label_style}">Geographic Scope</div>
            <div style="{val_style}">{scope_label}</div>
            <div style="{sub_style}"><span style="color: #059669;">Stratified Analysis</span></div>
        </div>
        """,
            unsafe_allow_html=True,
        )
    with k3:
        st.markdown(
            f"""
        <div style="{card_style}">
            <div style="{label_style}">Cohort Size</div>
            <div style="{val_style}">{tot_cohort:,} Samples</div>
            <div style="{sub_style}">{num_patterns} Distinct Patterns</div>
        </div>
        """,
            unsafe_allow_html=True,
        )
    with k4:
        st.markdown(
            f"""
        <div style="{card_style}">
            <div style="{label_style}">Predominant Genotype</div>
            <div style="{val_style}"><span style="color: #166534;">{top_pattern}</span></div>
            <div style="{sub_style}">{top_pct} of Active Cohort</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    st.markdown(
        "<div style='margin-top: 14px;'></div>", unsafe_allow_html=True
    )

    # -------------------------------------------------------------------------
    # 4. DUAL VISUAL CHARTS (DONUT / PIE CHART + COMPARATIVE BAR CHART)
    # -------------------------------------------------------------------------
    c_pie, c_bar = st.columns([1.0, 1.0])

    with c_pie:
        fig_donut = px.pie(
            combo_grp,
            names="Combination Pattern",
            values="Sample Count (N)",
            hole=0.48,
            title=f"Proportion Distribution ({scope_label})",
            color_discrete_sequence=[
                "#166534",
                "#1E40AF",
                "#D97706",
                "#991B1B",
                "#0D9488",
                "#6B21A8",
                "#475569",
            ],
        )
        fig_donut.update_layout(
            height=280,
            margin=dict(l=10, r=10, t=32, b=10),
            title_font=dict(
                size=13,
                color="#0F172A",
                family="Inter, -apple-system, sans-serif",
            ),
            legend=dict(font=dict(size=10), orientation="v"),
        )
        st.plotly_chart(fig_donut, use_container_width=True)

    with c_bar:
        fig_bars = px.bar(
            combo_grp,
            x="Sample Count (N)",
            y="Combination Pattern",
            orientation="h",
            text="Frequency (%)",
            title=f"Sample Counts & Frequencies ({scope_label})",
            color="Combination Pattern",
            color_discrete_sequence=[
                "#166534",
                "#1E40AF",
                "#D97706",
                "#991B1B",
                "#0D9488",
                "#6B21A8",
                "#475569",
            ],
        )
        fig_bars.update_layout(
            height=280,
            margin=dict(l=10, r=10, t=32, b=10),
            yaxis=dict(autorange="reversed", tickfont=dict(size=11)),
            xaxis=dict(tickfont=dict(size=11), gridcolor="#F1F5F9"),
            title_font=dict(
                size=13,
                color="#0F172A",
                family="Inter, -apple-system, sans-serif",
            ),
            showlegend=False,
        )
        fig_bars.update_traces(textposition="outside", textfont=dict(size=10))
        st.plotly_chart(fig_bars, use_container_width=True)

    # -------------------------------------------------------------------------
    # 5. GEOGRAPHICAL & STATE-WISE BREAKDOWN CHARTS & TABLES
    # -------------------------------------------------------------------------
    if "Geographic Groups" in geo_scope:
        st.markdown(f"#### Geographic Group Comparison for {combo_title}")
        if selected_region != "All Regions (Comparative Overview)":
            location_subset = (
                processed_df[
                    processed_df["Region"].astype(str) == str(selected_region)
                ]
                if "Region" in processed_df.columns
                else pd.DataFrame()
            )
            if state_col and state_col in location_subset.columns:
                location_counts = (
                    location_subset[state_col]
                    .value_counts()
                    .rename_axis("State / Locality")
                    .reset_index(name="Samples")
                )
                st.dataframe(
                    location_counts, use_container_width=True, hide_index=True
                )
        elif "Region" in processed_df.columns:
            pattern_series = processed_df[active_snps].apply(
                lambda row: " + ".join(row.astype(str)), axis=1
            )
            st.dataframe(
                pd.crosstab(
                    processed_df["Region"],
                    pattern_series,
                    margins=True,
                    margins_name="Total",
                ),
                use_container_width=True,
            )
            region_patterns = (
                processed_df.assign(_pattern=pattern_series)
                .groupby(["Region", "_pattern"])
                .size()
                .reset_index(name="Samples")
            )
            fig_region = px.bar(
                region_patterns,
                x="Region",
                y="Samples",
                color="_pattern",
                barmode="stack",
                title=f"Genotype Combinations by Geographic Group ({combo_title})",
                color_discrete_sequence=px.colors.qualitative.Safe,
            )
            st.plotly_chart(fig_region, use_container_width=True)
        else:
            st.info(
                "Map a geographic group column during upload to compare uploaded groups."
            )
    elif "State / Locality Groups" in geo_scope:
        st.markdown(f"#### State / Locality Analysis for {combo_title}")
        if (
            selected_state != "All States (Comparative Overview)"
            and state_col
            and state_col in active_df.columns
        ):
            group_col = (
                "Gender_Clean"
                if "Gender_Clean" in active_df.columns
                else ("Gender" if "Gender" in active_df.columns else None)
            )
            if group_col:
                group_counts = (
                    active_df[group_col]
                    .value_counts()
                    .rename_axis("Population Group")
                    .reset_index(name="Samples")
                )
                st.dataframe(
                    group_counts, use_container_width=True, hide_index=True
                )
        elif state_col and state_col in processed_df.columns:
            pattern_series = processed_df[active_snps].apply(
                lambda row: " + ".join(row.astype(str)), axis=1
            )
            st.dataframe(
                pd.crosstab(
                    processed_df[state_col],
                    pattern_series,
                    margins=True,
                    margins_name="Total",
                ),
                use_container_width=True,
            )
            state_counts = (
                processed_df[state_col]
                .value_counts()
                .rename_axis("State / Locality")
                .reset_index(name="Samples")
                .head(10)
            )
            fig_state = px.bar(
                state_counts,
                x="State / Locality",
                y="Samples",
                title="State / Locality Group Sizes",
                color="State / Locality",
                color_discrete_sequence=px.colors.qualitative.Safe,
            )
            st.plotly_chart(fig_state, use_container_width=True)
        else:
            st.info(
                "Map a state or locality column during upload to compare state-wise groups."
            )
    elif "Region" in processed_df.columns:
        region_counts = (
            processed_df["Region"]
            .value_counts()
            .rename_axis("Geographic Group")
            .reset_index(name="Samples")
        )
        fig_region = px.bar(
            region_counts,
            x="Geographic Group",
            y="Samples",
            title="Uploaded Geographic Group Distribution",
            color="Geographic Group",
            color_discrete_sequence=px.colors.qualitative.Safe,
        )
        st.plotly_chart(fig_region, use_container_width=True)

    # -------------------------------------------------------------------------
    # 6. STATISTICAL COMBINATION TABLE (STANDARDIZED SCIENTIFIC HEADERS)
    # -------------------------------------------------------------------------
    st.markdown(
        f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin: 14px 0 6px 0;'>📋 Statistical Combination Table — {combo_title} ({scope_label})</div>",
        unsafe_allow_html=True,
    )

    table_display = combo_grp.rename(
        columns={
            "Combination Pattern": "Genotype Combination",
            "Sample Count (N)": "Samples (N)",
            "CPIC Phenotype(s)": "Configured Phenotype",
        }
    )

    disp_cols = [c for c in table_display.columns if c != "pct_val"]
    st.dataframe(
        table_display[disp_cols], use_container_width=True, hide_index=True
    )


def render_diplotype_phenotype(
    active_gene: str,
    active_markers: List[str],
    active_snp_conf: Dict[str, Any],
    processed_df: pd.DataFrame,
    full_results: Dict[str, Any],
    render_subcohort_pills: Callable[[str], Tuple[str, Dict[str, Any]]],
) -> None:
    """Render the Diplotype and Phenotype translation and CPIC guidelines view."""
    st.markdown(f"### 💊 {active_gene} Calling & Translation")
    st.caption(
        f"Pharmacogenomic diplotype and metabolizer phenotype translation for target gene **{active_gene}**."
    )
    sel_key, sel_cohort = render_subcohort_pills("tab_pgx")
    pgx = (
        sel_cohort.get("pgx_summary")
        or sel_cohort.get("cyp2c9_summary")
        or sel_cohort.get("cyp2c19_summary", {})
    )
    dips = {
        k: v
        for k, v in pgx.get("diplotypes", {}).items()
        if k != "Not configured"
    }
    phenos = {
        k: v
        for k, v in pgx.get("phenotypes", {}).items()
        if k != "Not configured"
    }

    with st.expander(
        "⚙️ Configure diplotype and phenotype rules for this gene",
        expanded=(
            not bool(pgx.get("translation_configured"))
            or pgx.get("valid_pharmacogenomic_samples", 0) == 0
        ),
    ):
        st.caption(
            "Observed marker dosage patterns are calculated from the selected alleles and uploaded genotypes. Enter your validated diplotype and phenotype for each pattern; unconfigured patterns remain indeterminate."
        )
        dosage_maps = {}
        for marker in active_markers:
            marker_rule = active_snp_conf.get(marker, {}) or {}
            ref_allele = str(marker_rule.get("ref_allele", "")).strip().upper()
            var_allele = str(marker_rule.get("var_allele", "")).strip().upper()
            if not ref_allele or not var_allele or ref_allele == var_allele:
                continue
            dosage_maps[marker] = {
                f"{ref_allele}{ref_allele}": 0,
                f"{ref_allele}{var_allele}": 1,
                f"{var_allele}{ref_allele}": 1,
                f"{var_allele}{var_allele}": 2,
            }

        pattern_counts = Counter()
        if dosage_maps and not processed_df.empty:
            for _, sample in processed_df.iterrows():
                dosages = []
                for marker in dosage_maps:
                    raw_call = sample.get(marker)
                    if pd.isna(raw_call):
                        dosages = []
                        break
                    call = (
                        str(raw_call)
                        .strip()
                        .upper()
                        .replace("/", "")
                        .replace("|", "")
                        .replace(" ", "")
                    )
                    dosage = dosage_maps[marker].get(call)
                    if dosage is None and len(call) == 2:
                        dosage = dosage_maps[marker].get(call[::-1])
                    if dosage is None:
                        dosages = []
                        break
                    dosages.append(dosage)
                if dosages:
                    pattern_counts[",".join(map(str, dosages))] += 1

        if not dosage_maps:
            st.info(
                "Define distinct reference and variant alleles for at least one mapped marker before building translation patterns. Monomorphic markers remain available for population statistics."
            )
        elif not pattern_counts:
            st.info(
                "No complete, permitted genotype patterns were found for the selected markers."
            )
        else:
            st.write(
                "Marker order used in each dosage pattern:", ", ".join(dosage_maps)
            )
            pattern_rows = []
            existing_custom = st.session_state.get(
                "custom_translation_rules", {}
            ).get(str(active_gene).upper(), {})
            old_diplotype_rules = existing_custom.get("diplotype_rules", {})
            old_phenotype_map = existing_custom.get("phenotype_map", {})
            with st.form(
                f"translation_rules_form_{re.sub(r'[^A-Za-z0-9]+', '_', str(active_gene))}"
            ):
                st.markdown(
                    "**Observed dosage pattern · sample count · diplotype · phenotype**"
                )
                for pattern, count in sorted(pattern_counts.items()):
                    previous_diplotype = old_diplotype_rules.get(pattern, "")
                    previous_phenotype = old_phenotype_map.get(
                        previous_diplotype, ""
                    )
                    c_pattern, c_count, c_diplotype, c_phenotype = st.columns(
                        [1.25, 0.6, 1.1, 1.2]
                    )
                    c_pattern.write(pattern)
                    c_count.write(f"{count:,}")
                    gene_key = re.sub(r"[^A-Za-z0-9]+", "_", str(active_gene))
                    pattern_key = pattern.replace(",", "_")
                    diplotype_value = c_diplotype.text_input(
                        "Diplotype",
                        value=previous_diplotype,
                        key=f"rule_diplotype_{gene_key}_{pattern_key}",
                        label_visibility="collapsed",
                        placeholder="e.g. *1/*2",
                    )
                    phenotype_value = c_phenotype.text_input(
                        "Phenotype",
                        value=previous_phenotype,
                        key=f"rule_phenotype_{gene_key}_{pattern_key}",
                        label_visibility="collapsed",
                        placeholder="e.g. Intermediate",
                    )
                    pattern_rows.append(
                        (
                            pattern,
                            diplotype_value.strip(),
                            phenotype_value.strip(),
                        )
                    )
                save_rules = st.form_submit_button(
                    "Save translation rules and recalculate", type="primary"
                )
            if save_rules:
                diplotype_rules = {
                    pattern: diplotype
                    for pattern, diplotype, _ in pattern_rows
                    if diplotype
                }
                phenotype_map = {
                    diplotype: phenotype
                    for _, diplotype, phenotype in pattern_rows
                    if diplotype and phenotype
                }
                missing_phenotypes = [
                    pattern
                    for pattern, diplotype, phenotype in pattern_rows
                    if diplotype and not phenotype
                ]
                conflicting_phenotypes = any(
                    len(
                        {
                            phenotype
                            for _, dip, phenotype in pattern_rows
                            if dip == diplotype and phenotype
                        }
                    )
                    > 1
                    for diplotype in set(diplotype_rules.values())
                )
                if (
                    not diplotype_rules
                    or missing_phenotypes
                    or conflicting_phenotypes
                ):
                    st.error(
                        "Enter both a diplotype and phenotype for every dosage pattern you want classified."
                    )
                else:
                    st.session_state.setdefault(
                        "custom_translation_rules", {}
                    )[str(active_gene).upper()] = {
                        "markers": {
                            marker: {"genotype_dosage": dosage}
                            for marker, dosage in dosage_maps.items()
                        },
                        "diplotype_rules": diplotype_rules,
                        "phenotype_map": phenotype_map,
                        "phenotype_order": list(
                            dict.fromkeys(phenotype_map.values())
                        ),
                    }
                    save_active_session()
                    st.rerun()

    if not pgx.get("translation_configured"):
        st.info(
            "Diplotype and metabolizer classifications need a gene-specific, validated rule table. Marker and population statistics are available now; enter the interpretation rules above to calculate clinical categories."
        )
    elif (
        pgx.get("total_samples", 0)
        and pgx.get("valid_pharmacogenomic_samples", 0) == 0
    ):
        st.warning(
            "No determinate diplotype calls were produced. Check that the selected markers match the rule set and that every required genotype call is mapped, or configure rules for this uploaded marker panel above."
        )

    st.markdown("<br>", unsafe_allow_html=True)
    view_choice = st.radio(
        "🔎 Select View to Display:",
        [
            "👁️ Show All Views",
            "🧩 Combination",
            "📋 Diplotype distribution",
            "📊 Diplotype counts",
            "💊 Metabolizer / phenotype classification",
            "🗺️ Geographical & State-Wise Matrices",
        ],
        horizontal=True,
        key="tab_pgx_view_selector",
    )
    st.divider()

    # 1. Diplotype Distribution & Diplotype Counts
    if view_choice in [
        "👁️ Show All Views",
        "📋 Diplotype distribution",
        "📊 Diplotype counts",
    ]:
        if view_choice == "📋 Diplotype distribution":
            st.markdown("#### Diplotype distribution")
            dip_rows = [
                {
                    "DIPLOTYPE": k,
                    "COUNT": v["count"],
                    "FREQUENCY": f"{v['percentage']:.1f}%",
                }
                for k, v in dips.items()
                if v["count"] > 0
            ]
            st.table(pd.DataFrame(dip_rows))
        elif view_choice == "📊 Diplotype counts":
            st.markdown("#### Diplotype counts")
            dip_plot_rows = [
                {"Diplotype": str(k), "Count": int(v["count"])}
                for k, v in dips.items()
                if v.get("count", 0) > 0
            ]
            dip_fig_df = pd.DataFrame(
                dip_plot_rows, columns=["Diplotype", "Count"]
            )
            if not dip_fig_df.empty and len(dip_plot_rows) > 0:
                fig_dip = px.bar(
                    dip_fig_df,
                    x="Count",
                    y="Diplotype",
                    orientation="h",
                    title="Diplotype counts",
                    color_discrete_sequence=["#166534"],
                )
                fig_dip.update_layout(
                    height=280,
                    margin=dict(l=20, r=20, t=35, b=20),
                    yaxis=dict(autorange="reversed"),
                )
                st.plotly_chart(fig_dip, use_container_width=True)
            else:
                st.info("No diplotype counts to display.")
        else:
            c_dip, c_dip_chart = st.columns([1.1, 0.9])
            with c_dip:
                st.markdown("#### Diplotype distribution")
                dip_rows = [
                    {
                        "DIPLOTYPE": k,
                        "COUNT": v["count"],
                        "FREQUENCY": f"{v['percentage']:.1f}%",
                    }
                    for k, v in dips.items()
                    if v["count"] > 0
                ]
                st.table(pd.DataFrame(dip_rows))
            with c_dip_chart:
                dip_plot_rows = [
                    {"Diplotype": str(k), "Count": int(v["count"])}
                    for k, v in dips.items()
                    if v.get("count", 0) > 0
                ]
                dip_fig_df = pd.DataFrame(
                    dip_plot_rows, columns=["Diplotype", "Count"]
                )
                if not dip_fig_df.empty and len(dip_plot_rows) > 0:
                    fig_dip = px.bar(
                        dip_fig_df,
                        x="Count",
                        y="Diplotype",
                        orientation="h",
                        title="Diplotype counts",
                        color_discrete_sequence=["#166534"],
                    )
                    fig_dip.update_layout(
                        height=260,
                        margin=dict(l=20, r=20, t=35, b=20),
                        yaxis=dict(autorange="reversed"),
                    )
                    st.plotly_chart(fig_dip, use_container_width=True)
                else:
                    st.info("No diplotype counts to display.")

    # 2. Metabolizer / Phenotype Classification
    if view_choice in [
        "👁️ Show All Views",
        "💊 Metabolizer / phenotype classification",
    ]:
        if view_choice == "👁️ Show All Views":
            st.markdown("<br>", unsafe_allow_html=True)
        c_phe, c_phe_chart = st.columns([1.1, 0.9])
        with c_phe:
            st.markdown("#### Metabolizer / phenotype classification")
            phe_rows = [
                {
                    "PHENOTYPE": k,
                    "COUNT": v["count"],
                    "PERCENTAGE": f"{v['percentage']:.1f}%",
                }
                for k, v in phenos.items()
                if v["count"] > 0
            ]
            st.table(pd.DataFrame(phe_rows))
        with c_phe_chart:
            phe_labels = [
                str(k) for k, v in phenos.items() if v.get("count", 0) > 0
            ]
            phe_counts = [
                int(v["count"])
                for k, v in phenos.items()
                if v.get("count", 0) > 0
            ]
            if phe_labels and phe_counts and sum(phe_counts) > 0:
                fig_donut = px.pie(
                    names=phe_labels,
                    values=phe_counts,
                    hole=0.5,
                    title="Phenotype distribution",
                    color_discrete_sequence=[
                        "#166534",
                        "#1E40AF",
                        "#4B5563",
                        "#991B1B",
                        "#D97706",
                        "#0D9488",
                    ],
                )
                fig_donut.update_layout(
                    height=260, margin=dict(l=20, r=20, t=35, b=20)
                )
                st.plotly_chart(fig_donut, use_container_width=True)
            else:
                st.info("No phenotype calls to display.")

    # 3. Combination & Geographic Studio
    if view_choice in [
        "👁️ Show All Views",
        "🧩 Combination",
        "🗺️ Geographical & State-Wise Matrices",
    ]:
        render_cyp2c19_combination_studio(
            processed_df, full_results, view_mode=view_choice
        )

    # 4. Clinical CPIC Dosing Guidelines
    if active_gene == "CYP2C9":
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### 💊 CPIC® Clinical Dosing Guidelines for CYP2C9")
        st.caption(
            "Standardized CPIC level A/B clinical recommendations for narrow therapeutic index substrates."
        )
        c_g1, c_g2, c_g3 = st.columns(3)
        with c_g1:
            st.markdown(
                """
            <div style="background: #FEF3C7; border: 1px solid #F59E0B; border-radius: 8px; padding: 14px;">
                <div style="font-weight: 700; color: #92400E; font-size: 14px;">🩸 Warfarin</div>
                <div style="font-size: 12px; color: #78350F; margin-top: 6px; line-height: 1.5;">
                    <b>Normal (*1/*1):</b> Standard starting dose (titrate INR).<br>
                    <b>Intermediate (*1/*2, *1/*3, *2/*2):</b> Reduce initial dose by 20–50% to prevent supratherapeutic INR and bleeding.<br>
                    <b>Poor (*2/*3, *3/*3):</b> Reduce dose by 50–80% or consider alternative oral anticoagulant (e.g., DOAC).
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )
        with c_g2:
            st.markdown(
                """
            <div style="background: #EFF6FF; border: 1px solid #3B82F6; border-radius: 8px; padding: 14px;">
                <div style="font-weight: 700; color: #1E40AF; font-size: 14px;">🧠 Phenytoin</div>
                <div style="font-size: 12px; color: #1E3A8A; margin-top: 6px; line-height: 1.5;">
                    <b>Normal (*1/*1):</b> Standard starting dose (300 mg/day).<br>
                    <b>Intermediate (*1/*2, *1/*3, *2/*2):</b> Reduce initial maintenance dose by 25% to prevent neurotoxicity.<br>
                    <b>Poor (*2/*3, *3/*3):</b> Reduce initial maintenance dose by 50% and closely monitor serum trough.
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )
        with c_g3:
            st.markdown(
                """
            <div style="background: #ECFDF5; border: 1px solid #10B981; border-radius: 8px; padding: 14px;">
                <div style="font-weight: 700; color: #065F46; font-size: 14px;">🩹 Celecoxib</div>
                <div style="font-size: 12px; color: #064E3B; margin-top: 6px; line-height: 1.5;">
                    <b>Normal (*1/*1):</b> Standard starting dose (200 mg/day).<br>
                    <b>Intermediate (*1/*2, *1/*3, *2/*2):</b> Standard dose with cautious clinical titration.<br>
                    <b>Poor (*2/*3, *3/*3):</b> Reduce initial dose by 50% to prevent prolonged exposure and toxicity.
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )
    elif active_gene == "CYP2C19":
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### 💊 CPIC® Clinical Dosing Guidelines for CYP2C19")
        st.caption(
            "Standardized CPIC level A/B clinical recommendations for antiplatelet and psychiatric therapies."
        )
        c_g1, c_g2 = st.columns(2)
        with c_g1:
            st.markdown(
                """
            <div style="background: #EFF6FF; border: 1px solid #3B82F6; border-radius: 8px; padding: 14px;">
                <div style="font-weight: 700; color: #1E40AF; font-size: 14px;">🩸 Clopidogrel (Plavix)</div>
                <div style="font-size: 12px; color: #1E3A8A; margin-top: 6px; line-height: 1.5;">
                    <b>Normal (*1/*1) & Rapid (*1/*17):</b> Standard dose (75 mg/day).<br>
                    <b>Intermediate (*1/*2, *1/*3) & Poor (*2/*2, *2/*3, *3/*3):</b> Avoid clopidogrel; prescribe alternative antiplatelet (Prasugrel, Ticagrelor) due to reduced active metabolite and high stent thrombosis risk.
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )
        with c_g2:
            st.markdown(
                """
            <div style="background: #FDF2F8; border: 1px solid #EC4899; border-radius: 8px; padding: 14px;">
                <div style="font-weight: 700; color: #9D174D; font-size: 14px;">💊 Voriconazole & Escitalopram</div>
                <div style="font-size: 12px; color: #831843; margin-top: 6px; line-height: 1.5;">
                    <b>Voriconazole:</b> In Poor Metabolizers, select alternative antifungal (Isavuconazole) to prevent neurotoxicity.<br>
                    <b>Escitalopram:</b> In Poor Metabolizers, reduce starting dose by 50% to avoid QT prolongation.
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(
            f"### 📚 Pharmacogenomic Knowledgebase Guidelines & Reference for {active_gene}"
        )
        st.caption(
            f"Curated population markers and reference allele frequencies for {active_gene}."
        )
        catalog_df = load_marker_catalog()
        sub_catalog = (
            catalog_df[
                catalog_df["Gene"].str.upper() == str(active_gene).upper()
            ]
            if not catalog_df.empty
            else pd.DataFrame()
        )
        if not sub_catalog.empty:
            disp_cat = sub_catalog[
                [
                    "Gene",
                    "Marker",
                    "SA MAF",
                    "%",
                    "MAF",
                    "Level of evidence",
                    "PubMed reference",
                ]
            ].copy()
            st.dataframe(disp_cat, use_container_width=True, hide_index=True)
        else:
            st.info(
                f"Target gene {active_gene} analyzed successfully. Refer to CPIC / PharmGKB guidelines for drug dosing recommendations."
            )
