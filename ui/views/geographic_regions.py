"""Geographic Population Groups analytical view component.

Provides interactive selection options:
1. Multi-Region Overview & Comparison (with region filter)
2. Specific Region Deep-Dive (individual region breakdown)
3. Biomarker Parameter Matrix (detailed marker counts and freqs)
4. Phenotype & Diplotype Distribution (clinical metabolizer shift)
5. Regional Hardy-Weinberg Equilibrium (subpopulation HWE & Wahlund audit)
"""

from typing import Any, Dict, List
import pandas as pd
import plotly.express as px
import streamlit as st


def render_geographic_regions(
    active_gene: str,
    active_markers: List[str],
    active_snp_conf: Dict[str, Any],
    full_results: Dict[str, Any],
) -> None:
    """Render geographic population groups comparison with interactive selection options."""
    st.markdown(f"### 🌐 Geographic Population Groups Analysis ({active_gene})")
    reg = full_results.get("regional", {})
    ov = full_results.get("overall", {})
    regions = sorted(
        r for r, item in reg.items() if item.get("sample_count", 0) > 0
    )
    if not regions:
        st.info(
            "No geographic groups are available. Map a region or state column during upload to enable this analysis."
        )
        return

    phenotype_categories = list(
        dict.fromkeys(
            phenotype
            for group in reg.values()
            for phenotype in (
                (
                    group.get("pgx_summary")
                    or group.get(f"{active_gene.lower()}_summary")
                    or {}
                ).get("phenotypes", {})
            )
            if phenotype != "Not configured"
        )
    )

    # -------------------------------------------------------------------------
    # INTERACTIVE VIEW MODE SELECTOR
    # -------------------------------------------------------------------------
    geo_view_mode = st.radio(
        "Select Geographic Analysis View Option:",
        [
            "📊 Multi-Region Overview & Comparison",
            "🔍 Specific Region Deep-Dive",
            "🧬 Biomarker Parameter Matrix",
            "💊 Phenotype & Diplotype Distribution",
            "⚖️ Regional Hardy-Weinberg Equilibrium",
        ],
        horizontal=True,
        key="geo_view_mode_selector",
    )
    st.markdown("<div style='margin-bottom: 12px;'></div>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # OPTION 1: MULTI-REGION OVERVIEW & COMPARISON
    # -------------------------------------------------------------------------
    if geo_view_mode == "📊 Multi-Region Overview & Comparison":
        selected_regions = st.multiselect(
            "Filter Geographic Regions to Include:",
            regions,
            default=regions,
            key="geo_multiselect_regions",
        )
        if not selected_regions:
            st.warning("Please select at least one geographic region to display.")
            return

        # KPI Cards for Selected Regions
        cols = st.columns(min(5, len(selected_regions)))
        for idx, r in enumerate(selected_regions):
            r_data = reg.get(r, {})
            snps_r = r_data.get("snps", {})
            pgx_r = (
                r_data.get("pgx_summary")
                or r_data.get(f"{active_gene.lower()}_summary")
                or r_data.get("cyp2c19_summary", {})
            )
            phenotype_values = pgx_r.get("phenotypes", {})
            leading_phenotype = (
                max(
                    phenotype_values.items(),
                    key=lambda item: item[1].get("count", 0),
                )
                if phenotype_values
                else (None, {})
            )

            card_lines = []
            for m in active_markers:
                m_conf = active_snp_conf.get(m, {})
                var_al = m_conf.get("var_allele")
                al_freqs = snps_r.get(m, {}).get("allele_freqs", {})
                if not var_al:
                    ref_al = m_conf.get("ref_allele")
                    non_refs = [a for a in al_freqs.keys() if a != ref_al]
                    var_al = (
                        non_refs[0]
                        if non_refs
                        else (list(al_freqs.keys())[0] if al_freqs else "Var")
                    )
                maf_val = snps_r.get(m, {}).get("allele_freqs", {}).get(
                    var_al, snps_r.get(m, {}).get("maf", 0)
                )
                card_lines.append(
                    f"{m} ({var_al}): <strong>{maf_val}</strong>"
                )

            lines_html = "<br>".join(card_lines)
            if lines_html:
                lines_html += "<br>"

            with cols[idx % len(cols)]:
                st.markdown(
                    f"""
                <div class="kpi-card" style="background: #FFFFFF; border: 1.5px solid #E2E8F0; border-radius: 8px; padding: 12px; margin-bottom: 10px;">
                    <div class="kpi-lbl" style="font-size: 11px; font-weight: 700; color: #64748B; text-transform: uppercase;">{r}</div>
                    <div class="kpi-val" style="font-size: 18px; font-weight: 700; color: #0F172A; margin: 4px 0;">{r_data.get('sample_count', 0):,}</div>
                    <div style="font-size: 0.78rem; color: #475569; line-height: 1.4;">
                        {lines_html}{f"Top phenotype: {leading_phenotype[0]} ({leading_phenotype[1].get('percentage', 0):.1f}%)" if leading_phenotype[0] else "Phenotype translation not configured"}
                    </div>
                </div>
                """,
                    unsafe_allow_html=True,
                )

        st.markdown("<br>", unsafe_allow_html=True)

        c_geo1, c_geo2 = st.columns(2)
        with c_geo1:
            geo_marker_data = []
            for r in selected_regions:
                r_snps = reg.get(r, {}).get("snps", {})
                for m in active_markers:
                    m_conf = active_snp_conf.get(m, {})
                    var_al = m_conf.get("var_allele")
                    al_freqs = r_snps.get(m, {}).get("allele_freqs", {})
                    if not var_al:
                        ref_al = m_conf.get("ref_allele")
                        non_refs = [a for a in al_freqs.keys() if a != ref_al]
                        var_al = (
                            non_refs[0]
                            if non_refs
                            else (list(al_freqs.keys())[0] if al_freqs else "Var")
                        )
                    maf = float(
                        al_freqs.get(var_al, r_snps.get(m, {}).get("maf", 0.0))
                    )
                    geo_marker_data.append(
                        {
                            "Geographic Group": r,
                            "Biomarker": f"{m} (f_{var_al})",
                            "Minor Allele Frequency": maf,
                        }
                    )
            geo_df = pd.DataFrame(geo_marker_data)
            if not geo_df.empty:
                marker_colors = ["#2563EB", "#059669", "#D97706", "#9333EA", "#DC2626"]
                fig_geo_all = px.bar(
                    geo_df,
                    x="Geographic Group",
                    y="Minor Allele Frequency",
                    color="Biomarker",
                    barmode="group",
                    title=f"Minor Allele Frequencies Across Selected Regions",
                    color_discrete_sequence=marker_colors,
                )
                fig_geo_all.update_layout(height=300, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_geo_all, use_container_width=True)

        with c_geo2:
            geo_pheno_data = []
            for r in selected_regions:
                r_data = reg.get(r, {})
                p_dict = (
                    r_data.get("pgx_summary")
                    or r_data.get(f"{active_gene.lower()}_summary")
                    or r_data.get("cyp2c19_summary", {})
                ).get("phenotypes", {})
                for phenotype in phenotype_categories:
                    metric = p_dict.get(phenotype, {})
                    geo_pheno_data.append(
                        {
                            "Geographic Group": r,
                            "Phenotype": phenotype,
                            "Percentage": float(metric.get("percentage", 0.0)),
                        }
                    )
            geo_pm_df = pd.DataFrame(geo_pheno_data)
            if not geo_pm_df.empty:
                fig_geo_pm = px.bar(
                    geo_pm_df,
                    x="Geographic Group",
                    y="Percentage",
                    color="Phenotype",
                    barmode="group",
                    title=f"Phenotype Distribution Across Selected Regions ({active_gene})",
                    color_discrete_sequence=px.colors.qualitative.Safe,
                )
                fig_geo_pm.update_layout(height=300, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_geo_pm, use_container_width=True)

    # -------------------------------------------------------------------------
    # OPTION 2: SPECIFIC REGION DEEP-DIVE
    # -------------------------------------------------------------------------
    elif geo_view_mode == "🔍 Specific Region Deep-Dive":
        chosen_reg = st.selectbox(
            "Select Target Geographic Region to Inspect:",
            regions,
            key="geo_deepdive_reg_sel",
        )
        r_data = reg.get(chosen_reg, {})
        n_reg = r_data.get("sample_count", 0)
        snps_r = r_data.get("snps", {})
        pgx_r = (
            r_data.get("pgx_summary")
            or r_data.get(f"{active_gene.lower()}_summary")
            or {}
        )
        phenos_r = pgx_r.get("phenotypes", {})
        im_n = phenos_r.get("Intermediate Metabolizer", {}).get("count", 0)
        pm_n = phenos_r.get("Poor Metabolizer", {}).get("count", 0)
        nm_n = phenos_r.get("Normal Metabolizer", {}).get("count", 0)
        act_n = im_n + pm_n
        act_pct = (act_n / n_reg * 100.0) if n_reg > 0 else 0.0

        # KPI Row for Selected Region
        k1, k2, k3, k4 = st.columns(4)
        card_style = "background: #FFFFFF; border: 1.5px solid #E2E8F0; border-radius: 8px; padding: 12px; box-shadow: 0 1px 2px rgba(0,0,0,0.05);"
        with k1:
            st.markdown(
                f"""
            <div style="{card_style}">
                <div style="font-size: 11px; font-weight: 700; color: #64748B;">REGION</div>
                <div style="font-size: 18px; font-weight: 800; color: #1E40AF; margin-top: 4px;">{chosen_reg}</div>
                <div style="font-size: 12px; color: #475569;">{n_reg:,} Total Samples</div>
            </div>
            """,
                unsafe_allow_html=True,
            )
        with k2:
            st.markdown(
                f"""
            <div style="{card_style}">
                <div style="font-size: 11px; font-weight: 700; color: #64748B;">NORMAL METABOLIZERS (NM)</div>
                <div style="font-size: 18px; font-weight: 800; color: #059669; margin-top: 4px;">{nm_n}/{n_reg}</div>
                <div style="font-size: 12px; color: #475569;">{(nm_n/n_reg*100.0 if n_reg>0 else 0.0):.1f}% of cohort</div>
            </div>
            """,
                unsafe_allow_html=True,
            )
        with k3:
            st.markdown(
                f"""
            <div style="{card_style}">
                <div style="font-size: 11px; font-weight: 700; color: #64748B;">ACTIONABLE (IM + PM)</div>
                <div style="font-size: 18px; font-weight: 800; color: #D97706; margin-top: 4px;">{act_n}/{n_reg}</div>
                <div style="font-size: 12px; color: #475569;">{act_pct:.1f}% require adjustment</div>
            </div>
            """,
                unsafe_allow_html=True,
            )
        with k4:
            st.markdown(
                f"""
            <div style="{card_style}">
                <div style="font-size: 11px; font-weight: 700; color: #64748B;">POOR METABOLIZERS (PM)</div>
                <div style="font-size: 18px; font-weight: 800; color: #DC2626; margin-top: 4px;">{pm_n}/{n_reg}</div>
                <div style="font-size: 12px; color: #475569;">{(pm_n/n_reg*100.0 if n_reg>0 else 0.0):.1f}% high risk</div>
            </div>
            """,
                unsafe_allow_html=True,
            )

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(f"#### 🧬 Biomarker Allele & Genotype Breakdown in {chosen_reg}")

        m_deep_rows = []
        for m in active_markers:
            m_res = snps_r.get(m, {})
            m_conf = active_snp_conf.get(m, {})
            al_counts = m_res.get("allele_counts", {})
            gt_counts = m_res.get("genotype_counts", {})
            hw = m_res.get("hwe", {})

            ref_al = m_conf.get("ref_allele")
            var_al = m_conf.get("var_allele")
            if not ref_al or not var_al:
                sorted_als = sorted(al_counts.items(), key=lambda x: x[1], reverse=True)
                ref_al = ref_al or (sorted_als[0][0] if sorted_als else "Ref")
                var_al = var_al or (sorted_als[1][0] if len(sorted_als) > 1 else "Var")

            wt_g = m_conf.get("wildtype_genotype") or f"{ref_al}{ref_al}"
            het_g = m_conf.get("het_genotype") or f"{ref_al}{var_al}"
            var_g = m_conf.get("hom_var_genotype") or f"{var_al}{var_al}"

            c_wt = gt_counts.get(wt_g, 0)
            c_het = gt_counts.get(het_g, 0)
            c_var = gt_counts.get(var_g, 0)
            c_car = c_het + c_var
            car_pct = (c_car / n_reg * 100.0) if n_reg > 0 else 0.0

            ref_cnt = al_counts.get(ref_al, 0)
            var_cnt = al_counts.get(var_al, 0)
            tot_al = ref_cnt + var_cnt
            maf = (var_cnt / tot_al) if tot_al > 0 else 0.0

            exact_p = hw.get("exact_p_value", 1.0)
            f_is = hw.get("inbreeding_coeff", 0.0)

            m_deep_rows.append(
                {
                    "Biomarker": m,
                    "Ref / Var Allele": f"{ref_al} / {var_al}",
                    "Ref Count (2N)": f"{ref_cnt}/{tot_al}",
                    "Var Count (2N)": f"{var_cnt}/{tot_al}",
                    "Minor Allele Freq (MAF)": f"{maf:.4f}",
                    f"Wildtype {wt_g}": c_wt,
                    f"Heterozygote {het_g}": c_het,
                    f"Homozygous Var {var_g}": c_var,
                    "Carrier Count": f"{c_car}/{n_reg}",
                    "Carrier Rate (%)": f"{car_pct:.1f}%",
                    "HWE Exact P": f"{exact_p:.4f}",
                    "Inbreeding (F_IS)": f"{f_is:+.4f}",
                    "HWE Status": "✅ In HWE" if exact_p >= 0.05 else "⚠️ Departure (p < 0.05)",
                }
            )

        if m_deep_rows:
            st.dataframe(pd.DataFrame(m_deep_rows), use_container_width=True, hide_index=True)

        # Phenotype and Diplotype side-by-side charts for this region
        c_dp1, c_dp2 = st.columns(2)
        with c_dp1:
            p_rows = [
                {"Phenotype": p, "Count": data.get("count", 0), "Percentage": data.get("percentage", 0.0)}
                for p, data in phenos_r.items()
                if p != "Not configured"
            ]
            if p_rows:
                fig_p = px.bar(
                    pd.DataFrame(p_rows),
                    x="Phenotype",
                    y="Percentage",
                    text="Percentage",
                    title=f"Phenotype Distribution in {chosen_reg}",
                    color="Phenotype",
                    color_discrete_sequence=px.colors.qualitative.Safe,
                )
                fig_p.update_traces(texttemplate="%{text:.1f}%", textposition="outside")
                fig_p.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20), showlegend=False)
                st.plotly_chart(fig_p, use_container_width=True)

        with c_dp2:
            dips_r = pgx_r.get("diplotypes", {})
            d_rows = [
                {"Diplotype": d, "Count": data.get("count", 0), "Percentage": data.get("percentage", 0.0)}
                for d, data in dips_r.items()
                if d != "Not configured"
            ]
            if d_rows:
                fig_d = px.bar(
                    pd.DataFrame(d_rows),
                    x="Diplotype",
                    y="Percentage",
                    text="Percentage",
                    title=f"Diplotype Distribution in {chosen_reg}",
                    color="Diplotype",
                    color_discrete_sequence=px.colors.qualitative.Bold,
                )
                fig_d.update_traces(texttemplate="%{text:.1f}%", textposition="outside")
                fig_d.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20), showlegend=False)
                st.plotly_chart(fig_d, use_container_width=True)

    # -------------------------------------------------------------------------
    # OPTION 3: BIOMARKER PARAMETER MATRIX
    # -------------------------------------------------------------------------
    elif geo_view_mode == "🧬 Biomarker Parameter Matrix":
        st.markdown("#### Comprehensive Population Statistical Summary Matrix")
        st.caption("Displays complete observed allele, genotype, carrier, and HWE counts across all geographic zones.")

        def render_marker_geo_table(snp_key: str) -> None:
            snp_conf_sel = active_snp_conf.get(snp_key, {})
            ov_snp = ov.get("snps", {}).get(snp_key, {})
            al_counts = ov_snp.get("allele_counts", {})
            gt_counts = ov_snp.get("genotype_counts", {})

            ref_al = snp_conf_sel.get("ref_allele")
            var_al = snp_conf_sel.get("var_allele")
            if not ref_al or not var_al:
                sorted_als = sorted(al_counts.items(), key=lambda x: x[1], reverse=True)
                ref_al = ref_al or (sorted_als[0][0] if sorted_als else "Ref")
                var_al = var_al or (sorted_als[1][0] if len(sorted_als) > 1 else "Var")

            wt_g = snp_conf_sel.get("wildtype_genotype") or (ref_al + ref_al)
            var_g = snp_conf_sel.get("hom_var_genotype") or (var_al + var_al)
            het_g = snp_conf_sel.get("het_genotype") or f"{ref_al}{var_al}"

            rows_cfg = [
                (f"{var_al} Count (Variant Allele)", lambda s: s.get(snp_key, {}).get("allele_counts", {}).get(var_al, 0)),
                (f"{ref_al} Count (Reference Allele)", lambda s: s.get(snp_key, {}).get("allele_counts", {}).get(ref_al, 0)),
                ("Total Allele Count (2N)", lambda s: sum(s.get(snp_key, {}).get("allele_counts", {}).values())),
                (f"{het_g} COUNT (Heterozygote)", lambda s: s.get(snp_key, {}).get("genotype_counts", {}).get(het_g, 0)),
                (f"{var_g} COUNT (Homozygous Variant)", lambda s: s.get(snp_key, {}).get("genotype_counts", {}).get(var_g, 0)),
                (f"{wt_g} COUNT (Homozygous Wildtype)", lambda s: s.get(snp_key, {}).get("genotype_counts", {}).get(wt_g, 0)),
                ("Total Sample N", lambda s: sum(s.get(snp_key, {}).get("genotype_counts", {}).values())),
                (f"Allele Frequency F({var_al})", lambda s: f"{s.get(snp_key, {}).get('allele_freqs', {}).get(var_al, 0.0):.4f}"),
                (f"Allele Frequency F({ref_al})", lambda s: f"{s.get(snp_key, {}).get('allele_freqs', {}).get(ref_al, 0.0):.4f}"),
                ("Chi-Square (χ²)", lambda s: f"{s.get(snp_key, {}).get('hwe', {}).get('chi2_stat', 0.0):.3f}"),
                ("Exact P-Value", lambda s: f"{s.get(snp_key, {}).get('hwe', {}).get('exact_p_value', s.get(snp_key, {}).get('hwe', {}).get('p_value', 1.0)):.4f}"),
                ("Inbreeding (F_IS)", lambda s: f"{s.get(snp_key, {}).get('hwe', {}).get('inbreeding_coeff', 0.0):+.4f}"),
                (
                    f"Carrier Count ({het_g} + {var_g})",
                    lambda s: (
                        s.get(snp_key, {}).get("genotype_counts", {}).get(het_g, 0)
                        + s.get(snp_key, {}).get("genotype_counts", {}).get(var_g, 0)
                    ),
                ),
                (
                    "Carrier Frequency (%)",
                    lambda s: f"{(s.get(snp_key, {}).get('genotype_counts', {}).get(het_g, 0) + s.get(snp_key, {}).get('genotype_counts', {}).get(var_g, 0)) / max(1, sum(s.get(snp_key, {}).get('genotype_counts', {}).values())) * 100.0:.1f}%",
                ),
            ]

            m_matrix = []
            for lbl, fn in rows_cfg:
                row_item = {f"Statistical Parameter ({snp_key})": lbl}
                for r in regions:
                    row_item[r] = fn(reg.get(r, {}).get("snps", {}))
                row_item["Overall Population"] = fn(ov.get("snps", {}))
                m_matrix.append(row_item)
            st.dataframe(pd.DataFrame(m_matrix), use_container_width=True, hide_index=True)

        marker_tab_labels = [f"🧬 {m}" for m in active_markers]
        if marker_tab_labels:
            marker_tabs = st.tabs(marker_tab_labels)
            for m_idx, m_name in enumerate(active_markers):
                with marker_tabs[m_idx]:
                    render_marker_geo_table(m_name)

    # -------------------------------------------------------------------------
    # OPTION 4: PHENOTYPE & DIPLOTYPE DISTRIBUTION
    # -------------------------------------------------------------------------
    elif geo_view_mode == "💊 Phenotype & Diplotype Distribution":
        st.markdown(f"#### 🌐 Geographic Region Metabolizer Phenotype Distribution Matrix ({active_gene})")
        st.caption(f"Observed {active_gene} phenotype categories across the geographic groups present in this dataset.")

        pheno_order = phenotype_categories
        geo_pheno_matrix = []
        for r_name in regions:
            r_data = reg.get(r_name, {})
            p_dict = (
                r_data.get("pgx_summary")
                or r_data.get(f"{active_gene.lower()}_summary")
                or r_data.get("cyp2c19_summary", {})
            ).get("phenotypes", {})
            row = {
                "Geographic Region": r_name,
                "Total Sample N": r_data.get("sample_count", 0),
            }
            for p_cat in pheno_order:
                cnt = p_dict.get(p_cat, {}).get("count", 0)
                pct = p_dict.get(p_cat, {}).get("percentage", 0.0)
                p_label = p_cat.split(" (")[0]
                row[p_label] = f"{cnt} ({pct:.1f}%)"
            geo_pheno_matrix.append(row)

        tot_p = (
            ov.get("pgx_summary")
            or ov.get(f"{active_gene.lower()}_summary")
            or ov.get("cyp2c19_summary", {})
        ).get("phenotypes", {})
        tot_row = {
            "Geographic Region": "Total Cohort (All Regions Combined)",
            "Total Sample N": ov.get("sample_count", 0),
        }
        for p_cat in pheno_order:
            cnt = tot_p.get(p_cat, {}).get("count", 0)
            pct = tot_p.get(p_cat, {}).get("percentage", 0.0)
            p_label = p_cat.split(" (")[0]
            tot_row[p_label] = f"{cnt} ({pct:.1f}%)"
        geo_pheno_matrix.append(tot_row)

        st.dataframe(pd.DataFrame(geo_pheno_matrix), use_container_width=True, hide_index=True)

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(f"#### Diplotype Distribution by Geographic Group ({active_gene})")
        geo_diplotype_categories = list(
            dict.fromkeys(
                name
                for group in reg.values()
                for name in (
                    (group.get("pgx_summary") or group.get(f"{active_gene.lower()}_summary") or {}).get("diplotypes", {})
                )
                if name != "Not configured"
            )
        )
        geo_diplotype_matrix = []
        for r_name in regions:
            r_data = reg.get(r_name, {})
            dip_data = (
                r_data.get("pgx_summary") or r_data.get(f"{active_gene.lower()}_summary") or {}
            ).get("diplotypes", {})
            row = {
                "Geographic Group": r_name,
                "Total Sample N": r_data.get("sample_count", 0),
            }
            for category in geo_diplotype_categories:
                metric = dip_data.get(category, {})
                row[category] = f"{metric.get('count', 0)} ({metric.get('percentage', 0.0):.1f}%)"
            geo_diplotype_matrix.append(row)
        st.dataframe(pd.DataFrame(geo_diplotype_matrix), use_container_width=True, hide_index=True)

    # -------------------------------------------------------------------------
    # OPTION 5: REGIONAL HARDY-WEINBERG EQUILIBRIUM
    # -------------------------------------------------------------------------
    elif geo_view_mode == "⚖️ Regional Hardy-Weinberg Equilibrium":
        st.markdown("#### ⚖️ Regional Hardy-Weinberg Equilibrium & Wahlund Effect Matrix")
        st.caption(
            "Tests whether genotype frequencies in each individual geographic zone conform to HWE "
            "or whether population substructure causes a Wahlund effect (heterozygote deficiency)."
        )

        hwe_geo_rows = []
        for r_name in regions:
            r_data = reg.get(r_name, {})
            r_n = r_data.get("sample_count", 0)
            r_snps = r_data.get("snps", {})

            for m in active_markers:
                m_res = r_snps.get(m, {})
                hw = m_res.get("hwe", {})
                gc = m_res.get("genotype_counts", {})
                exp_c = hw.get("expected_counts", {})
                chi2_v = hw.get("chi2_stat", 0.0)
                exact_p = hw.get("exact_p_value", hw.get("p_value", 1.0))
                f_is = hw.get("inbreeding_coeff", 0.0)

                hwe_geo_rows.append(
                    {
                        "Geographic Region": r_name,
                        "Biomarker": m,
                        "Sample N": r_n,
                        "Observed Genotypes": " / ".join([f"{k}:{v}" for k, v in gc.items()]),
                        "Expected Genotypes": " / ".join([f"{k}:{exp_c.get(k, 0)}" for k in gc.keys()]),
                        "Chi-Square (χ²)": f"{chi2_v:.3f}",
                        "Exact P-Value": f"{exact_p:.4f}",
                        "Inbreeding (F_IS)": f"{f_is:+.4f}",
                        "HWE Status": "✅ In HWE" if exact_p >= 0.05 else "⚠️ Departure (p < 0.05)",
                    }
                )

        if hwe_geo_rows:
            st.dataframe(pd.DataFrame(hwe_geo_rows), use_container_width=True, hide_index=True)
