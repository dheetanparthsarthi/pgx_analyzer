"""Geographic Population Groups analytical view component."""

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
    """Render geographic population groups comparison, MAF breakdown, and phenotype distribution."""
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
    cols = st.columns(min(5, len(regions))) if regions else []
    for idx, r in enumerate(regions):
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
                f"{m} frequency ({var_al}): <strong>{maf_val}</strong>"
            )

        lines_html = "<br>".join(card_lines)
        if lines_html:
            lines_html += "<br>"

        with cols[idx % len(cols)]:
            st.markdown(
                f"""
            <div class="kpi-card">
                <div class="kpi-lbl">{r}</div>
                <div class="kpi-val">{r_data.get('sample_count', 0):,}</div>
                <div style="font-size: 0.78rem; color: #475569; margin-top: 0.3rem; line-height: 1.4;">
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
        for r in regions:
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
        geo_df = pd.DataFrame(
            geo_marker_data,
            columns=[
                "Geographic Group",
                "Biomarker",
                "Minor Allele Frequency",
            ],
        )
        if not geo_df.empty:
            marker_colors = [
                "#2563EB",
                "#059669",
                "#D97706",
                "#9333EA",
                "#DC2626",
                "#0891B2",
            ]
            fig_geo_all = px.bar(
                geo_df,
                x="Geographic Group",
                y="Minor Allele Frequency",
                color="Biomarker",
                barmode="group",
                title=f"All Biomarkers ({', '.join(active_markers)}) MAF Across Geographic Groups",
                color_discrete_sequence=marker_colors,
            )
            fig_geo_all.update_layout(
                height=280, margin=dict(l=20, r=20, t=35, b=20)
            )
            st.plotly_chart(fig_geo_all, use_container_width=True)

    with c_geo2:
        geo_pheno_data = []
        for r in regions:
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
        geo_pm_df = pd.DataFrame(
            geo_pheno_data,
            columns=["Geographic Group", "Phenotype", "Percentage"],
        )
        if not geo_pm_df.empty:
            fig_geo_pm = px.bar(
                geo_pm_df,
                x="Geographic Group",
                y="Percentage",
                color="Phenotype",
                barmode="group",
                title=f"Phenotype Distribution Across Geographic Groups ({active_gene})",
                color_discrete_sequence=px.colors.qualitative.Safe,
            )
            fig_geo_pm.update_layout(
                height=280, margin=dict(l=20, r=20, t=35, b=20)
            )
            st.plotly_chart(fig_geo_pm, use_container_width=True)

    st.markdown("#### Comprehensive Population Statistical Summary Matrix")

    def render_marker_geo_table(snp_key: str) -> None:
        snp_conf_sel = active_snp_conf.get(snp_key, {})
        ov_snp = ov.get("snps", {}).get(snp_key, {})
        al_counts = ov_snp.get("allele_counts", {})
        gt_counts = ov_snp.get("genotype_counts", {})

        ref_al = snp_conf_sel.get("ref_allele")
        var_al = snp_conf_sel.get("var_allele")
        if not ref_al or not var_al:
            sorted_als = sorted(
                al_counts.items(), key=lambda x: x[1], reverse=True
            )
            ref_al = ref_al or (sorted_als[0][0] if sorted_als else "Ref")
            var_al = var_al or (
                sorted_als[1][0] if len(sorted_als) > 1 else "Var"
            )

        wt_g = snp_conf_sel.get("wildtype_genotype") or (ref_al + ref_al)
        var_g = snp_conf_sel.get("hom_var_genotype") or (var_al + var_al)
        het_g = snp_conf_sel.get("het_genotype")
        if not het_g:
            if f"{ref_al}{var_al}" in gt_counts:
                het_g = f"{ref_al}{var_al}"
            elif f"{var_al}{ref_al}" in gt_counts:
                het_g = f"{var_al}{ref_al}"
            else:
                het_g = f"{ref_al}{var_al}"

        if ref_al == var_al:
            mono_genotype = f"{ref_al}{ref_al}"
            rows_cfg = [
                (
                    f"{ref_al} Count (Single Observed Allele)",
                    lambda s: s.get(snp_key, {})
                    .get("allele_counts", {})
                    .get(ref_al, 0),
                ),
                (
                    "Total Allele Count (2N)",
                    lambda s: sum(
                        s.get(snp_key, {}).get("allele_counts", {}).values()
                    ),
                ),
                (
                    f"{mono_genotype} Count",
                    lambda s: s.get(snp_key, {})
                    .get("genotype_counts", {})
                    .get(mono_genotype, 0),
                ),
                (
                    "Total Sample N",
                    lambda s: sum(
                        s.get(snp_key, {}).get("genotype_counts", {}).values()
                    ),
                ),
                (
                    f"Allele Frequency F({ref_al})",
                    lambda s: s.get(snp_key, {})
                    .get("allele_freqs", {})
                    .get(ref_al, 0),
                ),
                (
                    "Chi-Square (HWE not informative)",
                    lambda s: s.get(snp_key, {})
                    .get("hwe", {})
                    .get("chi2_stat", 0),
                ),
                (
                    "P VALUE (HWE not informative)",
                    lambda s: s.get(snp_key, {})
                    .get("hwe", {})
                    .get("p_value", 1.0),
                ),
            ]
        else:
            rows_cfg = [
                (
                    f"{var_al} Count (Variant Allele)",
                    lambda s: s.get(snp_key, {})
                    .get("allele_counts", {})
                    .get(var_al, 0),
                ),
                (
                    f"{ref_al} Count (Reference Allele)",
                    lambda s: s.get(snp_key, {})
                    .get("allele_counts", {})
                    .get(ref_al, 0),
                ),
                (
                    "Total Allele Count (2N)",
                    lambda s: sum(
                        s.get(snp_key, {}).get("allele_counts", {}).values()
                    ),
                ),
                (
                    f"{het_g} COUNT (Heterozygote)",
                    lambda s: s.get(snp_key, {})
                    .get("genotype_counts", {})
                    .get(het_g, 0)
                    if het_g
                    else 0,
                ),
                (
                    f"{var_g} COUNT (Homozygous Variant)",
                    lambda s: s.get(snp_key, {})
                    .get("genotype_counts", {})
                    .get(var_g, 0)
                    if var_g
                    else 0,
                ),
                (
                    f"{wt_g} COUNT (Homozygous Wildtype)",
                    lambda s: s.get(snp_key, {})
                    .get("genotype_counts", {})
                    .get(wt_g, 0),
                ),
                (
                    "Total Sample N",
                    lambda s: sum(
                        s.get(snp_key, {}).get("genotype_counts", {}).values()
                    ),
                ),
                (
                    f"Allele Frequency F({var_al})",
                    lambda s: s.get(snp_key, {})
                    .get("allele_freqs", {})
                    .get(var_al, 0),
                ),
                (
                    f"Allele Frequency F({ref_al})",
                    lambda s: s.get(snp_key, {})
                    .get("allele_freqs", {})
                    .get(ref_al, 0),
                ),
                (
                    "Chi-Square (χ²)",
                    lambda s: s.get(snp_key, {})
                    .get("hwe", {})
                    .get("chi2_stat", 0),
                ),
                (
                    "P VALUE",
                    lambda s: s.get(snp_key, {})
                    .get("hwe", {})
                    .get("p_value", 1.0),
                ),
            ]

        m_matrix = []
        for lbl, fn in rows_cfg:
            row_item = {f"Statistical Parameter ({snp_key})": lbl}
            for r in regions:
                row_item[r] = fn(reg.get(r, {}).get("snps", {}))
            m_matrix.append(row_item)
        st.dataframe(
            pd.DataFrame(m_matrix), use_container_width=True, hide_index=True
        )

    marker_tab_labels = []
    for m in active_markers:
        rs_id = active_snp_conf.get(m, {}).get("rs_id")
        label = f"🧬 {m}" + (f" ({rs_id})" if rs_id else "")
        marker_tab_labels.append(label)

    if marker_tab_labels:
        marker_tabs = st.tabs(marker_tab_labels)
        for m_idx, m_name in enumerate(active_markers):
            with marker_tabs[m_idx]:
                render_marker_geo_table(m_name)
    else:
        st.info("No active markers to display.")

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(
        f"#### 🌐 Geographic Region Metabolizer Phenotype Distribution Matrix ({active_gene})"
    )
    st.caption(
        f"Observed {active_gene} phenotype categories across the geographic groups present in this dataset."
    )

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

    st.dataframe(
        pd.DataFrame(geo_pheno_matrix),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown(
        f"#### Diplotype Distribution by Geographic Group ({active_gene})"
    )
    geo_diplotype_categories = list(
        dict.fromkeys(
            name
            for group in reg.values()
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
    geo_diplotype_matrix = []
    for r_name in regions:
        r_data = reg.get(r_name, {})
        dip_data = (
            r_data.get("pgx_summary")
            or r_data.get(f"{active_gene.lower()}_summary")
            or {}
        ).get("diplotypes", {})
        row = {
            "Geographic Group": r_name,
            "Total Sample N": r_data.get("sample_count", 0),
        }
        for category in geo_diplotype_categories:
            metric = dip_data.get(category, {})
            row[category] = (
                f"{metric.get('count', 0)} ({metric.get('percentage', 0.0):.1f}%)"
            )
        geo_diplotype_matrix.append(row)
    st.dataframe(
        pd.DataFrame(geo_diplotype_matrix),
        use_container_width=True,
        hide_index=True,
    )
