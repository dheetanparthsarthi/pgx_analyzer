"""
Genotype & Allele Analysis View.
Displays observed genotype counts, frequencies, zygosity classifications, and 2N allele breakdown charts.
"""

import pandas as pd
import plotly.express as px
import streamlit as st
from config import SNP_CONFIG, CYP2C9_SNP_CONFIG

def render_genotype_allele(
    active_gene: str = None,
    active_markers: list = None,
    active_snp_conf: dict = None,
    full_results: dict = None,
    render_subcohort_pills = None,
    **kwargs
):
    """Renders genotype frequencies, zygosity types, and allele charts for the selected cohort."""
    if render_subcohort_pills is not None:
        sel_key, sel_cohort = render_subcohort_pills("tab_gen")
    else:
        sel_key = kwargs.get('sel_key', 'Overall')
        sel_cohort = kwargs.get('sel_cohort', (full_results.get('overall', {}) if full_results else {}))

    snps_data = sel_cohort.get('snps', {})
    if active_snp_conf is None:
        active_snp_conf = {**SNP_CONFIG, **CYP2C9_SNP_CONFIG}


    def get_zygosity(gt_name, snp_name):
        snp_conf = active_snp_conf.get(snp_name, {})
        wt = snp_conf.get('wildtype_genotype', '')
        het = snp_conf.get('het_genotype', '')
        var = snp_conf.get('hom_var_genotype', '')
        if gt_name == wt:
            return "Homozygous Wildtype (Homogenous)"
        elif gt_name in [het, het[::-1]]:
            return "Heterozygote (Heterogenous)"
        elif gt_name == var:
            return "Homozygous Variant (Homogenous)"
        return "Genotype Call"

    def get_allele_type(al_name, snp_name):
        snp_conf = active_snp_conf.get(snp_name, {})
        ref = snp_conf.get('ref_allele', '')
        var = snp_conf.get('var_allele', '')
        if al_name == ref:
            return "Reference Allele"
        elif al_name == var:
            return "Variant Allele"
        return "Allele Call"

    snp_display_order = []
    for snp_code in snps_data.keys():
        m_conf = active_snp_conf.get(snp_code, {})
        rs_id = m_conf.get('rs_id') or m_conf.get('rsid')
        title = f"{snp_code}" + (f" ({rs_id})" if rs_id and rs_id not in snp_code else "")
        snp_display_order.append((snp_code, title))

    for snp_code, snp_title in snp_display_order:
        if snp_code in snps_data:
            s_res = snps_data[snp_code]
            v_cnt = s_res.get('valid_samples', 0)
            sub_total = sel_cohort.get('sample_count', v_cnt)
            sub_missing = max(0, sub_total - v_cnt)
            sub_m_pct = (sub_missing / sub_total * 100.0) if sub_total > 0 else 0.0
            
            st.markdown(f"#### {snp_title}")
            st.caption(f"Cohort Total (N): **{sub_total:,}** · Valid Genotyped: **{v_cnt:,}** · Missing Calls: **{sub_missing:,} ({sub_m_pct:.1f}%)**")
            
            c_tbl, c_fig = st.columns([1.1, 0.9])
            with c_tbl:
                g_c = s_res.get('genotype_counts', {})
                g_p = s_res.get('genotype_pcts', {})
                g_f = s_res.get('genotype_freqs', {})
                gt_rows = [{
                    'GENOTYPE': k, 
                    'ZYGOSITY': get_zygosity(k, snp_code),
                    'COUNT': v, 
                    'FREQUENCY': g_f.get(k, 0),
                    'PERCENTAGE': f"{g_p.get(k, 0):.1f}%"
                } for k, v in g_c.items() if v > 0]
                st.table(pd.DataFrame(gt_rows))
                
                a_c = s_res.get('allele_counts', {})
                a_p = s_res.get('allele_pcts', {})
                al_rows = [{
                    'ALLELE': k, 
                    'TYPE': get_allele_type(k, snp_code),
                    'COUNT': v, 
                    'FREQUENCY': f"{a_p.get(k, 0):.1f}%"
                } for k, v in a_c.items()]
                st.table(pd.DataFrame(al_rows))
                
            with c_fig:
                fig_rows = [{'Genotype': str(k), 'Count': int(v)} for k, v in g_c.items() if v > 0]
                fig_df = pd.DataFrame(fig_rows, columns=['Genotype', 'Count'])
                if not fig_df.empty and len(fig_rows) > 0:
                    fig_bar = px.bar(
                        fig_df, x='Genotype', y='Count', 
                        title=f"{snp_code} Genotype distribution ({sel_key})", 
                        color='Genotype',
                        color_discrete_sequence=['#15803D', '#1D4ED8', '#B91C1C']
                    )
                    fig_bar.update_layout(height=220, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_bar, use_container_width=True)
                else:
                    st.caption(f"No genotype calls to display for {snp_code}.")

                al_rows_plot = [{'Allele': str(k), 'Frequency (%)': float(a_p.get(k, 0))} for k, v in a_c.items()]
                al_df = pd.DataFrame(al_rows_plot, columns=['Allele', 'Frequency (%)'])
                if not al_df.empty and len(al_rows_plot) > 0:
                    fig_al = px.bar(
                        al_df, x='Allele', y='Frequency (%)',
                        title=f"{snp_code} Allele Frequencies ({sel_key})",
                        color='Allele',
                        color_discrete_sequence=['#0D9488', '#D97706']
                    )
                    fig_al.update_layout(height=200, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_al, use_container_width=True)
