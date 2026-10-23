"""
Executive Overview Analysis View.
Renders KPI metrics, cohort QC summary, phenotype overview, and regional sample distribution.
"""

import pandas as pd
import plotly.express as px
import streamlit as st

def render_executive_overview(full_results: dict, qc_report: dict, active_gene: str, active_markers: list = None):
    """Renders executive overview KPIs, phenotype summaries, and cohort quality breakdown."""
    ov = full_results.get('overall', {})
    tot_samples = ov.get('sample_count', qc_report.get('total_samples', 0))
    fully_valid = qc_report.get('fully_valid_samples', tot_samples)
    dup_count = qc_report.get('duplicate_sample_count', 0)
    dup_ids = qc_report.get('duplicate_sample_ids', [])
    invalid_count = qc_report.get('invalid_genotypes_count', 0)
    invalid_log = qc_report.get('invalid_genotypes_log', {})

    missing_count = qc_report.get('samples_with_missing_count', 0)
    missing_pct = qc_report.get('overall_missing_pct', 0.0)

    discarded_null = full_results.get('discarded_null_samples', qc_report.get('discarded_null_samples', 0))
    total_input = full_results.get('total_input_samples', qc_report.get('total_samples', tot_samples))

    st.markdown('<div class="metric-row-anchor"></div>', unsafe_allow_html=True)
    k1, k2, k3, k4, k5 = st.columns([1, 1, 1.35, 0.9, 1.05])
    if discarded_null > 0:
        k1.metric("ANALYZED SAMPLES", f"{tot_samples:,}", f"-{discarded_null} null discarded", delta_color="inverse")
    else:
        k1.metric("ANALYZED SAMPLES", f"{tot_samples:,}")
    k2.metric("FULLY VALID", f"{fully_valid:,}")
    k3.metric("MISSING CALLS", f"{missing_pct:.1f}%", f"{missing_count:,} samples", delta_color="off")
    k4.metric("DUPLICATE IDS", dup_count)
    k5.metric("INVALID GENOTYPES", invalid_count)

    if dup_ids:
        st.markdown(f"""
        <div style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 6px; padding: 0.8rem 1rem; margin-top: 0.5rem; margin-bottom: 1rem; font-family: monospace; font-size: 0.85rem; color: #334155;">
            <strong>Duplicate Sample IDs:</strong> {', '.join([str(x) for x in dup_ids[:10]])}{'...' if len(dup_ids) > 10 else ''}
        </div>
        """, unsafe_allow_html=True)
        
    if invalid_count > 0:
        inv_details = []
        for snp, samples in invalid_log.items():
            s_list = [f"{item['sample_id']} ('{item['raw_value']}')" for item in samples[:5]]
            inv_details.append(f"{snp}: {', '.join(s_list)}")
        st.warning(f"⚠️ **Warning: {invalid_count} Invalid Genotype Call(s) Detected!** Details: {'; '.join(inv_details)}")

    st.markdown("<br>", unsafe_allow_html=True)
    c_ov1, c_ov2 = st.columns(2)
    with c_ov1:
        st.markdown(f"### 🏆 {active_gene} Phenotype Summary")
        pgx = ov.get('pgx_summary') or ov.get('cyp2c9_summary') or ov.get('cyp2c19_summary', {})
        phe_dict = {k: v for k, v in pgx.get('phenotypes', {}).items() if k != 'Not configured'}
        phe_df = pd.DataFrame([{'Phenotype Category': k, 'Count': v['count'], 'Percentage': f"{v['percentage']:.1f}%"} for k, v in phe_dict.items() if v['count'] > 0])
        if not phe_df.empty:
            st.table(phe_df)
            fig_ov_phe = px.pie(
                names=[k for k, v in phe_dict.items() if v['count'] > 0],
                values=[v['count'] for k, v in phe_dict.items() if v['count'] > 0],
                hole=0.45,
                title="Overall Phenotype Distribution",
                color_discrete_sequence=px.colors.qualitative.Safe
            )
            fig_ov_phe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_ov_phe, use_container_width=True)
        else:
            st.info("Phenotype classification is unavailable until validated translation rules are configured for this gene and marker panel.")

    with c_ov2:
        st.markdown("### 🌐 Regional Sample Distribution")
        reg = full_results.get('regional', {})
        reg_df = pd.DataFrame([{'Region': r, 'Sample Count (N)': r_data.get('sample_count', 0), '% of Cohort': f"{(r_data.get('sample_count', 0)/tot_samples*100):.1f}%" if tot_samples>0 else "0%"} for r, r_data in reg.items()])
        st.table(reg_df)

        reg_chart_df = pd.DataFrame([{'Region': r, 'Samples (N)': r_data.get('sample_count', 0)} for r, r_data in reg.items()])
        fig_ov_reg = px.bar(
            reg_chart_df, x='Region', y='Samples (N)',
            title="Regional Population Group Sizes",
            color='Region',
            color_discrete_sequence=px.colors.qualitative.Set2
        )
        fig_ov_reg.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
        st.plotly_chart(fig_ov_reg, use_container_width=True)
