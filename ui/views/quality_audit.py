"""Data Quality Control and Audit analytical view component."""

from typing import Any, Dict
import pandas as pd
import plotly.express as px
import streamlit as st


def render_quality_audit(qc_report: Dict[str, Any]) -> None:
    """Render Data Quality Control (QC) & Missing Values Audit view."""
    st.markdown("### 🛡️ Data Quality Control (QC) & Missing Values Audit")
    tot_up = qc_report.get("total_samples", 0)
    tot_analyzed = qc_report.get("total_analyzed_samples", tot_up - qc_report.get("discarded_null_samples", 0))
    discarded_cnt = qc_report.get("discarded_null_samples", qc_report.get("discarded_samples_count", 0))

    q1, q2, q3, q4, q5 = st.columns(5)
    q1.metric("Total Cohort Uploaded", f"{tot_up:,}")
    q2.metric("Clean Analyzed", f"{tot_analyzed:,}", f"-{discarded_cnt} discarded" if discarded_cnt else None, delta_color="inverse")
    q3.metric("Discarded (Any Null)", f"{discarded_cnt:,}")
    q4.metric("Duplicate Sample IDs", qc_report.get("duplicate_sample_count", 0))
    q5.metric("Invalid Genotype Calls", qc_report.get("invalid_genotypes_count", 0))

    st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 1. 4-FIELD DATA COMPLETENESS AUDIT TABLE
    # -------------------------------------------------------------------------
    st.markdown("#### 📋 4-Field Completeness & Missing Value Audit")
    st.caption("Audit across Sample ID, Patient Name, State, and Genetic Marker calls. Rows with missing values in any of these fields are excluded from calculations.")

    completeness_data = qc_report.get("completeness_audit", [])
    if completeness_data:
        comp_df = pd.DataFrame([
            {
                "Field Name": item["field"],
                "Field Role": item["description"],
                "Uploaded (N)": f"{item['total']:,}",
                "Missing (N)": f"{item['missing']:,}",
                "Valid (N)": f"{item['valid']:,}",
                "Completion Rate": f"{item['completion_pct']:.1f}%",
                "Audit Status": "✅ PASS" if item["missing"] == 0 else f"⚠️ FLAGGED ({item['missing']} missing)"
            }
            for item in completeness_data
        ])
        st.dataframe(comp_df, use_container_width=True, hide_index=True)

    # -------------------------------------------------------------------------
    # 2. DISCARDED / EXCLUDED SAMPLES LOG & DOWNLOAD
    # -------------------------------------------------------------------------
    discarded_log = qc_report.get("discarded_samples_log", [])
    if discarded_log:
        st.markdown(f"#### 🚫 Discarded / Excluded Samples ({len(discarded_log):,} Rows)")
        st.warning(
            f"⚠️ **Strict Complete-Record Rule Applied:** {len(discarded_log):,} row(s) were excluded from all statistical, HWE, and CPIC phenotype calculations "
            f"because one or more required columns (Sample ID, Name, State, or Marker Genotypes) were blank or uncalled."
        )

        disc_rows = []
        for item in discarded_log:
            marker_str = ", ".join([f"{k}: {v}" for k, v in item.get("marker_calls", {}).items()]) if item.get("marker_calls") else "None"
            disc_rows.append({
                "Sample ID": item.get("sample_id", "N/A"),
                "Patient Name": item.get("name", "N/A"),
                "State": item.get("state", "Missing"),
                "Reason(s) Discarded": item.get("reasons", "Missing required field"),
                "Observed Markers": marker_str
            })

        disc_df = pd.DataFrame(disc_rows)
        c_disc1, c_disc2 = st.columns([1.2, 0.8])
        with c_disc1:
            st.dataframe(disc_df, use_container_width=True, hide_index=True)
        with c_disc2:
            st.markdown("<div style='margin-top: 4px;'></div>", unsafe_allow_html=True)
            csv_disc = disc_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                "⬇️ Download Discarded Records (.csv)",
                data=csv_disc,
                file_name="Discarded_Incomplete_Samples.csv",
                mime="text/csv",
                use_container_width=True,
                key="btn_download_disc_csv"
            )
            st.info(
                f"ℹ️ **Audit Summary:** {len(disc_df):,} sample(s) had missing values. "
                "Download the CSV to inspect or correct these patient records in your spreadsheet."
            )

    st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)

    if qc_report.get("invalid_genotypes_count", 0) > 0:
        st.error(
            f"⚠️ **Warning: {qc_report['invalid_genotypes_count']} Invalid Genotype Call(s) Detected!** Invalid genotype entries were flagged and quarantined."
        )
        inv_log_rows = []
        for snp, samples in qc_report.get("invalid_genotypes_log", {}).items():
            for item in samples:
                inv_log_rows.append(
                    {
                        "SNP Variant": snp,
                        "Sample ID": item["sample_id"],
                        "Raw Invalid Genotype Value": item["raw_value"],
                        "Audit Action": "Quarantined & Excluded",
                    }
                )
        st.table(pd.DataFrame(inv_log_rows))

    c_qc_tbl, c_qc_fig = st.columns([1.1, 0.9])
    with c_qc_tbl:
        st.markdown("#### Per-SNP Quality Control & Missing Values Table")
        qc_rows = []
        for snp, sdata in qc_report["snp_qc_stats"].items():
            qc_rows.append(
                {
                    "SNP Variant": snp,
                    "Valid Samples": sdata["valid_samples"],
                    "Missing Genotypes Count": sdata["missing_samples"],
                    "Invalid Genotypes Count": sdata.get("invalid_samples", 0),
                    "Missing Data Rate (%)": f"{sdata['missing_pct']:.1f}%",
                    "Quality Audit Status": "PASS"
                    if sdata.get("invalid_samples", 0) == 0
                    else "WARNING",
                }
            )
        st.table(pd.DataFrame(qc_rows))

    with c_qc_fig:
        qc_plot_rows = [
            {
                "SNP": str(snp),
                "Valid": int(sdata.get("valid_samples", 0)),
                "Missing": int(sdata.get("missing_samples", 0)),
                "Invalid": int(sdata.get("invalid_samples", 0)),
            }
            for snp, sdata in qc_report.get("snp_qc_stats", {}).items()
        ]
        qc_fig_df = pd.DataFrame(
            qc_plot_rows, columns=["SNP", "Valid", "Missing", "Invalid"]
        )
        if not qc_fig_df.empty and len(qc_plot_rows) > 0:
            fig_qc = px.bar(
                qc_fig_df,
                x="SNP",
                y=["Valid", "Missing", "Invalid"],
                title="Sample Quality Audit Breakdown per SNP",
                barmode="stack",
                color_discrete_sequence=["#16A34A", "#EAB308", "#DC2626"],
            )
            fig_qc.update_layout(
                height=260, margin=dict(l=20, r=20, t=35, b=20)
            )
            st.plotly_chart(fig_qc, use_container_width=True)
        else:
            st.info("No QC chart data available.")
