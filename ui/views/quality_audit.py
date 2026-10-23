"""Data Quality Control and Audit analytical view component."""

from typing import Any, Dict
import pandas as pd
import plotly.express as px
import streamlit as st


def render_quality_audit(qc_report: Dict[str, Any]) -> None:
    """Render Data Quality Control (QC) & Missing Values Audit view."""
    st.markdown("### 🛡️ Data Quality Control (QC) & Missing Values Audit")
    q1, q2, q3, q4, q5 = st.columns(5)
    q1.metric("Total Cohort Uploaded", qc_report["total_samples"])
    q2.metric("Fully Valid Samples", qc_report["fully_valid_samples"])
    q3.metric("Duplicate Sample IDs", qc_report["duplicate_sample_count"])
    q4.metric(
        "Missing Genotypes",
        f"{qc_report.get('samples_with_missing_count', 0)} ({qc_report.get('overall_missing_pct', 0.0):.1f}%)",
    )
    q5.metric("Invalid Genotype Calls", qc_report["invalid_genotypes_count"])

    if qc_report["invalid_genotypes_count"] > 0:
        st.error(
            f"⚠️ **Warning: {qc_report['invalid_genotypes_count']} Invalid Genotype Call(s) Detected!** Invalid genotype entries were flagged and quarantined."
        )
        inv_log_rows = []
        for snp, samples in qc_report["invalid_genotypes_log"].items():
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
