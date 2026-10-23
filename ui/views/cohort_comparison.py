"""Geographic and State-Wise Cohort Comparison & Pharmacogenomic Prediction analytical view component."""

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
try:
    import plotly.express as px
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except Exception:
    HAS_PLOTLY = False
    px = None
    go = None
from scipy import stats
import streamlit as st


def _get_cohort_data(
    cohort_name: str,
    cohort_type: str,
    full_results: Dict[str, Any],
    processed_df: pd.DataFrame,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """Retrieve statistics and filtered sample rows for a given cohort."""
    if cohort_type == "Region":
        cohort_stats = full_results.get("regional", {}).get(cohort_name, {})
        subset_df = (
            processed_df[
                processed_df["Region"].astype(str) == str(cohort_name)
            ]
            if "Region" in processed_df.columns
            else pd.DataFrame()
        )
    elif cohort_type == "State":
        cohort_stats = full_results.get("state_wise", {}).get(cohort_name, {})
        state_col = (
            "State_Clean"
            if "State_Clean" in processed_df.columns
            else ("State" if "State" in processed_df.columns else None)
        )
        if cohort_name == "Unclassified / Other State":
            subset_df = (
                processed_df[
                    processed_df[state_col].isna()
                    | processed_df[state_col]
                    .astype(str)
                    .str.strip()
                    .isin(["", "nan", "null", "none"])
                ]
                if state_col and state_col in processed_df.columns
                else pd.DataFrame()
            )
        else:
            subset_df = (
                processed_df[
                    processed_df[state_col].astype(str) == str(cohort_name)
                ]
                if state_col and state_col in processed_df.columns
                else pd.DataFrame()
            )
    else:  # Overall
        cohort_stats = full_results.get("overall", {})
        subset_df = processed_df.copy()

    return cohort_stats, subset_df


def _generate_prediction_insights(
    active_gene: str,
    name_a: str,
    stats_a: Dict[str, Any],
    name_b: str,
    stats_b: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """Generate clinical drug response risk predictions and outcome forecasts based on phenotype proportions."""
    insights = []

    pgx_a = (
        stats_a.get("pgx_summary")
        or stats_a.get(f"{active_gene.lower()}_summary")
        or stats_a.get("cyp2c19_summary")
        or stats_a.get("cyp2c9_summary")
        or {}
    )
    pgx_b = (
        stats_b.get("pgx_summary")
        or stats_b.get(f"{active_gene.lower()}_summary")
        or stats_b.get("cyp2c19_summary")
        or stats_b.get("cyp2c9_summary")
        or {}
    )

    phenos_a = pgx_a.get("phenotypes", {})
    phenos_b = pgx_b.get("phenotypes", {})

    pm_a = float(phenos_a.get("Poor Metabolizer", {}).get("percentage", 0.0))
    pm_b = float(phenos_b.get("Poor Metabolizer", {}).get("percentage", 0.0))

    im_a = float(
        phenos_a.get("Intermediate Metabolizer", {}).get("percentage", 0.0)
    )
    im_b = float(
        phenos_b.get("Intermediate Metabolizer", {}).get("percentage", 0.0)
    )

    rm_um_a = float(
        phenos_a.get("Rapid Metabolizer", {}).get("percentage", 0.0)
    ) + float(
        phenos_a.get("Ultrarapid Metabolizer", {}).get("percentage", 0.0)
    )
    rm_um_b = float(
        phenos_b.get("Rapid Metabolizer", {}).get("percentage", 0.0)
    ) + float(
        phenos_b.get("Ultrarapid Metabolizer", {}).get("percentage", 0.0)
    )

    gene_upper = active_gene.upper()

    if gene_upper == "CYP2C19":
        # Clopidogrel (Plavix) Antiplatelet Efficacy
        risk_a = pm_a + im_a
        risk_b = pm_b + im_b
        higher_risk_cohort = name_a if risk_a > risk_b else name_b
        diff_risk = abs(risk_a - risk_b)

        insights.append(
            {
                "drug": "🩸 Clopidogrel (Antiplatelet Efficacy & Stent Thrombosis)",
                "risk_title": f"High-Risk Non-Responders (PM + IM): {name_a} ({risk_a:.1f}%) vs {name_b} ({risk_b:.1f}%)",
                "divergence": f"{diff_risk:.1f}% absolute divergence",
                "prediction": (
                    f"Patients from **{higher_risk_cohort}** exhibit a higher rate of CYP2C19 loss-of-function (*2/*3), "
                    f"resulting in reduced formation of active clopidogrel metabolite. "
                    f"In **{higher_risk_cohort}**, approximately {max(risk_a, risk_b):.1f}% of individuals are predicted to have elevated risk of post-PCI ischemic events."
                ),
                "clinical_action": (
                    "**Recommendation:** Favor alternative P2Y12 inhibitors (Prasugrel or Ticagrelor) "
                    f"in high-risk patients from {higher_risk_cohort} unless contraindications apply (CPIC Level A)."
                ),
                "severity": "high" if max(risk_a, risk_b) > 40 else "medium",
            }
        )

        # Voriconazole & Escitalopram
        insights.append(
            {
                "drug": "💊 Voriconazole & SSRIs (Escitalopram / Citalopram)",
                "risk_title": f"Toxicity Risk (PM): {name_a} ({pm_a:.1f}%) vs {name_b} ({pm_b:.1f}%)",
                "divergence": f"{abs(pm_a - pm_b):.1f}% divergence in Poor Metabolizers",
                "prediction": (
                    f"Poor metabolizers experience marked supratherapeutic plasma concentrations. "
                    f"{name_a} has {pm_a:.1f}% PMs, whereas {name_b} has {pm_b:.1f}% PMs. "
                    f"High PM rates correlate with higher incidence of QT prolongation and voriconazole neurotoxicity."
                ),
                "clinical_action": (
                    "**Recommendation:** For Poor Metabolizers, reduce starting SSRI dose by 50% or choose alternative antifungal agents (e.g. Isavuconazole)."
                ),
                "severity": "high" if max(pm_a, pm_b) > 10 else "medium",
            }
        )

    elif gene_upper == "CYP2C9":
        # Warfarin Sensitivity
        sub_normal_a = pm_a + im_a
        sub_normal_b = pm_b + im_b
        insights.append(
            {
                "drug": "🩸 Warfarin (Oral Anticoagulation Sensitivity)",
                "risk_title": f"Dose Reduction Substrates (PM + IM): {name_a} ({sub_normal_a:.1f}%) vs {name_b} ({sub_normal_b:.1f}%)",
                "divergence": f"{abs(sub_normal_a - sub_normal_b):.1f}% divergence",
                "prediction": (
                    f"Allelic variants *2 and *3 significantly impair clearance of S-warfarin. "
                    f"The combined frequency of reduced metabolizers is {sub_normal_a:.1f}% in {name_a} and {sub_normal_b:.1f}% in {name_b}. "
                    f"Patients carrying these alleles have prolonged INR and 2-3x higher major bleeding risk on standard 5mg starting doses."
                ),
                "clinical_action": (
                    "**Recommendation:** Initiate with 20% to 50% dose reductions for Intermediate Metabolizers and up to 80% reduction (or DOAC alternatives) for Poor Metabolizers."
                ),
                "severity": "high" if max(sub_normal_a, sub_normal_b) > 20 else "medium",
            }
        )

        # Phenytoin & Celecoxib
        insights.append(
            {
                "drug": "🧠 Phenytoin & NSAIDs (Celecoxib)",
                "risk_title": f"Impaired Clearance Risk: {name_a} vs {name_b}",
                "divergence": f"PM rates: {pm_a:.1f}% vs {pm_b:.1f}%",
                "prediction": (
                    f"Phenytoin maintenance dose requirements are markedly lower in CYP2C9 intermediate and poor metabolizers. "
                    f"Toxicity risk is disproportionately elevated in cohorts with higher carriage of reduced-function alleles."
                ),
                "clinical_action": (
                    "**Recommendation:** Monitor serum trough concentrations closely and reduce initial maintenance dosing by 25-50% in predicted IM/PM patients."
                ),
                "severity": "medium",
            }
        )

    else:
        # Generic Target Gene Insights
        insights.append(
            {
                "drug": f"🧬 {active_gene} Substrate Forecast",
                "risk_title": f"Observed Cohort Variation: {name_a} vs {name_b}",
                "divergence": f"Comparative Phenotype / Allele Variance",
                "prediction": (
                    f"Comparing {name_a} and {name_b} reveals distinct population frequency structures for target gene {active_gene}. "
                    f"Sub-population genetic drift and regional allele differences should be taken into account when interpreting dosing recommendations."
                ),
                "clinical_action": (
                    "**Recommendation:** Refer to CPIC / PharmGKB guidelines for specific substrate adjustments corresponding to the observed diplotypes."
                ),
                "severity": "low",
            }
        )

    return insights


def render_cohort_comparison(
    active_gene: str,
    active_markers: List[str],
    full_results: Dict[str, Any],
    processed_df: pd.DataFrame,
) -> None:
    """Render the Geographic & State Comparison and Pharmacogenomic Prediction view."""
    st.markdown(f"### ⚖️ Geographic & State-Wise Cohort Comparison ({active_gene})")
    st.caption(
        "Directly compare genetic marker frequencies, Hardy-Weinberg equilibrium, and clinical metabolizer phenotypes "
        "between two selected geographical regions or individual states, complete with automated outcome predictions."
    )

    reg_dict = full_results.get("regional", {})
    state_dict = full_results.get("state_wise", {})

    available_regions = sorted(
        [r for r, d in reg_dict.items() if d.get("sample_count", 0) > 0]
    )
    available_states = sorted(
        [s for s, d in state_dict.items() if d.get("sample_count", 0) > 0]
    )

    if not available_regions and not available_states:
        st.info(
            "ℹ️ Cohort comparison requires mapped geographic regions or state columns in your dataset. "
            "Please ensure Region or State columns are mapped in Step 3."
        )
        return

    # Comparison Mode Selector
    modes = []
    if len(available_regions) >= 2:
        modes.append("🗺️ Region vs Region")
    if len(available_states) >= 2:
        modes.append("📍 State vs State")
    if available_regions and available_states:
        modes.append("🌐 Region vs State")
    modes.append("📊 Cohort vs Overall Population")

    c_mode, _ = st.columns([1.2, 0.8])
    with c_mode:
        selected_mode = st.radio(
            "Select Comparison Scope:",
            modes,
            horizontal=True,
            key="cohort_cmp_mode",
        )

    st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)

    # Pick Cohort A and Cohort B based on mode
    col_a, col_b = st.columns(2)

    if "Region vs Region" in selected_mode:
        with col_a:
            name_a = st.selectbox(
                "Select Cohort A (First Region):",
                available_regions,
                index=0,
                key="cmp_cohort_a_reg",
            )
            type_a = "Region"
        with col_b:
            b_default_idx = 1 if len(available_regions) > 1 else 0
            name_b = st.selectbox(
                "Select Cohort B (Second Region):",
                available_regions,
                index=b_default_idx,
                key="cmp_cohort_b_reg",
            )
            type_b = "Region"

    elif "State vs State" in selected_mode:
        with col_a:
            name_a = st.selectbox(
                "Select Cohort A (First State):",
                available_states,
                index=0,
                key="cmp_cohort_a_st",
            )
            type_a = "State"
        with col_b:
            b_default_idx = 1 if len(available_states) > 1 else 0
            name_b = st.selectbox(
                "Select Cohort B (Second State):",
                available_states,
                index=b_default_idx,
                key="cmp_cohort_b_st",
            )
            type_b = "State"

    elif "Region vs State" in selected_mode:
        with col_a:
            name_a = st.selectbox(
                "Select Cohort A (Region):",
                available_regions,
                index=0,
                key="cmp_cohort_a_reg_st",
            )
            type_a = "Region"
        with col_b:
            name_b = st.selectbox(
                "Select Cohort B (State):",
                available_states,
                index=0,
                key="cmp_cohort_b_reg_st",
            )
            type_b = "State"

    else:  # Cohort vs Overall
        all_options = [f"Region: {r}" for r in available_regions] + [
            f"State: {s}" for s in available_states
        ]
        with col_a:
            chosen_sub = st.selectbox(
                "Select Sub-Cohort:",
                all_options,
                index=0,
                key="cmp_cohort_a_ov",
            )
            if chosen_sub.startswith("Region: "):
                name_a = chosen_sub.replace("Region: ", "")
                type_a = "Region"
            else:
                name_a = chosen_sub.replace("State: ", "")
                type_a = "State"
        with col_b:
            st.selectbox(
                "Reference Cohort:",
                ["Total Overall Population"],
                disabled=True,
                key="cmp_cohort_b_ov",
            )
            name_b = "Total Overall Population"
            type_b = "Overall"

    # Fetch Cohort Data
    stats_a, df_a = _get_cohort_data(name_a, type_a, full_results, processed_df)
    stats_b, df_b = _get_cohort_data(name_b, type_b, full_results, processed_df)

    n_a = stats_a.get("sample_count", len(df_a))
    n_b = stats_b.get("sample_count", len(df_b))

    # Top KPI Comparison Header Cards
    st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)
    st.markdown('<div class="metric-row-anchor"></div>', unsafe_allow_html=True)
    k1, k2, k3, k4 = st.columns(4)

    card_style = "background: #FFFFFF; border: 1.5px solid #E2E8F0; border-radius: 8px; padding: 12px 14px; box-shadow: 0 1px 3px rgba(0,0,0,0.05);"
    lbl_style = "font-size: 11px; font-weight: 700; text-transform: uppercase; color: #64748B; letter-spacing: 0.05em;"
    val_style = "font-size: 16px; font-weight: 700; color: #0F172A; margin-top: 4px;"

    with k1:
        st.markdown(
            f"""
        <div class="kpi-card cohort-kpi-card" style="{card_style}">
            <div class="kpi-lbl" style="{lbl_style}">Cohort A ({type_a})</div>
            <div class="kpi-val" style="{val_style}"><span style="color: #2563EB;">{name_a}</span></div>
            <div style="font-size: 12px; color: #475569; margin-top: 2px;">{n_a:,} Patients</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    with k2:
        st.markdown(
            f"""
        <div class="kpi-card cohort-kpi-card" style="{card_style}">
            <div class="kpi-lbl" style="{lbl_style}">Cohort B ({type_b})</div>
            <div class="kpi-val" style="{val_style}"><span style="color: #059669;">{name_b}</span></div>
            <div style="font-size: 12px; color: #475569; margin-top: 2px;">{n_b:,} Patients</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    pgx_a = (
        stats_a.get("pgx_summary")
        or stats_a.get(f"{active_gene.lower()}_summary")
        or stats_a.get("cyp2c19_summary")
        or stats_a.get("cyp2c9_summary")
        or {}
    )
    pgx_b = (
        stats_b.get("pgx_summary")
        or stats_b.get(f"{active_gene.lower()}_summary")
        or stats_b.get("cyp2c19_summary")
        or stats_b.get("cyp2c9_summary")
        or {}
    )

    pm_pct_a = float(
        pgx_a.get("phenotypes", {})
        .get("Poor Metabolizer", {})
        .get("percentage", 0.0)
    )
    pm_pct_b = float(
        pgx_b.get("phenotypes", {})
        .get("Poor Metabolizer", {})
        .get("percentage", 0.0)
    )
    pm_count_a = int(
        pgx_a.get("phenotypes", {})
        .get("Poor Metabolizer", {})
        .get("count", 0)
    )
    pm_count_b = int(
        pgx_b.get("phenotypes", {})
        .get("Poor Metabolizer", {})
        .get("count", 0)
    )
    pm_diff = abs(pm_pct_a - pm_pct_b)
    pm_diff_color = "#64748B" if pm_diff == 0.0 else "#DC2626"

    with k3:
        st.markdown(
            f"""
        <div class="kpi-card cohort-kpi-card" style="{card_style}">
            <div class="kpi-lbl" style="{lbl_style}">Poor Metabolizers (PM)</div>
            <div class="kpi-val" style="{val_style}">{pm_pct_a:.1f}% vs {pm_pct_b:.1f}%</div>
            <div style="font-size: 11px; color: #64748B; margin-top: 2px;">({pm_count_a}/{n_a} vs {pm_count_b}/{n_b} samples)</div>
            <div style="font-size: 12px; color: {pm_diff_color}; font-weight: 600; margin-top: 2px;">Δ {pm_diff:.1f}% difference</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    im_pct_a = float(
        pgx_a.get("phenotypes", {})
        .get("Intermediate Metabolizer", {})
        .get("percentage", 0.0)
    )
    im_pct_b = float(
        pgx_b.get("phenotypes", {})
        .get("Intermediate Metabolizer", {})
        .get("percentage", 0.0)
    )
    im_count_a = int(
        pgx_a.get("phenotypes", {})
        .get("Intermediate Metabolizer", {})
        .get("count", 0)
    )
    im_count_b = int(
        pgx_b.get("phenotypes", {})
        .get("Intermediate Metabolizer", {})
        .get("count", 0)
    )
    im_diff = abs(im_pct_a - im_pct_b)
    im_diff_color = "#64748B" if im_diff == 0.0 else "#D97706"

    with k4:
        st.markdown(
            f"""
        <div class="kpi-card cohort-kpi-card" style="{card_style}">
            <div class="kpi-lbl" style="{lbl_style}">Intermediate Metabolizers (IM)</div>
            <div class="kpi-val" style="{val_style}">{im_pct_a:.1f}% vs {im_pct_b:.1f}%</div>
            <div style="font-size: 11px; color: #64748B; margin-top: 2px;">({im_count_a}/{n_a} vs {im_count_b}/{n_b} samples)</div>
            <div style="font-size: 12px; color: {im_diff_color}; font-weight: 600; margin-top: 2px;">Δ {im_diff:.1f}% difference</div>
        </div>
        """,
            unsafe_allow_html=True,
        )

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 1. BIOMARKER ALLELE FREQUENCY COMPARISON & CHI-SQUARE DIVERGENCE TEST
    # -------------------------------------------------------------------------
    st.markdown(f"#### 🧬 Marker Allele Frequency Divergence ({name_a} vs {name_b})")
    st.caption(
        "Compares Minor Allele Frequencies (MAF) and performs a 2x2 contingency Chi-Square test of allele distribution independence."
    )

    marker_rows = []
    chart_marker_rows = []

    snps_a = stats_a.get("snps", {})
    snps_b = stats_b.get("snps", {})

    for marker in active_markers:
        m_a = snps_a.get(marker, {})
        m_b = snps_b.get(marker, {})

        al_counts_a = m_a.get("allele_counts", {})
        al_counts_b = m_b.get("allele_counts", {})

        var_al = None
        for al in al_counts_a:
            if al_counts_a.get(al, 0) < sum(al_counts_a.values()) / 2:
                var_al = al
                break
        if not var_al and al_counts_a:
            sorted_a = sorted(al_counts_a.items(), key=lambda x: x[1])
            var_al = sorted_a[0][0]
        if not var_al:
            var_al = "Var"

        maf_a = float(m_a.get("allele_freqs", {}).get(var_al, m_a.get("maf", 0.0)))
        maf_b = float(m_b.get("allele_freqs", {}).get(var_al, m_b.get("maf", 0.0)))

        delta_maf = maf_a - maf_b

        # 2x2 contingency test on alleles if counts available
        p_val_divergence = 1.0
        ref_count_a = sum(al_counts_a.values()) - al_counts_a.get(var_al, 0)
        var_count_a = al_counts_a.get(var_al, 0)
        ref_count_b = sum(al_counts_b.values()) - al_counts_b.get(var_al, 0)
        var_count_b = al_counts_b.get(var_al, 0)

        if (ref_count_a + var_count_a > 0) and (ref_count_b + var_count_b > 0):
            try:
                table_2x2 = [[var_count_a, ref_count_a], [var_count_b, ref_count_b]]
                _, p_val_divergence, _, _ = stats.chi2_contingency(table_2x2)
            except Exception:
                p_val_divergence = 1.0

        sig_label = (
            "⚠️ Significant (p < 0.05)"
            if p_val_divergence < 0.05
            else "Not Significant"
        )

        marker_rows.append(
            {
                "Biomarker": marker,
                "Allele Analyzed": var_al,
                f"MAF in {name_a}": f"{maf_a:.4f}",
                f"MAF in {name_b}": f"{maf_b:.4f}",
                "Frequency Delta (Δ)": f"{delta_maf:+.4f}",
                "Divergence P-Value": f"{p_val_divergence:.4f}"
                if p_val_divergence < 0.0001
                else f"{p_val_divergence:.4f}",
                "Significance": sig_label,
            }
        )

        chart_marker_rows.append(
            {"Biomarker": marker, "Cohort": name_a, "Frequency": maf_a}
        )
        chart_marker_rows.append(
            {"Biomarker": marker, "Cohort": name_b, "Frequency": maf_b}
        )

    c_tbl_m, c_fig_m = st.columns([1.2, 0.8])
    with c_tbl_m:
        st.dataframe(pd.DataFrame(marker_rows), use_container_width=True, hide_index=True)

    with c_fig_m:
        if chart_marker_rows:
            if HAS_PLOTLY and px is not None:
                fig_cmp_maf = px.bar(
                    pd.DataFrame(chart_marker_rows),
                    x="Biomarker",
                    y="Frequency",
                    color="Cohort",
                    barmode="group",
                    title="Allele Frequency Comparison",
                    color_discrete_sequence=["#2563EB", "#059669"],
                )
                fig_cmp_maf.update_layout(
                    height=260, margin=dict(l=10, r=10, t=32, b=10)
                )
                st.plotly_chart(fig_cmp_maf, use_container_width=True)
            else:
                m_df = pd.DataFrame(chart_marker_rows)
                st.bar_chart(m_df, x="Biomarker", y="Frequency", color="Cohort")

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 2. METABOLIZER PHENOTYPE DISTRIBUTION COMPARISON
    # -------------------------------------------------------------------------
    st.markdown(f"#### 💊 Metabolizer Phenotype Distribution ({name_a} vs {name_b})")
    st.caption("Distribution of clinical metabolizer categories across both cohorts.")

    phenos_a = pgx_a.get("phenotypes", {})
    phenos_b = pgx_b.get("phenotypes", {})

    all_pheno_keys = list(
        dict.fromkeys(
            list(phenos_a.keys()) + list(phenos_b.keys())
        )
    )
    all_pheno_keys = [k for k in all_pheno_keys if k != "Not configured"]

    pheno_rows = []
    pheno_chart_data = []

    for phe in all_pheno_keys:
        cnt_a = int(phenos_a.get(phe, {}).get("count", 0))
        pct_a = float(phenos_a.get(phe, {}).get("percentage", 0.0))

        cnt_b = int(phenos_b.get(phe, {}).get("count", 0))
        pct_b = float(phenos_b.get(phe, {}).get("percentage", 0.0))

        pheno_rows.append(
            {
                "Phenotype Category": phe,
                f"{name_a} (N)": cnt_a,
                f"{name_a} (%)": f"{pct_a:.1f}%",
                f"{name_b} (N)": cnt_b,
                f"{name_b} (%)": f"{pct_b:.1f}%",
                "Percentage Shift (Δ)": f"{(pct_a - pct_b):+.1f}%",
            }
        )

        pheno_chart_data.append(
            {"Phenotype": phe.split(" (")[0], "Cohort": name_a, "Percentage": pct_a}
        )
        pheno_chart_data.append(
            {"Phenotype": phe.split(" (")[0], "Cohort": name_b, "Percentage": pct_b}
        )

    c_phe_t, c_phe_f = st.columns([1.1, 0.9])
    with c_phe_t:
        if pheno_rows:
            st.dataframe(pd.DataFrame(pheno_rows), use_container_width=True, hide_index=True)
        else:
            st.info("No phenotype categories configured for this gene.")

    with c_phe_f:
        if pheno_chart_data:
            if HAS_PLOTLY and px is not None:
                fig_cmp_phe = px.bar(
                    pd.DataFrame(pheno_chart_data),
                    x="Phenotype",
                    y="Percentage",
                    color="Cohort",
                    barmode="group",
                    title="Phenotype Shift (% Total Cohort)",
                    color_discrete_sequence=["#2563EB", "#059669"],
                )
                fig_cmp_phe.update_layout(
                    height=260, margin=dict(l=10, r=10, t=32, b=10)
                )
                st.plotly_chart(fig_cmp_phe, use_container_width=True)
            else:
                p_df = pd.DataFrame(pheno_chart_data)
                st.bar_chart(p_df, x="Phenotype", y="Percentage", color="Cohort")

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 3. CLINICAL OUTCOME & DRUG SENSITIVITY PREDICTIONS
    # -------------------------------------------------------------------------
    st.markdown("### 🔮 Automated Clinical Predictions & Actionability Forecast")
    st.caption(
        "Predictive clinical modeling based on observed genetic frequencies and CPIC guideline thresholds."
    )

    insights = _generate_prediction_insights(
        active_gene, name_a, stats_a, name_b, stats_b
    )

    for item in insights:
        border_color = (
            "#EF4444"
            if item["severity"] == "high"
            else ("#F59E0B" if item["severity"] == "medium" else "#3B82F6")
        )
        bg_color = (
            "#FEF2F2"
            if item["severity"] == "high"
            else ("#FFFBEB" if item["severity"] == "medium" else "#EFF6FF")
        )
        badge_color = (
            "#DC2626"
            if item["severity"] == "high"
            else ("#D97706" if item["severity"] == "medium" else "#2563EB")
        )

        st.markdown(
            f"""
        <div style="background: {bg_color}; border: 1.5px solid {border_color}; border-radius: 10px; padding: 16px 18px; margin-bottom: 14px;">
            <div style="display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 8px; margin-bottom: 6px;">
                <div style="font-size: 15px; font-weight: 700; color: #0F172A;">
                    {item['drug']}
                </div>
                <div style="font-size: 11px; font-weight: 700; background: {badge_color}; color: #FFFFFF; padding: 3px 8px; border-radius: 4px;">
                    {item['divergence']}
                </div>
            </div>
            <div style="font-weight: 600; font-size: 13px; color: #1E293B; margin-bottom: 4px;">
                {item['risk_title']}
            </div>
            <div style="font-size: 13px; color: #334155; line-height: 1.5; margin-bottom: 8px;">
                {item['prediction']}
            </div>
            <div style="background: #FFFFFF; border-left: 3px solid {border_color}; padding: 8px 12px; border-radius: 4px; font-size: 12px; color: #1E293B; line-height: 1.4;">
                {item['clinical_action']}
            </div>
        </div>
        """,
            unsafe_allow_html=True,
        )
