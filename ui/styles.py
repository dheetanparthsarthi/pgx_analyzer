"""
Custom CSS styling and HTML design elements.
Includes Georgia serif headers, forest green buttons, KPI cards, and responsive table CSS.
"""

import streamlit as st

CUSTOM_CSS = """
<style>
    /* Main Background & Clean Typography */
    .stApp {
        background-color: #F8FAFC !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
        color: #0F172A !important;
    }
    
    /* Georgia Serif Header Titles */
    h1, h2, h3, .serif-header {
        font-family: Georgia, "Times New Roman", Times, serif !important;
        color: #0F172A !important;
        font-weight: 700 !important;
    }
    
    /* Primary Buttons Styling (Dark Forest Green #166534) */
    button[kind="primary"] {
        background-color: #166534 !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
    }
    button[kind="primary"]:hover {
        background-color: #15803D !important;
    }
    
    /* Secondary Action Buttons */
    button[kind="secondary"] {
        background-color: #475569 !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
    }
    button[kind="secondary"]:hover {
        background-color: #334155 !important;
    }

    /* Clean Card Box Containers */
    .clean-card-box {
        background-color: #FFFFFF;
        border: 1px solid #CBD5E1;
        border-radius: 8px;
        padding: 1.25rem;
        margin-bottom: 1.2rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.03);
    }
    
    /* Upload Box Styling */
    .upload-container-box {
        background-color: #FAFAFA;
        border: 2px dashed #94A3B8;
        border-radius: 8px;
        padding: 2.2rem 1.5rem;
        text-align: center;
        margin-bottom: 1rem;
    }
    .upload-icon-style {
        font-size: 2.5rem;
        color: #475569;
        margin-bottom: 0.3rem;
    }
    .upload-title-style {
        font-size: 1.1rem;
        font-weight: 700;
        color: #0F172A;
    }
    .upload-sub-style {
        font-size: 0.85rem;
        color: #64748B;
        margin-top: 0.2rem;
    }

    /* Privacy Banner */
    .privacy-info-banner {
        background-color: #E2E8F0;
        border-radius: 4px;
        padding: 0.75rem 1rem;
        font-size: 0.85rem;
        color: #334155;
        margin-bottom: 1.2rem;
    }

    /* Metric Cards */
    .kpi-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 6px;
        padding: 0.9rem;
        text-align: left;
        box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    }
    .kpi-val {
        font-family: Georgia, serif;
        font-size: 1.9rem;
        font-weight: 700;
        color: #0F172A;
        line-height: 1.1;
    }
    .kpi-lbl {
        font-size: 0.68rem;
        font-weight: 700;
        text-transform: uppercase;
        color: #64748B;
        letter-spacing: 0.05em;
        margin-top: 0.3rem;
    }

    /* Prevent Metric Text Truncation */
    div[data-testid="stMetricValue"] {
        white-space: nowrap !important;
        overflow: visible !important;
        text-overflow: unset !important;
    }
    div[data-testid="stMetricValue"] > div {
        white-space: nowrap !important;
        overflow: visible !important;
        text-overflow: unset !important;
    }
    div[data-testid="stMetricDelta"] > div {
        white-space: nowrap !important;
        overflow: visible !important;
        text-overflow: unset !important;
    }

    /* Table Clean Border Header Styling */
    div[data-testid="stTable"] table {
        border-collapse: collapse !important;
        width: 100% !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    }
    div[data-testid="stTable"] th {
        border-top: 1.5px solid #0F172A !important;
        border-bottom: 1.5px solid #0F172A !important;
        background-color: transparent !important;
        color: #0F172A !important;
        font-size: 0.75rem !important;
        font-weight: 700 !important;
        letter-spacing: 0.05em !important;
        padding: 0.6rem 0.8rem !important;
    }
    div[data-testid="stTable"] td {
        border-bottom: 1px solid #E2E8F0 !important;
        padding: 0.55rem 0.8rem !important;
        font-size: 0.88rem !important;
        color: #1E293B !important;
    }

    /* Sub-cohort Horizontal Radio Pills Styling */
    div[role="radiogroup"] {
        gap: 0.5rem !important;
    }
    div[role="radiogroup"] label {
        background-color: #E2E8F0 !important;
        border-radius: 20px !important;
        padding: 0.35rem 0.9rem !important;
        font-size: 0.82rem !important;
        font-weight: 600 !important;
        color: #1E293B !important;
        border: none !important;
    }
    div[role="radiogroup"] label[aria-checked="true"], div[role="radiogroup"] label:has(input:checked) {
        background-color: #166534 !important;
        color: #FFFFFF !important;
    }
    
    /* Ensure all tabs are visible, horizontally scrollable, and never cut off */
    div[data-testid="stTabs"] [data-baseweb="tab-list"] {
        gap: 6px !important;
        overflow-x: auto !important;
        white-space: nowrap !important;
        padding-bottom: 4px !important;
        scrollbar-width: thin !important;
    }
    div[data-testid="stTabs"] [data-baseweb="tab"] {
        white-space: nowrap !important;
        padding: 8px 14px !important;
        font-size: 0.88rem !important;
    }

    /* Footer Disclaimer */
    .footer-disclaimer-text {
        text-align: center;
        font-size: 0.8rem;
        color: #64748B;
        margin-top: 2.5rem;
        margin-bottom: 1rem;
    }

    /* Header container transparent with click passthrough */
    [data-testid="stHeader"] {
        background-color: transparent !important;
        z-index: 99999 !important;
        display: flex !important;
        visibility: visible !important;
        pointer-events: none !important;
    }
    
    [data-testid="stToolbar"] {
        visibility: visible !important;
        display: flex !important;
        background-color: transparent !important;
        pointer-events: none !important;
    }

    [data-testid="stHeader"] button,
    [data-testid="stToolbar"] button {
        pointer-events: auto !important;
    }

    /* Ensure Sidebar Expand and Collapse Buttons are visible */
    [data-testid="stExpandSidebarButton"], 
    button[data-testid="stExpandSidebarButton"],
    [data-testid="stSidebarCollapsedControl"],
    [data-testid="collapsedControl"],
    button[aria-label="Expand sidebar"] {
        visibility: visible !important;
        display: flex !important;
        opacity: 1 !important;
        color: #0F172A !important;
        background-color: #FFFFFF !important;
        border: 1.5px solid #94A3B8 !important;
        border-radius: 8px !important;
        padding: 5px 8px !important;
        margin: 0.4rem !important;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12) !important;
        cursor: pointer !important;
        pointer-events: auto !important;
        z-index: 1000000 !important;
    }

    [data-testid="stExpandSidebarButton"]:hover,
    button[data-testid="stExpandSidebarButton"]:hover {
        background-color: #F1F5F9 !important;
        border-color: #334155 !important;
        transform: scale(1.05);
    }

    [data-testid="stSidebarCollapseButton"],
    button[data-testid="stSidebarCollapseButton"],
    button[aria-label="Collapse sidebar"],
    button[aria-label="Close sidebar"] {
        visibility: visible !important;
        display: flex !important;
        color: #0F172A !important;
        background-color: #F8FAFC !important;
        border: 1px solid #CBD5E1 !important;
        border-radius: 6px !important;
        padding: 4px 6px !important;
        cursor: pointer !important;
        pointer-events: auto !important;
    }

    [data-testid="stSidebarCollapseButton"]:hover,
    button[data-testid="stSidebarCollapseButton"]:hover {
        background-color: #E2E8F0 !important;
    }

    /* Hide Streamlit Footer, MainMenu & Badges */
    #MainMenu {visibility: hidden !important; display: none !important;}
    footer {visibility: hidden !important; display: none !important;}
    [data-testid="stStatusWidget"] {visibility: hidden !important; display: none !important;}
    [data-testid="stDecoration"] {visibility: hidden !important; display: none !important;}
    [data-testid="stReportViewer"] footer {display: none !important;}
    .viewerBadge_container__1QSVn {display: none !important;}
    .viewerBadge_link__1S137 {display: none !important;}
    div[class*="viewerBadge"] {display: none !important;}
    div[data-testid="user-avatar"] {display: none !important;}
    button[title*="Streamlit"] {display: none !important;}
    button[title*="View app in Streamlit Community Cloud"] {display: none !important;}
</style>
"""

def inject_custom_styles():
    """Injects custom CSS stylesheet into the Streamlit app."""
    st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


def apply_custom_styles():
    """Alias for inject_custom_styles."""
    inject_custom_styles()


def render_sidebar_controls(df_raw, reset_all_data_fn):
    """Renders persistent sidebar branding and dataset management controls."""
    with st.sidebar:
        st.image("https://cdn-icons-png.flaticon.com/512/3004/3004458.png", width=52)
        st.markdown("<h2 style='font-size: 1.25rem; margin: 0; padding: 0;'>PGx Analytics Pro</h2>", unsafe_allow_html=True)
        st.caption("Automated Population Pharmacogenomics Platform")

        st.divider()

        # Workspace Controls (Reset / Delete / Download)
        st.markdown("<div style='font-size: 11px; font-weight: 700; text-transform: uppercase; color: #64748B; letter-spacing: 0.05em; margin-bottom: 6px;'>⚙️ Workspace Actions</div>", unsafe_allow_html=True)
        if st.session_state.get('show_delete_confirm', False):
            st.warning("⚠️ **Delete Active Dataset?**\nResets workspace to 0 samples.")
            c_del1, c_del2 = st.columns(2)
            with c_del1:
                if st.button("🔴 Delete", type="primary", use_container_width=True, key="btn_sb_confirm_del"):
                    reset_all_data_fn()
                    st.success("Dataset deleted.")
                    st.rerun()
            with c_del2:
                if st.button("❌ Cancel", type="secondary", use_container_width=True, key="btn_sb_cancel_del"):
                    st.session_state['show_delete_confirm'] = False
                    st.rerun()
        else:
            if st.button("🗑️ Delete Active Dataset", type="secondary", use_container_width=True, key="btn_sb_init_del"):
                st.session_state['show_delete_confirm'] = True
                st.rerun()

        if df_raw is not None and not df_raw.empty:
            st.markdown("<div style='font-size: 11px; font-weight: 700; text-transform: uppercase; color: #64748B; letter-spacing: 0.05em; margin-top: 14px; margin-bottom: 6px;'>📥 Download Complete Dataset</div>", unsafe_allow_html=True)
            from ui.views.results_dashboard import get_df_download_excel, get_df_download_csv
            c_sb_dl1, c_sb_dl2 = st.columns(2)
            with c_sb_dl1:
                st.download_button(
                    label="⬇️ Excel (.xlsx)",
                    data=get_df_download_excel(df_raw),
                    file_name="Complete_PGx_Dataset.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True,
                    key="btn_sb_download_xlsx"
                )
            with c_sb_dl2:
                st.download_button(
                    label="⬇️ CSV (.csv)",
                    data=get_df_download_csv(df_raw),
                    file_name="Complete_PGx_Dataset.csv",
                    mime="text/csv",
                    use_container_width=True,
                    key="btn_sb_download_csv"
                )

