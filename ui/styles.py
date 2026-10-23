"""
Custom CSS styling and HTML design elements.
Includes Georgia serif headers, forest green buttons, KPI cards, and responsive table CSS.
"""

import streamlit as st

CUSTOM_CSS = """
<style>
    /* =========================================================================
       1. GLOBAL RESET, TYPOGRAPHY & DESKTOP BASELINE
       ========================================================================= */
    html, body {
        overflow-x: hidden !important;
        max-width: 100vw !important;
        -webkit-font-smoothing: antialiased !important;
        -moz-osx-font-smoothing: grayscale !important;
    }

    *, *::before, *::after {
        box-sizing: border-box !important;
    }

    .stApp {
        background-color: #F8FAFC !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif !important;
        color: #0F172A !important;
        overflow-x: hidden !important;
        max-width: 100vw !important;
    }
    
    /* Georgia Serif Header Titles */
    h1, h2, h3, .serif-header {
        font-family: Georgia, "Times New Roman", Times, serif !important;
        color: #0F172A !important;
        font-weight: 700 !important;
        letter-spacing: -0.01em !important;
    }

    /* Primary Action Buttons (#166534 Dark Forest Green) */
    button[kind="primary"],
    button[data-testid="baseButton-primary"] {
        background-color: #166534 !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
        transition: all 0.15s ease-in-out !important;
    }
    button[kind="primary"]:hover,
    button[data-testid="baseButton-primary"]:hover {
        background-color: #15803D !important;
        box-shadow: 0 2px 6px rgba(22, 101, 52, 0.25) !important;
    }
    
    /* Secondary Action Buttons */
    button[kind="secondary"],
    button[data-testid="baseButton-secondary"] {
        background-color: #475569 !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
        transition: all 0.15s ease-in-out !important;
    }
    button[kind="secondary"]:hover,
    button[data-testid="baseButton-secondary"]:hover {
        background-color: #334155 !important;
        box-shadow: 0 2px 6px rgba(51, 65, 85, 0.25) !important;
    }

    /* Container Card Boxes */
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
        border-radius: 6px;
        padding: 0.75rem 1rem;
        font-size: 0.85rem;
        color: #334155;
        margin-bottom: 1.2rem;
    }

    /* KPI Cards */
    .kpi-card, .cohort-kpi-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 0.9rem;
        text-align: left;
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
        transition: transform 0.15s ease, box-shadow 0.15s ease;
    }
    .kpi-val {
        font-family: Georgia, serif;
        font-size: 1.85rem;
        font-weight: 700;
        color: #0F172A;
        line-height: 1.15;
    }
    .kpi-lbl {
        font-size: 0.68rem;
        font-weight: 700;
        text-transform: uppercase;
        color: #64748B;
        letter-spacing: 0.05em;
        margin-top: 0.3rem;
    }

    /* Clean Metric Value & Label Wrapping (Eliminates text overlap when shrinking page) */
    div[data-testid="stMetric"] {
        min-width: 0 !important;
        overflow: hidden !important;
    }
    div[data-testid="stMetricValue"],
    div[data-testid="stMetricValue"] > div {
        font-family: Georgia, serif !important;
        font-size: clamp(1.15rem, 2.2vw, 1.85rem) !important;
        line-height: 1.2 !important;
        white-space: normal !important;
        overflow-wrap: break-word !important;
        word-break: break-word !important;
    }
    div[data-testid="stMetricLabel"],
    div[data-testid="stMetricLabel"] > div {
        white-space: normal !important;
        overflow-wrap: break-word !important;
        word-break: break-word !important;
        line-height: 1.3 !important;
    }
    div[data-testid="stMetricDelta"],
    div[data-testid="stMetricDelta"] > div {
        white-space: normal !important;
        overflow-wrap: break-word !important;
        word-break: break-word !important;
        line-height: 1.2 !important;
    }

    /* Data Tables & DataFrames */
    div[data-testid="stDataFrame"],
    div[data-testid="stTable"],
    div[data-testid="stTable"] > div {
        width: 100% !important;
        max-width: 100% !important;
        overflow-x: auto !important;
        -webkit-overflow-scrolling: touch !important;
        margin-bottom: 0.75rem !important;
    }

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
        white-space: nowrap !important;
    }
    div[data-testid="stTable"] td {
        border-bottom: 1px solid #E2E8F0 !important;
        padding: 0.55rem 0.8rem !important;
        font-size: 0.88rem !important;
        color: #1E293B !important;
    }

    /* Sub-cohort Horizontal Radio Pills Styling */
    div[role="radiogroup"] {
        display: flex !important;
        flex-wrap: wrap !important;
        gap: 0.45rem !important;
        align-items: center !important;
    }
    div[role="radiogroup"] label {
        background-color: #E2E8F0 !important;
        border-radius: 20px !important;
        padding: 0.35rem 0.9rem !important;
        font-size: 0.82rem !important;
        font-weight: 600 !important;
        color: #1E293B !important;
        border: none !important;
        cursor: pointer !important;
        white-space: nowrap !important;
        transition: all 0.15s ease !important;
    }
    div[role="radiogroup"] label[aria-checked="true"], 
    div[role="radiogroup"] label:has(input:checked) {
        background-color: #166534 !important;
        color: #FFFFFF !important;
    }

    /* =========================================================================
       2. TABS BAR - FULL HORIZONTAL TOUCH SCROLLING & SLEEK TABS
       ========================================================================= */
    div[data-testid="stTabs"] {
        width: 100% !important;
        max-width: 100% !important;
        overflow-x: hidden !important;
    }
    div[data-testid="stTabs"] [role="tablist"],
    div[data-testid="stTabs"] [data-baseweb="tab-list"],
    div[role="tablist"] {
        display: flex !important;
        flex-wrap: nowrap !important;
        overflow-x: auto !important;
        -webkit-overflow-scrolling: touch !important;
        gap: 6px !important;
        padding: 4px 4px 8px 4px !important;
        white-space: nowrap !important;
        scrollbar-width: thin !important;
        scrollbar-color: #94A3B8 transparent !important;
        border-bottom: 2px solid #E2E8F0 !important;
    }
    div[data-testid="stTabs"] [role="tablist"]::-webkit-scrollbar,
    div[data-testid="stTabs"] [data-baseweb="tab-list"]::-webkit-scrollbar,
    div[role="tablist"]::-webkit-scrollbar {
        height: 5px !important;
    }
    div[data-testid="stTabs"] [role="tablist"]::-webkit-scrollbar-track,
    div[data-testid="stTabs"] [data-baseweb="tab-list"]::-webkit-scrollbar-track,
    div[role="tablist"]::-webkit-scrollbar-track {
        background: transparent !important;
    }
    div[data-testid="stTabs"] [role="tablist"]::-webkit-scrollbar-thumb,
    div[data-testid="stTabs"] [data-baseweb="tab-list"]::-webkit-scrollbar-thumb,
    div[role="tablist"]::-webkit-scrollbar-thumb {
        background-color: #CBD5E1 !important;
        border-radius: 4px !important;
    }
    div[data-testid="stTabs"] [role="tablist"]::-webkit-scrollbar-thumb:hover,
    div[data-testid="stTabs"] [data-baseweb="tab-list"]::-webkit-scrollbar-thumb:hover,
    div[role="tablist"]::-webkit-scrollbar-thumb:hover {
        background-color: #94A3B8 !important;
    }
    div[data-testid="stTabs"] [role="tab"],
    div[data-testid="stTabs"] [data-baseweb="tab"],
    div[role="tab"] {
        flex-shrink: 0 !important;
        white-space: nowrap !important;
        padding: 8px 15px !important;
        font-size: 0.88rem !important;
        font-weight: 600 !important;
        border-radius: 6px !important;
        min-height: 40px !important;
        transition: all 0.15s ease !important;
    }
    div[data-testid="stTabs"] [role="tab"]:hover,
    div[data-testid="stTabs"] [data-baseweb="tab"]:hover,
    div[role="tab"]:hover {
        background-color: #F1F5F9 !important;
        color: #166534 !important;
    }
        color: #166534 !important;
    }

    /* =========================================================================
       3. PLOTLY CHARTS & VISUALIZATIONS RESPONSIVENESS
       ========================================================================= */
    .js-plotly-plot, .plot-container, .user-select-none {
        width: 100% !important;
        max-width: 100% !important;
    }

    /* =========================================================================
       4. MOBILE & TABLET RESPONSIVE MEDIA QUERIES (DEVICE COMPATIBILITY)
       ========================================================================= */

    /* Tablet Landscape, Laptops & Shrunk Windows (<= 1024px) */
    @media (max-width: 1024px) {
        .block-container,
        div[data-testid="stMainBlockContainer"] {
            padding-left: 1.25rem !important;
            padding-right: 1.25rem !important;
            padding-top: 1.8rem !important;
            padding-bottom: 3.5rem !important;
            max-width: 100% !important;
        }

        .kpi-val {
            font-size: 1.5rem !important;
        }

        /* Responsive Column Wrapping: Prevents columns and cards from squishing and overlapping text */
        div[data-testid="stHorizontalBlock"] {
            display: flex !important;
            flex-wrap: wrap !important;
            gap: 8px !important;
            width: 100% !important;
        }

        /* Default columns wrap cleanly into 50% width on shrunk viewports */
        div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"] > div[data-testid="column"] {
            min-width: calc(50% - 8px) !important;
            flex: 1 1 calc(50% - 8px) !important;
            max-width: 100% !important;
            box-sizing: border-box !important;
        }

        /* Tables & charts stack into 100% full width to prevent column header overlap */
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stDataFrame"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stDataFrame"]) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stTable"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stTable"]) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(table) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(table) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(.js-plotly-plot) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(.js-plotly-plot) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stPlotlyChart"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stPlotlyChart"]) > div[data-testid="column"] {
            min-width: 100% !important;
            flex: 1 1 100% !important;
            margin-bottom: 0.6rem !important;
        }
    }

    /* Tablets Portrait & Mobile Devices (<= 768px) */
    @media (max-width: 768px) {
        /* Container Spacing */
        .block-container,
        div[data-testid="stMainBlockContainer"],
        div[data-testid="stAppViewBlockContainer"] {
            padding-left: 0.75rem !important;
            padding-right: 0.75rem !important;
            padding-top: 1.5rem !important;
            padding-bottom: 3rem !important;
            max-width: 100% !important;
            overflow-x: hidden !important;
        }

        /* Typography Scaling */
        h1, .serif-header {
            font-size: clamp(1.35rem, 5vw, 1.8rem) !important;
            line-height: 1.25 !important;
        }
        h2 {
            font-size: clamp(1.15rem, 4.2vw, 1.45rem) !important;
        }
        h3 {
            font-size: clamp(1.05rem, 3.6vw, 1.25rem) !important;
        }

        /* Content blocks with tables, dataframes, or plotly charts stack into 100% full width */
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stDataFrame"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stDataFrame"]) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stTable"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stTable"]) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(table) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(table) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(.js-plotly-plot) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(.js-plotly-plot) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stPlotlyChart"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stPlotlyChart"]) > div[data-testid="column"] {
            min-width: 100% !important;
            flex: 1 1 100% !important;
            margin-bottom: 0.6rem !important;
        }

        /* Top 4 Step Navigation Buttons Styling */
        .st-key-nav_btn_step1 button,
        .st-key-nav_btn_step2 button,
        .st-key-nav_btn_step3 button,
        .st-key-nav_btn_step4 button,
        button[key^="nav_btn_step"] {
            width: 100% !important;
            font-size: 0.82rem !important;
            padding: 0.5rem 0.35rem !important;
            min-height: 42px !important;
            white-space: nowrap !important;
            overflow: hidden !important;
            text-overflow: ellipsis !important;
            display: flex !important;
            align-items: center !important;
            justify-content: center !important;
        }

        /* Metric & KPI Cards Styling */
        div[data-testid="stMetric"] {
            background: #FFFFFF !important;
            border: 1px solid #E2E8F0 !important;
            border-radius: 8px !important;
            padding: 0.65rem 0.75rem !important;
            box-shadow: 0 1px 2px rgba(0,0,0,0.03) !important;
            width: 100% !important;
            box-sizing: border-box !important;
        }
        div[data-testid="stMetricValue"] {
            font-size: clamp(1.2rem, 4vw, 1.55rem) !important;
            line-height: 1.2 !important;
        }
        div[data-testid="stMetricLabel"] {
            font-size: 0.66rem !important;
            letter-spacing: 0.04em !important;
            font-weight: 700 !important;
            text-transform: uppercase !important;
            color: #64748B !important;
        }
        .kpi-card, .cohort-kpi-card {
            padding: 0.65rem 0.75rem !important;
            width: 100% !important;
            box-sizing: border-box !important;
        }
        .kpi-val {
            font-size: clamp(1.2rem, 4vw, 1.55rem) !important;
        }
        .kpi-lbl {
            font-size: 0.64rem !important;
        }

        /* Export Center Action Buttons: 2-Column Responsive Grid */
        .export-nav-anchor + div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"],
        .export-nav-anchor + div[data-testid="stHorizontalBlock"] > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(.st-key-btn_toggle_pdf_prev) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(.st-key-btn_toggle_pdf_prev) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(button[key="btn_toggle_pdf_prev"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(button[key="btn_toggle_pdf_prev"]) > div[data-testid="column"],
        div[data-testid="stHorizontalBlock"]:has(button[key="btn_exp_excel"]) > div[data-testid="stColumn"],
        div[data-testid="stHorizontalBlock"]:has(button[key="btn_exp_excel"]) > div[data-testid="column"] {
            min-width: calc(50% - 8px) !important;
            flex: 1 1 calc(50% - 8px) !important;
        }

        /* Tabs adjustments for mobile */
        div[data-testid="stTabs"] [role="tab"],
        div[data-testid="stTabs"] [data-baseweb="tab"],
        div[role="tab"] {
            padding: 7px 12px !important;
            font-size: 0.82rem !important;
            min-height: 38px !important;
        }

        /* Box containers on mobile */
        .clean-card-box {
            padding: 0.9rem !important;
            margin-bottom: 0.9rem !important;
        }
        .upload-container-box {
            padding: 1.4rem 0.9rem !important;
        }
        .upload-icon-style {
            font-size: 2rem !important;
        }

        /* Tables font and padding */
        div[data-testid="stTable"] th {
            padding: 0.45rem 0.6rem !important;
            font-size: 0.7rem !important;
        }
        div[data-testid="stTable"] td {
            padding: 0.45rem 0.6rem !important;
            font-size: 0.8rem !important;
        }

        /* Radio Pills */
        div[role="radiogroup"] label {
            padding: 0.3rem 0.75rem !important;
            font-size: 0.78rem !important;
        }
    }

    /* Small Mobile Phones (e.g. iPhone SE, <= 480px) */
    @media (max-width: 480px) {
        .block-container,
        div[data-testid="stMainBlockContainer"] {
            padding-left: 0.5rem !important;
            padding-right: 0.5rem !important;
            padding-top: 1rem !important;
            padding-bottom: 2.5rem !important;
        }

        button[key^="nav_btn_step"] {
            font-size: 0.74rem !important;
            padding: 0.4rem 0.25rem !important;
        }

        .metric-row-anchor + div[data-testid="stHorizontalBlock"],
        div[data-testid="stHorizontalBlock"]:has(div[data-testid="stMetric"]),
        div[data-testid="stHorizontalBlock"]:has(.kpi-card),
        div[data-testid="stHorizontalBlock"]:has(.cohort-kpi-card) {
            grid-template-columns: 1fr 1fr !important;
            gap: 6px !important;
        }

        .export-nav-anchor + div[data-testid="stHorizontalBlock"],
        div[data-testid="stHorizontalBlock"]:has(button[key="btn_toggle_pdf_prev"]) {
            grid-template-columns: 1fr !important;
        }
    }

    /* =========================================================================
       5. TOUCH-FRIENDLY SIDEBAR CONTROLS & HEADER
       ========================================================================= */
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

    /* Sidebar Expand & Collapse Buttons - Accessible Touch Targets */
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
        padding: 6px 10px !important;
        min-height: 40px !important;
        min-width: 40px !important;
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
        padding: 5px 8px !important;
        min-height: 36px !important;
        cursor: pointer !important;
        pointer-events: auto !important;
    }

    [data-testid="stSidebarCollapseButton"]:hover,
    button[data-testid="stSidebarCollapseButton"]:hover {
        background-color: #E2E8F0 !important;
    }

    /* Clean UI: Hide Streamlit Badges, Footer & Status */
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

