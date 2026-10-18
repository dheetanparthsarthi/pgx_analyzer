"""
Automated Population Pharmacogenomics Analysis Platform
Comprehensive Multi-Page Streamlit SaaS Application

Full implementation of the 12-Section Specification:
1. Title & Aim
2. Input Data Normalization (CYP2C19*2 rs4244285, CYP2C19*3 rs4986893, CYP2C19*17 rs12248560)
3. Visual Welcome Landing Screen with Interactive Manual Data Entry Form & File Upload
4. Data Quality Control (QC) & Missing Values Audit with Invalid Genotype Call Detection & Warnings
5. Genotype Counts & Allele Frequency Engine (p & q calculation, Homogenous vs Heterogenous Zygosity)
6. Hardy-Weinberg Equilibrium Engine (Chi2, P-Value, Haldane Exact Test, Table 3 Reference)
7. Geographic Population Groups (South India, North India, East India, West India, Central India)
8. Dedicated State-Wise Analysis Pages (Individual State Profiles)
9. CYP2C19 Star Allele, Diplotype Calling & CPIC Phenotype Classification (UM, RM, NM, IM, PM)
10. Live Report Preview & Multi-Format Exporter (Excel, CSV, PDF)
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
import io
import re
from collections import Counter
from scipy import stats

import importlib
import core.qc
try:
    importlib.reload(core.qc)
except Exception:
    pass
from core.qc import DataQCEngine
from core.stratification import DemographicStratifier
from core.stats import GeneticStatsEngine
from core.cyp2c19 import CYP2C19Translator
from core.cyp2c9 import CYP2C9Translator
from database.db_manager import DatabaseManager
from reports.report_generator import ReportGenerator
from config import SNP_CONFIG, CYP2C9_SNP_CONFIG, PHENOTYPE_ORDER, CYP2C9_PHENOTYPE_ORDER, STATE_TO_REGION

# -----------------------------------------------------------------------------
# 1. PAGE CONFIG & CUSTOM CSS
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Population Pharmacogenomics Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

CUSTOM_CSS = """
<style>
    /* Main Background & Clean Typography */
    .stApp {
        background-color: #F8FAFC !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
        color: #0F172A !important;
    }
    
    /* Georgia Serif Header Titles matching Screenshots 1-5 */
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
    
    /* Upload Box Styling (Matching Screenshot 1) */
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

    /* Privacy Banner (Matching Screenshot 1) */
    .privacy-info-banner {
        background-color: #E2E8F0;
        border-radius: 4px;
        padding: 0.75rem 1rem;
        font-size: 0.85rem;
        color: #334155;
        margin-bottom: 1.2rem;
    }

    /* Metric Cards (Matching Screenshots 4 & 5) */
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

    /* Table Clean Border Header Styling matching Screenshots 1, 2, 3 */
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

    /* Sub-cohort Horizontal Radio Pills Styling (Matching Screenshot 4 & 5) */
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
    
    /* Footer Disclaimer (Matching Screenshot 1 & 3) */
    .footer-disclaimer-text {
        text-align: center;
        font-size: 0.8rem;
        color: #64748B;
        margin-top: 2.5rem;
        margin-bottom: 1rem;
    }

    /* Header container transparent while keeping sidebar toggle button visible */
    /* Header container transparent with click passthrough so it never blocks page buttons */
    [data-testid="stHeader"] {
        background-color: transparent !important;
        z-index: 99999 !important;
        display: flex !important;
        visibility: visible !important;
        pointer-events: none !important;
    }
    
    /* Keep Toolbar visible so the expand button can render cleanly with click passthrough */
    [data-testid="stToolbar"] {
        visibility: visible !important;
        display: flex !important;
        background-color: transparent !important;
        pointer-events: none !important;
    }

    /* Interactive controls inside header/toolbar have pointer-events auto */
    [data-testid="stHeader"] button,
    [data-testid="stToolbar"] button {
        pointer-events: auto !important;
    }

    /* Ensure BOTH Sidebar Expand and Collapse Buttons are ALWAYS prominently visible & clickable */
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

    /* Hide Streamlit Footer, MainMenu, Status Widget, and Community Cloud Badges (without hiding expand button) */
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
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# 2. SESSION STATE MANAGEMENT & FAST PIPELINE CACHING
# -----------------------------------------------------------------------------
DELETED_FLAG_FILE = ".data_deleted.flag"
ACTIVE_DATASET_CACHE = ".active_dataset.pkl"

def save_active_dataset(df: pd.DataFrame):
    """Saves active dataset to disk so it persists across page reloads and browser refreshes."""
    if df is not None and len(df) > 0:
        try:
            df.to_pickle(ACTIVE_DATASET_CACHE)
        except Exception:
            pass

def load_persisted_dataset() -> pd.DataFrame:
    """Loads dataset from persistent disk cache file if present."""
    if os.path.exists(ACTIVE_DATASET_CACHE):
        try:
            return pd.read_pickle(ACTIVE_DATASET_CACHE)
        except Exception:
            pass
    return None

def clear_active_dataset_cache():
    """Removes persistent dataset cache file."""
    if os.path.exists(ACTIVE_DATASET_CACHE):
        try:
            os.remove(ACTIVE_DATASET_CACHE)
        except Exception:
            pass

if os.path.exists(DELETED_FLAG_FILE):
    st.session_state['is_deleted'] = True

if 'is_deleted' not in st.session_state:
    st.session_state['is_deleted'] = False

if 'show_delete_confirm' not in st.session_state:
    st.session_state['show_delete_confirm'] = False

if 'uploaded_df' not in st.session_state:
    st.session_state['uploaded_df'] = None

if 'mapped_cols' not in st.session_state:
    st.session_state['mapped_cols'] = {}

if 'demographic_mapping' not in st.session_state:
    st.session_state['demographic_mapping'] = {}

if 'datasets_dict' not in st.session_state:
    st.session_state['datasets_dict'] = {}

if 'selected_gene' not in st.session_state:
    st.session_state['selected_gene'] = 'CYP2C19'

if 'selected_markers' not in st.session_state:
    st.session_state['selected_markers'] = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']

if 'workflow_step' not in st.session_state:
    st.session_state['workflow_step'] = '1 · Select Target'

if 'custom_genes' not in st.session_state:
    st.session_state['custom_genes'] = []

if 'custom_markers' not in st.session_state:
    st.session_state['custom_markers'] = {}

if 'extra_mapping_rows' not in st.session_state:
    st.session_state['extra_mapping_rows'] = []

if 'mapping_view_mode' not in st.session_state:
    st.session_state['mapping_view_mode'] = '🗂️ Interactive Marker Cards (with Live QC)'

if 'custom_snp_rules' not in st.session_state:
    st.session_state['custom_snp_rules'] = {}

if 'mapped_marker_details' not in st.session_state:
    st.session_state['mapped_marker_details'] = {}

EMPTY_COLS = [
    'sample ID', 'Gender', 'Date of Birth', 'Native place ', 'State',
    'Test requested', 'Is their family lived at Native place for past 3 generations?',
    'CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 ( rs12248560)'
]

# Auto-restore active dataset from disk cache if session state was cleared by a page reload/refresh
if not st.session_state.get('is_deleted', False) and (st.session_state['uploaded_df'] is None or len(st.session_state['uploaded_df']) == 0):
    persisted_df = load_persisted_dataset()
    if persisted_df is not None and len(persisted_df) > 0:
        st.session_state['uploaded_df'] = persisted_df
        if 'workflow_step' not in st.session_state:
            st.session_state['workflow_step'] = '1 · Select Target'

if st.session_state['uploaded_df'] is None:
    st.session_state['uploaded_df'] = pd.DataFrame(columns=EMPTY_COLS)

db_manager = DatabaseManager()

@st.cache_data
def load_marker_catalog() -> pd.DataFrame:
    """Loads and caches the 100 genes and 147 biomarkers from List of markers_v1.xlsx."""
    excel_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "List of markers_v1.xlsx")
    if os.path.exists(excel_path):
        try:
            df = pd.read_excel(excel_path)
            df.columns = [str(c).strip() for c in df.columns]
            df['Gene'] = df['Gene'].astype(str).str.strip()
            df['Marker'] = df['Marker'].astype(str).str.strip()
            return df
        except Exception:
            pass
    return pd.DataFrame()

def parse_maf_value(val) -> tuple[float, str]:
    """Parses numeric MAF frequency and formatted string from raw excel value."""
    if val is None or pd.isna(val) or str(val).strip().lower() in ['nan', 'none', '']:
        return 0.0, "N/A"
    s = str(val).strip()
    is_pct = '%' in s
    if '=' in s:
        part = s.split('=')[-1].strip()
        try:
            return float(part), s
        except Exception:
            pass
    clean_s = s.replace('~', '').replace('%', '').strip()
    if '-' in clean_s:
        try:
            low, high = clean_s.split('-')
            avg = (float(low) + float(high)) / 2.0
            return avg / 100.0, s
        except Exception:
            pass
    try:
        f = float(clean_s)
        if is_pct or f > 1.0:
            return f / 100.0, s
        return f, s
    except Exception:
        return 0.0, s

def analyze_generic_marker(series: pd.Series, marker_name: str, ref_al: str = None, var_al: str = None) -> dict:
    """Computes genotype counts, allele frequencies (p & q), and HWE for any marker."""
    clean = series.dropna().astype(str).str.strip().str.upper()
    clean = clean[~clean.isin(['', 'NAN', 'NONE', 'MISSING', 'UNKNOWN', 'NULL'])]
    n = len(clean)
    if n == 0:
        return {
            'valid_samples': 0,
            'genotype_counts': {},
            'allele_counts': {},
            'allele_freqs': {},
            'p': 0.0, 'q': 0.0,
            'hwe_chi2': 0.0, 'hwe_p': 1.0,
            'interpretation': 'No Data'
        }
    
    gt_counts = dict(clean.value_counts())
    all_alleles = Counter()
    for gt, count in gt_counts.items():
        clean_gt = gt.replace('/', '').replace('|', '').replace(' ', '')
        if len(clean_gt) == 2:
            all_alleles[clean_gt[0]] += count
            all_alleles[clean_gt[1]] += count
        elif len(clean_gt) == 1:
            all_alleles[clean_gt[0]] += count * 2
            
    sorted_als = [al for al, _ in all_alleles.most_common()]
    if not ref_al and len(sorted_als) > 0:
        ref_al = sorted_als[0]
    if not var_al and len(sorted_als) > 1:
        var_al = sorted_als[1]
    elif not var_al and len(sorted_als) == 1:
        var_al = 'Var'
        
    tot_al = sum(all_alleles.values())
    al_freqs = {al: round(c / tot_al, 4) for al, c in all_alleles.items()} if tot_al > 0 else {}
    p = al_freqs.get(ref_al, 1.0)
    q = al_freqs.get(var_al, 0.0)
    
    exp_wt = n * (p ** 2)
    exp_het = n * (2 * p * q)
    exp_var = n * (q ** 2)
    
    wt_gt = (ref_al + ref_al) if ref_al else ""
    het_gt1 = (ref_al + var_al) if ref_al and var_al else ""
    het_gt2 = (var_al + ref_al) if ref_al and var_al else ""
    var_gt = (var_al + var_al) if var_al else ""
    
    obs_wt = gt_counts.get(wt_gt, 0)
    obs_het = gt_counts.get(het_gt1, 0) + (gt_counts.get(het_gt2, 0) if het_gt1 != het_gt2 else 0)
    obs_var = gt_counts.get(var_gt, 0)
    
    chi2 = 0.0
    if exp_wt > 0: chi2 += ((obs_wt - exp_wt)**2) / exp_wt
    if exp_het > 0: chi2 += ((obs_het - exp_het)**2) / exp_het
    if exp_var > 0: chi2 += ((obs_var - exp_var)**2) / exp_var
    
    p_val = 1 - stats.chi2.cdf(chi2, df=1) if chi2 > 0 else 1.0
    interp = "Departure from HWE (χ² ≥ 3.4)" if chi2 >= 3.4 else "In HWE"
    
    return {
        'valid_samples': n,
        'genotype_counts': gt_counts,
        'allele_counts': dict(all_alleles),
        'allele_freqs': al_freqs,
        'ref_allele': ref_al,
        'var_allele': var_al,
        'p': p, 'q': q,
        'hwe_chi2': round(chi2, 4),
        'hwe_p': round(p_val, 4),
        'interpretation': interp
    }

def inspect_uploaded_file(df: pd.DataFrame) -> dict:
    """Checks for Sample ID (compulsory), Gender, and State in uploaded dataframe."""
    sid_col = find_matching_gene_column('Sample_ID', df)
    gen_col = find_matching_gene_column('Gender', df)
    state_col = find_matching_gene_column('State', df)
    return {
        'has_sample_id': sid_col is not None,
        'sample_id_col': sid_col,
        'has_gender': gen_col is not None,
        'gender_col': gen_col,
        'has_state': state_col is not None,
        'state_col': state_col
    }

def find_matching_gene_column(target_key: str, df: pd.DataFrame) -> str:
    """Finds best matching column for a given gene/SNP/demographic in a dataframe."""
    if df is None or df.empty:
        return None
    cols = list(df.columns)
    
    patterns = {
        'CYP2C19*2': ['cyp2c19*2', 'rs4244285', 'c.681g>a', 'c.681', 'cyp2c19_2', 'cyp2'],
        'CYP2C19*3': ['cyp2c19*3', 'rs4986893', 'c.636g>a', 'c.636', 'cyp2c19_3', 'cyp3'],
        'CYP2C19*17': ['*17', 'cyp2c19*17', 'rs12248560', 'c.-806c>t', 'c.-806', 'cyp2c19_17', 'cyp17'],
        'CYP2C9*2': ['cyp2c9*2', 'rs1799853', 'cyp2c9_2', 'cyp2c9-2', '430c>t'],
        'CYP2C9*3': ['cyp2c9*3', 'rs1057910', 'cyp2c9_3', 'cyp2c9-3', '1075a>c'],
        'Sample_ID': ['sample id', 'sample_id', 'sample', 'patient id', 'patient_id', 'lab id', 'id'],
        'Gender': ['gender', 'sex'],
        'State': ['state / country', 'state/country', 'state', 'region', 'native place', 'native', 'location'],
        'DOB': ['date of birth', 'dob', 'birth']
    }
    
    target_patterns = patterns.get(target_key, [])
    for p in target_patterns:
        for c in cols:
            c_clean = str(c).lower().replace(" ", "").replace("_", "").replace("-", "")
            p_clean = p.lower().replace(" ", "").replace("_", "").replace("-", "")
            if p_clean in c_clean:
                return c
                
    # Generic fallback: match target_key directly in column names
    tk_clean = str(target_key).lower().replace(" ", "").replace("_", "").replace("-", "")
    for c in cols:
        c_clean = str(c).lower().replace(" ", "").replace("_", "").replace("-", "")
        if tk_clean == c_clean or tk_clean in c_clean or c_clean in tk_clean:
            return c
    return None

KNOWN_PGX_GENES = [
    'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP3A4', 'VKORC1', 
    'DPYD', 'TPMT', 'NUDT15', 'SLCO1B1', 'UGT1A1', 'HLA-B', 'HLA-A', 'G6PD', 'CYP1A2', 'CYP2B6'
]

def get_column_gene(col_name: str) -> str:
    """Returns the gene symbol found in a column header, if any."""
    c_clean = str(col_name).upper().replace("_", "").replace("-", "").replace(" ", "")
    for g in sorted(KNOWN_PGX_GENES, key=len, reverse=True):
        g_clean = g.upper().replace("_", "").replace("-", "")
        if g_clean in c_clean:
            return g
    return ""

def detect_genes_in_columns(cols: list) -> list:
    """Detects any known pharmacogenetic genes mentioned across column headers."""
    found_genes = []
    for c in cols:
        g = get_column_gene(c)
        if g and g not in found_genes:
            found_genes.append(g)
    return found_genes

def smart_merge_datasets(base_df: pd.DataFrame, new_df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Intelligently aligns gene and demographic columns between base_df and new_df
    before concatenating to ensure genotypes merge into the exact correct gene columns.
    """
    if base_df is None or base_df.empty:
        return new_df.copy(), ["Base dataset empty; adopted newly uploaded dataset."]
    if new_df is None or new_df.empty:
        return base_df.copy(), ["New dataset empty; retained current base dataset."]

    new_aligned = new_df.copy()
    alignment_log = []
    rename_dict = {}

    target_keys = ['Sample_ID', 'CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17', 'Gender', 'State', 'DOB']

    for key in target_keys:
        base_match = find_matching_gene_column(key, base_df)
        new_match = find_matching_gene_column(key, new_aligned)

        if base_match and new_match and base_match != new_match:
            rename_dict[new_match] = base_match
            alignment_log.append(f"🧬 **{key}**: Aligned new column `{new_match}` ➔ existing column `{base_match}`")
        elif base_match and new_match and base_match == new_match:
            alignment_log.append(f"✅ **{key}**: Identical column `{base_match}` matched")

    if rename_dict:
        new_aligned = new_aligned.rename(columns=rename_dict)

    merged_df = pd.concat([base_df, new_aligned], ignore_index=True)
    return merged_df, alignment_log

def clear_deletion_flag():
    """Clears the disk-based permanent deletion marker file and purges cache."""
    st.session_state['is_deleted'] = False
    if os.path.exists(DELETED_FLAG_FILE):
        try:
            os.remove(DELETED_FLAG_FILE)
        except Exception:
            pass
    try:
        st.cache_data.clear()
    except Exception:
        pass

def reset_all_data():
    """Clears session dataset, resets custom column mappings, creates permanent deletion marker, and purges database & cache."""
    st.session_state['is_deleted'] = True
    st.session_state['uploaded_df'] = pd.DataFrame(columns=EMPTY_COLS)
    st.session_state['mapped_cols'] = {}
    st.session_state['demographic_mapping'] = {}
    st.session_state['datasets_dict'] = {}
    st.session_state['last_uploaded_file_id'] = None
    st.session_state['pending_new_df'] = None
    st.session_state['pending_file_id'] = None
    st.session_state['pending_file_name'] = None
    st.session_state['merge_alignment_log'] = None
    st.session_state['show_delete_confirm'] = False
    st.session_state['workflow_step'] = '1 · Select Target'
    
    clear_active_dataset_cache()

    # Write persistent disk marker so deletion survives app restarts, browser refreshes, and new tabs
    try:
        with open(DELETED_FLAG_FILE, "w") as f:
            f.write("deleted")
    except Exception:
        pass

    try:
        db_manager.clear_database()
    except Exception:
        pass

    try:
        st.cache_data.clear()
    except Exception:
        pass

# High performance in-memory cached analysis execution pipeline
@st.cache_data(show_spinner=False)
def run_fast_pipeline(df: pd.DataFrame, mapped_cols_tuple: tuple, snp_rules_tuple: tuple = ()):
    mapped_dict = dict(mapped_cols_tuple)
    snp_rules = {k: dict(v) for k, v in snp_rules_tuple} if snp_rules_tuple else None
    if df is None or len(df) == 0:
        clean_df = pd.DataFrame()
        qc_report = {
            'total_samples': 0,
            'fully_valid_samples': 0,
            'duplicate_sample_count': 0,
            'duplicate_sample_ids': [],
            'invalid_genotypes_count': 0,
            'invalid_genotypes_log': {},
            'snp_qc_stats': {}
        }
        full_results = {
            'overall': {'sample_count': 0, 'snps': {}, 'cyp2c19_summary': {'phenotypes': {}, 'diplotypes': {}}},
            'regional': {},
            'state_wise': {},
            'gender': {},
            'processed_dataframe': pd.DataFrame()
        }
        return clean_df, qc_report, full_results

    try:
        qc_engine = DataQCEngine(df, custom_col_mapping=mapped_dict, custom_snp_rules=snp_rules)
    except TypeError:
        qc_engine = DataQCEngine(df, custom_col_mapping=mapped_dict)
        if snp_rules:
            try:
                if hasattr(qc_engine, 'snp_rules'):
                    for k, v in snp_rules.items():
                        if k in qc_engine.snp_rules:
                            qc_engine.snp_rules[k] = {**qc_engine.snp_rules[k], **v}
                        else:
                            qc_engine.snp_rules[k] = v
            except Exception:
                pass
    clean_df, qc_report = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    return clean_df, qc_report, full_results

# Execute Pipeline
df_raw = st.session_state['uploaded_df']
mapped_tuple = tuple(sorted(st.session_state['mapped_cols'].items()))

def serialize_rules(rules_dict):
    if not rules_dict:
        return ()
    serialized = []
    for k, v in sorted(rules_dict.items()):
        inner_items = []
        for sk, sv in sorted(v.items()):
            if isinstance(sv, (set, list, tuple)):
                inner_items.append((sk, tuple(sorted(sv))))
            else:
                inner_items.append((sk, sv))
        serialized.append((k, tuple(inner_items)))
    return tuple(serialized)

snp_rules_items = serialize_rules(st.session_state.get('custom_snp_rules', {}))
clean_df, qc_report, full_results = run_fast_pipeline(df_raw, mapped_tuple, snp_rules_items)
processed_df = full_results.get('processed_dataframe', pd.DataFrame())

# Helper to check active dataset state
def has_active_data() -> bool:
    if full_results is None or 'overall' not in full_results:
        return False
    return full_results['overall'].get('sample_count', 0) > 0

@st.cache_data(show_spinner=False)
def get_df_download_excel(df: pd.DataFrame) -> bytes:
    if df is None or df.empty:
        return b""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='Complete Dataset', index=False)
    return output.getvalue()

@st.cache_data(show_spinner=False)
def get_df_download_csv(df: pd.DataFrame) -> bytes:
    if df is None or df.empty:
        return b""
    return df.to_csv(index=False).encode('utf-8')

# -----------------------------------------------------------------------------
# 3. SIDEBAR NAVIGATION MENU (WITH ACTIVE WORKFLOW & STATUS)
# -----------------------------------------------------------------------------
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
                reset_all_data()
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

# -----------------------------------------------------------------------------
# CYP2C19 COMBINATION & GEOGRAPHIC STUDIO (DOWNSIDE OF TAB 5)
# -----------------------------------------------------------------------------
def render_cyp2c19_combination_studio(processed_df: pd.DataFrame, full_results: dict, view_mode: str = "All Views"):
    """
    Renders interactive CYP2C19 Multi-SNP Combination, Diplotype, Phenotype,
    and Geographic/State-Wise matrices with active, working scientific rsID nomenclature.
    """
    if processed_df is None or processed_df.empty:
        st.info("No pharmacogenomic data available for combination analysis.")
        return

    # -------------------------------------------------------------------------
    # STUDIO CONTAINER & HEADER (CLEAN, COMPACT TYPOGRAPHY)
    # -------------------------------------------------------------------------
    target_gene = full_results.get('target_gene') or st.session_state.get('selected_gene', 'CYP2C19')
    st.markdown(f"""
    <div style="border-top: 1px solid #E2E8F0; margin-top: 1.5rem; padding-top: 1rem; margin-bottom: 0.75rem;">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div>
                <span style="font-size: 16px; font-weight: 700; color: #0F172A;">🧩 {target_gene} Multi-Variant & Geographic Studio</span>
                <span style="margin-left: 8px; font-size: 11px; font-weight: 600; background: #EEF2F6; color: #334155; padding: 2px 7px; border-radius: 4px; border: 1px solid #CBD5E1;">CPIC / PharmGKB</span>
            </div>
        </div>
        <p style="margin: 3px 0 0 0; font-size: 12px; color: #64748B;">
            Standardized multi-locus genotype stratification, scientific nomenclature mapping, and sub-population distribution analysis.
        </p>
    </div>
    """, unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 1. VARIANT COMBINATION FOCUS SELECTOR (DYNAMIC ACROSS GENES)
    # -------------------------------------------------------------------------
    mapped_snps_in_df = [c for c in full_results.get('snps', {}).keys() if c in processed_df.columns]
    if not mapped_snps_in_df:
        mapped_snps_in_df = [c for c in processed_df.columns if c not in [
            'Sample_ID', 'Gender_Clean', 'Gender', 'State_Clean', 'State', 
            'Native_Place_Clean', 'Native_Place', 'Region', 'Diplotype', 
            'Phenotype', 'Activity_Score', 'Test_Requested', 'Date_of_Birth'
        ]]

    cyp2c19_markers_present = [c for c in ['CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17'] if c in processed_df.columns]
    cyp2c9_markers_present = [c for c in ['CYP2C9*2', 'CYP2C9*3'] if c in processed_df.columns]

    if len(cyp2c19_markers_present) >= 2:
        st.markdown("<div style='font-size: 12px; font-weight: 600; color: #334155; margin-bottom: 4px;'>1️⃣ Select Variant Combination to Analyze:</div>", unsafe_allow_html=True)
        avail_options = []
        if 'CYP2C19*2' in cyp2c19_markers_present and 'CYP2C19*3' in cyp2c19_markers_present:
            avail_options.append("CYP2C19 *2 + *3 (c.681G>A / c.636G>A · Loss of Function)")
        if 'CYP2C19*2' in cyp2c19_markers_present and 'CYP2C19*17' in cyp2c19_markers_present:
            avail_options.append("CYP2C19 *2 + *17 (c.681G>A / c.-806C>T · LOF + Increased Function)")
        if 'CYP2C19*3' in cyp2c19_markers_present and 'CYP2C19*17' in cyp2c19_markers_present:
            avail_options.append("CYP2C19 *3 + *17 (c.636G>A / c.-806C>T · LOF + Increased Function)")
        if len(cyp2c19_markers_present) == 3:
            avail_options.append("CYP2C19 Core Panel (*2 + *3 + *17 · Full 3-Locus Trio)")
        if not avail_options:
            avail_options = [f"CYP2C19 Panel ({' + '.join(cyp2c19_markers_present)})"]

        variant_combo_choice = st.radio(
            "Variant Combination:",
            avail_options,
            horizontal=True,
            key="studio_variant_combo_choice",
            label_visibility="collapsed"
        )
        if "*2 + *3" in variant_combo_choice:
            active_snps = ['CYP2C19*2', 'CYP2C19*3']
            combo_title = "CYP2C19 *2 + *3"
            std_snps_label = "c.681G>A & c.636G>A"
            rs_label = "rs4244285 + rs4986893"
        elif "*2 + *17" in variant_combo_choice:
            active_snps = ['CYP2C19*2', 'CYP2C19*17']
            combo_title = "CYP2C19 *2 + *17"
            std_snps_label = "c.681G>A & c.-806C>T"
            rs_label = "rs4244285 + rs12248560"
        elif "*3 + *17" in variant_combo_choice:
            active_snps = ['CYP2C19*3', 'CYP2C19*17']
            combo_title = "CYP2C19 *3 + *17"
            std_snps_label = "c.636G>A & c.-806C>T"
            rs_label = "rs4986893 + rs12248560"
        elif "Core Panel" in variant_combo_choice:
            active_snps = ['CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17']
            combo_title = "CYP2C19 *2 + *3 + *17"
            std_snps_label = "c.681G>A, c.636G>A & c.-806C>T"
            rs_label = "rs4244285 + rs4986893 + rs12248560"
        else:
            active_snps = cyp2c19_markers_present
            combo_title = " + ".join(active_snps)
            std_snps_label = f"{len(active_snps)}-Locus Panel"
            rs_label = combo_title
    elif len(cyp2c9_markers_present) >= 2:
        st.markdown("<div style='font-size: 12px; font-weight: 600; color: #334155; margin-bottom: 4px;'>1️⃣ Select Variant Combination to Analyze:</div>", unsafe_allow_html=True)
        avail_options = ["CYP2C9 *2 + *3 (c.430C>T / c.1075A>C · Decreased Function Panel)"]
        variant_combo_choice = st.radio(
            "Variant Combination:",
            avail_options,
            horizontal=True,
            key="studio_variant_combo_choice_cyp2c9",
            label_visibility="collapsed"
        )
        active_snps = ['CYP2C9*2', 'CYP2C9*3']
        combo_title = "CYP2C9 *2 + *3"
        std_snps_label = "c.430C>T & c.1075A>C"
        rs_label = "rs1799853 + rs1057910"
    elif len(mapped_snps_in_df) >= 1:
        active_snps = mapped_snps_in_df
        combo_title = f"{target_gene} ({' + '.join(active_snps)})"
        std_snps_label = f"{len(active_snps)}-Marker Panel"
        rs_label = combo_title
        st.markdown(f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin-bottom: 6px;'>1️⃣ Analyzed Multi-Locus Panel: <b>{combo_title}</b></div>", unsafe_allow_html=True)
    else:
        st.info("ℹ️ Multi-Locus Combination Studio requires at least one mapped marker column with valid genotype calls.")
        return

    active_snps = [c for c in active_snps if c in processed_df.columns]
    if not active_snps:
        st.info("ℹ️ No valid marker columns found in the active dataset to build a combination studio matrix.")
        return

    state_col = 'State_Clean' if 'State_Clean' in processed_df.columns else ('State' if 'State' in processed_df.columns else None)

    # -------------------------------------------------------------------------
    # 2. GEOGRAPHIC LEVEL, SPECIFIC DRILL-DOWN, AND NOMENCLATURE FORMAT
    # -------------------------------------------------------------------------
    c_geo, c_specific, c_nom = st.columns([1.1, 1.4, 1.1])

    with c_geo:
        default_scope_idx = 1 if "Geographical" in view_mode else 0
        geo_scope = st.selectbox(
            "🌍 Geographic Level:",
            ["🌐 Overall Cohort (Pan-India)", "🗺️ Geographical Regions", "📍 State-Wise Profiles"],
            index=default_scope_idx,
            key="studio_geo_scope"
        )

    selected_region = "All Regions (Comparative Overview)"
    selected_state = "All States (Comparative Overview)"

    with c_specific:
        if "Geographical Regions" in geo_scope:
            region_list = ["All Regions (Comparative Overview)"]
            if 'Region' in processed_df.columns:
                unique_regs = sorted([str(r) for r in processed_df['Region'].dropna().unique() if str(r).strip() != ''])
                region_list.extend(unique_regs)
            selected_region = st.selectbox(
                "🗺️ Specific Geographical Region:",
                region_list,
                key="studio_sel_region"
            )
        elif "State-Wise Profiles" in geo_scope:
            state_list = ["All States (Comparative Overview)"]
            if state_col and state_col in processed_df.columns:
                unique_states = sorted([str(s) for s in processed_df[state_col].dropna().unique() if str(s).strip() != ''])
                state_list.extend(unique_states)
            selected_state = st.selectbox(
                "📍 Specific State Profile:",
                state_list,
                key="studio_sel_state"
            )
        else:
            st.selectbox(
                "📍 Stratification Status:",
                ["Pan-Cohort (All 28 States & UTs Included)"],
                disabled=True,
                key="studio_disabled_geo"
            )

    with c_nom:
        nomenclature = st.selectbox(
            "🔬 Nomenclature Format:",
            ["🧬 Star-Allele (*2, *3, *17)", "🔬 Scientific rsID (NCBI / HGVS)"],
            key="studio_nomenclature"
        )

    is_sci = "Scientific" in nomenclature

    # -------------------------------------------------------------------------
    # 3. DYNAMIC DATA FILTERING BASED ON GEOGRAPHY / STATE
    # -------------------------------------------------------------------------
    active_df = processed_df.copy()
    scope_label = "Overall Cohort"

    if "Geographical Regions" in geo_scope:
        if selected_region != "All Regions (Comparative Overview)":
            active_df = active_df[active_df['Region'].astype(str) == str(selected_region)]
            scope_label = f"Region: {selected_region}"
        else:
            scope_label = "All Regions (Comparative)"
    elif "State-Wise Profiles" in geo_scope:
        if selected_state != "All States (Comparative Overview)" and state_col:
            active_df = active_df[active_df[state_col].astype(str) == str(selected_state)]
            scope_label = f"State: {selected_state}"
        else:
            scope_label = "All States (Comparative)"

    tot_cohort = len(active_df)
    if tot_cohort == 0:
        st.warning(f"⚠️ No patient records found for the selected geographic filter: **{scope_label}**.")
        return

    # Compute combination breakdown across the active filtered dataset
    combo_grp = active_df.groupby(active_snps, dropna=False).size().reset_index(name='Sample Count (N)')
    combo_grp = combo_grp.sort_values(by='Sample Count (N)', ascending=False)
    combo_grp['Frequency (%)'] = combo_grp['Sample Count (N)'].apply(lambda x: f"{(x / tot_cohort * 100):.1f}%" if tot_cohort > 0 else "0.0%")
    combo_grp['pct_val'] = combo_grp['Sample Count (N)'].apply(lambda x: round(x / tot_cohort * 100, 1) if tot_cohort > 0 else 0.0)

    # Add Combination Pattern string
    if len(active_snps) == 1:
        combo_grp['Combination Pattern'] = combo_grp[active_snps[0]].astype(str)
    elif len(active_snps) == 2:
        combo_grp['Combination Pattern'] = combo_grp.apply(lambda r: f"{r[active_snps[0]]} + {r[active_snps[1]]}", axis=1)
    else:
        combo_grp['Combination Pattern'] = combo_grp.apply(lambda r: " + ".join([str(r[c]) for c in active_snps]), axis=1)

    # Associated Diplotypes and Phenotypes for each combination in active_df
    def get_associated_dips(row):
        if 'Diplotype' not in active_df.columns:
            return "N/A"
        cond = True
        for col in active_snps:
            cond = cond & (active_df[col] == row[col])
        sub = active_df[cond]
        dips = ", ".join(sorted([str(d) for d in sub['Diplotype'].dropna().unique()]))
        return dips if dips else "N/A"

    def get_associated_phes(row):
        if 'Phenotype' not in active_df.columns:
            return "N/A"
        cond = True
        for col in active_snps:
            cond = cond & (active_df[col] == row[col])
        sub = active_df[cond]
        phes = ", ".join(sorted([str(p) for p in sub['Phenotype'].dropna().unique()]))
        return phes if phes else "N/A"

    combo_grp['Associated Diplotypes'] = combo_grp.apply(get_associated_dips, axis=1)
    combo_grp['CPIC Phenotype(s)'] = combo_grp.apply(get_associated_phes, axis=1)

    # Top KPI Metrics reflecting selected geography
    top_pattern = combo_grp.iloc[0]['Combination Pattern'] if not combo_grp.empty else 'N/A'
    top_pct = combo_grp.iloc[0]['Frequency (%)'] if not combo_grp.empty else '0%'
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
        st.markdown(f"""
        <div style="{card_style}">
            <div style="{label_style}">Variant Target</div>
            <div style="{val_style}">{combo_title}</div>
            <div style="{sub_style}"><span style="color: #2563EB;">{std_snps_label}</span></div>
        </div>
        """, unsafe_allow_html=True)
    with k2:
        st.markdown(f"""
        <div style="{card_style}">
            <div style="{label_style}">Geographic Scope</div>
            <div style="{val_style}">{scope_label}</div>
            <div style="{sub_style}"><span style="color: #059669;">Stratified Analysis</span></div>
        </div>
        """, unsafe_allow_html=True)
    with k3:
        st.markdown(f"""
        <div style="{card_style}">
            <div style="{label_style}">Cohort Size</div>
            <div style="{val_style}">{tot_cohort:,} Samples</div>
            <div style="{sub_style}">{num_patterns} Distinct Patterns</div>
        </div>
        """, unsafe_allow_html=True)
    with k4:
        st.markdown(f"""
        <div style="{card_style}">
            <div style="{label_style}">Predominant Genotype</div>
            <div style="{val_style}"><span style="color: #166534;">{top_pattern}</span></div>
            <div style="{sub_style}">{top_pct} of Active Cohort</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # 4. DUAL VISUAL CHARTS (DONUT / PIE CHART + COMPARATIVE BAR CHART)
    # -------------------------------------------------------------------------
    c_pie, c_bar = st.columns([1.0, 1.0])
    
    with c_pie:
        fig_donut = px.pie(
            combo_grp,
            names='Combination Pattern',
            values='Sample Count (N)',
            hole=0.48,
            title=f"Proportion Distribution ({scope_label})",
            color_discrete_sequence=['#166534', '#1E40AF', '#D97706', '#991B1B', '#0D9488', '#6B21A8', '#475569']
        )
        fig_donut.update_layout(
            height=280,
            margin=dict(l=10, r=10, t=32, b=10),
            title_font=dict(size=13, color='#0F172A', family='Inter, -apple-system, sans-serif'),
            legend=dict(font=dict(size=10), orientation='v')
        )
        st.plotly_chart(fig_donut, use_container_width=True)

    with c_bar:
        fig_bars = px.bar(
            combo_grp,
            x='Sample Count (N)',
            y='Combination Pattern',
            orientation='h',
            text='Frequency (%)',
            title=f"Sample Counts & Frequencies ({scope_label})",
            color='Combination Pattern',
            color_discrete_sequence=['#166534', '#1E40AF', '#D97706', '#991B1B', '#0D9488', '#6B21A8', '#475569']
        )
        fig_bars.update_layout(
            height=280,
            margin=dict(l=10, r=10, t=32, b=10),
            yaxis=dict(autorange="reversed", tickfont=dict(size=11)),
            xaxis=dict(tickfont=dict(size=11), gridcolor='#F1F5F9'),
            title_font=dict(size=13, color='#0F172A', family='Inter, -apple-system, sans-serif'),
            showlegend=False
        )
        fig_bars.update_traces(textposition='outside', textfont=dict(size=10))
        st.plotly_chart(fig_bars, use_container_width=True)

    # -------------------------------------------------------------------------
    # 5. GEOGRAPHICAL & STATE-WISE BREAKDOWN CHARTS & TABLES
    # -------------------------------------------------------------------------
    if "Geographical Regions" in geo_scope:
        if selected_region != "All Regions (Comparative Overview)":
            st.markdown(f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin: 12px 0 6px 0;'>🗺️ State Breakdown within {selected_region} (N = {tot_cohort})</div>", unsafe_allow_html=True)
            if state_col and state_col in active_df.columns:
                st_counts = active_df[state_col].value_counts().reset_index()
                st_counts.columns = ['State', 'Patients (N)']
                if len(st_counts) > 1:
                    fig_st_reg = px.bar(
                        st_counts, x='State', y='Patients (N)',
                        title=f"State Distribution in {selected_region}",
                        color='Patients (N)',
                        color_continuous_scale='Greens'
                    )
                    fig_st_reg.update_layout(
                        height=260, margin=dict(l=10, r=10, t=30, b=10),
                        title_font=dict(size=12, color='#0F172A'),
                        xaxis=dict(tickfont=dict(size=11)), yaxis=dict(tickfont=dict(size=11))
                    )
                    st.plotly_chart(fig_st_reg, use_container_width=True)
                else:
                    st.caption(f"All {tot_cohort} samples in `{selected_region}` are from State: **{st_counts.iloc[0]['State']}**.")
        else:
            st.markdown(f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin: 12px 0 6px 0;'>🗺️ Regional Cross-Comparison for {combo_title} Across India</div>", unsafe_allow_html=True)
            if 'Region' in processed_df.columns:
                pattern_series = processed_df[active_snps].apply(lambda r: " + ".join(r.astype(str)), axis=1)
                pivot_reg = pd.crosstab(processed_df['Region'], pattern_series, margins=True, margins_name="Total")
                st.dataframe(pivot_reg, use_container_width=True)

                reg_agg = processed_df.groupby(['Region'] + active_snps).size().reset_index(name='Patients (N)')
                reg_agg['Pattern'] = reg_agg[active_snps].apply(lambda r: " + ".join(r.astype(str)), axis=1)
                fig_reg_bar = px.bar(
                    reg_agg, x='Region', y='Patients (N)', color='Pattern',
                    barmode='stack', title=f"Regional Composition of {combo_title} (Stacked)",
                    color_discrete_sequence=['#166534', '#1E40AF', '#D97706', '#991B1B', '#0D9488', '#6B21A8']
                )
                fig_reg_bar.update_layout(
                    height=280, margin=dict(l=10, r=10, t=30, b=10),
                    title_font=dict(size=12, color='#0F172A'),
                    xaxis=dict(tickfont=dict(size=11)), yaxis=dict(tickfont=dict(size=11))
                )
                st.plotly_chart(fig_reg_bar, use_container_width=True)

    elif "State-Wise Profiles" in geo_scope:
        if selected_state != "All States (Comparative Overview)":
            st.markdown(f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin: 12px 0 6px 0;'>📍 State Profile: {selected_state} (N = {tot_cohort})</div>", unsafe_allow_html=True)
            if 'Gender_Clean' in active_df.columns or 'Gender' in active_df.columns:
                g_col = 'Gender_Clean' if 'Gender_Clean' in active_df.columns else 'Gender'
                g_counts = active_df[g_col].value_counts().reset_index()
                g_counts.columns = ['Gender', 'Count']
                fig_g = px.pie(g_counts, names='Gender', values='Count', hole=0.45, title=f"Gender Distribution in {selected_state}")
                fig_g.update_layout(height=240, margin=dict(l=10, r=10, t=30, b=10), title_font=dict(size=12))
                st.plotly_chart(fig_g, use_container_width=True)
        else:
            st.markdown(f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin: 12px 0 6px 0;'>📍 State Cross-Tabulation & Top States by Cohort Size</div>", unsafe_allow_html=True)
            if state_col and state_col in processed_df.columns:
                pattern_series = processed_df[active_snps].apply(lambda r: " + ".join(r.astype(str)), axis=1)
                pivot_st = pd.crosstab(processed_df[state_col], pattern_series, margins=True, margins_name="Total")
                st.dataframe(pivot_st, use_container_width=True)

                st_counts = processed_df[state_col].value_counts().head(10).reset_index()
                st_counts.columns = ['State', 'Total Patients']
                fig_st = px.bar(
                    st_counts, x='State', y='Total Patients',
                    title="Top 10 Indian States by Cohort Size",
                    color='Total Patients',
                    color_continuous_scale='Greens'
                )
                fig_st.update_layout(
                    height=280, margin=dict(l=10, r=10, t=30, b=10),
                    title_font=dict(size=12, color='#0F172A'),
                    xaxis=dict(tickfont=dict(size=11)), yaxis=dict(tickfont=dict(size=11))
                )
                st.plotly_chart(fig_st, use_container_width=True)

    else:
        if 'Region' in processed_df.columns:
            reg_counts = processed_df['Region'].value_counts().reset_index()
            reg_counts.columns = ['Region', 'Patients (N)']
            fig_reg_ov = px.bar(
                reg_counts, x='Region', y='Patients (N)',
                title="Pan-Indian Regional Population Distribution",
                color='Region',
                color_discrete_sequence=px.colors.qualitative.Set2
            )
            fig_reg_ov.update_layout(
                height=250, margin=dict(l=10, r=10, t=30, b=10),
                title_font=dict(size=12, color='#0F172A'),
                xaxis=dict(tickfont=dict(size=11)), yaxis=dict(tickfont=dict(size=11))
            )
            st.plotly_chart(fig_reg_ov, use_container_width=True)

    # -------------------------------------------------------------------------
    # 6. STATISTICAL COMBINATION TABLE (STANDARDIZED SCIENTIFIC HEADERS)
    # -------------------------------------------------------------------------
    st.markdown(f"<div style='font-size: 13px; font-weight: 600; color: #1E293B; margin: 14px 0 6px 0;'>📋 Statistical Combination Table — {combo_title} ({scope_label})</div>", unsafe_allow_html=True)

    table_display = combo_grp.copy()
    if is_sci:
        rename_dict = {
            'CYP2C19*2': 'rs4244285 (*2, c.681G>A)',
            'CYP2C19*3': 'rs4986893 (*3, c.636G>A)',
            'CYP2C19*17': 'rs12248560 (*17, c.-806C>T)',
            'Combination Pattern': 'Multi-Locus Pattern',
            'Sample Count (N)': 'Samples (N)',
            'CPIC Phenotype(s)': 'CPIC Clinical Phenotype'
        }
        table_display = table_display.rename(columns=rename_dict)
    else:
        rename_dict = {
            'CYP2C19*2': 'CYP2C19*2 (c.681G>A)',
            'CYP2C19*3': 'CYP2C19*3 (c.636G>A)',
            'CYP2C19*17': 'CYP2C19*17 (c.-806C>T)',
            'Combination Pattern': 'Genotype Combination',
            'Sample Count (N)': 'Samples (N)',
            'CPIC Phenotype(s)': 'CPIC Clinical Phenotype'
        }
        table_display = table_display.rename(columns=rename_dict)

    disp_cols = [c for c in table_display.columns if c != 'pct_val']
    st.dataframe(table_display[disp_cols], use_container_width=True, hide_index=True)

# -----------------------------------------------------------------------------
# REUSABLE UNIFIED RESULTS RENDERER (MATCHING SCREENSHOTS 1, 2, 3, 4)
# -----------------------------------------------------------------------------
def render_unified_results(full_results: dict, qc_report: dict):
    if not has_active_data() or qc_report is None:
        st.info("ℹ️ Workspace is empty (0 samples). All existing data was deleted. Upload a new Excel/CSV dataset file or add sample records on Page 1.")
        return

    if len(st.session_state.get('datasets_dict', {})) > 1:
        st.markdown("### 📂 Multi-Dataset Profile View")
        ds_list = list(st.session_state['datasets_dict'].keys())
        sel_ds_name = st.selectbox("Select Active Dataset Profile for Separate Calculation:", ds_list, key="sel_multi_ds_switcher")
        if sel_ds_name and st.session_state.get('active_ds_name') != sel_ds_name:
            st.session_state['active_ds_name'] = sel_ds_name
            st.session_state['uploaded_df'] = st.session_state['datasets_dict'][sel_ds_name]
            st.rerun()
        st.divider()

    active_gene = full_results.get('active_gene') or full_results.get('target_gene') or st.session_state.get('selected_gene', 'CYP2C19')
    ov = full_results['overall']
    tot_samples = qc_report.get('total_samples', ov['sample_count'])
    fully_valid = qc_report.get('fully_valid_samples', ov['sample_count'])
    dup_count = qc_report.get('duplicate_sample_count', 0)
    dup_ids = qc_report.get('duplicate_sample_ids', [])
    invalid_count = qc_report.get('invalid_genotypes_count', 0)
    invalid_log = qc_report.get('invalid_genotypes_log', {})

    # 9 Main Tabs Consolidating All Analysis Views
    tab_overview, tab_gen, tab_hwe, tab_pop, tab_pgx, tab_qc, tab_geo, tab_gender, tab_state = st.tabs([
        "📊 Executive Overview",
        "🧬 Genotype & Allele", 
        "⚖️ Hardy-Weinberg", 
        "🌐 Population Stratification", 
        f"💊 {active_gene} Calling",
        "🛡️ Data Quality Audit",
        "🗺️ Geographic Groups",
        "👫 Gender-Wise",
        "📍 Dedicated State Profiles"
    ])

    # Cohort Dictionary for Sub-cohort Filtering
    cohort_dict = {
        'Overall': full_results.get('overall', {}),
        'South India': full_results.get('regional', {}).get('South India', {}),
        'North India': full_results.get('regional', {}).get('North India', {}),
        'East India': full_results.get('regional', {}).get('East India', {}),
        'West India': full_results.get('regional', {}).get('West India', {}),
        'Central India': full_results.get('regional', {}).get('Central India', {}),
        'Female': full_results.get('gender', {}).get('Female', {}),
        'Male': full_results.get('gender', {}).get('Male', {})
    }

    # Helper for Sub-cohort & State-Wise Filter UI
    def render_subcohort_pills(key_prefix):
        st_dict = full_results.get('state_wise', {})
        pill_options = [f"{k} (n={v.get('sample_count', 0)})" for k, v in cohort_dict.items() if v]
        
        c_pill, c_state = st.columns([1.1, 0.9])
        with c_pill:
            selected_pill = st.radio("Select Sub-Cohort Profile:", pill_options, horizontal=True, key=f"{key_prefix}_subcohort_pills")
        
        state_items = sorted([(s, data.get('sample_count', 0)) for s, data in st_dict.items() if data.get('sample_count', 0) > 0])
        state_options = ["(All States / Default)"] + [f"{s} (n={cnt})" for s, cnt in state_items]
        
        with c_state:
            selected_state_str = st.selectbox("🗺️ Select State-Wise Profile Data:", state_options, key=f"{key_prefix}_state_select")
            
        if selected_state_str != "(All States / Default)":
            state_name = selected_state_str.split(" (n=")[0]
            sel_key = f"State: {state_name}"
            sel_cohort = st_dict.get(state_name, full_results['overall'])
        else:
            sel_key = selected_pill.split(' (')[0]
            sel_cohort = cohort_dict.get(sel_key, full_results['overall'])
            
        return sel_key, sel_cohort

    # -------------------------------------------------------------------------
    # TAB 1: EXECUTIVE OVERVIEW
    # -------------------------------------------------------------------------
    with tab_overview:
        missing_count = qc_report.get('samples_with_missing_count', 0)
        missing_pct = qc_report.get('overall_missing_pct', 0.0)

        k1, k2, k3, k4, k5 = st.columns(5)
        k1.metric("TOTAL SAMPLES", f"{tot_samples:,}")
        k2.metric("FULLY VALID", f"{fully_valid:,}")
        k3.metric("MISSING GENOTYPES", f"{missing_pct:.1f}% ({missing_count} samples)")
        k4.metric("DUPLICATE IDS", dup_count)
        k5.metric("INVALID GENOTYPE CALLS", invalid_count)

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
            phe_dict = pgx.get('phenotypes', {})
            phe_df = pd.DataFrame([{'Phenotype Category': k, 'Count': v['count'], 'Percentage': f"{v['percentage']:.1f}%"} for k, v in phe_dict.items() if v['count'] > 0])
            st.table(phe_df)

            fig_ov_phe = px.pie(
                names=[k for k, v in phe_dict.items() if v['count'] > 0],
                values=[v['count'] for k, v in phe_dict.items() if v['count'] > 0],
                hole=0.45,
                title="Overall Phenotype Distribution",
                color_discrete_sequence=['#166534', '#1E40AF', '#4B5563', '#991B1B', '#D97706', '#0D9488']
            )
            fig_ov_phe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_ov_phe, use_container_width=True)

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

    # -------------------------------------------------------------------------
    # TAB 2: GENOTYPE & ALLELE
    # -------------------------------------------------------------------------
    with tab_gen:
        sel_key, sel_cohort = render_subcohort_pills("tab_gen")
        snps_data = sel_cohort.get('snps', {})
        
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

        snp_display_order = [
            ('CYP2C19*2', 'rs4244285 (CYP2C19*2, c.681G>A)'),
            ('CYP2C19*3', 'rs4986893 (CYP2C19*3, c.636G>A)'),
            ('CYP2C19*17', 'rs12248560 (CYP2C19*17, c.-806C>T)'),
            ('CYP2C9*2', 'rs1799853 (CYP2C9*2, c.430C>T)'),
            ('CYP2C9*3', 'rs1057910 (CYP2C9*3, c.1075A>C)')
        ]
        for snp_code in snps_data.keys():
            if snp_code not in [s[0] for s in snp_display_order]:
                snp_display_order.append((snp_code, f"{snp_code} Marker"))

        for snp_code, snp_title in snp_display_order:
            if snp_code in snps_data:
                s_res = snps_data[snp_code]
                v_cnt = s_res.get('valid_samples', 0)
                inv_cnt = qc_report['snp_qc_stats'].get(snp_code, {}).get('invalid_samples', 0)
                m_cnt = qc_report['snp_qc_stats'].get(snp_code, {}).get('missing_samples', 0)
                
                st.markdown(f"#### {snp_title}")
                st.caption(f"Valid: {v_cnt:,} · Missing: {m_cnt} · Invalid calls: {inv_cnt}")
                
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

    # -------------------------------------------------------------------------
    # TAB 3: HARDY-WEINBERG
    # -------------------------------------------------------------------------
    with tab_hwe:
        sel_key, sel_cohort = render_subcohort_pills("tab_hwe")
        snps_data = sel_cohort.get('snps', {})
        
        st.markdown(f"### Hardy-Weinberg Equilibrium Audit ({sel_key} Cohort)")
        hwe_table_rows = []
        for snp_name, s_res in snps_data.items():
            hw = s_res.get('hwe', {})
            gc = s_res.get('genotype_counts', {})
            exp_c = hw.get('expected_counts', {})
            chi2_val = hw.get('chi2_stat', 0.0)
            hwe_table_rows.append({
                'SNP Variant': snp_name,
                'Valid Samples (N)': s_res.get('valid_samples', 0),
                'Observed Genotypes': " / ".join([f"{k}:{v}" for k, v in gc.items()]),
                'Expected Genotypes': " / ".join([f"{k}:{exp_c.get(k, 0)}" for k in gc.keys()]),
                '_raw_chi2': chi2_val,
                'Degrees of Freedom': hw.get('df', 1),
                'P-Value': f"{hw.get('p_value', 1.0):.4f}",
                'Exact Test P-Value': f"{hw.get('exact_p_value', 1.0):.4f}"
            })
            
        c_hwe_tbl, c_hwe_fig = st.columns([1.1, 0.9])
        with c_hwe_tbl:
            HWE_HEADERS = [
                'SNP Variant', 'Valid Samples (N)', 'Observed Genotypes', 
                'Expected Genotypes', 'Chi-Square (χ²)', 'df', 'P-Value', 'Exact P-Val'
            ]
            if not hwe_table_rows:
                hwe_df = pd.DataFrame(columns=HWE_HEADERS)
            else:
                hwe_df = pd.DataFrame([{
                    'SNP Variant': r['SNP Variant'],
                    'Valid Samples (N)': r['Valid Samples (N)'],
                    'Observed Genotypes': r['Observed Genotypes'],
                    'Expected Genotypes': r['Expected Genotypes'],
                    'Chi-Square (χ²)': f"⚠️ {r['_raw_chi2']:.4f} (Chi² ≥ 3.4)" if r['_raw_chi2'] >= 3.4 else f"{r['_raw_chi2']:.4f}",
                    'df': r['Degrees of Freedom'],
                    'P-Value': r['P-Value'],
                    'Exact P-Val': r['Exact Test P-Value']
                } for r in hwe_table_rows])

            def style_hwe_chi2(val):
                if '⚠️' in str(val) or 'Chi²' in str(val):
                    return 'background-color: #FEE2E2; color: #DC2626; font-weight: bold; border-radius: 4px;'
                return 'color: #166534; font-weight: 600;'

            if not hwe_df.empty and 'Chi-Square (χ²)' in hwe_df.columns:
                styled_hwe = hwe_df.style.map(style_hwe_chi2, subset=['Chi-Square (χ²)'])
                st.dataframe(styled_hwe, use_container_width=True, hide_index=True)
            else:
                st.dataframe(hwe_df, use_container_width=True, hide_index=True)

        with c_hwe_fig:
            chi2_rows = []
            for snp_name, s_res in snps_data.items():
                raw_c = s_res.get('hwe', {}).get('chi2_stat', 0.0)
                try:
                    c_num = float(raw_c) if (raw_c is not None and str(raw_c).strip() not in ['', 'nan', 'none', '-']) else 0.0
                except (ValueError, TypeError):
                    c_num = 0.0
                chi2_rows.append({'SNP': str(snp_name), 'Chi2 Stat': c_num})

            chi2_chart_df = pd.DataFrame(chi2_rows, columns=['SNP', 'Chi2 Stat'])
            if not chi2_chart_df.empty and len(chi2_rows) > 0:
                fig_hwe = px.bar(
                    chi2_chart_df, x='SNP', y='Chi2 Stat',
                    title=f"HWE Chi-Square (χ²) Stats vs Threshold 3.4 ({sel_key})",
                    color='SNP',
                    color_discrete_sequence=['#2563EB', '#7C3AED', '#DB2777']
                )
                fig_hwe.add_hline(y=3.4, line_dash="dash", line_color="red", annotation_text="Critical Threshold (3.4)")
                fig_hwe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_hwe, use_container_width=True)
            else:
                st.info("No Chi-Square data available to plot for this cohort.")
        
        st.markdown("##### Table 3. Chi-Square Distribution Reference Table (df = 1, α = 0.05, Critical = 3.4)")
        chi2_ref_df = pd.DataFrame([{
            'Level of Significance (α)': '0.50',
            '0.10': '2.706',
            '0.05 (Critical Threshold)': '3.4',
            '0.02': '5.412',
            '0.01': '6.635',
            '0.001': '10.827'
        }])
        st.table(chi2_ref_df)
        st.info("📌 **HWE Decision Rule:** If $\\chi^2 \\ge 3.4$, population deviates from Hardy-Weinberg Equilibrium.")

    # -------------------------------------------------------------------------
    # TAB 4: POPULATION STRATIFICATION
    # -------------------------------------------------------------------------
    with tab_pop:
        st.markdown("### Group sizes")
        reg_dict = full_results.get('regional', {})
        gen_dict = full_results.get('gender', {})
        tot_n = ov['sample_count']
        
        group_rows = [
            {'GROUP': 'Overall', 'N': tot_n, '% OF TOTAL': '100.0%'},
            {'GROUP': 'Region — South India', 'N': reg_dict.get('South India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('South India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — North India', 'N': reg_dict.get('North India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('North India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — East India', 'N': reg_dict.get('East India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('East India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — West India', 'N': reg_dict.get('West India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('West India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — Central India', 'N': reg_dict.get('Central India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('Central India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Gender — Female', 'N': gen_dict.get('Female', {}).get('sample_count', 0), '% OF TOTAL': f"{(gen_dict.get('Female', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Gender — Male', 'N': gen_dict.get('Male', {}).get('sample_count', 0), '% OF TOTAL': f"{(gen_dict.get('Male', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
        ]
        
        c_pop1, c_pop2 = st.columns([1.0, 1.0])
        with c_pop1:
            st.table(pd.DataFrame(group_rows))

        with c_pop2:
            pop_chart_df = pd.DataFrame([
                {'Group': r['GROUP'], 'N': r['N']}
                for r in group_rows if r['GROUP'] != 'Overall'
            ])
            fig_pop = px.bar(
                pop_chart_df, x='N', y='Group', orientation='h',
                title="Population Subgroup Sizes (N)",
                color='Group',
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            fig_pop.update_layout(height=320, margin=dict(l=20, r=20, t=35, b=20), yaxis=dict(autorange="reversed"))
            st.plotly_chart(fig_pop, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 5: CYP2C19 CALLING
    # -------------------------------------------------------------------------
    with tab_pgx:
        sel_key, sel_cohort = render_subcohort_pills("tab_pgx")
        pgx = sel_cohort.get('pgx_summary') or sel_cohort.get('cyp2c9_summary') or sel_cohort.get('cyp2c19_summary', {})
        dips = pgx.get('diplotypes', {})
        phenos = pgx.get('phenotypes', {})

        st.markdown("<br>", unsafe_allow_html=True)
        view_choice = st.radio(
            "🔎 Select View to Display:",
            [
                "👁️ Show All Views",
                "🧩 Combination",
                "📋 Diplotype distribution",
                "📊 Diplotype counts",
                "💊 Metabolizer / phenotype classification",
                "🗺️ Geographical & State-Wise Matrices"
            ],
            horizontal=True,
            key="tab_pgx_view_selector"
        )
        st.divider()

        # 1. Diplotype Distribution & Diplotype Counts
        if view_choice in ["👁️ Show All Views", "📋 Diplotype distribution", "📊 Diplotype counts"]:
            if view_choice == "📋 Diplotype distribution":
                st.markdown("#### Diplotype distribution")
                dip_rows = [{'DIPLOTYPE': k, 'COUNT': v['count'], 'FREQUENCY': f"{v['percentage']:.1f}%"} for k, v in dips.items() if v['count'] > 0]
                st.table(pd.DataFrame(dip_rows))
            elif view_choice == "📊 Diplotype counts":
                st.markdown("#### Diplotype counts")
                dip_plot_rows = [{'Diplotype': str(k), 'Count': int(v['count'])} for k, v in dips.items() if v.get('count', 0) > 0]
                dip_fig_df = pd.DataFrame(dip_plot_rows, columns=['Diplotype', 'Count'])
                if not dip_fig_df.empty and len(dip_plot_rows) > 0:
                    fig_dip = px.bar(dip_fig_df, x='Count', y='Diplotype', orientation='h', title="Diplotype counts", color_discrete_sequence=['#166534'])
                    fig_dip.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20), yaxis=dict(autorange="reversed"))
                    st.plotly_chart(fig_dip, use_container_width=True)
                else:
                    st.info("No diplotype counts to display.")
            else:
                c_dip, c_dip_chart = st.columns([1.1, 0.9])
                with c_dip:
                    st.markdown("#### Diplotype distribution")
                    dip_rows = [{'DIPLOTYPE': k, 'COUNT': v['count'], 'FREQUENCY': f"{v['percentage']:.1f}%"} for k, v in dips.items() if v['count'] > 0]
                    st.table(pd.DataFrame(dip_rows))
                with c_dip_chart:
                    dip_plot_rows = [{'Diplotype': str(k), 'Count': int(v['count'])} for k, v in dips.items() if v.get('count', 0) > 0]
                    dip_fig_df = pd.DataFrame(dip_plot_rows, columns=['Diplotype', 'Count'])
                    if not dip_fig_df.empty and len(dip_plot_rows) > 0:
                        fig_dip = px.bar(dip_fig_df, x='Count', y='Diplotype', orientation='h', title="Diplotype counts", color_discrete_sequence=['#166534'])
                        fig_dip.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20), yaxis=dict(autorange="reversed"))
                        st.plotly_chart(fig_dip, use_container_width=True)
                    else:
                        st.info("No diplotype counts to display.")

        # 2. Metabolizer / Phenotype Classification
        if view_choice in ["👁️ Show All Views", "💊 Metabolizer / phenotype classification"]:
            if view_choice == "👁️ Show All Views":
                st.markdown("<br>", unsafe_allow_html=True)
            c_phe, c_phe_chart = st.columns([1.1, 0.9])
            with c_phe:
                st.markdown("#### Metabolizer / phenotype classification")
                phe_rows = [{'PHENOTYPE': k, 'COUNT': v['count'], 'PERCENTAGE': f"{v['percentage']:.1f}%"} for k, v in phenos.items() if v['count'] > 0]
                st.table(pd.DataFrame(phe_rows))
            with c_phe_chart:
                phe_labels = [str(k) for k, v in phenos.items() if v.get('count', 0) > 0]
                phe_counts = [int(v['count']) for k, v in phenos.items() if v.get('count', 0) > 0]
                if phe_labels and phe_counts and sum(phe_counts) > 0:
                    fig_donut = px.pie(
                        names=phe_labels, 
                        values=phe_counts, 
                        hole=0.5, 
                        title="Phenotype distribution", 
                        color_discrete_sequence=['#166534', '#1E40AF', '#4B5563', '#991B1B', '#D97706', '#0D9488']
                    )
                    fig_donut.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_donut, use_container_width=True)
                else:
                    st.info("No phenotype calls to display.")

        # 3. Combination & Geographic Studio
        if view_choice in ["👁️ Show All Views", "🧩 Combination", "🗺️ Geographical & State-Wise Matrices"]:
            render_cyp2c19_combination_studio(processed_df, full_results, view_mode=view_choice)

        # 4. Clinical CPIC Dosing Guidelines
        if active_gene == 'CYP2C9':
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown("### 💊 CPIC® Clinical Dosing Guidelines for CYP2C9")
            st.caption("Standardized CPIC level A/B clinical recommendations for narrow therapeutic index substrates.")
            c_g1, c_g2, c_g3 = st.columns(3)
            with c_g1:
                st.markdown("""
                <div style="background: #FEF3C7; border: 1px solid #F59E0B; border-radius: 8px; padding: 14px;">
                    <div style="font-weight: 700; color: #92400E; font-size: 14px;">🩸 Warfarin</div>
                    <div style="font-size: 12px; color: #78350F; margin-top: 6px; line-height: 1.5;">
                        <b>Normal (*1/*1):</b> Standard starting dose (titrate INR).<br>
                        <b>Intermediate (*1/*2, *1/*3, *2/*2):</b> Reduce initial dose by 20–50% to prevent supratherapeutic INR and bleeding.<br>
                        <b>Poor (*2/*3, *3/*3):</b> Reduce dose by 50–80% or consider alternative oral anticoagulant (e.g., DOAC).
                    </div>
                </div>
                """, unsafe_allow_html=True)
            with c_g2:
                st.markdown("""
                <div style="background: #EFF6FF; border: 1px solid #3B82F6; border-radius: 8px; padding: 14px;">
                    <div style="font-weight: 700; color: #1E40AF; font-size: 14px;">🧠 Phenytoin</div>
                    <div style="font-size: 12px; color: #1E3A8A; margin-top: 6px; line-height: 1.5;">
                        <b>Normal (*1/*1):</b> Standard starting dose (300 mg/day).<br>
                        <b>Intermediate (*1/*2, *1/*3, *2/*2):</b> Reduce initial maintenance dose by 25% to prevent neurotoxicity.<br>
                        <b>Poor (*2/*3, *3/*3):</b> Reduce initial maintenance dose by 50% and closely monitor serum trough.
                    </div>
                </div>
                """, unsafe_allow_html=True)
            with c_g3:
                st.markdown("""
                <div style="background: #ECFDF5; border: 1px solid #10B981; border-radius: 8px; padding: 14px;">
                    <div style="font-weight: 700; color: #065F46; font-size: 14px;">🩹 Celecoxib</div>
                    <div style="font-size: 12px; color: #064E3B; margin-top: 6px; line-height: 1.5;">
                        <b>Normal (*1/*1):</b> Standard starting dose (200 mg/day).<br>
                        <b>Intermediate (*1/*2, *1/*3, *2/*2):</b> Standard dose with cautious clinical titration.<br>
                        <b>Poor (*2/*3, *3/*3):</b> Reduce initial dose by 50% to prevent prolonged exposure and toxicity.
                    </div>
                </div>
                """, unsafe_allow_html=True)
        elif active_gene == 'CYP2C19':
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown("### 💊 CPIC® Clinical Dosing Guidelines for CYP2C19")
            st.caption("Standardized CPIC level A/B clinical recommendations for antiplatelet and psychiatric therapies.")
            c_g1, c_g2 = st.columns(2)
            with c_g1:
                st.markdown("""
                <div style="background: #EFF6FF; border: 1px solid #3B82F6; border-radius: 8px; padding: 14px;">
                    <div style="font-weight: 700; color: #1E40AF; font-size: 14px;">🩸 Clopidogrel (Plavix)</div>
                    <div style="font-size: 12px; color: #1E3A8A; margin-top: 6px; line-height: 1.5;">
                        <b>Normal (*1/*1) & Rapid (*1/*17):</b> Standard dose (75 mg/day).<br>
                        <b>Intermediate (*1/*2, *1/*3) & Poor (*2/*2, *2/*3, *3/*3):</b> Avoid clopidogrel; prescribe alternative antiplatelet (Prasugrel, Ticagrelor) due to reduced active metabolite and high stent thrombosis risk.
                    </div>
                </div>
                """, unsafe_allow_html=True)
            with c_g2:
                st.markdown("""
                <div style="background: #FDF2F8; border: 1px solid #EC4899; border-radius: 8px; padding: 14px;">
                    <div style="font-weight: 700; color: #9D174D; font-size: 14px;">💊 Voriconazole & Escitalopram</div>
                    <div style="font-size: 12px; color: #831843; margin-top: 6px; line-height: 1.5;">
                        <b>Voriconazole:</b> In Poor Metabolizers, select alternative antifungal (Isavuconazole) to prevent neurotoxicity.<br>
                        <b>Escitalopram:</b> In Poor Metabolizers, reduce starting dose by 50% to avoid QT prolongation.
                    </div>
                </div>
                """, unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # TAB 6: DATA QUALITY AUDIT
    # -------------------------------------------------------------------------
    with tab_qc:
        st.markdown("### 🛡️ Data Quality Control (QC) & Missing Values Audit")
        q1, q2, q3, q4, q5 = st.columns(5)
        q1.metric("Total Cohort Uploaded", qc_report['total_samples'])
        q2.metric("Fully Valid Samples", qc_report['fully_valid_samples'])
        q3.metric("Duplicate Sample IDs", qc_report['duplicate_sample_count'])
        q4.metric("Missing Genotypes", f"{qc_report.get('samples_with_missing_count', 0)} ({qc_report.get('overall_missing_pct', 0.0):.1f}%)")
        q5.metric("Invalid Genotype Calls", qc_report['invalid_genotypes_count'])
        
        if qc_report['invalid_genotypes_count'] > 0:
            st.error(f"⚠️ **Warning: {qc_report['invalid_genotypes_count']} Invalid Genotype Call(s) Detected!** Invalid genotype entries were flagged and quarantined.")
            inv_log_rows = []
            for snp, samples in qc_report['invalid_genotypes_log'].items():
                for item in samples:
                    inv_log_rows.append({'SNP Variant': snp, 'Sample ID': item['sample_id'], 'Raw Invalid Genotype Value': item['raw_value'], 'Audit Action': 'Quarantined & Excluded'})
            st.table(pd.DataFrame(inv_log_rows))

        c_qc_tbl, c_qc_fig = st.columns([1.1, 0.9])
        with c_qc_tbl:
            st.markdown("#### Per-SNP Quality Control & Missing Values Table")
            qc_rows = []
            for snp, sdata in qc_report['snp_qc_stats'].items():
                qc_rows.append({
                    'SNP Variant': snp,
                    'Valid Samples': sdata['valid_samples'],
                    'Missing Genotypes Count': sdata['missing_samples'],
                    'Invalid Genotypes Count': sdata.get('invalid_samples', 0),
                    'Missing Data Rate (%)': f"{sdata['missing_pct']:.1f}%",
                    'Quality Audit Status': 'PASS' if sdata.get('invalid_samples', 0) == 0 else 'WARNING'
                })
            st.table(pd.DataFrame(qc_rows))

        with c_qc_fig:
            qc_plot_rows = [
                {'SNP': str(snp), 'Valid': int(sdata.get('valid_samples', 0)), 'Missing': int(sdata.get('missing_samples', 0)), 'Invalid': int(sdata.get('invalid_samples', 0))}
                for snp, sdata in qc_report.get('snp_qc_stats', {}).items()
            ]
            qc_fig_df = pd.DataFrame(qc_plot_rows, columns=['SNP', 'Valid', 'Missing', 'Invalid'])
            if not qc_fig_df.empty and len(qc_plot_rows) > 0:
                fig_qc = px.bar(
                    qc_fig_df, x='SNP', y=['Valid', 'Missing', 'Invalid'],
                    title="Sample Quality Audit Breakdown per SNP",
                    barmode='stack',
                    color_discrete_sequence=['#16A34A', '#EAB308', '#DC2626']
                )
                fig_qc.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_qc, use_container_width=True)
            else:
                st.info("No QC chart data available.")

    # -------------------------------------------------------------------------
    # TAB 7: GEOGRAPHIC POPULATION GROUPS
    # -------------------------------------------------------------------------
    with tab_geo:
        st.markdown("### 🌐 Geographic Population Groups Analysis")
        reg = full_results['regional']
        regions = ['South India', 'North India', 'East India', 'West India', 'Central India']
        
        cols = st.columns(5)
        for idx, r in enumerate(regions):
            r_data = reg.get(r, {})
            snps_r = r_data.get('snps', {})
            c2_maf = snps_r.get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
            c3_maf = snps_r.get('CYP2C19*3', {}).get('allele_freqs', {}).get('A', 0)
            c17_maf = snps_r.get('CYP2C19*17', {}).get('allele_freqs', {}).get('T', 0)
            pm_val = r_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer', r_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {})).get('percentage', 0.0)

            with cols[idx]:
                st.markdown(f"""
                <div class="kpi-card">
                    <div class="kpi-lbl">{r}</div>
                    <div class="kpi-val">{r_data.get('sample_count', 0):,}</div>
                    <div style="font-size: 0.78rem; color: #475569; margin-top: 0.3rem; line-height: 1.4;">
                        *2 MAF (A): <strong>{c2_maf}</strong><br>
                        *3 MAF (A): <strong>{c3_maf}</strong><br>
                        *17 MAF (T): <strong>{c17_maf}</strong><br>
                        PM: <strong style="color: #EF4444;">{pm_val:.1f}%</strong>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
        st.markdown("<br>", unsafe_allow_html=True)

        c_geo1, c_geo2 = st.columns(2)
        with c_geo1:
            geo_marker_data = []
            for r in regions:
                r_snps = reg.get(r, {}).get('snps', {})
                geo_marker_data.append({'Geographic Group': r, 'Biomarker': 'CYP2C19*2 (f_A)', 'Minor Allele Frequency': float(r_snps.get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0))})
                geo_marker_data.append({'Geographic Group': r, 'Biomarker': 'CYP2C19*3 (f_A)', 'Minor Allele Frequency': float(r_snps.get('CYP2C19*3', {}).get('allele_freqs', {}).get('A', 0))})
                geo_marker_data.append({'Geographic Group': r, 'Biomarker': 'CYP2C19*17 (f_T)', 'Minor Allele Frequency': float(r_snps.get('CYP2C19*17', {}).get('allele_freqs', {}).get('T', 0))})
            geo_df = pd.DataFrame(geo_marker_data, columns=['Geographic Group', 'Biomarker', 'Minor Allele Frequency'])
            if not geo_df.empty:
                fig_geo_all = px.bar(
                    geo_df,
                    x='Geographic Group', y='Minor Allele Frequency',
                    color='Biomarker', barmode='group',
                    title="All Biomarkers (CYP2C19*2, *3, *17) MAF Across Geographic Groups",
                    color_discrete_map={'CYP2C19*2 (f_A)': '#2563EB', 'CYP2C19*3 (f_A)': '#059669', 'CYP2C19*17 (f_T)': '#D97706'}
                )
                fig_geo_all.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_geo_all, use_container_width=True)

        with c_geo2:
            geo_pheno_data = []
            for r in regions:
                p_dict = reg.get(r, {}).get('cyp2c19_summary', {}).get('phenotypes', {})
                pm_val = p_dict.get('Poor Metabolizer', p_dict.get('Poor Metabolizer (PM)', {})).get('percentage', 0.0)
                im_val = p_dict.get('Intermediate Metabolizer', p_dict.get('Intermediate Metabolizer (IM)', {})).get('percentage', 0.0)
                try:
                    pm_pct = float(pm_val)
                except (ValueError, TypeError):
                    pm_pct = 0.0
                try:
                    im_pct = float(im_val)
                except (ValueError, TypeError):
                    im_pct = 0.0
                geo_pheno_data.append({'Geographic Group': r, 'Phenotype': 'Poor Metabolizer (PM %)', 'Percentage': pm_pct})
                geo_pheno_data.append({'Geographic Group': r, 'Phenotype': 'Intermediate Metabolizer (IM %)', 'Percentage': im_pct})
            geo_pm_df = pd.DataFrame(geo_pheno_data, columns=['Geographic Group', 'Phenotype', 'Percentage'])
            if not geo_pm_df.empty:
                fig_geo_pm = px.bar(
                    geo_pm_df,
                    x='Geographic Group', y='Percentage',
                    color='Phenotype', barmode='group',
                    title="Clinical Phenotypes (PM & IM %) Across Geographic Groups",
                    color_discrete_map={'Poor Metabolizer (PM %)': '#DC2626', 'Intermediate Metabolizer (IM %)': '#D97706'}
                )
                fig_geo_pm.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_geo_pm, use_container_width=True)

        st.markdown("#### Comprehensive Population Statistical Summary Matrix")
        
        t_cyp2, t_cyp3, t_cyp17 = st.tabs([
            "🧬 CYP2C19*2 (rs4244285)", 
            "🧬 CYP2C19*3 (rs4986893)", 
            "🧬 CYP2C19*17 (rs12248560)"
        ])

        def render_marker_geo_table(snp_key):
            snp_conf_sel = SNP_CONFIG.get(snp_key, {})
            ref_al = snp_conf_sel.get('ref_allele', 'G')
            var_al = snp_conf_sel.get('var_allele', 'A')
            wt_g = snp_conf_sel.get('wildtype_genotype', 'GG')
            het_g = snp_conf_sel.get('het_genotype', 'GA')
            var_g = snp_conf_sel.get('hom_var_genotype', 'AA')

            rows_cfg = [
                (f'{var_al} Count (Variant Allele)', lambda s: s.get(snp_key, {}).get('allele_counts', {}).get(var_al, 0)),
                (f'{ref_al} Count (Reference Allele)', lambda s: s.get(snp_key, {}).get('allele_counts', {}).get(ref_al, 0)),
                ('Total Allele Count (2N)', lambda s: sum(s.get(snp_key, {}).get('allele_counts', {}).values())),
                (f'{het_g} COUNT (Heterozygote)', lambda s: s.get(snp_key, {}).get('genotype_counts', {}).get(het_g, 0) if het_g else 0),
                (f'{var_g} COUNT (Homozygous Variant)', lambda s: s.get(snp_key, {}).get('genotype_counts', {}).get(var_g, 0) if var_g else 0),
                (f'{wt_g} COUNT (Homozygous Wildtype)', lambda s: s.get(snp_key, {}).get('genotype_counts', {}).get(wt_g, 0)),
                ('Total Sample N', lambda s: sum(s.get(snp_key, {}).get('genotype_counts', {}).values())),
                (f'Allele Frequency F({var_al})', lambda s: s.get(snp_key, {}).get('allele_freqs', {}).get(var_al, 0)),
                (f'Allele Frequency F({ref_al})', lambda s: s.get(snp_key, {}).get('allele_freqs', {}).get(ref_al, 0)),
                ('Chi-Square (χ²)', lambda s: s.get(snp_key, {}).get('hwe', {}).get('chi2_stat', 0)),
                ('P VALUE', lambda s: s.get(snp_key, {}).get('hwe', {}).get('p_value', 1.0))
            ]
            
            m_matrix = []
            for lbl, fn in rows_cfg:
                row_item = {f'Statistical Parameter ({snp_key})': lbl}
                for r in regions:
                    row_item[r] = fn(reg.get(r, {}).get('snps', {}))
                m_matrix.append(row_item)
            st.dataframe(pd.DataFrame(m_matrix), use_container_width=True, hide_index=True)

        with t_cyp2:
            render_marker_geo_table('CYP2C19*2')
        with t_cyp3:
            render_marker_geo_table('CYP2C19*3')
        with t_cyp17:
            render_marker_geo_table('CYP2C19*17')

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("#### 💊 State-Wise Metabolizer Phenotype Distribution Matrix")
        st.caption("Distribution of CPIC CYP2C19 Metabolizer Phenotypes across individual states in your cohort.")

        st_dict = full_results.get('state_wise', {})
        state_pheno_matrix = []
        for s_name, s_data in sorted(st_dict.items(), key=lambda x: x[1].get('sample_count', 0), reverse=True):
            if s_data.get('sample_count', 0) > 0:
                p_dict = s_data.get('cyp2c19_summary', {}).get('phenotypes', {})
                row = {'State': s_name, 'Total Sample N': s_data.get('sample_count', 0)}
                for p_cat in PHENOTYPE_ORDER:
                    cnt = p_dict.get(p_cat, {}).get('count', 0)
                    pct = p_dict.get(p_cat, {}).get('percentage', 0.0)
                    p_label = p_cat.split(' (')[0]
                    row[p_label] = f"{cnt} ({pct:.1f}%)"
                state_pheno_matrix.append(row)

        # Include summary row for total cohort
        tot_p = ov.get('cyp2c19_summary', {}).get('phenotypes', {})
        tot_row = {'State': 'Total Cohort (All States Combined)', 'Total Sample N': ov.get('sample_count', 0)}
        for p_cat in PHENOTYPE_ORDER:
            cnt = tot_p.get(p_cat, {}).get('count', 0)
            pct = tot_p.get(p_cat, {}).get('percentage', 0.0)
            p_label = p_cat.split(' (')[0]
            tot_row[p_label] = f"{cnt} ({pct:.1f}%)"
        state_pheno_matrix.append(tot_row)

        st.dataframe(pd.DataFrame(state_pheno_matrix), use_container_width=True, hide_index=True)

    # -------------------------------------------------------------------------
    # TAB 8: GENDER-WISE ANALYSIS
    # -------------------------------------------------------------------------
    with tab_gender:
        st.markdown("### 👫 Gender-Wise Population Pharmacogenomics Analysis")
        gender_data = full_results.get('gender', {})
        male = gender_data.get('Male', {})
        female = gender_data.get('Female', {})
        
        m_count = male.get('sample_count', 0)
        f_count = female.get('sample_count', 0)
        tot_count = ov.get('sample_count', 0)
        
        gk1, gk2, gk3 = st.columns(3)
        gk1.metric("Male Samples (N)", f"{m_count:,}")
        gk2.metric("Female Samples (N)", f"{f_count:,}")
        gk3.metric("Total Cohort N", f"{tot_count:,}")
        
        m_phenos = male.get('cyp2c19_summary', {}).get('phenotypes', {})
        f_phenos = female.get('cyp2c19_summary', {}).get('phenotypes', {})
        
        pheno_chart_df = []
        for p_name in PHENOTYPE_ORDER:
            if p_name in m_phenos:
                pheno_chart_df.append({'Gender': 'Male', 'Phenotype': p_name.split(' (')[0], 'Count': m_phenos[p_name]['count']})
            if p_name in f_phenos:
                pheno_chart_df.append({'Gender': 'Female', 'Phenotype': p_name.split(' (')[0], 'Count': f_phenos[p_name]['count']})
                
        st.markdown("<br>", unsafe_allow_html=True)
        cg1, cg2 = st.columns(2)
        with cg1:
            if pheno_chart_df:
                fig_g2 = px.bar(
                    pd.DataFrame(pheno_chart_df), 
                    x='Phenotype', y='Count', 
                    color='Gender', barmode='group', 
                    title="Metabolizer Phenotypes Comparison (Male vs Female)",
                    color_discrete_sequence=['#2563EB', '#EC4899']
                )
                fig_g2.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_g2, use_container_width=True)

        with cg2:
            m_snps = male.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
            f_snps = female.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
            gen_gt_df = []
            for gt in ['GG', 'GA', 'AA']:
                if gt in m_snps:
                    gen_gt_df.append({'Gender': 'Male', 'Genotype': gt, 'Count': m_snps[gt]})
                if gt in f_snps:
                    gen_gt_df.append({'Gender': 'Female', 'Genotype': gt, 'Count': f_snps[gt]})
            if gen_gt_df:
                fig_gen_gt = px.bar(
                    pd.DataFrame(gen_gt_df), x='Genotype', y='Count', color='Gender', barmode='group',
                    title="CYP2C19*2 Genotypes Comparison (Male vs Female)",
                    color_discrete_sequence=['#2563EB', '#EC4899']
                )
                fig_gen_gt.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_gen_gt, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 9: DEDICATED STATE-WISE PAGES
    # -------------------------------------------------------------------------
    with tab_state:
        st.markdown("### 📍 Dedicated Individual State Profiles")
        state_dict = full_results.get('state_wise', {})
        state_names = sorted(list(state_dict.keys()))
        
        if state_names:
            selected_st = st.selectbox("Select Dedicated State Profile:", state_names, index=0, key="st_tab_selectbox")
            st_data = state_dict[selected_st]
            
            sk1, sk2, sk3, sk4 = st.columns(4)
            sk1.metric("State Sample N", st_data.get('sample_count', 0))
            sk2.metric("Missing Data Rate", f"{st_data.get('missing_pct', 0)}%")
            cyp2_freq = st_data.get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
            sk3.metric("CYP2C19*2 Var Freq", cyp2_freq)
            pm_pct = st_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer', st_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {})).get('percentage', 0)
            sk4.metric("Poor Metabolizers %", f"{pm_pct:.1f}%")

            st_phenos = st_data.get('cyp2c19_summary', {}).get('phenotypes', {})
            st_snps = st_data.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
            
            cs1, cs2 = st.columns(2)
            with cs1:
                st_gt_df = pd.DataFrame([{'Genotype': k, 'Count': v} for k, v in st_snps.items()])
                if not st_gt_df.empty:
                    fig_st_gt = px.bar(
                        st_gt_df, x='Genotype', y='Count',
                        title=f"{selected_st} — CYP2C19*2 Genotypes",
                        color='Genotype',
                        color_discrete_sequence=['#16A34A', '#2563EB', '#DC2626']
                    )
                    fig_st_gt.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_st_gt, use_container_width=True)

            with cs2:
                st_phe_df = pd.DataFrame([{'Phenotype': k.split(' (')[0], 'Count': v['count']} for k, v in st_phenos.items() if v['count'] > 0])
                if not st_phe_df.empty:
                    fig_st_phe = px.pie(
                        st_phe_df, names='Phenotype', values='Count', hole=0.45,
                        title=f"{selected_st} — Metabolizer Phenotypes",
                        color_discrete_sequence=['#166534', '#1E40AF', '#4B5563', '#991B1B', '#D97706', '#0D9488']
                    )
                    fig_st_phe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_st_phe, use_container_width=True)

# -----------------------------------------------------------------------------
# 3-STEP WORKFLOW PIPELINE RENDERER (EXACT MATCH FOR USER SCREENSHOTS 1, 2, 3, 4, 5)
# -----------------------------------------------------------------------------
def render_workflow_pipeline(df_raw, full_results, qc_report):
    VALID_STEPS = ['1 · Select Target', '2 · Upload', '3 · Map columns', '4 · Results']
    curr_step = st.session_state.get('workflow_step', '')
    if curr_step not in VALID_STEPS:
        if 'Target' in curr_step:
            st.session_state['workflow_step'] = '1 · Select Target'
        elif 'Upload' in curr_step:
            st.session_state['workflow_step'] = '2 · Upload'
        elif 'Map' in curr_step:
            st.session_state['workflow_step'] = '3 · Map columns'
        elif 'Result' in curr_step:
            st.session_state['workflow_step'] = '4 · Results'
        else:
            st.session_state['workflow_step'] = '1 · Select Target'

    main_tab_home, main_tab_about = st.tabs(["🏠 Home Workspace", "ℹ️ About Platform"])

    with main_tab_about:
        # Premium Executive Header Card with Glowing Title Pill Badge
        st.markdown("""
        <div style="background: linear-gradient(135deg, #0F172A 0%, #1E1B4B 45%, #2E1065 100%); border-radius: 16px; padding: 2.5rem 2.8rem; color: #FFFFFF; margin-top: 0.5rem; margin-bottom: 2rem; box-shadow: 0 12px 30px rgba(15, 23, 42, 0.4); border: 1px solid #4338CA;">
            <div style="background: linear-gradient(90deg, #4F46E5 0%, #7C3AED 100%); border: 1px solid #A78BFA; border-radius: 12px; padding: 0.8rem 1.5rem; display: inline-block; margin-bottom: 1.2rem; box-shadow: 0 4px 14px rgba(124, 58, 237, 0.35);">
                <h1 style="color: #FFFFFF; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; font-weight: 800; font-size: 1.85rem; margin: 0; letter-spacing: -0.01em;">
                    🧬 Automated Population Pharmacogenomics Analysis Platform
                </h1>
            </div>
            <p style="color: #F1F5F9; font-size: 1.1rem; line-height: 1.7; font-weight: 400; max-width: 900px; margin-bottom: 1.8rem;">
                <strong>PGx Analytics Pro</strong> is a standardized, high-performance pharmacogenomics engine engineered for instant quality control, Hardy–Weinberg Equilibrium (HWE) evaluation, regional demographic stratification, and CPIC metabolizer phenotype translation.
            </p>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 1rem;">
                <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #34D399; backdrop-filter: blur(8px);">
                    <div style="font-size: 0.78rem; text-transform: uppercase; color: #A7F3D0; font-weight: 800; letter-spacing: 0.06em;">Target Gene & SNPs</div>
                    <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">CYP2C19 (*2, *3, *17)</div>
                    <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">rs4244285, rs4986893, rs12248560</div>
                </div>
                <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #60A5FA; backdrop-filter: blur(8px);">
                    <div style="font-size: 0.78rem; text-transform: uppercase; color: #93C5FD; font-weight: 800; letter-spacing: 0.06em;">Statistical & HWE Engine</div>
                    <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">Chi-Square (χ²) & Haldane</div>
                    <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">Threshold χ² ≥ 3.4 (p ≤ 0.05)</div>
                </div>
                <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #FBBF24; backdrop-filter: blur(8px);">
                    <div style="font-size: 0.78rem; text-transform: uppercase; color: #FDE68A; font-weight: 800; letter-spacing: 0.06em;">Demographic Matrix</div>
                    <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">5 Regional Macro-Zones</div>
                    <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">South, North, East, West, Central</div>
                </div>
                <div style="background: rgba(255, 255, 255, 0.08); border-radius: 10px; padding: 1rem 1.2rem; border-left: 4px solid #F472B6; backdrop-filter: blur(8px);">
                    <div style="font-size: 0.78rem; text-transform: uppercase; color: #F472B6; font-weight: 800; letter-spacing: 0.06em;">CPIC Clinical Calling</div>
                    <div style="font-size: 1.15rem; color: #FFFFFF; font-weight: 700; margin-top: 0.3rem;">5 Metabolizer Classes</div>
                    <div style="font-size: 0.82rem; color: #CBD5E1; margin-top: 0.2rem;">UM, RM, NM, IM, PM Phenotypes</div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("<h3 style='color:#1E293B; font-weight:700; margin-bottom:1.2rem;'>⚡ Platform Functionality Overview</h3>", unsafe_allow_html=True)

        # 6 Sleek Visual Cards in 2-Column Grid
        st.markdown("""
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 1.2rem; margin-bottom: 2rem;">
            <div style="background: #FFFFFF; border-radius: 12px; padding: 1.4rem; border: 1px solid #E2E8F0; box-shadow: 0 4px 12px rgba(0,0,0,0.04); border-top: 4px solid #4F46E5;">
                <div style="display: flex; align-items: center; gap: 0.8rem; margin-bottom: 0.6rem;">
                    <span style="font-size: 1.6rem; background: #EEF2FF; padding: 0.4rem; border-radius: 8px;">📁</span>
                    <h4 style="margin:0; color:#1E293B; font-weight:700; font-size:1.1rem;">Dataset Import & Smart Column Mapping</h4>
                </div>
                <p style="color:#64748B; font-size:0.92rem; line-height:1.5; margin:0;">
                    Accepts raw Excel (<code>.xlsx</code>, <code>.xls</code>) and CSV files with intelligent auto-detection of Sample IDs, location, and genotype columns.
                </p>
            </div>
            <div style="background: #FFFFFF; border-radius: 12px; padding: 1.4rem; border: 1px solid #E2E8F0; box-shadow: 0 4px 12px rgba(0,0,0,0.04); border-top: 4px solid #059669;">
                <div style="display: flex; align-items: center; gap: 0.8rem; margin-bottom: 0.6rem;">
                    <span style="font-size: 1.6rem; background: #ECFDF5; padding: 0.4rem; border-radius: 8px;">🛡️</span>
                    <h4 style="margin:0; color:#1E293B; font-weight:700; font-size:1.1rem;">Automated Quality Control & Deduplication</h4>
                </div>
                <p style="color:#64748B; font-size:0.92rem; line-height:1.5; margin:0;">
                    Instantly screens uploaded files to evaluate valid genotype call rates, flag missing genotypes, and detect duplicate Sample IDs within the dataset.
                </p>
            </div>
            <div style="background: #FFFFFF; border-radius: 12px; padding: 1.4rem; border: 1px solid #E2E8F0; box-shadow: 0 4px 12px rgba(0,0,0,0.04); border-top: 4px solid #7C3AED;">
                <div style="display: flex; align-items: center; gap: 0.8rem; margin-bottom: 0.6rem;">
                    <span style="font-size: 1.6rem; background: #F3E8FF; padding: 0.4rem; border-radius: 8px;">🧬</span>
                    <h4 style="margin:0; color:#1E293B; font-weight:700; font-size:1.1rem;">Target CYP2C19 SNP Panel Genotyping</h4>
                </div>
                <p style="color:#64748B; font-size:0.92rem; line-height:1.5; margin:0;">
                    Evaluates key pharmacogenetic star-allele variants: <code>CYP2C19*2</code> (rs4244285), <code>CYP2C19*3</code> (rs4986893), and <code>CYP2C19*17</code> (rs12248560).
                </p>
            </div>
            <div style="background: #FFFFFF; border-radius: 12px; padding: 1.4rem; border: 1px solid #E2E8F0; box-shadow: 0 4px 12px rgba(0,0,0,0.04); border-top: 4px solid #DC2626;">
                <div style="display: flex; align-items: center; gap: 0.8rem; margin-bottom: 0.6rem;">
                    <span style="font-size: 1.6rem; background: #FEE2E2; padding: 0.4rem; border-radius: 8px;">⚠️</span>
                    <h4 style="margin:0; color:#1E293B; font-weight:700; font-size:1.1rem;">Hardy–Weinberg Equilibrium Red Warning Alerts</h4>
                </div>
                <p style="color:#64748B; font-size:0.92rem; line-height:1.5; margin:0;">
                    Calculates observed vs expected genotype counts and automatically highlights statistically significant HWE deviations (<code>χ² ≥ 3.4</code>) in bright red.
                </p>
            </div>
            <div style="background: #FFFFFF; border-radius: 12px; padding: 1.4rem; border: 1px solid #E2E8F0; box-shadow: 0 4px 12px rgba(0,0,0,0.04); border-top: 4px solid #D97706;">
                <div style="display: flex; align-items: center; gap: 0.8rem; margin-bottom: 0.6rem;">
                    <span style="font-size: 1.6rem; background: #FEF3C7; padding: 0.4rem; border-radius: 8px;">🗺️</span>
                    <h4 style="margin:0; color:#1E293B; font-weight:700; font-size:1.1rem;">Demographic & Regional Stratification</h4>
                </div>
                <p style="color:#64748B; font-size:0.92rem; line-height:1.5; margin:0;">
                    Groups state sample locations into 5 Indian macro-zones (South, North, East, West, Central) to generate regional allele distribution matrices.
                </p>
            </div>
            <div style="background: #FFFFFF; border-radius: 12px; padding: 1.4rem; border: 1px solid #E2E8F0; box-shadow: 0 4px 12px rgba(0,0,0,0.04); border-top: 4px solid #0284C7;">
                <div style="display: flex; align-items: center; gap: 0.8rem; margin-bottom: 0.6rem;">
                    <span style="font-size: 1.6rem; background: #E0F2FE; padding: 0.4rem; border-radius: 8px;">💊</span>
                    <h4 style="margin:0; color:#1E293B; font-weight:700; font-size:1.1rem;">CPIC Metabolizer Phenotype Translation</h4>
                </div>
                <p style="color:#64748B; font-size:0.92rem; line-height:1.5; margin:0;">
                    Translates star-allele diplotypes into 5 standardized clinical phenotype classes: Ultrarapid (UM), Rapid (RM), Normal (NM), Intermediate (IM), and Poor (PM).
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)

        # Bottom Privacy & Export Highlights
        st.info("🔒 **100% Local Browser Privacy**: All dataset parsing, statistical math, and SQLite indexing (`pgx_analysis.db`) run locally on your device with zero cloud data transmission.")
        st.success("📊 **Multi-Format Publication Export**: One-click download of ready-to-present multi-tab Excel (`.xlsx`), raw summary CSV (`.csv`), and formatted PDF (`.pdf`) reports.")

        st.markdown("<br>", unsafe_allow_html=True)

    with main_tab_home:
        step_cols = st.columns(4)
        t1 = "primary" if st.session_state['workflow_step'] == '1 · Select Target' else "secondary"
        t2 = "primary" if st.session_state['workflow_step'] == '2 · Upload' else "secondary"
        t3 = "primary" if st.session_state['workflow_step'] == '3 · Map columns' else "secondary"
        t4 = "primary" if st.session_state['workflow_step'] == '4 · Results' else "secondary"

        if step_cols[0].button("1 · Select Target", use_container_width=True, type=t1, key="nav_btn_step1"):
            st.session_state['workflow_step'] = '1 · Select Target'
            st.rerun()

        if step_cols[1].button("2 · Upload", use_container_width=True, type=t2, key="nav_btn_step2"):
            st.session_state['workflow_step'] = '2 · Upload'
            st.rerun()

        if step_cols[2].button("3 · Map columns", use_container_width=True, type=t3, key="nav_btn_step3"):
            st.session_state['workflow_step'] = '3 · Map columns'
            st.rerun()

        if step_cols[3].button("4 · Results", use_container_width=True, type=t4, key="nav_btn_step4"):
            st.session_state['workflow_step'] = '4 · Results'
            st.rerun()

        st.markdown("<br>", unsafe_allow_html=True)

        catalog_df = load_marker_catalog()

        # -------------------------------------------------------------------------
        # STEP 1: SELECT TARGET (GENE & BIOMARKERS)
        # -------------------------------------------------------------------------
        if st.session_state['workflow_step'] == '1 · Select Target':
            st.markdown("## 1️⃣ Select Pharmacogenetic Target (Gene & Biomarkers)")
            st.caption("Choose which gene and specific biomarker(s) to analyze. Only your selected markers will be uploaded, mapped, and visualized.")

            catalog_genes = sorted(catalog_df['Gene'].unique()) if not catalog_df.empty else ['CYP2C19']
            custom_genes = st.session_state.get('custom_genes', [])
            gene_list = sorted(list(set(catalog_genes + custom_genes)))
            curr_gene = st.session_state.get('selected_gene', 'CYP2C19')
            if curr_gene not in gene_list:
                curr_gene = 'CYP2C19' if 'CYP2C19' in gene_list else gene_list[0]
            g_idx = gene_list.index(curr_gene)

            c_g1, c_g2 = st.columns([1, 1.8])
            with c_g1:
                sel_gene = st.selectbox("1️⃣ Select Target Gene:", gene_list, index=g_idx, key="step1_sel_gene")
                if sel_gene != st.session_state.get('selected_gene'):
                    st.session_state['selected_gene'] = sel_gene
                    st.session_state['selected_markers'] = []
                    st.rerun()

            with c_g2:
                if sel_gene == 'CYP2C19':
                    base_markers = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                elif sel_gene == 'CYP2C9':
                    base_markers = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
                else:
                    raw_ms = catalog_df[catalog_df['Gene'] == sel_gene]['Marker'].tolist() if not catalog_df.empty else []
                    base_markers = raw_ms if raw_ms else [f"{sel_gene}_marker1"]

                custom_m_for_gene = st.session_state.get('custom_markers', {}).get(sel_gene, [])
                available_markers = []
                for m in base_markers + custom_m_for_gene:
                    if m not in available_markers:
                        available_markers.append(m)

                curr_selected = st.session_state.get('selected_markers')
                if curr_selected is not None:
                    valid_defaults = [m for m in curr_selected if m in available_markers]
                else:
                    valid_defaults = available_markers

                sel_markers = st.multiselect(
                    f"2️⃣ Select Biomarker(s) for {sel_gene}:",
                    available_markers,
                    default=valid_defaults,
                    key="step1_sel_markers"
                )
                st.session_state['selected_markers'] = sel_markers

            # Expandable Section: Add Custom Gene / Biomarker Manually
            with st.expander("➕ Can't find your Gene or Biomarker? Add Custom Targets Manually", expanded=False):
                st.markdown("Add custom target genes or biomarkers not present in the default catalog:")
                t_mode = st.radio(
                    "Target Addition Type:",
                    [f"Add Biomarker(s) to Selected Gene ({sel_gene})", "Add a Completely New Gene & Biomarker(s)"],
                    horizontal=True,
                    key="step1_custom_type_radio"
                )
                if "Completely New" in t_mode:
                    c_cg1, c_cg2 = st.columns([1, 1.8])
                    with c_cg1:
                        in_new_gene = st.text_input("New Gene Symbol / Name:", placeholder="e.g. CYP2D6, DPYD, TPMT", key="in_new_custom_gene_step1").strip().upper()
                    with c_cg2:
                        in_new_markers = st.text_input("Biomarker Name(s) (comma-separated):", placeholder="e.g. *4 (rs3892097), *10 (rs1065852)", key="in_new_custom_markers_step1").strip()
                else:
                    in_new_gene = sel_gene
                    in_new_markers = st.text_input(f"New Biomarker Name(s) for {sel_gene} (comma-separated):", placeholder="e.g. *4 (rs28399504), *6", key="in_existing_custom_markers_step1").strip()

                c_btn_add, c_btn_reset = st.columns([0.35, 0.65])
                with c_btn_add:
                    if st.button("➕ Add to Target List", type="secondary", key="btn_add_custom_target_step1"):
                        if not in_new_gene:
                            st.error("Please provide a valid Gene name.")
                        else:
                            c_genes = st.session_state.get('custom_genes', [])
                            if in_new_gene not in c_genes and (catalog_df.empty or in_new_gene not in catalog_df['Gene'].unique()):
                                c_genes.append(in_new_gene)
                                st.session_state['custom_genes'] = c_genes

                            c_markers_dict = st.session_state.get('custom_markers', {})
                            current_m = c_markers_dict.get(in_new_gene, [])
                            if in_new_markers:
                                parsed_m = [m.strip() for m in in_new_markers.split(',') if m.strip()]
                                for m in parsed_m:
                                    if m not in current_m:
                                        current_m.append(m)
                            else:
                                if not current_m:
                                    current_m.append(f"{in_new_gene}_marker1")
                            c_markers_dict[in_new_gene] = current_m
                            st.session_state['custom_markers'] = c_markers_dict

                            st.session_state['selected_gene'] = in_new_gene
                            st.session_state['selected_markers'] = current_m
                            st.success(f"Successfully added **{in_new_gene}** ({', '.join(current_m)}) to your active targets!")
                            st.rerun()

                with c_btn_reset:
                    custom_genes_list = st.session_state.get('custom_genes', [])
                    custom_markers_dict = st.session_state.get('custom_markers', {})
                    if custom_genes_list or custom_markers_dict:
                        if st.button("🗑️ Clear Custom Targets", key="btn_clear_custom_targets_step1"):
                            st.session_state['custom_genes'] = []
                            st.session_state['custom_markers'] = {}
                            st.session_state['selected_gene'] = 'CYP2C19'
                            st.session_state['selected_markers'] = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                            st.rerun()

            # Active Target Summary Card
            st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
            active_gene = st.session_state.get('selected_gene', 'CYP2C19')
            active_markers = st.session_state.get('selected_markers', [])
            num_m = len(active_markers)

            st.markdown(f"""
            <div style="background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 12px 16px;">
                <div style="font-size: 11px; font-weight: 700; text-transform: uppercase; color: #64748B;">🎯 Active Configuration</div>
                <div style="font-size: 14px; font-weight: 600; color: #0F172A; margin-top: 4px;">
                    Gene: <span style="color: #2563EB;">{active_gene}</span> &nbsp;|&nbsp; 
                    Selected Markers ({num_m}): <span style="color: #059669;">{', '.join(active_markers) if active_markers else 'None'}</span>
                </div>
            </div>
            """, unsafe_allow_html=True)

            st.markdown("<br>", unsafe_allow_html=True)
            c_sp, c_next1 = st.columns([0.7, 0.3])
            with c_next1:
                if st.button("Continue to Upload Dataset →", type="primary", use_container_width=True, key="btn_step0_cont", disabled=num_m == 0):
                    st.session_state['workflow_step'] = '2 · Upload'
                    st.rerun()

        # -------------------------------------------------------------------------
        # STEP 2: UPLOAD DATASET (STRICT GATEKEEPER & VALIDATION CHECK)
        # -------------------------------------------------------------------------
        elif st.session_state['workflow_step'] == '2 · Upload':
            st.markdown("## 2️⃣ Upload Genotype Dataset")
            st.caption("Upload your spreadsheet (.xlsx, .xls, .csv). Compulsory fields (Sample ID, Gender, State) will be inspected immediately.")

            st.markdown("""
            <div class="upload-container-box">
                <div class="upload-icon-style">⇧</div>
                <div class="upload-title-style">Click to choose a file, or drag it here</div>
                <div class="upload-sub-style">.xlsx, .xls, or .csv — first row must be headers</div>
            </div>
            """, unsafe_allow_html=True)


            # Display active workspace status if data is loaded
            if has_active_data():
                c_status, c_act1, c_down_top, c_act2 = st.columns([1.1, 0.7, 0.7, 0.5])
                with c_status:
                    st.info(f"📂 **Active Workspace Data:** `{len(df_raw):,} samples` loaded.")
                with c_act1:
                    if st.button("📊 View Results →", use_container_width=True, key="btn_step1_quick_res"):
                        st.session_state['workflow_step'] = '4 · Results'
                        st.rerun()
                with c_down_top:
                    st.download_button(
                        label="⬇️ Download (.xlsx)",
                        data=get_df_download_excel(df_raw),
                        file_name=f"complete_dataset_{len(df_raw)}_samples.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        key="btn_step2_top_download"
                    )
                with c_act2:
                    if st.button("🗑️ Delete", use_container_width=True, key="btn_step1_quick_del"):
                        reset_all_data()
                        st.rerun()

            file_upload = st.file_uploader("Select dataset file", type=["xlsx", "xls", "csv"], key="w_step1_uploader", label_visibility="collapsed")
            if file_upload is not None:
                uploaded_file_id = f"{file_upload.name}_{file_upload.size}"
                if st.session_state.get('last_uploaded_file_id') != uploaded_file_id:
                    try:
                        file_bytes = io.BytesIO(file_upload.getvalue())
                        if file_upload.name.lower().endswith('.csv'):
                            df_load = pd.read_csv(file_bytes)
                        else:
                            excel_file = pd.ExcelFile(file_bytes)
                            sheet_names = excel_file.sheet_names
                            if len(sheet_names) > 1:
                                best_sheet = sheet_names[0]
                                max_density = -1
                                for s in sheet_names:
                                    s_df = pd.read_excel(excel_file, sheet_name=s)
                                    gt_cols = [c for c in s_df.columns if any(k in str(c).lower() for k in ['*2', '*3', '*17', 'cyp', 'rs', 'snp', 'marker', 'gene'])]
                                    score = int(s_df[gt_cols].notna().sum().sum()) if gt_cols else int(s_df.notna().sum().sum())
                                    if score > max_density:
                                        max_density = score
                                        best_sheet = s
                                df_load = pd.read_excel(excel_file, sheet_name=best_sheet)
                                st.session_state['uploaded_sheet_names'] = sheet_names
                                st.session_state['active_sheet_name'] = best_sheet
                                st.session_state['uploaded_excel_bytes'] = file_upload.getvalue()
                            else:
                                df_load = pd.read_excel(file_bytes)
                        
                        # POST-UPLOAD INSPECTION FOR COMPULSORY FIELDS
                        chk_res = inspect_uploaded_file(df_load)
                        
                        if not chk_res['has_sample_id']:
                            st.session_state['upload_error'] = (
                                "🚫 **Upload Blocked: Mandatory 'Sample ID' column is missing!**\n\n"
                                "A unique Sample ID (e.g. `Sample ID`, `sample_id`, `Patient ID`, `ID`) is strictly compulsory to track samples and calculate frequencies.\n\n"
                                "**You cannot proceed without a Sample ID column.** Please add a Sample ID column to your spreadsheet and re-upload."
                            )
                            st.session_state['pending_new_df'] = None
                        else:
                            st.session_state['upload_error'] = None
                            st.session_state['upload_inspection'] = chk_res
                            st.session_state['demographic_mapping'] = {
                                'Sample_ID': chk_res['sample_id_col'],
                                'Gender': chk_res['gender_col'],
                                'State': chk_res['state_col']
                            }

                            if has_active_data():
                                st.session_state['pending_new_df'] = df_load
                                st.session_state['pending_file_id'] = uploaded_file_id
                                st.session_state['pending_file_name'] = file_upload.name
                            else:
                                clear_deletion_flag()
                                st.session_state['last_uploaded_file_id'] = uploaded_file_id
                                st.session_state['uploaded_df'] = df_load
                                st.session_state['datasets_dict'] = {f"Active Dataset ({file_upload.name})": df_load}
                                st.session_state['mapped_cols'] = {}
                                save_active_dataset(df_load)
                                st.success(f"✅ Loaded **{file_upload.name}** ({len(df_load):,} samples, {len(df_load.columns)} columns).")
                                st.rerun()
                    except Exception as e:
                        st.error(f"Error loading file: {e}")

            # Multi-Sheet Switcher if uploaded workbook has multiple sheets
            if len(st.session_state.get('uploaded_sheet_names', [])) > 1 and st.session_state.get('uploaded_excel_bytes'):
                all_sheets = st.session_state['uploaded_sheet_names']
                cur_s = st.session_state.get('active_sheet_name', all_sheets[0])
                s_idx = all_sheets.index(cur_s) if cur_s in all_sheets else 0
                c_sh1, c_sh2 = st.columns([1.2, 1.8])
                with c_sh1:
                    sel_s = st.selectbox("📑 Multi-Sheet Detected — Switch Active Sheet:", all_sheets, index=s_idx, key="sel_active_excel_sheet")
                    if sel_s != cur_s:
                        st.session_state['active_sheet_name'] = sel_s
                        new_df = pd.read_excel(io.BytesIO(st.session_state['uploaded_excel_bytes']), sheet_name=sel_s)
                        chk_res = inspect_uploaded_file(new_df)
                        st.session_state['uploaded_df'] = new_df
                        st.session_state['upload_inspection'] = chk_res
                        st.session_state['demographic_mapping'] = {
                            'Sample_ID': chk_res['sample_id_col'],
                            'Gender': chk_res['gender_col'],
                            'State': chk_res['state_col']
                        }
                        save_active_dataset(new_df)
                        st.rerun()

            # Show Error Banner if Sample ID Missing
            if st.session_state.get('upload_error'):
                st.error(st.session_state['upload_error'])

            # Required Metadata Verification (Sample ID, Gender, State)
            st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)
            st.markdown("#### 📋 Required Metadata Fields (Compulsory)")
            st.caption("These 3 fields are required for patient tracking and demographic stratification.")

            cur_insp = st.session_state.get('upload_inspection')
            if not cur_insp and df_raw is not None and not df_raw.empty:
                cur_insp = inspect_uploaded_file(df_raw)
                st.session_state['upload_inspection'] = cur_insp

            raw_cols = list(df_raw.columns) if df_raw is not None and not df_raw.empty else []
            demo_map = st.session_state.get('demographic_mapping', {})
            if not demo_map and cur_insp:
                demo_map = {
                    'Sample_ID': cur_insp.get('sample_id_col'),
                    'Gender': cur_insp.get('gender_col'),
                    'State': cur_insp.get('state_col')
                }
                st.session_state['demographic_mapping'] = demo_map

            sid_mapped = demo_map.get('Sample_ID') or (cur_insp['sample_id_col'] if cur_insp else None)
            gen_mapped = demo_map.get('Gender') or (cur_insp['gender_col'] if cur_insp else None)
            state_mapped = demo_map.get('State') or (cur_insp['state_col'] if cur_insp else None)

            c_aud1, c_aud2, c_aud3 = st.columns(3)
            with c_aud1:
                if sid_mapped and sid_mapped in raw_cols:
                    st.success(f"✅ **Sample ID (Required):** `{sid_mapped}`")
                else:
                    st.error("❌ **Sample ID (Required):** Missing")
            with c_aud2:
                if gen_mapped and gen_mapped in raw_cols:
                    st.success(f"✅ **Gender (Required):** `{gen_mapped}`")
                else:
                    st.warning("⚠️ **Gender (Required):** Not detected in columns")
            with c_aud3:
                if state_mapped and state_mapped in raw_cols:
                    st.success(f"✅ **State / Region (Required):** `{state_mapped}`")
                else:
                    st.warning("⚠️ **State (Required):** Not detected in columns")

            if raw_cols:
                with st.expander("🔧 Confirm / Adjust Required Metadata Columns", expanded=False):
                    c_adj1, c_adj2, c_adj3 = st.columns(3)
                    sid_opt = raw_cols
                    gen_opt = ["(None)"] + raw_cols
                    state_opt = ["(None)"] + raw_cols

                    idx_sid = sid_opt.index(sid_mapped) if sid_mapped in sid_opt else 0
                    idx_gen = gen_opt.index(gen_mapped) if gen_mapped in gen_opt else 0
                    idx_state = state_opt.index(state_mapped) if state_mapped in state_opt else 0

                    chosen_sid = c_adj1.selectbox("Sample ID Column *", sid_opt, index=idx_sid, key="u_map_sid")
                    chosen_gen = c_adj2.selectbox("Gender Column *", gen_opt, index=idx_gen, key="u_map_gen")
                    chosen_state = c_adj3.selectbox("State / Region Column *", state_opt, index=idx_state, key="u_map_state")

                    st.session_state['demographic_mapping'] = {
                        'Sample_ID': chosen_sid,
                        'Gender': chosen_gen if chosen_gen != "(None)" else None,
                        'State': chosen_state if chosen_state != "(None)" else None
                    }

            # -------------------------------------------------------------------------
            # TARGET BIOMARKER & COLUMN HEADER AUDIT (STEP 1 TARGET VS UPLOADED HEADERS)
            # -------------------------------------------------------------------------
            if raw_cols:
                target_gene_active = st.session_state.get('selected_gene', 'CYP2C19')
                target_markers_active = st.session_state.get('selected_markers', [])
                if not target_markers_active:
                    if target_gene_active == 'CYP2C19':
                        target_markers_active = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                    elif target_gene_active == 'CYP2C9':
                        target_markers_active = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']

                st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)
                st.markdown("#### 🎯 Target Biomarker & Header Audit")
                st.caption(f"Verifying uploaded column headers against your Step 1 selection: Target Gene **{target_gene_active}** ({len(target_markers_active)} biomarker(s)).")

                detected_file_genes = detect_genes_in_columns(raw_cols)
                gene_mismatch = bool(detected_file_genes and target_gene_active.upper() not in [g.upper() for g in detected_file_genes])

                if gene_mismatch:
                    found_genes_str = ", ".join(detected_file_genes)
                    st.markdown(f"""
                    <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 14px; margin-bottom: 12px;">
                        <div style="display: flex; align-items: center; gap: 8px;">
                            <span style="font-size: 20px;">🚨</span>
                            <div style="font-weight: 800; color: #991B1B; font-size: 15px;">TARGET GENE & COLUMN HEADER MISMATCH DETECTED</div>
                        </div>
                        <div style="color: #7F1D1D; font-size: 13px; margin-top: 6px; line-height: 1.5;">
                            In Step 1 you selected target gene: <b>{target_gene_active}</b>.<br>
                            However, the column headers in your uploaded spreadsheet belong to: <b>{found_genes_str}</b>.<br>
                            If you intended to analyze <b>{found_genes_str}</b>, click the button below to switch your target gene to match this dataset.
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                    c_sw_btn, _ = st.columns([0.4, 0.6])
                    with c_sw_btn:
                        primary_detected = detected_file_genes[0]
                        if st.button(f"👉 Switch Target Gene to {primary_detected}", type="primary", use_container_width=True, key="btn_switch_gene_step2"):
                            st.session_state['selected_gene'] = primary_detected
                            if primary_detected == 'CYP2C9':
                                st.session_state['selected_markers'] = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
                            elif primary_detected == 'CYP2C19':
                                st.session_state['selected_markers'] = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                            else:
                                st.session_state['selected_markers'] = [f"{primary_detected}_marker1"]
                            for k in list(st.session_state.keys()):
                                if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', 'ex_')):
                                    del st.session_state[k]
                            st.rerun()

                # Audit each selected marker against uploaded columns
                audit_cols_display = []
                for m in target_markers_active:
                    m_base = m.split(' ')[0] if ' ' in m else m
                    m_star = re.search(r'\*\d+[a-zA-Z]?', m)
                    star_val = m_star.group(0).lower() if m_star else ""
                    rs_m = re.search(r'rs\d+', m.lower())
                    rs_val = rs_m.group(0) if rs_m else ""

                    matched_col = None
                    for c in raw_cols:
                        c_l = str(c).lower()
                        if target_gene_active.lower() in c_l and star_val and star_val in c_l:
                            matched_col = c
                            break
                        if rs_val and rs_val in c_l:
                            matched_col = c
                            break
                        if m_base.lower() in c_l:
                            matched_col = c
                            break

                    audit_cols_display.append({
                        'Biomarker (Step 1)': m,
                        'Header in File': matched_col if matched_col else "❌ Not found in file",
                        'Match Status': "✅ Header Match" if matched_col else "⚠️ Different / Missing Header"
                    })

                c_m_aud1, c_m_aud2 = st.columns([1.6, 1.4])
                with c_m_aud1:
                    st.dataframe(pd.DataFrame(audit_cols_display), use_container_width=True, hide_index=True)
                with c_m_aud2:
                    n_matched = sum(1 for a in audit_cols_display if "✅" in a['Match Status'])
                    if n_matched == len(target_markers_active) and not gene_mismatch:
                        st.success(f"✅ All {n_matched} selected biomarker header(s) matched in uploaded file.")
                    elif gene_mismatch:
                        st.error(f"🚨 **Target Gene Mismatch:** Uploaded file headers contain `{found_genes_str}`, but target gene is `{target_gene_active}`.")
                    else:
                        st.warning(f"⚠️ **Header Notice:** {n_matched} of {len(target_markers_active)} biomarker header(s) match file columns. You can manually map different headers in Step 3.")

            if st.session_state.get('pending_new_df') is not None:
                pending_df = st.session_state['pending_new_df']
                pending_fname = st.session_state.get('pending_file_name', 'Uploaded File')
                st.warning(f"⚠️ **An active dataset is already loaded ({len(st.session_state['uploaded_df']):,} samples).**\n\nHow would you like to process new file **{pending_fname}** ({len(pending_df):,} samples)?")
                
                c_mrg1, c_mrg2, c_mrg3 = st.columns(3)
                with c_mrg1:
                    if st.button("🔄 Merge with Active Dataset", type="primary", use_container_width=True, key="btn_merge_dataset"):
                        base_df = st.session_state['uploaded_df']
                        merged_df, align_log = smart_merge_datasets(base_df, pending_df)
                        clear_deletion_flag()
                        st.session_state['last_uploaded_file_id'] = st.session_state['pending_file_id']
                        st.session_state['uploaded_df'] = merged_df
                        st.session_state.setdefault('datasets_dict', {})['Master Merged Dataset'] = merged_df
                        st.session_state['merge_alignment_log'] = align_log
                        st.session_state['pending_new_df'] = None
                        st.session_state['pending_file_id'] = None
                        st.session_state['pending_file_name'] = None
                        save_active_dataset(merged_df)
                        st.success(f"✅ Successfully merged {len(pending_df):,} samples! Total dataset now contains {len(merged_df):,} samples.")
                        st.rerun()
                with c_mrg2:
                    if st.button("📁 Keep as Separate Dataset", type="primary", use_container_width=True, key="btn_separate_dataset"):
                        clear_deletion_flag()
                        st.session_state.setdefault('datasets_dict', {})[f"Dataset 1 (Active)"] = st.session_state['uploaded_df']
                        st.session_state['datasets_dict'][f"Dataset 2 ({pending_fname})"] = pending_df
                        st.session_state['last_uploaded_file_id'] = st.session_state['pending_file_id']
                        st.session_state['uploaded_df'] = pending_df
                        st.session_state['mapped_cols'] = {}
                        st.session_state['pending_new_df'] = None
                        st.session_state['pending_file_id'] = None
                        st.session_state['pending_file_name'] = None
                        save_active_dataset(pending_df)
                        st.success(f"📁 Added '{pending_fname}' as a separate dataset for independent analysis!")
                        st.rerun()
                with c_mrg3:
                    if st.button("🆕 Replace Active Dataset", type="secondary", use_container_width=True, key="btn_replace_dataset"):
                        clear_deletion_flag()
                        st.session_state['last_uploaded_file_id'] = st.session_state['pending_file_id']
                        st.session_state['uploaded_df'] = pending_df
                        st.session_state['datasets_dict'] = {f"Active Dataset ({pending_fname})": pending_df}
                        st.session_state['mapped_cols'] = {}
                        st.session_state['pending_new_df'] = None
                        st.session_state['pending_file_id'] = None
                        st.session_state['pending_file_name'] = None
                        save_active_dataset(pending_df)
                        st.success(f"🆕 Replaced active workspace with '{pending_fname}' ({len(pending_df):,} samples).")
                        st.rerun()

            # Show Gene Alignment Log if merge occurred
            if st.session_state.get('merge_alignment_log'):
                with st.expander("🧬 Gene & Column Alignment Verification Log", expanded=True):
                    for log_entry in st.session_state['merge_alignment_log']:
                        st.markdown(f"• {log_entry}")

            # Show active dataset preview and download options if loaded
            if has_active_data():
                st.markdown("<div style='margin-top: 14px;'></div>", unsafe_allow_html=True)
                c_prev_head, c_dl_xl, c_dl_csv = st.columns([1.2, 0.9, 0.9])
                with c_prev_head:
                    st.markdown(f"##### 📋 Active Dataset Preview ({len(df_raw):,} samples, {len(df_raw.columns)} columns):")
                with c_dl_xl:
                    st.download_button(
                        label="⬇️ Download Complete Data (.xlsx)",
                        data=get_df_download_excel(df_raw),
                        file_name=f"complete_pgx_dataset_{len(df_raw)}_samples.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        key="btn_step2_prev_download_xlsx"
                    )
                with c_dl_csv:
                    st.download_button(
                        label="⬇️ Download Complete Data (.csv)",
                        data=get_df_download_csv(df_raw),
                        file_name=f"complete_pgx_dataset_{len(df_raw)}_samples.csv",
                        mime="text/csv",
                        use_container_width=True,
                        key="btn_step2_prev_download_csv"
                    )
                st.dataframe(df_raw.head(5), use_container_width=True)

            st.markdown("""
            <div class="privacy-info-banner">
                No data leaves your browser. Parsing and all statistics run locally.
            </div>
            """, unsafe_allow_html=True)

            c_back0, c_space, c_next = st.columns([0.25, 0.45, 0.3])
            with c_back0:
                if st.button("← Back to Select Target", use_container_width=True, key="btn_step1_back"):
                    st.session_state['workflow_step'] = '1 · Select Target'
                    st.rerun()
            with c_next:
                can_proceed = has_active_data() and not st.session_state.get('upload_error')
                if st.button("Continue to Map Columns →", type="primary", use_container_width=True, key="btn_step1_cont", disabled=not can_proceed):
                    st.session_state['workflow_step'] = '3 · Map columns'
                    st.rerun()

            st.divider()

            with st.expander("🛠️ Advanced Workspace Actions & Manual Record Entry", expanded=False):
                if st.session_state.get('show_delete_confirm', False):
                    st.warning("⚠️ **Warning: Delete Active Dataset?**\nThis will permanently delete the current dataset and reset workspace to 0 samples. This action cannot be undone.")
                    c_del_exp1, c_del_exp2 = st.columns(2)
                    with c_del_exp1:
                        if st.button("🔴 Permanently Delete Dataset", type="primary", use_container_width=True, key="btn_exp_confirm_del"):
                            reset_all_data()
                            st.success("Active dataset permanently deleted.")
                            st.rerun()
                    with c_del_exp2:
                        if st.button("❌ Cancel", type="secondary", use_container_width=True, key="btn_exp_cancel_del"):
                            st.session_state['show_delete_confirm'] = False
                            st.rerun()
                else:
                    if st.button("🗑️ Delete Active Dataset", type="secondary", use_container_width=True, key="btn_exp_reset"):
                        st.session_state['show_delete_confirm'] = True
                        st.rerun()

                st.markdown("<br>", unsafe_allow_html=True)
                st.markdown("#### ➕ Manual Sample Data Entry Form")
                st.caption("Insert individual sample records.")
                
                with st.form("manual_sample_form_single", clear_on_submit=True):
                    f_c1, f_c2 = st.columns(2)
                    in_sid = f_c1.text_input("Sample ID *", value=f"PATIENT_{len(df_raw)+1 if df_raw is not None else 1:04d}")
                    in_gender = f_c2.selectbox("Gender *", ["Male", "Female", "Unknown"])
                    
                    f_c3, f_c4 = st.columns(2)
                    in_dob = f_c3.text_input("Date of Birth (YYYY-MM-DD)", value="1990-01-01")
                    in_state = f_c4.selectbox("State / Region *", list(STATE_TO_REGION.keys()))
                    
                    in_native = st.text_input("Native Place (City / District)", value=in_state)
                    in_3gen = st.selectbox("Is their family lived at Native place for past 3 generations?", ["Yes", "No", "Unknown"])
                    
                    st.markdown("##### CYP2C19 Target Genotypes")
                    fg1, fg2, fg3 = st.columns(3)
                    in_cyp2 = fg1.selectbox("CYP2C19*2 (rs4244285) *", ["GA", "AA", "GG", "Missing", "CC (Invalid Call)"])
                    in_cyp3 = fg2.selectbox("CYP2C19*3 (rs4986893) *", ["GG", "Missing"])
                    in_cyp17 = fg3.selectbox("CYP2C19*17 (rs12248560) *", ["CC", "CT", "TT", "Missing", "GG (Invalid Call)"])
                    
                    submit_sample = st.form_submit_button("➕ Add Sample & Recalculate Statistics", type="primary", use_container_width=True)
                    
                    if submit_sample:
                        cyp2_val = in_cyp2.split(' ')[0]
                        cyp3_val = in_cyp3.split(' ')[0]
                        cyp17_val = in_cyp17.split(' ')[0]

                        new_row = {
                            'sample ID': in_sid,
                            'Gender': in_gender,
                            'Date of Birth': in_dob,
                            'Native place ': in_native,
                            'State': in_state,
                            'Test requested': 'CYP2C19 Genotyping',
                            'Is their family lived at Native place for past 3 generations?': in_3gen,
                            'CYP2C19*2 (rs4244285)': np.nan if "Missing" in in_cyp2 else cyp2_val,
                            'CYP2C19*3 (rs4986893)': np.nan if "Missing" in in_cyp3 else cyp3_val,
                            'CYP2C19*17 ( rs12248560)': np.nan if "Missing" in in_cyp17 else cyp17_val
                        }
                        
                        clear_deletion_flag()
                        if st.session_state['uploaded_df'] is not None and len(st.session_state['uploaded_df']) > 0:
                            st.session_state['uploaded_df'] = pd.concat([st.session_state['uploaded_df'], pd.DataFrame([new_row])], ignore_index=True)
                        else:
                            st.session_state['uploaded_df'] = pd.DataFrame([new_row])
                        
                        save_active_dataset(st.session_state['uploaded_df'])
                        st.session_state['workflow_step'] = '4 · Results'
                        st.success(f"Sample {in_sid} added successfully!")
                        st.rerun()

        # -------------------------------------------------------------------------
        # STEP 3: MAP COLUMNS (MANUAL ALLELE & GENOTYPE MAPPING, MULTIPLE ROW ENTRY)
        # -------------------------------------------------------------------------
        elif st.session_state['workflow_step'] == '3 · Map columns':
            st.markdown("## 3️⃣ Map Columns: Gene & Marker Selection")
            st.caption("First select or add your target gene and biomarker(s), then map each marker to your spreadsheet columns. You can manually customize Reference and Variant alleles, define permitted genotypes, and enter multiple rows.")

            # -------------------------------------------------------------------------
            # 1. CASCADE DROPDOWNS: SELECT GENE FIRST, THEN SELECT MARKER(S)
            # -------------------------------------------------------------------------
            catalog_genes = sorted(catalog_df['Gene'].unique()) if not catalog_df.empty else ['CYP2C19']
            custom_genes = st.session_state.get('custom_genes', [])
            all_catalog_genes = sorted(list(set(catalog_genes + custom_genes)))

            curr_gene = st.session_state.get('selected_gene', 'CYP2C19')
            if curr_gene not in all_catalog_genes:
                curr_gene = 'CYP2C19' if 'CYP2C19' in all_catalog_genes else all_catalog_genes[0]
            g_idx = all_catalog_genes.index(curr_gene)

            c_g1, c_g2 = st.columns([1, 1.8])
            with c_g1:
                sel_map_gene = st.selectbox("1️⃣ Select Target Gene:", all_catalog_genes, index=g_idx, key="map_step_gene_sel")
                if sel_map_gene != st.session_state.get('selected_gene'):
                    st.session_state['selected_gene'] = sel_map_gene
                    st.session_state['selected_markers'] = []
                    st.rerun()

            with c_g2:
                # Dynamic list of markers for the chosen gene
                if sel_map_gene == 'CYP2C19':
                    base_markers = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                elif sel_map_gene == 'CYP2C9':
                    base_markers = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
                else:
                    raw_ms = catalog_df[catalog_df['Gene'] == sel_map_gene]['Marker'].tolist() if not catalog_df.empty else []
                    base_markers = raw_ms if raw_ms else [f"{sel_map_gene}_marker1"]

                custom_m_for_gene = st.session_state.get('custom_markers', {}).get(sel_map_gene, [])
                available_markers = []
                for m in base_markers + custom_m_for_gene:
                    if m not in available_markers:
                        available_markers.append(m)

                # Ensure default selections
                curr_selected = st.session_state.get('selected_markers')
                if curr_selected is not None:
                    valid_defaults = [m for m in curr_selected if m in available_markers]
                else:
                    valid_defaults = available_markers

                sel_map_markers = st.multiselect(
                    f"2️⃣ Select Biomarker(s) for {sel_map_gene}:",
                    available_markers,
                    default=valid_defaults,
                    key="map_step_marker_sel"
                )
                st.session_state['selected_markers'] = sel_map_markers

            # Expandable Section: Add Custom Gene / Biomarker Manually directly in Step 3
            with st.expander("➕ Can't find your Gene or Biomarker? Add Custom Targets Manually", expanded=False):
                st.markdown("Add custom target genes or biomarkers not present in the default catalog:")
                t_mode_step3 = st.radio(
                    "Target Addition Type:",
                    [f"Add Biomarker(s) to Selected Gene ({sel_map_gene})", "Add a Completely New Gene & Biomarker(s)"],
                    horizontal=True,
                    key="step3_custom_type_radio"
                )
                if "Completely New" in t_mode_step3:
                    c_cg1_3, c_cg2_3 = st.columns([1, 1.8])
                    with c_cg1_3:
                        in_new_gene_3 = st.text_input("New Gene Symbol / Name:", placeholder="e.g. CYP2D6, DPYD, TPMT", key="in_new_custom_gene_step3").strip().upper()
                    with c_cg2_3:
                        in_new_markers_3 = st.text_input("Biomarker Name(s) (comma-separated):", placeholder="e.g. *4 (rs3892097), *10 (rs1065852)", key="in_new_custom_markers_step3").strip()
                else:
                    in_new_gene_3 = sel_map_gene
                    in_new_markers_3 = st.text_input(f"New Biomarker Name(s) for {sel_map_gene} (comma-separated):", placeholder="e.g. *4 (rs28399504), *6", key="in_existing_custom_markers_step3").strip()

                c_btn_add_3, c_btn_reset_3 = st.columns([0.35, 0.65])
                with c_btn_add_3:
                    if st.button("➕ Add to Target List", type="secondary", key="btn_add_custom_target_step3"):
                        if not in_new_gene_3:
                            st.error("Please provide a valid Gene name.")
                        else:
                            c_genes = st.session_state.get('custom_genes', [])
                            if in_new_gene_3 not in c_genes and (catalog_df.empty or in_new_gene_3 not in catalog_df['Gene'].unique()):
                                c_genes.append(in_new_gene_3)
                                st.session_state['custom_genes'] = c_genes

                            c_markers_dict = st.session_state.get('custom_markers', {})
                            current_m = c_markers_dict.get(in_new_gene_3, [])
                            if in_new_markers_3:
                                parsed_m = [m.strip() for m in in_new_markers_3.split(',') if m.strip()]
                                for m in parsed_m:
                                    if m not in current_m:
                                        current_m.append(m)
                            else:
                                if not current_m:
                                    current_m.append(f"{in_new_gene_3}_marker1")
                            c_markers_dict[in_new_gene_3] = current_m
                            st.session_state['custom_markers'] = c_markers_dict

                            st.session_state['selected_gene'] = in_new_gene_3
                            st.session_state['selected_markers'] = current_m
                            st.success(f"Successfully added **{in_new_gene_3}** ({', '.join(current_m)}) to your active targets!")
                            st.rerun()

                with c_btn_reset_3:
                    custom_genes_list_3 = st.session_state.get('custom_genes', [])
                    custom_markers_dict_3 = st.session_state.get('custom_markers', {})
                    if custom_genes_list_3 or custom_markers_dict_3:
                        if st.button("🗑️ Clear Custom Targets", key="btn_clear_custom_targets_step3"):
                            st.session_state['custom_genes'] = []
                            st.session_state['custom_markers'] = {}
                            st.session_state['selected_gene'] = 'CYP2C19'
                            st.session_state['selected_markers'] = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                            st.rerun()

            if not sel_map_markers:
                st.warning(f"⚠️ Please select at least one biomarker for **{sel_map_gene}** from the dropdown above to map columns.")
            elif df_raw is not None and len(df_raw.columns) > 0:
                all_cols = list(df_raw.columns)
                demo_mapping = st.session_state.get('demographic_mapping', {})
                confirmed_demo_cols = {str(v).strip().lower() for k, v in demo_mapping.items() if v}
                cur_insp = st.session_state.get('upload_inspection', {})
                if cur_insp:
                    for k in ['sample_id_col', 'gender_col', 'state_col']:
                        if cur_insp.get(k):
                            confirmed_demo_cols.add(str(cur_insp[k]).strip().lower())

                NON_GENOTYPE_KEYWORDS = [
                    'sample', 'id', 'gender', 'sex', 'state', 'region', 'place', 
                    'date', 'birth', 'dob', 'test', 'family', 'generation', 'native', 
                    'unnamed', 'patient', 'name', 'age', 'center', 'hospital', 'clinic', 'india'
                ]

                def is_demographic_or_metadata(col_name: str) -> bool:
                    clean = str(col_name).strip().lower()
                    if clean.startswith('unnamed:'):
                        return True
                    if clean in confirmed_demo_cols:
                        return True
                    has_marker_symbol = any(sym in clean for sym in ['*', 'rs', 'cyp', 'dpyd', 'tpmt', 'vkorc', 'nudt', 'ugt', 'slco', 'hla', 'g6pd'])
                    if not has_marker_symbol:
                        for kw in NON_GENOTYPE_KEYWORDS:
                            if kw in clean:
                                return True
                    return False

                # Clean candidate genotype columns
                genotype_cols = [c for c in all_cols if not is_demographic_or_metadata(c)]
                if not genotype_cols:
                    genotype_cols = [c for c in all_cols if not str(c).strip().startswith('Unnamed:')]
                if not genotype_cols:
                    genotype_cols = all_cols

                # Check if uploaded columns suggest a specific target gene (e.g. CYP2C9, CYP2D6)
                detected_gene_matches = []
                for c in all_cols:
                    c_clean = str(c).upper().replace(' ', '').replace('_', '').replace('-', '')
                    for cat_gene in (catalog_df['Gene'].unique() if not catalog_df.empty else []):
                        g_clean = str(cat_gene).upper().replace(' ', '').replace('_', '').replace('-', '')
                        if g_clean != sel_map_gene.upper().replace(' ', '') and len(g_clean) >= 4 and g_clean in c_clean:
                            if cat_gene not in detected_gene_matches:
                                detected_gene_matches.append(cat_gene)

                # Define helper functions and callbacks for flexible marker mapping
                def extract_alleles_from_calls(calls):
                    al_counts = Counter()
                    for call in calls:
                        clean_c = str(call).replace('/', '').replace('|', '').replace(' ', '').upper()
                        for ch in clean_c:
                            if ch.isalpha():
                                al_counts[ch] += 1
                    most = [al for al, _ in al_counts.most_common()]
                    r = most[0] if len(most) > 0 else 'A'
                    v = most[1] if len(most) > 1 else 'None'
                    return r, v

                def find_default(patterns, cols):
                    for pattern in patterns:
                        for c in cols:
                            if str(pattern).strip().lower() == str(c).strip().lower():
                                return c
                    for pattern in patterns:
                        p_norm = str(pattern).lower().replace(" ", "").replace("*", "").replace("(", "").replace(")", "").replace("-", "")
                        if not p_norm: continue
                        for c in cols:
                            c_norm = str(c).lower().replace(" ", "").replace("*", "").replace("(", "").replace(")", "").replace("-", "")
                            if p_norm == c_norm or p_norm in c_norm or c_norm in p_norm:
                                return c
                    for pattern in patterns:
                        for c in cols:
                            if str(pattern).lower() in str(c).lower():
                                return c
                    return cols[0] if cols else ""

                def find_best_column_match(marker_label: str, std_name: str, sel_gene: str, cols: list) -> str:
                    if not cols:
                        return ""
                    # 1. Exact match with label, name, or gene_marker
                    for target in [marker_label, std_name, f"{sel_gene}_{std_name}", f"{sel_gene} {std_name}", f"{sel_gene}_{marker_label}"]:
                        for c in cols:
                            if str(target).strip().lower() == str(c).strip().lower():
                                return c

                    # 2. Extract rsID if present in marker_label or std_name
                    import re
                    rs_match = re.search(r'rs\d+', f"{marker_label} {std_name}".lower())
                    if rs_match:
                        rs_id = rs_match.group(0)
                        for c in cols:
                            if rs_id in str(c).lower():
                                return c

                    # 3. Check for gene + star allele e.g. CYP2C9 and *2
                    star_match = re.search(r'\*\d+[a-zA-Z]?', f"{marker_label} {std_name}")
                    if star_match:
                        star_str = star_match.group(0).lower()
                        for c in cols:
                            c_low = str(c).lower()
                            if sel_gene.lower() in c_low and star_str in c_low:
                                return c
                        # Only fallback to star_str if column header does NOT belong to another gene
                        for c in cols:
                            cg = get_column_gene(c)
                            if star_str in str(c).lower() and (not cg or cg.upper() == sel_gene.upper()):
                                return c

                    # 4. Fallback: column containing sel_gene
                    for c in cols:
                        if sel_gene.lower() in str(c).lower():
                            return c

                    # 5. Fallback to find_default
                    return find_default([marker_label, std_name], cols)

                def get_marker_defaults(marker_label: str, sel_gene: str, detected_calls: list = None, col_name: str = ""):
                    """
                    Generates clean, user-driven biomarker configuration without hardcoding.
                    Permitted genotypes are NOT pre-filled with catalog assumptions,
                    leaving full flexibility to the user to type or click observed genotypes.
                    """
                    clean_label = str(marker_label).strip()
                    if ' ' in clean_label and '(' in clean_label:
                        base_name = clean_label.split(' ')[0]
                    else:
                        base_name = clean_label

                    auto_r, auto_v = 'A', 'G'
                    if detected_calls and len(detected_calls) > 0:
                        cand_r, cand_v = extract_alleles_from_calls(detected_calls)
                        if cand_r:
                            auto_r = cand_r
                        if cand_v:
                            auto_v = cand_v

                    badge_short = f"Column: {col_name}" if col_name else "Spreadsheet Mapping"
                    if detected_calls:
                        badge_detail = f"Observed calls in data: <b>{', '.join(detected_calls)}</b>"
                    else:
                        badge_detail = f"Configure permitted genotypes and alleles for <b>{base_name}</b>"

                    # DO NOT pre-fill with hardcoded catalog strings - let the user enter or click chips!
                    def_permitted = ""
                    return base_name, auto_r, auto_v, def_permitted, badge_short, badge_detail

                # Define callbacks to update widget session state safely before widgets instantiate
                def cb_switch_target_gene(new_g, new_markers):
                    st.session_state['selected_gene'] = new_g
                    st.session_state['selected_markers'] = new_markers[:3] if new_markers else [f"{new_g}_marker1"]
                    for k in list(st.session_state.keys()):
                        if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', 'ex_')):
                            del st.session_state[k]

                def cb_remove_marker_card(idx_to_remove):
                    cur_markers = list(st.session_state.get('selected_markers', []))
                    if 0 <= idx_to_remove < len(cur_markers):
                        cur_markers.pop(idx_to_remove)
                        st.session_state['selected_markers'] = cur_markers
                        prefixes = ['w_map_name_', 'w_map_col_', 'w_map_ref_', 'w_map_var_', 'w_map_perm_', 'sel_ref_choice_', 'sel_var_choice_', 'custom_ref_in_', 'custom_var_in_']
                        for i in range(idx_to_remove, len(cur_markers) + 1):
                            for prefix in prefixes:
                                cur_k = f"{prefix}{i}"
                                next_k = f"{prefix}{i+1}"
                                if next_k in st.session_state:
                                    st.session_state[cur_k] = st.session_state[next_k]
                                    del st.session_state[next_k]
                                elif cur_k in st.session_state:
                                    del st.session_state[cur_k]

                def cb_add_new_marker_card(gene_prefix):
                    cur_markers = list(st.session_state.get('selected_markers', []))
                    new_num = len(cur_markers) + 1
                    new_marker_lbl = f"{gene_prefix}_Marker_{new_num}"
                    while new_marker_lbl in cur_markers:
                        new_num += 1
                        new_marker_lbl = f"{gene_prefix}_Marker_{new_num}"
                    cur_markers.append(new_marker_lbl)
                    st.session_state['selected_markers'] = cur_markers

                def cb_toggle_permitted_call(p_k, call, r_k, v_k, sel_r_k, sel_v_k):
                    cur_p = str(st.session_state.get(p_k, "")).strip()
                    calls = [g.strip().upper() for g in cur_p.split(',') if g.strip()]
                    call_u = str(call).strip().upper()
                    if call_u in calls:
                        calls.remove(call_u)
                    else:
                        calls.append(call_u)
                    new_p = ", ".join(calls)
                    st.session_state[p_k] = new_p
                    if calls:
                        r_c, v_c = extract_alleles_from_calls(calls)
                        if r_c:
                            st.session_state[r_k] = r_c
                            if sel_r_k:
                                st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                        if v_c:
                            st.session_state[v_k] = v_c
                            if sel_v_k:
                                st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'

                def cb_allow_all_detected(p_k, calls_list, r_k, v_k, sel_r_k, sel_v_k):
                    st.session_state[p_k] = ", ".join(calls_list)
                    r_c, v_c = extract_alleles_from_calls(calls_list)
                    if r_c:
                        st.session_state[r_k] = r_c
                        if sel_r_k:
                            st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                    if v_c:
                        st.session_state[v_k] = v_c
                        if sel_v_k:
                            st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'

                def cb_clear_permitted(p_k):
                    st.session_state[p_k] = ""

                def cb_automap_detected(p_k, calls_list, r_k, v_k, sel_r_k, sel_v_k):
                    st.session_state[p_k] = ", ".join(calls_list)
                    r_c, v_c = extract_alleles_from_calls(calls_list)
                    if r_c:
                        st.session_state[r_k] = r_c
                        if sel_r_k:
                            st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                            if r_c not in ['A', 'C', 'G', 'T']:
                                idx_str = sel_r_k.split('_')[-1]
                                st.session_state[f"custom_ref_in_{idx_str}"] = r_c
                    if v_c:
                        st.session_state[v_k] = v_c
                        if sel_v_k:
                            st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'
                            if v_c not in ['A', 'C', 'G', 'T', 'None']:
                                idx_str = sel_v_k.split('_')[-1]
                                st.session_state[f"custom_var_in_{idx_str}"] = v_c

                def cb_add_unmapped(p_k, unexp_list, r_k, v_k, sel_r_k, sel_v_k):
                    cur_p = str(st.session_state.get(p_k, "")).strip()
                    add_str = ", ".join(unexp_list)
                    if cur_p:
                        if not cur_p.endswith(','):
                            cur_p += ", "
                        new_p = cur_p + add_str
                    else:
                        new_p = add_str
                    st.session_state[p_k] = new_p
                    all_c = [g.strip().upper() for g in new_p.split(',') if g.strip()]
                    r_c, v_c = extract_alleles_from_calls(all_c)
                    if r_c:
                        st.session_state[r_k] = r_c
                        if sel_r_k:
                            st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                            if r_c not in ['A', 'C', 'G', 'T']:
                                idx_str = sel_r_k.split('_')[-1]
                                st.session_state[f"custom_ref_in_{idx_str}"] = r_c
                    if v_c:
                        st.session_state[v_k] = v_c
                        if sel_v_k:
                            st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'
                            if v_c not in ['A', 'C', 'G', 'T', 'None']:
                                idx_str = sel_v_k.split('_')[-1]
                                st.session_state[f"custom_var_in_{idx_str}"] = v_c

                def cb_align_alleles(r_k, v_k, r_c, v_c, sel_r_k, sel_v_k):
                    st.session_state[r_k] = r_c
                    st.session_state[v_k] = v_c
                    if sel_r_k:
                        st.session_state[sel_r_k] = r_c if r_c in ['A', 'C', 'G', 'T'] else 'Other (Manual)'
                        if r_c not in ['A', 'C', 'G', 'T']:
                            idx_str = sel_r_k.split('_')[-1]
                            st.session_state[f"custom_ref_in_{idx_str}"] = r_c
                    if sel_v_k:
                        st.session_state[sel_v_k] = v_c if v_c in ['A', 'C', 'G', 'T', 'None'] else 'Other (Manual)'
                        if v_c not in ['A', 'C', 'G', 'T', 'None']:
                            idx_str = sel_v_k.split('_')[-1]
                            st.session_state[f"custom_var_in_{idx_str}"] = v_c

                def cb_extra_automap(p_k, r_k, v_k, calls_list):
                    st.session_state[p_k] = ", ".join(calls_list)
                    r_c, v_c = extract_alleles_from_calls(calls_list)
                    if r_c:
                        st.session_state[r_k] = r_c
                    if v_c:
                        st.session_state[v_k] = v_c

                def cb_extra_add_unmapped(p_k, unexp_list, r_k, v_k):
                    cur_p = str(st.session_state.get(p_k, "")).strip()
                    add_str = ", ".join(unexp_list)
                    if cur_p:
                        if not cur_p.endswith(','):
                            cur_p += ", "
                        new_p = cur_p + add_str
                    else:
                        new_p = add_str
                    st.session_state[p_k] = new_p
                    all_c = [g.strip().upper() for g in new_p.split(',') if g.strip()]
                    r_c, v_c = extract_alleles_from_calls(all_c)
                    if r_c:
                        st.session_state[r_k] = r_c
                    if v_c:
                        st.session_state[v_k] = v_c


                c_opt_filter, c_mode_sel = st.columns([0.45, 0.55])
                with c_opt_filter:
                    show_all_cols = st.checkbox("Show all spreadsheet columns in options", value=False, key="cb_show_all_cols_step3", help="Enable if your genotype column was auto-classified as metadata")
                col_options = all_cols if show_all_cols else genotype_cols

                with c_mode_sel:
                    mapping_view_mode = st.radio(
                        "Mapping Interface:",
                        ["🗂️ Interactive Marker Cards", "📋 Multi-Row Table Editor"],
                        horizontal=True,
                        key="step3_mapping_view_mode"
                    )

                mapped_marker_dict = {}
                mapping_has_errors = False
                extra_rows = st.session_state.get('extra_mapping_rows', [])

                # -------------------------------------------------------------------------
                # -------------------------------------------------------------------------
                # MODE A: INTERACTIVE MARKER CARDS (WITH LIVE DATA QC & MULTI-VALUE COMBINATIONS)
                # -------------------------------------------------------------------------
                if mapping_view_mode == "🗂️ Interactive Marker Cards":
                    st.markdown("#### 🧬 Biomarker Mapping & Genotype Combination Rules")
                    st.caption("Map spreadsheet columns and configure any combination of permitted genotypes (e.g. GA, GG, AG or CC, CT, TT). Check live data breakdown and verify cohort calls instantly.")

                    detected_file_genes = detect_genes_in_columns(all_cols)
                    if detected_file_genes and sel_map_gene.upper() not in [g.upper() for g in detected_file_genes]:
                        found_genes_str = ", ".join(detected_file_genes)
                        st.markdown(f"""
                        <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 12px 16px; margin-bottom: 14px;">
                            <div style="display: flex; align-items: center; gap: 8px;">
                                <span style="font-size: 18px;">🚨</span>
                                <div style="font-weight: 800; color: #991B1B; font-size: 14px;">TARGET GENE & COLUMN HEADER DISCREPANCY DETECTED</div>
                            </div>
                            <div style="color: #7F1D1D; font-size: 13px; margin-top: 4px;">
                                You are mapping markers for <b>{sel_map_gene}</b>, but spreadsheet column headers contain <b>{found_genes_str}</b>.<br>
                                If you uploaded data for <b>{detected_file_genes[0]}</b>, click the button below to switch your target gene:
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
                        c_sw_map, _ = st.columns([0.4, 0.6])
                        with c_sw_map:
                            if st.button(f"👉 Switch Target Gene to {detected_file_genes[0]}", key="btn_switch_gene_step3_banner", type="primary"):
                                st.session_state['selected_gene'] = detected_file_genes[0]
                                if detected_file_genes[0] == 'CYP2C9':
                                    st.session_state['selected_markers'] = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
                                elif detected_file_genes[0] == 'CYP2C19':
                                    st.session_state['selected_markers'] = ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
                                else:
                                    st.session_state['selected_markers'] = [f"{detected_file_genes[0]}_marker1"]
                                for k in list(st.session_state.keys()):
                                    if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', 'ex_')):
                                        del st.session_state[k]
                                st.rerun()

                    for idx, marker_label in enumerate(sel_map_markers):
                        col_key = f"w_map_col_{idx}"
                        name_key = f"w_map_name_{idx}"
                        ref_key = f"w_map_ref_{idx}"
                        var_key = f"w_map_var_{idx}"
                        perm_key = f"w_map_perm_{idx}"

                        # Determine best matching column
                        best_match = find_best_column_match(marker_label, marker_label.split(' ')[0], sel_map_gene, col_options)
                        if col_key not in st.session_state or st.session_state[col_key] not in col_options:
                            st.session_state[col_key] = best_match if (best_match and best_match in col_options) else col_options[0]

                        sel_col = st.session_state[col_key]

                        # Live Data Inspection for this column
                        observed_counts = {}
                        detected_calls_list = []
                        norm_calls = pd.Series(dtype=str)
                        if sel_col in df_raw.columns:
                            raw_s = df_raw[sel_col].dropna().astype(str).str.strip().str.upper()
                            valid_s = raw_s[~raw_s.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])]
                            norm_calls = valid_s.map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                            observed_counts = dict(norm_calls.value_counts())
                            detected_calls_list = list(observed_counts.keys())
                        tot_samples_col = len(norm_calls)

                        # Determine defaults prioritizing detected calls in the data
                        std_name, def_ref, def_var, def_permitted, badge_short, badge_detail = get_marker_defaults(marker_label, sel_map_gene, detected_calls_list, col_name=sel_col)

                        # Initialize session state for this card
                        if name_key not in st.session_state:
                            st.session_state[name_key] = std_name
                        if ref_key not in st.session_state:
                            st.session_state[ref_key] = def_ref
                        if var_key not in st.session_state:
                            st.session_state[var_key] = def_var
                        if perm_key not in st.session_state:
                            st.session_state[perm_key] = def_permitted

                        cur_marker_name = st.session_state.get(name_key, std_name)
                        obs_str = ', '.join([f"<b>{c}</b> ({cnt:,})" for c, cnt in observed_counts.items()]) if observed_counts else "<i>No calls detected in this column</i>"
                        st.markdown(f"""
                        <div class="clean-card-box" style="margin-bottom: 8px; border-left: 4px solid #3B82F6; padding: 12px 16px; background: #F8FAFC; border-radius: 8px; border: 1px solid #E2E8F0;">
                            <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;">
                                <h5 style="margin:0; color:#0F172A; font-weight:700;">🧬 Biomarker #{idx+1}: {cur_marker_name}</h5>
                                <span style="background-color: #EFF6FF; color: #1D4ED8; font-size: 0.8rem; font-weight: 600; padding: 3px 10px; border-radius: 9999px; border: 1px solid #BFDBFE;">Column: <code>{sel_col}</code></span>
                            </div>
                            <div style="margin-top: 6px; font-size: 0.85rem; color: #475569;">
                                Observed Genotypes in column: {obs_str}
                            </div>
                        </div>
                        """, unsafe_allow_html=True)

                        # Check for header discrepancies (different gene, star allele, or rsID)
                        col_gene = get_column_gene(sel_col)
                        is_diff_gene = bool(col_gene and col_gene.upper() != sel_map_gene.upper())

                        m_star = re.search(r'\*\d+[a-zA-Z]?', cur_marker_name)
                        c_star = re.search(r'\*\d+[a-zA-Z]?', str(sel_col))
                        is_diff_star = bool(m_star and c_star and m_star.group(0).lower() != c_star.group(0).lower())

                        m_rs = re.search(r'rs\d+', cur_marker_name.lower())
                        c_rs = re.search(r'rs\d+', str(sel_col).lower())
                        is_diff_rs = bool(m_rs and c_rs and m_rs.group(0) != c_rs.group(0))

                        if is_diff_gene:
                            st.markdown(f"""
                            <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 10px 14px; margin-top: 6px; margin-bottom: 8px;">
                                <div style="font-weight: 800; color: #991B1B; font-size: 13px;">🚨 GENE HEADER MISMATCH: Column <code>{sel_col}</code> belongs to {col_gene}, but active gene is {sel_map_gene}!</div>
                                <div style="color: #7F1D1D; font-size: 12px; margin-top: 3px;">
                                    Target marker is <b>{cur_marker_name}</b> ({sel_map_gene}), but selected column <code>{sel_col}</code> represents gene <b>{col_gene}</b>.<br>
                                    Please verify if you selected the wrong column from the dropdown above.
                                </div>
                            </div>
                            """, unsafe_allow_html=True)
                        elif is_diff_star:
                            st.markdown(f"""
                            <div style="background: #FFFBEB; border: 2px solid #F59E0B; border-radius: 8px; padding: 10px 14px; margin-top: 6px; margin-bottom: 8px;">
                                <div style="font-weight: 800; color: #B45309; font-size: 13px;">⚠️ STAR ALLELE MISMATCH: Target is {m_star.group(0)}, but column header is {c_star.group(0)}!</div>
                                <div style="color: #92400E; font-size: 12px; margin-top: 3px;">
                                    Target marker <b>{cur_marker_name}</b> is for <b>{m_star.group(0)}</b>, but column header <code>{sel_col}</code> indicates <b>{c_star.group(0)}</b>.
                                </div>
                            </div>
                            """, unsafe_allow_html=True)
                        elif is_diff_rs:
                            st.markdown(f"""
                            <div style="background: #FFFBEB; border: 2px solid #F59E0B; border-radius: 8px; padding: 10px 14px; margin-top: 6px; margin-bottom: 8px;">
                                <div style="font-weight: 800; color: #B45309; font-size: 13px;">⚠️ RSID MISMATCH: Target is {m_rs.group(0)}, but column header has {c_rs.group(0)}!</div>
                                <div style="color: #92400E; font-size: 12px; margin-top: 3px;">
                                    Please verify if column <code>{sel_col}</code> corresponds to marker <b>{cur_marker_name}</b>.
                                </div>
                            </div>
                            """, unsafe_allow_html=True)

                        c_top1, c_top2, c_top3 = st.columns([1.1, 1.5, 0.4])
                        with c_top1:
                            user_marker_name = st.text_input(
                                "Biomarker Name:",
                                value=st.session_state.get(name_key, std_name),
                                key=name_key,
                                help="Edit marker name (e.g. CYP2C9*2, CYP2C19*2, rs1799853)"
                            ).strip()
                            if not user_marker_name:
                                user_marker_name = std_name
                        with c_top2:
                            def_col_idx = col_options.index(st.session_state[col_key]) if st.session_state[col_key] in col_options else 0
                            new_sel_col = st.selectbox(
                                "Spreadsheet Column:",
                                col_options,
                                index=def_col_idx,
                                key=col_key,
                                help=f"Spreadsheet column for {user_marker_name}"
                            )
                            if new_sel_col != sel_col:
                                sel_col = new_sel_col
                                if sel_col in df_raw.columns:
                                    raw_s = df_raw[sel_col].dropna().astype(str).str.strip().str.upper()
                                    valid_s = raw_s[~raw_s.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])]
                                    norm_calls = valid_s.map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                                    observed_counts = dict(norm_calls.value_counts())
                                    detected_calls_list = list(observed_counts.keys())
                                    tot_samples_col = len(norm_calls)
                        with c_top3:
                            st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
                            st.button("🗑️ Remove", key=f"btn_remove_marker_{idx}", on_click=cb_remove_marker_card, args=(idx,), help=f"Remove marker from mapping")

                        # Ref and Var Alleles Row
                        c2_m, c3_m = st.columns(2)
                        with c2_m:
                            ref_choices = ['A', 'C', 'G', 'T', 'Other (Manual)']
                            cur_ref = st.session_state.get(ref_key, def_ref)
                            if cur_ref not in ref_choices and cur_ref != 'Other (Manual)':
                                ref_choices.insert(0, cur_ref)
                            ref_idx = ref_choices.index(cur_ref) if cur_ref in ref_choices else (ref_choices.index(def_ref) if def_ref in ref_choices else 0)
                            sel_ref_choice = st.selectbox("Reference allele:", ref_choices, index=ref_idx, key=f"sel_ref_choice_{idx}", help="Reference allele")
                            if sel_ref_choice == 'Other (Manual)':
                                ref_val = st.text_input("Custom Ref Allele:", value=cur_ref if cur_ref != 'Other (Manual)' else '', key=f"custom_ref_in_{idx}").strip().upper()
                            else:
                                ref_val = sel_ref_choice
                            st.session_state[ref_key] = ref_val

                        with c3_m:
                            var_choices = ['A', 'C', 'G', 'T', 'None', 'Other (Manual)']
                            cur_var = st.session_state.get(var_key, def_var)
                            if cur_var not in var_choices and cur_var != 'Other (Manual)':
                                var_choices.insert(0, cur_var)
                            var_idx = var_choices.index(cur_var) if cur_var in var_choices else (var_choices.index(def_var) if def_var in var_choices else 0)
                            sel_var_choice = st.selectbox("Variant allele:", var_choices, index=var_idx, key=f"sel_var_choice_{idx}", help="Variant allele")
                            if sel_var_choice == 'Other (Manual)':
                                var_val = st.text_input("Custom Var Allele:", value=cur_var if cur_var != 'Other (Manual)' else '', key=f"custom_var_in_{idx}").strip().upper()
                            else:
                                var_val = sel_var_choice
                            st.session_state[var_key] = var_val

                        # Permitted Genotypes Input Row
                        c_perm_in, c_perm_info = st.columns([1.6, 1.4])
                        with c_perm_in:
                            perm_val = st.text_input(
                                "Permitted Genotypes (comma-separated):",
                                value=st.session_state.get(perm_key, def_permitted),
                                placeholder="e.g. CC, CT, TT or GA, GG, AA (or click chips below)",
                                key=perm_key,
                                help="Type any combinations e.g. CC, CT or GG, GA or AA or TT, or click the observed genotype chips below"
                            )
                        with c_perm_info:
                            hom_ref_str = f"{ref_val}{ref_val}" if ref_val else "N/A"
                            het_str = f"{ref_val}{var_val}" if (ref_val and var_val and var_val != 'None') else "N/A"
                            hom_var_str = f"{var_val}{var_val}" if (var_val and var_val != 'None') else "N/A"
                            st.markdown(f"""
                            <div style="background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px; margin-top: 5px;">
                                <div style="font-size: 11px; font-weight: 700; color: #64748B; text-transform: uppercase;">🧬 Genotype-to-Allele Mapping:</div>
                                <div style="font-size: 13px; color: #0F172A; margin-top: 4px;">
                                    Hom Ref = <b>{hom_ref_str}</b> &nbsp;|&nbsp; 
                                    Het = <b>{het_str}</b> &nbsp;|&nbsp; 
                                    Hom Var = <b>{hom_var_str}</b>
                                </div>
                            </div>
                            """, unsafe_allow_html=True)

                        # Interactive Call Toggle Chips (Manual Selection directly from Observed Data)
                        if detected_calls_list:
                            st.markdown("<div style='font-size: 11px; font-weight: 700; color: #475569; margin-top: 6px; margin-bottom: 4px; text-transform: uppercase;'>Observed Calls in this Column — Click chip to toggle in/out of permitted rules:</div>", unsafe_allow_html=True)
                            cur_perm_tokens_set = {x.strip().upper().replace('/', '').replace('|', '').replace(' ', '') for x in str(st.session_state.get(perm_key, '')).split(',') if x.strip()}

                            n_items = len(detected_calls_list)
                            chunk_size = 5
                            for c_start in range(0, n_items, chunk_size):
                                c_chunk = detected_calls_list[c_start:c_start + chunk_size]
                                is_first_row = (c_start == 0)

                                col_widths = [1] * len(c_chunk)
                                if is_first_row:
                                    col_widths.extend([1.3, 0.9])
                                    chip_cols = st.columns(col_widths)
                                    for c_i, call_item in enumerate(c_chunk):
                                        cnt = observed_counts.get(call_item, 0)
                                        is_active = call_item in cur_perm_tokens_set
                                        btn_lbl = f"✅ {call_item} ({cnt:,})" if is_active else f"➕ {call_item} ({cnt:,})"
                                        chip_cols[c_i].button(
                                            btn_lbl,
                                            key=f"btn_chip_{idx}_{call_item}",
                                            type="primary" if is_active else "secondary",
                                            use_container_width=True,
                                            on_click=cb_toggle_permitted_call,
                                            args=(perm_key, call_item, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}"),
                                            help=f"Click to {'remove' if is_active else 'permit'} {call_item}"
                                        )
                                    chip_cols[-2].button("⚡ Allow All", key=f"btn_chip_all_{idx}", use_container_width=True,
                                                        on_click=cb_allow_all_detected, args=(perm_key, detected_calls_list, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}"), help="Add all observed calls to permitted list")
                                    chip_cols[-1].button("🧹 Clear", key=f"btn_chip_clr_{idx}", use_container_width=True,
                                                        on_click=cb_clear_permitted, args=(perm_key,), help="Clear permitted list")
                                else:
                                    chip_cols = st.columns(col_widths)
                                    for c_i, call_item in enumerate(c_chunk):
                                        cnt = observed_counts.get(call_item, 0)
                                        is_active = call_item in cur_perm_tokens_set
                                        btn_lbl = f"✅ {call_item} ({cnt:,})" if is_active else f"➕ {call_item} ({cnt:,})"
                                        chip_cols[c_i].button(
                                            btn_lbl,
                                            key=f"btn_chip_{idx}_{call_item}",
                                            type="primary" if is_active else "secondary",
                                            use_container_width=True,
                                            on_click=cb_toggle_permitted_call,
                                            args=(perm_key, call_item, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}"),
                                            help=f"Click to {'remove' if is_active else 'permit'} {call_item}"
                                        )

                        # Parse permitted genotypes & expanded set
                        raw_perm_list = [g.strip().upper() for g in perm_val.split(',') if g.strip()]
                        valid_clean_set = {g.replace('/', '').replace('|', '').replace(' ', '') for g in raw_perm_list}
                        expanded_set = set(valid_clean_set)
                        for g in valid_clean_set:
                            if len(g) == 2:
                                expanded_set.add(g[1] + g[0])

                        # LIVE DATA AUDIT & STRICT DETECTION / FLAGGING ENGINE
                        # CHECK 1: Detect Non-Genomic / Demographic / Text Column
                        METADATA_WORDS = {
                            'MALE', 'FEMALE', 'MAN', 'WOMAN', 'UNKNOWN', 'NULL', 'NONE', 'NAN',
                            'CASE', 'CONTROL', 'YES', 'NO', 'TRUE', 'FALSE', 'POS', 'NEG',
                            'POSITIVE', 'NEGATIVE', 'NORMAL', 'MUTANT', 'HETEROZYGOUS', 'HOMOZYGOUS',
                            'INDIA', 'DELHI', 'PUNJAB', 'MAHARASHTRA', 'TAMIL', 'GUJARAT', 'KERALA',
                            'PATIENT', 'SAMPLE', 'CONTROL_GROUP', 'CASE_GROUP'
                        }
                        suspicious_calls = []
                        for c in detected_calls_list:
                            c_u = str(c).strip().upper()
                            if c_u in METADATA_WORDS or len(c_u) > 6 or any(ch.isdigit() for ch in c_u) or any(ch not in 'ACGTID/*|- ' for ch in c_u):
                                suspicious_calls.append(c)

                        is_metadata_col = bool(suspicious_calls and (len(suspicious_calls) >= len(detected_calls_list) * 0.5 or any(c.strip().upper() in METADATA_WORDS for c in detected_calls_list)))
                        if is_metadata_col:
                            mapping_has_errors = True
                            st.markdown(f"""
                            <div style="background: #FEF2F2; border: 2px solid #EF4444; border-radius: 8px; padding: 12px; margin-top: 10px; margin-bottom: 8px;">
                                <div style="font-weight: 800; color: #991B1B; font-size: 14px;">🚨 WRONG COLUMN DETECTED: Non-Genotype Data in `{sel_col}`</div>
                                <div style="color: #7F1D1D; font-size: 13px; margin-top: 4px;">
                                    Column <code>{sel_col}</code> appears to contain demographic or text values (e.g. <b>{suspicious_calls[:4]}</b>) rather than genetic marker calls. Please select a valid genotype column above.
                                </div>
                            </div>
                            """, unsafe_allow_html=True)

                        # CHECK 2: Severe Allele Disconnect / Wrong Marker Column
                        if raw_perm_list:
                            col_nucleotides = {ch for call in detected_calls_list for ch in str(call).upper() if ch in 'ACGT'}
                            configured_alleles = {ch for ch in [ref_val, var_val] if ch in 'ACGT'}
                            is_severe_disconnect = bool(col_nucleotides and configured_alleles and not col_nucleotides.intersection(configured_alleles))
                            if is_severe_disconnect and not is_metadata_col:
                                st.markdown(f"""
                                <div style="background: #FFF7ED; border: 2px solid #F97316; border-radius: 8px; padding: 12px; margin-top: 10px; margin-bottom: 8px;">
                                    <div style="font-weight: 800; color: #9A3412; font-size: 14px;">🚨 SEVERE ALLELE MISMATCH: Possible Wrong Column Selected!</div>
                                    <div style="color: #7C2D12; font-size: 13px; margin-top: 4px;">
                                        Observed calls in <code>{sel_col}</code> are <b>{detected_calls_list}</b> (nucleotides: <b>{sorted(col_nucleotides)}</b>).<br>
                                        These have <b>ZERO overlap</b> with configured alleles for <b>{user_marker_name}</b> (Ref: <b>{ref_val}</b>, Var: <b>{var_val}</b>).<br>
                                        Please verify if you selected the wrong column, or click below to align your alleles to this column if intended.
                                    </div>
                                </div>
                                """, unsafe_allow_html=True)
                                cand_col_r, cand_col_v = extract_alleles_from_calls(detected_calls_list)
                                st.button(
                                    f"⚡ Align Alleles to this Column (Ref: {cand_col_r}, Var: {cand_col_v})",
                                    key=f"btn_align_severe_{idx}",
                                    on_click=cb_align_alleles,
                                    args=(ref_key, var_key, cand_col_r, cand_col_v, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}")
                                )

                        # CHECK 3: Extra / Unmapped Calls in Dataset
                        unexpected_calls = [k for k in observed_counts.keys() if k not in expanded_set] if (observed_counts and raw_perm_list) else []
                        tot_unmapped_cnt = sum(observed_counts[u] for u in unexpected_calls) if unexpected_calls else 0
                        unmapped_pct = (tot_unmapped_cnt / tot_samples_col * 100) if tot_samples_col > 0 else 0

                        cand_r, cand_v = extract_alleles_from_calls(raw_perm_list)
                        perm_chars = {ch for g in raw_perm_list for ch in g if ch.isalpha()}
                        ref_mismatch = bool(ref_val and ref_val not in perm_chars) if raw_perm_list else False
                        var_mismatch = bool(var_val and var_val != 'None' and var_val not in perm_chars) if raw_perm_list else False

                        identical_alleles = bool(ref_val == var_val and ref_val != '' and var_val != 'None')
                        if identical_alleles:
                            mapping_has_errors = True

                        st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)
                        c_flg1, c_flg2 = st.columns([1.5, 1.5])
                        with c_flg1:
                            st.markdown(f"**Data Audit for `{sel_col}` ({tot_samples_col:,} valid samples):**")
                            if observed_counts:
                                check_rows = []
                                for call, cnt in sorted(observed_counts.items(), key=lambda x: x[1], reverse=True):
                                    pct = (cnt / tot_samples_col * 100) if tot_samples_col > 0 else 0
                                    is_covered = call in expanded_set
                                    if is_covered:
                                        status_str = "✅ Permitted"
                                    elif not raw_perm_list:
                                        status_str = "⚪ Pending User Entry"
                                    else:
                                        status_str = "🚨 Flagged: Unmapped (Extra)"
                                    check_rows.append({
                                        'Observed Genotype': call,
                                        'Sample Count': cnt,
                                        'Cohort %': f"{pct:.1f}%",
                                        'Mapping Status': status_str
                                    })
                                st.dataframe(pd.DataFrame(check_rows), use_container_width=True, hide_index=True)
                            else:
                                st.info(f"ℹ️ No valid genotype calls found in column `{sel_col}`.")

                        with c_flg2:
                            if not raw_perm_list:
                                st.info(f"ℹ️ **Action Needed:** Please enter permitted genotypes for **{user_marker_name}** above (type calls or click **[⚡ Allow All]** / chip buttons).")
                                mapping_has_errors = True
                            elif unexpected_calls:
                                st.markdown(f"""
                                <div style="background: #FEF2F2; border: 1.5px solid #F87171; border-radius: 8px; padding: 12px; margin-bottom: 8px;">
                                    <div style="font-weight: 700; color: #991B1B; font-size: 13px;">🚨 FLAGGED: Extra / Unmapped Entry in Dataset</div>
                                    <div style="color: #7F1D1D; font-size: 12px; margin-top: 4px;">
                                        Found <b>{len(unexpected_calls)} unmapped call(s)</b>: <code>{unexpected_calls}</code> ({tot_unmapped_cnt:,} sample(s), {unmapped_pct:.1f}% of cohort).<br>
                                        Permitted genotypes are: <code>{perm_val}</code>.<br>
                                        <i>These samples will be treated as invalid/missing during QC unless you explicitly allow them.</i>
                                    </div>
                                </div>
                                """, unsafe_allow_html=True)
                                unexp_str = ", ".join(unexpected_calls)
                                st.button(
                                    f"➕ Allow Flagged Call(s) ({unexp_str}) into Permitted Genotypes",
                                    key=f"btn_add_unexp_{idx}",
                                    on_click=cb_add_unmapped,
                                    args=(perm_key, unexpected_calls, ref_key, var_key, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}")
                                )
                            elif tot_samples_col > 0:
                                st.markdown(f"""
                                <div style="background: #F0FDF4; border: 1.5px solid #86EFAC; border-radius: 8px; padding: 12px; margin-bottom: 8px;">
                                    <div style="font-weight: 700; color: #166534; font-size: 13px;">✅ 100% Cohort Verified: All {tot_samples_col:,} Samples Permitted</div>
                                    <div style="color: #14532D; font-size: 12px; margin-top: 3px;">
                                        All observed calls in <code>{sel_col}</code> match your configured genotype rules (<code>{', '.join(raw_perm_list)}</code>).
                                    </div>
                                </div>
                                """, unsafe_allow_html=True)

                            # Flag: Allele Inconsistency
                            if (ref_mismatch or var_mismatch) and cand_r:
                                st.markdown(f"""
                                <div style="background: #FFFBEB; border: 1.5px solid #FCD34D; border-radius: 8px; padding: 10px; margin-bottom: 8px;">
                                    <div style="font-weight: 700; color: #92400E; font-size: 13px;">⚠️ FLAG: Allele Rule Inconsistency</div>
                                    <div style="color: #78350F; font-size: 12px; margin-top: 3px;">
                                        Configured Ref (<b>{ref_val}</b>) or Var (<b>{var_val}</b>) allele is absent from permitted genotypes. Expected alleles: Ref = <b>{cand_r}</b>, Var = <b>{cand_v}</b>.
                                    </div>
                                </div>
                                """, unsafe_allow_html=True)
                                st.button(
                                    f"👉 Align Alleles to Ref: {cand_r}, Var: {cand_v}",
                                    key=f"btn_align_alleles_{idx}",
                                    on_click=cb_align_alleles,
                                    args=(ref_key, var_key, cand_r, cand_v, f"sel_ref_choice_{idx}", f"sel_var_choice_{idx}")
                                )

                            # Flag: Identical Alleles
                            if identical_alleles:
                                st.error(f"⛔ FLAG: Reference allele ({ref_val}) and Variant allele ({var_val}) cannot be identical!")

                        v_al_stored = 'A' if var_val == 'None' else var_val.strip().upper()
                        mapped_marker_dict[sel_col] = {
                            'std_name': user_marker_name,
                            'ref_allele': ref_val.strip().upper(),
                            'var_allele': v_al_stored,
                            'marker_label': user_marker_name,
                            'permitted_genotypes': raw_perm_list
                        }

                    # Button to add new marker card row
                    st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                    c_add_card, _ = st.columns([0.4, 0.6])
                    with c_add_card:
                        st.button("➕ Add Marker Mapping Row", key="btn_add_marker_card_row", on_click=cb_add_new_marker_card, args=(sel_map_gene,))

                    # Render Extra Custom Mapping Rows
                    if extra_rows:
                        st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                        st.markdown("#### ➕ Custom Marker Mapping Rows")
                        rows_to_remove = []
                        for r_idx, row_item in enumerate(extra_rows):
                            r_id = row_item.get('id', f"r_{r_idx}")
                            st.markdown(f"""
                            <div class="clean-card-box" style="margin-bottom: 8px; border-left: 4px solid #8B5CF6;">
                                <div style="display: flex; justify-content: space-between; align-items: center;">
                                    <h5 style="margin:0; color:#0F172A; font-weight:700;">🧬 Custom Row #{r_idx+1}: {row_item.get('marker_name', 'Custom Marker')}</h5>
                                    <span style="background-color: #F3E8FF; color: #7E22CE; font-size: 0.8rem; font-weight: 600; padding: 2px 8px; border-radius: 9999px;">Manual Entry</span>
                                </div>
                            </div>
                            """, unsafe_allow_html=True)

                            c_r1, c_r2, c_r3, c_r4, c_del = st.columns([1.4, 1.4, 0.6, 0.6, 0.4])
                            with c_r1:
                                c_name = st.text_input("Marker Name", value=row_item.get('marker_name', f"Marker_{r_idx+1}"), key=f"ex_name_{r_id}")
                                row_item['marker_name'] = c_name
                            with c_r2:
                                def_col_ex = row_item.get('col', col_options[0] if col_options else '')
                                col_idx_ex = col_options.index(def_col_ex) if def_col_ex in col_options else 0
                                c_col = st.selectbox("Genotype Column", col_options, index=col_idx_ex, key=f"ex_col_{r_id}")
                                row_item['col'] = c_col
                            with c_r3:
                                c_ref = st.text_input("Ref Allele", value=row_item.get('ref_allele', 'A'), max_chars=4, key=f"ex_ref_{r_id}").strip().upper()
                                row_item['ref_allele'] = c_ref
                            with c_r4:
                                c_var = st.text_input("Var Allele", value=row_item.get('var_allele', 'G'), max_chars=4, key=f"ex_var_{r_id}").strip().upper()
                                row_item['var_allele'] = c_var
                            with c_del:
                                st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
                                if st.button("🗑️", key=f"ex_del_{r_id}", help="Delete this row"):
                                    rows_to_remove.append(r_idx)

                            # Live column detection for custom row
                            ex_observed_counts = {}
                            ex_detected_list = []
                            ex_norm = pd.Series(dtype=str)
                            if c_col in df_raw.columns:
                                raw_ex = df_raw[c_col].dropna().astype(str).str.strip().str.upper()
                                valid_ex = raw_ex[~raw_ex.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])]
                                ex_norm = valid_ex.map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                                ex_observed_counts = dict(ex_norm.value_counts())
                                ex_detected_list = list(ex_observed_counts.keys())

                            c_p_in, c_p_desc = st.columns([1.6, 1.4])
                            with c_p_in:
                                def_p_ex = row_item.get('permitted', f"{c_ref}{c_ref}, {c_ref}{c_var}, {c_var}{c_var}")
                                c_perm = st.text_input("Permitted Combinations (comma-separated):", value=def_p_ex, key=f"ex_perm_{r_id}")
                                row_item['permitted'] = c_perm

                                if ex_detected_list:
                                    st.button(
                                        f"⚡ Map Detected Calls ({', '.join(ex_detected_list)})",
                                        key=f"btn_automap_ex_{r_id}",
                                        on_click=cb_extra_automap,
                                        args=(f"ex_perm_{r_id}", f"ex_ref_{r_id}", f"ex_var_{r_id}", ex_detected_list)
                                    )

                            with c_p_desc:
                                st.markdown(f"""
                                <div style="font-size: 12px; color: #475569; padding-top: 28px;">
                                    🧬 <b>Genotype Mapping:</b> Hom Ref: <code>{c_ref}{c_ref}</code> | Het: <code>{c_ref}{c_var}</code> | Hom Var: <code>{c_var}{c_var}</code>
                                </div>
                                """, unsafe_allow_html=True)

                            # Custom row live check table
                            if c_col in df_raw.columns and ex_observed_counts:
                                allowed_c = {g.strip().upper().replace('/', '').replace('|', '').replace(' ', '') for g in c_perm.split(',') if g.strip()}
                                exp_c = set(allowed_c)
                                for g in allowed_c:
                                    if len(g) == 2:
                                        exp_c.add(g[1] + g[0])
                                c_unexp = [k for k in ex_observed_counts.keys() if k not in exp_c]

                                ex_chk_rows = []
                                for call, cnt in sorted(ex_observed_counts.items(), key=lambda x: x[1], reverse=True):
                                    pct = (cnt / len(ex_norm) * 100) if len(ex_norm) > 0 else 0
                                    ex_chk_rows.append({
                                        'Observed Genotype': call,
                                        'Sample Count': cnt,
                                        'Cohort %': f"{pct:.1f}%",
                                        'Mapping Status': "✅ Permitted" if call in exp_c else "🚨 Flagged: Unmapped"
                                    })
                                st.markdown(f"**Data Check for `{c_col}` ({len(ex_norm):,} valid samples):**")
                                st.dataframe(pd.DataFrame(ex_chk_rows), use_container_width=True, hide_index=True)

                                if c_unexp:
                                    st.warning(f"🚨 Detected unmapped calls: `{c_unexp}`.")
                                    st.button(
                                        f"➕ Add `{', '.join(c_unexp)}` to Permitted List",
                                        key=f"btn_add_unexp_ex_{r_id}",
                                        on_click=cb_extra_add_unmapped,
                                        args=(f"ex_perm_{r_id}", c_unexp, f"ex_ref_{r_id}", f"ex_var_{r_id}")
                                    )
                                else:
                                    st.success(f"✅ All {len(ex_norm):,} samples verified for `{c_name}`.")

                            mapped_marker_dict[c_col] = {
                                'std_name': c_name,
                                'ref_allele': c_ref,
                                'var_allele': c_var,
                                'marker_label': c_name,
                                'permitted_genotypes': [g.strip().upper() for g in c_perm.split(',') if g.strip()]
                            }

                        if rows_to_remove:
                            for idx_to_del in sorted(rows_to_remove, reverse=True):
                                extra_rows.pop(idx_to_del)
                            st.session_state['extra_mapping_rows'] = extra_rows
                            st.rerun()

                    # Button to add extra mapping row
                    st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)
                    if st.button("➕ Add Custom Marker Mapping Row", key="btn_add_extra_mapping_row"):
                        new_id = len(extra_rows) + 1
                        extra_rows.append({
                            'id': f"extra_{len(extra_rows)}_{pd.Timestamp.now().strftime('%H%M%S')}",
                            'marker_name': f"Custom_Marker_{new_id}",
                            'col': col_options[0] if col_options else '',
                            'ref_allele': 'A',
                            'var_allele': 'G',
                            'permitted': ''
                        })
                        st.session_state['extra_mapping_rows'] = extra_rows
                        st.rerun()

                # -------------------------------------------------------------------------
                # MODE B: MULTI-ROW TABLE EDITOR (FAST BATCH ENTRY WITH LIVE VERIFICATION)
                # -------------------------------------------------------------------------
                else:
                    st.markdown("#### 📋 Batch Multi-Row Mapping Table")
                    st.caption("Enter and edit all mapping rows in an interactive table. Directly edit cells, select columns, configure alleles, and add rows with '➕ Add row'.")

                    editor_data = []
                    for m_idx, marker_label in enumerate(sel_map_markers):
                        s_name, d_ref, d_var, d_perm, _, _ = get_marker_defaults(marker_label, sel_map_gene)
                        cur_c = st.session_state.get(f"w_map_col_{m_idx}", find_best_column_match(marker_label, s_name, sel_map_gene, col_options))
                        cur_r = st.session_state.get(f"w_map_ref_{m_idx}", d_ref)
                        cur_v = st.session_state.get(f"w_map_var_{m_idx}", d_var)
                        cur_p = st.session_state.get(f"w_map_perm_{m_idx}", d_perm)
                        cur_n = st.session_state.get(f"w_map_name_{m_idx}", s_name)
                        editor_data.append({
                            'Marker Name': cur_n,
                            'Spreadsheet Column': cur_c if cur_c in col_options else col_options[0],
                            'Reference Allele': cur_r,
                            'Variant Allele': cur_v,
                            'Permitted Genotypes': cur_p
                        })

                    for r_idx, row_item in enumerate(extra_rows):
                        cur_c_ex = row_item.get('col', col_options[0])
                        editor_data.append({
                            'Marker Name': row_item.get('marker_name', f"Custom_{r_idx+1}"),
                            'Spreadsheet Column': cur_c_ex if cur_c_ex in col_options else col_options[0],
                            'Reference Allele': row_item.get('ref_allele', 'A'),
                            'Variant Allele': row_item.get('var_allele', 'G'),
                            'Permitted Genotypes': row_item.get('permitted', '')
                        })

                    init_table_df = pd.DataFrame(editor_data)

                    edited_table = st.data_editor(
                        init_table_df,
                        num_rows="dynamic",
                        use_container_width=True,
                        column_config={
                            'Marker Name': st.column_config.TextColumn("Marker Name", required=True),
                            'Spreadsheet Column': st.column_config.SelectboxColumn("Spreadsheet Column", options=col_options, required=True),
                            'Reference Allele': st.column_config.TextColumn("Reference Allele", required=True),
                            'Variant Allele': st.column_config.TextColumn("Variant Allele", required=True),
                            'Permitted Genotypes': st.column_config.TextColumn("Permitted Genotypes (comma-separated)", required=True)
                        },
                        key="batch_data_editor_step3"
                    )

                    c_app1, c_app2 = st.columns([0.3, 0.7])
                    with c_app1:
                        if st.button("💾 Apply Table Mappings", type="secondary", key="btn_apply_table_step3"):
                            new_extra = []
                            for t_idx, row_val in edited_table.iterrows():
                                m_name = str(row_val.get('Marker Name', '')).strip()
                                m_col = str(row_val.get('Spreadsheet Column', '')).strip()
                                m_ref = str(row_val.get('Reference Allele', 'A')).strip().upper()
                                m_var = str(row_val.get('Variant Allele', 'G')).strip().upper()
                                m_perm = str(row_val.get('Permitted Genotypes', '')).strip().upper()

                                if t_idx < len(sel_map_markers):
                                    st.session_state[f"w_map_name_{t_idx}"] = m_name
                                    st.session_state[f"w_map_col_{t_idx}"] = m_col
                                    st.session_state[f"w_map_ref_{t_idx}"] = m_ref
                                    st.session_state[f"w_map_var_{t_idx}"] = m_var
                                    st.session_state[f"w_map_perm_{t_idx}"] = m_perm
                                else:
                                    new_extra.append({
                                        'id': f"extra_{t_idx}_{pd.Timestamp.now().strftime('%H%M%S')}",
                                        'marker_name': m_name,
                                        'col': m_col,
                                        'ref_allele': m_ref,
                                        'var_allele': m_var,
                                        'permitted': m_perm
                                    })
                            st.session_state['extra_mapping_rows'] = new_extra
                            st.success("All table mappings saved and synchronized!")
                            st.rerun()

                    for _, row_val in edited_table.iterrows():
                        m_name = str(row_val.get('Marker Name', '')).strip()
                        m_col = str(row_val.get('Spreadsheet Column', '')).strip()
                        m_ref = str(row_val.get('Reference Allele', 'A')).strip().upper()
                        m_var = str(row_val.get('Variant Allele', 'G')).strip().upper()
                        m_perm = str(row_val.get('Permitted Genotypes', '')).strip().upper()
                        if m_col and m_col in df_raw.columns:
                            mapped_marker_dict[m_col] = {
                                'std_name': m_name,
                                'ref_allele': m_ref,
                                'var_allele': 'A' if m_var == 'None' else m_var,
                                'marker_label': m_name,
                                'permitted_genotypes': [g.strip().upper() for g in m_perm.split(',') if g.strip()]
                            }

                    # Live Batch Data Verification Table
                    st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                    st.markdown("##### 📊 Batch Live Data Verification Across All Mapped Columns")
                    batch_check_rows = []
                    for _, r_val in edited_table.iterrows():
                        m_nm = str(r_val.get('Marker Name', '')).strip()
                        m_cl = str(r_val.get('Spreadsheet Column', '')).strip()
                        m_rf = str(r_val.get('Reference Allele', 'A')).strip().upper()
                        m_vr = str(r_val.get('Variant Allele', 'G')).strip().upper()
                        m_pm = str(r_val.get('Permitted Genotypes', '')).strip().upper()

                        if m_cl in df_raw.columns:
                            raw_s = df_raw[m_cl].dropna().astype(str).str.strip().str.upper()
                            valid_s = raw_s[~raw_s.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])]
                            norm_c = valid_s.map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                            c_cnts = dict(norm_c.value_counts())
                            tot_c = len(norm_c)

                            pm_list = [g.strip().upper() for g in m_pm.split(',') if g.strip()]
                            pm_set = {g.replace('/', '').replace('|', '').replace(' ', '') for g in pm_list}
                            for g in list(pm_set):
                                if len(g) == 2:
                                    pm_set.add(g[1] + g[0])

                            unexp = [k for k in c_cnts.keys() if k not in pm_set]
                            unexp_cnt = sum(c_cnts[k] for k in unexp)
                            dist_str = ", ".join([f"{k} ({v})" for k, v in sorted(c_cnts.items())])

                            status_str = f"✅ Verified ({tot_c:,} samples)" if not unexp else f"🚨 Flagged: {len(unexp)} Unmapped ({unexp_cnt} samples)"
                            batch_check_rows.append({
                                'Marker': m_nm,
                                'Column': m_cl,
                                'Ref / Var': f"{m_rf} / {m_vr}",
                                'Permitted Combinations': m_pm,
                                'Valid Samples': tot_c,
                                'Observed Call Breakdown': dist_str,
                                'QC Status': status_str
                            })

                    if batch_check_rows:
                        st.dataframe(pd.DataFrame(batch_check_rows), use_container_width=True, hide_index=True)

                # Check for Duplicate Column Mapping Conflicts
                mapped_cols_list = list(mapped_marker_dict.keys())
                dup_cols = {c for c in mapped_cols_list if mapped_cols_list.count(c) > 1}
                if dup_cols:
                    st.error(f"⛔ **Duplicate Column Conflict:** Spreadsheet column(s) `{dup_cols}` are mapped to multiple markers! Each biomarker must be mapped to a distinct column.")
                    mapping_has_errors = True

                # Check for Identical Reference and Variant Alleles
                for col_name, m_info in mapped_marker_dict.items():
                    r_al = m_info.get('ref_allele', '')
                    v_al = m_info.get('var_allele', '')
                    if r_al and v_al and r_al == v_al and v_al != 'None':
                        st.error(f"⛔ For marker **{m_info.get('marker_label', col_name)}**, Reference allele ({r_al}) and Variant allele ({v_al}) cannot be identical.")
                        mapping_has_errors = True

                # Check for unmapped calls across all mapped columns to display summary advisory
                total_unmapped_samples = 0
                for col_name, m_info in mapped_marker_dict.items():
                    if col_name in df_raw.columns:
                        p_list = m_info.get('permitted_genotypes', [])
                        p_set = {g.replace('/', '').replace('|', '').replace(' ', '') for g in p_list if g.strip()}
                        for g in list(p_set):
                            if len(g) == 2: p_set.add(g[1] + g[0])
                        r_s = df_raw[col_name].dropna().astype(str).str.strip().str.upper()
                        v_s = r_s[~r_s.isin(['', 'NAN', 'NONE', 'NULL', 'NA', '?', 'MISSING', 'UNKNOWN'])].map(lambda v: v.replace('/', '').replace('|', '').replace(' ', ''))
                        unm_cnt = sum(cnt for call, cnt in v_s.value_counts().items() if call not in p_set)
                        total_unmapped_samples += unm_cnt

                if total_unmapped_samples > 0:
                    st.warning(f"⚠️ **Advisory:** Detected {total_unmapped_samples:,} sample call(s) across your dataset that are not permitted by your genotype rules. Any unpermitted calls will be flagged during QC audit and treated as missing.")

                st.markdown("<br>", unsafe_allow_html=True)
                c_b1, c_b2, c_b3 = st.columns([0.25, 0.5, 0.25])
                with c_b1:
                    if st.button("← Back to Upload", use_container_width=True, key="btn_map_back"):
                        st.session_state['workflow_step'] = '2 · Upload'
                        st.rerun()
                with c_b3:
                    if mapping_has_errors:
                        st.button("Continue to Results →", type="primary", disabled=True, use_container_width=True, key="btn_map_cont_dis", help="Resolve mapping errors above to proceed")
                    elif st.button("Continue to Results →", type="primary", use_container_width=True, key="btn_map_cont"):
                        # BUILD FINAL MAPPING WITH CONFIRMED METADATA & ALL MARKERS
                        final_mapping = {}
                        demo = st.session_state.get('demographic_mapping', {})
                        if demo.get('Sample_ID') and demo['Sample_ID'] in all_cols:
                            final_mapping[demo['Sample_ID']] = 'Sample_ID'
                        if demo.get('Gender') and demo['Gender'] in all_cols:
                            final_mapping[demo['Gender']] = 'Gender'
                        if demo.get('State') and demo['State'] in all_cols:
                            final_mapping[demo['State']] = 'State'

                        custom_rules_for_qc = {}
                        for m_col, m_info in mapped_marker_dict.items():
                            final_mapping[m_col] = m_info['std_name']
                            perm_list = m_info.get('permitted_genotypes', [])
                            if perm_list:
                                clean_perms = set(g.replace('/', '').replace('|', '').replace(' ', '') for g in perm_list if g.strip())
                                augmented_p = set(clean_perms)
                                for g in clean_perms:
                                    if len(g) == 2:
                                        augmented_p.add(g[1] + g[0])
                                custom_rules_for_qc[m_info['std_name']] = {
                                    'ref_allele': m_info['ref_allele'],
                                    'var_allele': m_info['var_allele'],
                                    'valid_genotypes': tuple(sorted(augmented_p))
                                }
                            else:
                                r_a = m_info['ref_allele']
                                v_a = m_info['var_allele']
                                dyn_perms = {r_a + r_a, r_a + v_a, v_a + r_a, v_a + v_a}
                                custom_rules_for_qc[m_info['std_name']] = {
                                    'ref_allele': r_a,
                                    'var_allele': v_a,
                                    'valid_genotypes': tuple(sorted(dyn_perms))
                                }

                        st.session_state['mapped_cols'] = final_mapping
                        st.session_state['mapped_marker_details'] = mapped_marker_dict
                        st.session_state['custom_snp_rules'] = custom_rules_for_qc

                        # RESTRICT DATASET: RETAIN ONLY MAPPED COLUMNS & ENFORCE PERMITTED GENOTYPES
                        if df_raw is not None and not df_raw.empty:
                            cols_to_keep = [c for c in df_raw.columns if c in final_mapping]
                            if cols_to_keep:
                                filtered_df = df_raw[cols_to_keep].copy()

                                for m_col, m_info in mapped_marker_dict.items():
                                    if m_col in filtered_df.columns:
                                        user_permitted = m_info.get('permitted_genotypes', [])
                                        if user_permitted:
                                            valid_set = {str(g).strip().upper().replace('/', '').replace('|', '').replace(' ', '') for g in user_permitted if str(g).strip()}
                                            augmented = set(valid_set)
                                            for g in valid_set:
                                                if len(g) == 2:
                                                    augmented.add(g[1] + g[0])
                                            def clean_user_calls(v):
                                                if pd.isna(v): return np.nan
                                                clean = str(v).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
                                                if clean in valid_set: return clean
                                                if clean in augmented:
                                                    rev = clean[1] + clean[0] if len(clean) == 2 else clean
                                                    if rev in valid_set: return rev
                                                    return clean
                                                return np.nan
                                            filtered_df[m_col] = filtered_df[m_col].apply(clean_user_calls)
                                        else:
                                            r_a = m_info['ref_allele']
                                            v_a = m_info['var_allele']
                                            valid_dyn = {r_a + r_a, r_a + v_a, v_a + r_a, v_a + v_a}
                                            def clean_dyn(v):
                                                if pd.isna(v): return np.nan
                                                clean = str(v).strip().upper().replace('/', '').replace('|', '').replace(' ', '')
                                                if clean in valid_dyn: return clean
                                                return np.nan
                                            filtered_df[m_col] = filtered_df[m_col].apply(clean_dyn)

                                st.session_state['uploaded_df'] = filtered_df
                                save_active_dataset(filtered_df)

                        st.session_state['workflow_step'] = '4 · Results'
                        st.rerun()

        # -------------------------------------------------------------------------
        # STEP 4: RESULTS (GENE & BIOMARKER INTERACTIVE VISUALIZATION STUDIO)
        # -------------------------------------------------------------------------
        elif st.session_state['workflow_step'] == '4 · Results':
            active_gene = full_results.get('active_gene') or st.session_state.get('selected_gene', 'CYP2C19')

            # Render Full 9-Tab Studio if CYP2C19 or CYP2C9 was selected and analyzed
            if active_gene in ('CYP2C19', 'CYP2C9'):
                render_unified_results(full_results, qc_report)
            else:
                st.markdown(f"## 📊 Genomic Studio & Biomarker Analysis for {active_gene}")
                mapped_details = st.session_state.get('mapped_marker_details', {})
                if df_raw is not None and not df_raw.empty and mapped_details:
                    st.caption(f"Interactive genotype analysis, allele frequencies, and Hardy-Weinberg equilibrium test for {len(mapped_details)} mapped biomarker(s).")

                    for m_col, m_info in mapped_details.items():
                        m_label = m_info.get('marker_label', m_info.get('std_name', m_col))
                        r_al = m_info.get('ref_allele', 'A')
                        v_al = m_info.get('var_allele', 'G')

                        res = analyze_generic_marker(df_raw[m_col], m_label, ref_al=r_al, var_al=v_al)

                        st.markdown(f"""
                        <div class="clean-card-box" style="margin-top: 15px; margin-bottom: 10px; border-left: 4px solid #2563EB;">
                            <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap;">
                                <h4 style="margin: 0; color: #0F172A; font-weight: 700;">🧬 {m_label}</h4>
                                <span style="background-color: #EFF6FF; color: #1D4ED8; font-size: 0.85rem; font-weight: 600; padding: 3px 12px; border-radius: 9999px;">
                                    {res.get('interpretation', 'N/A')}
                                </span>
                            </div>
                        </div>
                        """, unsafe_allow_html=True)

                        k1, k2, k3, k4 = st.columns(4)
                        k1.metric("Valid Samples", f"{res.get('valid_samples', 0):,}")
                        k2.metric(f"Ref Allele ({r_al}) Freq", f"{res.get('p', 0.0):.4f}")
                        k3.metric(f"Var Allele ({v_al}) Freq", f"{res.get('q', 0.0):.4f}")
                        k4.metric("HWE χ² Stat (p-val)", f"{res.get('hwe_chi2', 0.0):.2f} (p={res.get('hwe_p', 1.0):.4f})")

                        gt_counts = res.get('genotype_counts', {})
                        if gt_counts:
                            c_ch1, c_ch2 = st.columns([1.2, 1])
                            with c_ch1:
                                fig_bar = px.bar(
                                    x=list(gt_counts.keys()),
                                    y=list(gt_counts.values()),
                                    labels={'x': 'Genotype Call', 'y': 'Sample Count'},
                                    title=f"Genotype Distribution ({m_label})",
                                    color_discrete_sequence=['#2563EB']
                                )
                                fig_bar.update_layout(height=280, margin=dict(l=20, r=20, t=40, b=20))
                                st.plotly_chart(fig_bar, use_container_width=True)
                            with c_ch2:
                                gt_df = pd.DataFrame([
                                    {'Genotype': k, 'Count': v, 'Percentage': f"{(v / res['valid_samples'] * 100):.1f}%"}
                                    for k, v in gt_counts.items()
                                ])
                                st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                                st.dataframe(gt_df, use_container_width=True, hide_index=True)

                st.divider()
                st.markdown(f"### 📚 Genomic Knowledgebase Markers for {active_gene}")
                sub_catalog = catalog_df[catalog_df['Gene'].str.upper() == active_gene.upper()] if not catalog_df.empty else pd.DataFrame()
                if not sub_catalog.empty:
                    disp_cat = sub_catalog[['Gene', 'Marker', 'SA MAF', '%', 'MAF', 'Level of evidence', 'PubMed reference']].copy()
                    st.dataframe(disp_cat, use_container_width=True, hide_index=True)
                else:
                    st.info(f"No additional markers recorded for {active_gene}.")

            st.divider()
            st.markdown("### Export")
            st.caption("Download the full result set for the current group selection.")

            e1, e2, e3, e_remap = st.columns([0.25, 0.25, 0.25, 0.25])

            rep_gen = ReportGenerator(full_results, qc_report) if has_active_data() else None

            with e1:
                if rep_gen:
                    excel_bytes = rep_gen.generate_excel()
                    st.download_button("⬇ Excel (.xlsx, all sheets)", data=excel_bytes, file_name="PGx_Population_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True, key="btn_exp_excel")
                else:
                    st.button("⬇ Excel (.xlsx, all sheets)", disabled=True, use_container_width=True)

            with e2:
                if rep_gen:
                    csv_bytes = rep_gen.generate_csv()
                    st.download_button("⬇ CSV (summary)", data=csv_bytes, file_name="PGx_Population_Summary.csv", mime="text/csv", use_container_width=True, key="btn_exp_csv")
                else:
                    st.button("⬇ CSV (summary)", disabled=True, use_container_width=True)

            with e3:
                if rep_gen:
                    pdf_bytes = rep_gen.generate_pdf()
                    st.download_button("⬇ PDF (print view)", data=pdf_bytes, file_name="PGx_Population_Report.pdf", mime="application/pdf", use_container_width=True, key="btn_exp_pdf")
                else:
                    st.button("⬇ PDF (print view)", disabled=True, use_container_width=True)

            with e_remap:
                if st.button("← Re-map columns", use_container_width=True, key="btn_res_remap"):
                    st.session_state['workflow_step'] = '3 · Map columns'
                    st.rerun()

            if df_raw is not None and not df_raw.empty:
                st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)
                target_full_df = processed_df if (processed_df is not None and not processed_df.empty) else df_raw
                c_full_dl1, c_full_dl2 = st.columns(2)
                with c_full_dl1:
                    label_desc = "⬇️ Download Complete Classified Dataset (.xlsx)" if (processed_df is not None and not processed_df.empty) else "⬇️ Download Complete Dataset (.xlsx)"
                    st.download_button(
                        label=label_desc,
                        data=get_df_download_excel(target_full_df),
                        file_name=f"complete_pgx_dataset_{len(target_full_df)}_samples.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        key="btn_res_full_dl_xlsx"
                    )
                with c_full_dl2:
                    st.download_button(
                        label="⬇️ Download Complete Dataset (.csv)",
                        data=get_df_download_csv(target_full_df),
                        file_name=f"complete_pgx_dataset_{len(target_full_df)}_samples.csv",
                        mime="text/csv",
                        use_container_width=True,
                        key="btn_res_full_dl_csv"
                    )

            st.markdown("""
            <div class="footer-disclaimer-text">
                PGx Pop MVP · CYP2C19 calling uses a simplified unphased-genotype heuristic — see notes in the CYP2C19 tab · Not for clinical use without independent validation.
            </div>
            """, unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# MAIN WORKFLOW EXECUTION
# -----------------------------------------------------------------------------
render_workflow_pipeline(df_raw, full_results, qc_report)

# -----------------------------------------------------------------------------
# FOOTER NOTE
# -----------------------------------------------------------------------------
st.markdown("""
<br><hr>
<div style="font-size: 0.8rem; color: #94A3B8; text-align: center;">
    Automated Population Pharmacogenomics Analysis Platform — Built for Reproducible Genotype, Allele, HWE, and CPIC Phenotype Analysis.
</div>
""", unsafe_allow_html=True)
