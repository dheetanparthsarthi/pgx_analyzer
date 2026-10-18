"""
Gene & Biomarker Detection and Utilities
Provides standardized discovery, parsing, and defaults for Pharmacogenetic targets.
"""
import os
import re
from typing import List, Optional
import pandas as pd

KNOWN_PGX_GENES = [
    'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP3A4', 'VKORC1', 
    'DPYD', 'TPMT', 'NUDT15', 'SLCO1B1', 'UGT1A1', 'HLA-B', 'HLA-A', 'G6PD', 'CYP1A2', 'CYP2B6'
]

def load_marker_catalog_file() -> pd.DataFrame:
    """Loads markers from List of markers_v1.xlsx if present."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    excel_path = os.path.join(base_dir, "List of markers_v1.xlsx")
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

def get_column_gene(col_name: str) -> str:
    """Returns the gene symbol found in a column header, if any."""
    c_clean = str(col_name).upper().replace("_", "").replace("-", "").replace(" ", "")
    for g in sorted(KNOWN_PGX_GENES, key=len, reverse=True):
        g_clean = g.upper().replace("_", "").replace("-", "")
        if g_clean in c_clean:
            return g
    return ""

def detect_genes_in_columns(cols: list) -> List[str]:
    """Detects any known pharmacogenetic genes mentioned across column headers."""
    found_genes = []
    for c in cols:
        g = get_column_gene(c)
        if g and g not in found_genes:
            found_genes.append(g)
    return found_genes

def get_default_markers_for_gene(gene_symbol: str, catalog_df: Optional[pd.DataFrame] = None) -> List[str]:
    """Returns standard biomarker labels for a given gene."""
    g = str(gene_symbol).upper().strip() if gene_symbol else ""
    if g == 'CYP2C19':
        return ['CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 (rs12248560)']
    elif g == 'CYP2C9':
        return ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']

    if catalog_df is None or catalog_df.empty:
        catalog_df = load_marker_catalog_file()

    if catalog_df is not None and not catalog_df.empty:
        ms = catalog_df[catalog_df['Gene'].str.upper() == g]['Marker'].tolist()
        if ms:
            return ms
    return [f"{gene_symbol}_marker1"] if gene_symbol else []

def find_best_column_match(marker_label: str, std_name: str, sel_gene: str, cols: list) -> str:
    """Finds the best matching spreadsheet column for a given biomarker label without false fallbacks."""
    if not cols:
        return ""

    real_cols = [c for c in cols if c and str(c).strip() != "(Select Column)"]
    if not real_cols:
        return ""

    # 1. Exact or stripped match with label, name, or gene_marker
    targets = [marker_label, std_name, f"{sel_gene}_{std_name}", f"{sel_gene} {std_name}", f"{sel_gene}_{marker_label}"]
    for target in targets:
        t_clean = str(target).strip().lower().replace(" ", "").replace("_", "").replace("-", "")
        for c in real_cols:
            c_clean = str(c).strip().lower().replace(" ", "").replace("_", "").replace("-", "")
            if t_clean and (t_clean == c_clean):
                return c

    # 2. Extract rsID if present in marker_label or std_name
    rs_match = re.search(r'rs\d+', f"{marker_label} {std_name}".lower())
    if rs_match:
        rs_id = rs_match.group(0)
        for c in real_cols:
            if rs_id in str(c).lower():
                return c

    # 3. Check for gene + star allele e.g. CYP2C9 and *2
    star_match = re.search(r'\*\d+[a-zA-Z]?', f"{marker_label} {std_name}")
    if star_match:
        star_str = star_match.group(0).lower()
        for c in real_cols:
            c_low = str(c).lower()
            if sel_gene and sel_gene.lower() in c_low and star_str in c_low:
                return c
        for c in real_cols:
            cg = get_column_gene(c)
            if star_str in str(c).lower() and (not cg or (sel_gene and cg.upper() == sel_gene.upper())):
                return c

    # 4. Normalized substring match (only if length >= 3 and not conflicting with another gene)
    clean_targets = [str(t).lower().replace(" ", "").replace("_", "").replace("-", "").replace("*", "").replace("(", "").replace(")", "") for t in [marker_label, std_name] if t]
    for c in real_cols:
        c_clean = str(c).lower().replace(" ", "").replace("_", "").replace("-", "").replace("*", "").replace("(", "").replace(")", "")
        for t_clean in clean_targets:
            if t_clean and len(t_clean) >= 3 and len(c_clean) >= 3 and (t_clean in c_clean or c_clean in t_clean):
                cg = get_column_gene(c)
                if not cg or (sel_gene and cg.upper() == sel_gene.upper()):
                    return c

    return ""
