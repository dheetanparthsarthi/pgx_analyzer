"""
Gene & Biomarker Detection and Utilities
Provides standardized discovery, parsing, and defaults for Pharmacogenetic targets.
"""
import os
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
    g = str(gene_symbol).upper().strip()
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
    return [f"{gene_symbol}_marker1"]
