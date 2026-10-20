"""
Gene & Biomarker Detection and Utilities
Provides standardized discovery, parsing, and defaults for Pharmacogenetic targets.
"""
import os
import re
from functools import lru_cache
from typing import List, Optional
import pandas as pd

KNOWN_PGX_GENES = [
    'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP3A4', 'VKORC1', 
    'DPYD', 'TPMT', 'NUDT15', 'SLCO1B1', 'UGT1A1', 'HLA-B', 'HLA-A', 'G6PD', 'CYP1A2', 'CYP2B6'
]

@lru_cache(maxsize=1)
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
    genes = set(KNOWN_PGX_GENES)
    catalog = load_marker_catalog_file()
    if not catalog.empty and 'Gene' in catalog.columns:
        genes.update(catalog['Gene'].dropna().astype(str).str.strip().tolist())
    for g in sorted(genes, key=len, reverse=True):
        g_clean = g.upper().replace("_", "").replace("-", "")
        if g_clean in c_clean:
            return g

    # rsIDs provide a safer dynamic association than generic labels such as *2.
    rsids = set(re.findall(r'rs\d+', str(col_name).lower()))
    if rsids and not catalog.empty and {'Gene', 'Marker'}.issubset(catalog.columns):
        for _, item in catalog.iterrows():
            marker_rsids = set(re.findall(r'rs\d+', str(item.get('Marker', '')).lower()))
            if rsids.intersection(marker_rsids):
                gene = str(item.get('Gene', '')).strip()
                if gene:
                    return gene
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

def smart_merge_datasets(base_df: pd.DataFrame, new_df: pd.DataFrame, active_gene: Optional[str] = None) -> tuple[pd.DataFrame, list[str]]:
    """
    Intelligently aligns gene, marker, and demographic columns between base_df and new_df:
    - Demographic columns (Sample_ID, Gender, State, DOB) are aligned.
    - If datasets share the same target gene:
        - Common markers are merged into the same column.
        - Markers unique to new_df are preserved as separate distinct columns (outer merge).
    - If datasets belong to completely different genes (e.g. CYP2C9 vs CYP2C19),
      returns an advisory error log explaining that merging different genes is blocked.
    """
    if base_df is None or base_df.empty:
        return new_df.copy() if new_df is not None else pd.DataFrame(), ["Base dataset empty; adopted newly uploaded dataset."]
    if new_df is None or new_df.empty:
        return base_df.copy(), ["New dataset empty; retained current base dataset."]

    base_genes = detect_genes_in_columns(list(base_df.columns))
    new_genes = detect_genes_in_columns(list(new_df.columns))

    # Verify gene compatibility
    if base_genes and new_genes:
        shared_genes = [g for g in new_genes if any(g.upper() == bg.upper() for bg in base_genes)]
        if not shared_genes and (not active_gene or (active_gene.upper() not in [g.upper() for g in new_genes])):
            err_msg = f"⛔ Cannot merge: Base dataset represents **{', '.join(base_genes)}** but newly uploaded dataset represents **{', '.join(new_genes)}**. Merging is only permitted when datasets share the same gene. Please keep them as separate datasets."
            return base_df.copy(), [err_msg]

    new_aligned = new_df.copy()
    alignment_log = []
    rename_dict = {}

    # Helper for demographic matching
    demo_patterns = {
        'Sample_ID': ['sample id', 'sample_id', 'sample', 'patient id', 'patient_id', 'lab id', 'id'],
        'Gender': ['gender', 'sex'],
        'State': ['state / country', 'state/country', 'state', 'region', 'native place', 'native', 'location'],
        'DOB': ['date of birth', 'dob', 'birth']
    }

    def _find_demo_col(key: str, df: pd.DataFrame) -> Optional[str]:
        p_list = demo_patterns.get(key, [])
        for p in p_list:
            for c in df.columns:
                c_clean = str(c).lower().replace(" ", "").replace("_", "").replace("-", "")
                p_clean = p.lower().replace(" ", "").replace("_", "").replace("-", "")
                if p_clean in c_clean:
                    return c
        return None

    # 1. Align core metadata fields
    for demo_key in ['Sample_ID', 'Gender', 'State', 'DOB']:
        b_col = _find_demo_col(demo_key, base_df)
        n_col = _find_demo_col(demo_key, new_aligned)
        if b_col and n_col and b_col != n_col:
            rename_dict[n_col] = b_col
            alignment_log.append(f"📋 **{demo_key}**: Aligned new column `{n_col}` ➔ existing `{b_col}`")
        elif b_col and n_col and b_col == n_col:
            alignment_log.append(f"✅ **{demo_key}**: Identical column `{b_col}` matched")

    # 2. Check each marker column in new dataset against base dataset
    b_cols = list(base_df.columns)
    n_cols = list(new_df.columns)

    for n_c in n_cols:
        if n_c in rename_dict:
            continue
        g_cand = get_column_gene(n_c) or active_gene or ""
        b_match = find_best_column_match(n_c, n_c.split(' ')[0], g_cand, b_cols)
        if b_match and b_match != n_c:
            rename_dict[n_c] = b_match
            alignment_log.append(f"🧬 **Same Marker**: Aligned `{n_c}` ➔ `{b_match}` (merged into same column)")
        elif b_match and b_match == n_c:
            alignment_log.append(f"✅ **Same Marker**: Identical column `{n_c}` matched")
        else:
            # Different marker: will create separate column in merged DataFrame!
            alignment_log.append(f"➕ **Distinct Marker**: `{n_c}` added as a separate column")

    if rename_dict:
        new_aligned = new_aligned.rename(columns=rename_dict)

    merged_df = pd.concat([base_df, new_aligned], ignore_index=True)
    return merged_df, alignment_log
