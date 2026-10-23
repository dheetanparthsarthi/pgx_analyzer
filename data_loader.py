"""Dataset loading, validation, inspection, and merging utilities."""

import re
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd
from core.gene_utils import (
    KNOWN_PGX_GENES,
    detect_genes_in_columns,
    find_best_column_match,
    get_column_gene,
    get_default_markers_for_gene,
    smart_merge_datasets,
)


def find_matching_gene_column(target_key: str, df: pd.DataFrame) -> Optional[str]:
    """Finds best matching column for a given gene/SNP/demographic in a dataframe."""
    if df is None or df.empty:
        return None
    cols = list(df.columns)

    patterns = {
        "CYP2C19*2": ["cyp2c19*2", "rs4244285", "c.681g>a", "c.681", "cyp2c19_2", "cyp2"],
        "CYP2C19*3": ["cyp2c19*3", "rs4986893", "c.636g>a", "c.636", "cyp2c19_3", "cyp3"],
        "CYP2C19*17": ["*17", "cyp2c19*17", "rs12248560", "c.-806c>t", "c.-806", "cyp2c19_17", "cyp17"],
        "CYP2C9*2": ["cyp2c9*2", "rs1799853", "cyp2c9_2", "cyp2c9-2", "430c>t"],
        "CYP2C9*3": ["cyp2c9*3", "rs1057910", "cyp2c9_3", "cyp2c9-3", "1075a>c"],
        "Sample_ID": ["sample id", "sample_id", "sample", "patient id", "patient_id", "lab id", "id"],
        "Gender": ["gender", "sex"],
        "State": ["state / country", "state/country", "state", "region", "native place", "native", "location"],
        "DOB": ["date of birth", "dob", "birth"],
    }

    target_patterns = patterns.get(target_key, [])
    for p in target_patterns:
        p_clean = p.lower().replace(" ", "").replace("_", "").replace("-", "")
        for c in cols:
            c_clean = str(c).lower().replace(" ", "").replace("_", "").replace("-", "")
            if p_clean in c_clean:
                return c

    # Generic fallback: match target_key directly in column names
    tk_clean = str(target_key).lower().replace(" ", "").replace("_", "").replace("-", "")
    for c in cols:
        c_clean = str(c).lower().replace(" ", "").replace("_", "").replace("-", "")
        if tk_clean == c_clean or tk_clean in c_clean or c_clean in tk_clean:
            return c
    return None


def inspect_uploaded_file(df: pd.DataFrame) -> Dict[str, Any]:
    """Checks for Sample ID (compulsory), Gender, and State in uploaded dataframe."""
    sid_col = find_matching_gene_column("Sample_ID", df)
    gen_col = find_matching_gene_column("Gender", df)
    state_col = find_matching_gene_column("State", df)
    return {
        "has_sample_id": sid_col is not None,
        "sample_id_col": sid_col,
        "has_gender": gen_col is not None,
        "gender_col": gen_col,
        "has_state": state_col is not None,
        "state_col": state_col,
    }
