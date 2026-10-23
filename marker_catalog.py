"""Backward-compatible marker catalog and gene utilities adapter.

Delegates to core.gene_utils to provide detect_genes_in_columns,
get_default_markers_for_gene, find_best_column_match, and load_marker_catalog.
"""

from typing import List, Optional
import pandas as pd
from core.gene_utils import (
    KNOWN_PGX_GENES,
    detect_genes_in_columns,
    find_best_column_match,
    get_column_gene,
    get_default_markers_for_gene,
    load_marker_catalog_file,
    smart_merge_datasets,
)


def load_marker_catalog() -> pd.DataFrame:
    """Loads and caches the pharmacogenetic marker catalog."""
    return load_marker_catalog_file()
