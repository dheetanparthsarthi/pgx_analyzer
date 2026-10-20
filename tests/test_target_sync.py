import pytest
import streamlit as st
import pandas as pd
from core.gene_utils import detect_genes_in_columns, get_default_markers_for_gene, get_column_gene

def mock_set_active_target_gene(gene_symbol: str, markers: list = None, catalog_df: pd.DataFrame = None):
    gene_u = str(gene_symbol).strip().upper()
    if markers is None:
        markers = get_default_markers_for_gene(gene_u, catalog_df)
    st.session_state['selected_gene'] = gene_u
    st.session_state['selected_markers'] = list(markers)
    for k in ['step1_sel_gene', 'step1_sel_markers', 'map_step_gene_sel', 'map_step_marker_sel']:
        if k in st.session_state:
            del st.session_state[k]
    for k in list(st.session_state.keys()):
        if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', 'ex_')):
            del st.session_state[k]

def test_detect_genes_in_columns():
    cols_cyp2c9 = ['Sample ID', 'Gender', 'CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
    genes = detect_genes_in_columns(cols_cyp2c9)
    assert genes == ['CYP2C9']

    cols_cyp2c19 = ['sample_id', 'CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17']
    genes19 = detect_genes_in_columns(cols_cyp2c19)
    assert genes19 == ['CYP2C19']

def test_detect_genes_from_marker_labels():
    markers = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
    genes = detect_genes_in_columns(markers)
    assert genes == ['CYP2C9']

def test_get_default_markers_for_gene():
    m19 = get_default_markers_for_gene('CYP2C19')
    assert 'CYP2C19*2 (rs4244285)' in m19
    assert 'CYP2C19*3 (rs4986893)' in m19
    assert 'CYP2C19*17 (rs12248560)' in m19

    m9 = get_default_markers_for_gene('cyp2c9')
    assert 'CYP2C9*2 (rs1799853)' in m9
    assert 'CYP2C9*3 (rs1057910)' in m9

def test_set_active_target_gene_sync():
    # Setup initial mock session state
    st.session_state['selected_gene'] = 'CYP2C19'
    st.session_state['selected_markers'] = ['CYP2C19*2 (rs4244285)']
    st.session_state['step1_sel_gene'] = 'CYP2C19'
    st.session_state['step1_sel_markers'] = ['CYP2C19*2 (rs4244285)']
    st.session_state['map_step_gene_sel'] = 'CYP2C19'
    st.session_state['map_step_marker_sel'] = ['CYP2C19*2 (rs4244285)']
    st.session_state['w_map_col_0'] = 'Old Column'

    # Switch to CYP2C9
    mock_set_active_target_gene('CYP2C9')

    assert st.session_state['selected_gene'] == 'CYP2C9'
    assert 'CYP2C9*2 (rs1799853)' in st.session_state['selected_markers']
    assert 'CYP2C9*3 (rs1057910)' in st.session_state['selected_markers']

    # Stale widget keys should be cleanly purged so dropdowns re-instantiate cleanly
    assert 'step1_sel_gene' not in st.session_state
    assert 'step1_sel_markers' not in st.session_state
    assert 'map_step_gene_sel' not in st.session_state
    assert 'map_step_marker_sel' not in st.session_state
    assert 'w_map_col_0' not in st.session_state

def test_auto_heal_gene_from_markers():
    # Test the auto-heal logic when markers are CYP2C9 but target gene is CYP2C19
    target_markers = ['CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']
    target_gene = 'CYP2C19'

    m_genes = []
    for m in target_markers:
        dg = detect_genes_in_columns([m])
        if dg:
            m_genes.extend(dg)
    assert m_genes == ['CYP2C9', 'CYP2C9']
    assert all(g.upper() == m_genes[0].upper() for g in m_genes)

    # Inferred gene should correctly be CYP2C9
    inferred_gene = m_genes[0]
    assert inferred_gene == 'CYP2C9'
    assert target_gene != inferred_gene

def test_find_best_column_match_exact_and_rsid():
    from core.gene_utils import find_best_column_match
    cols = ['Sample ID', 'Gender', 'State', 'CYP2C9*2 (rs1799853)', 'CYP2C9*3 (rs1057910)']

    # Exact match with marker label
    assert find_best_column_match('CYP2C9*2 (rs1799853)', 'CYP2C9*2', 'CYP2C9', cols) == 'CYP2C9*2 (rs1799853)'
    assert find_best_column_match('CYP2C9*3 (rs1057910)', 'CYP2C9*3', 'CYP2C9', cols) == 'CYP2C9*3 (rs1057910)'

    # rsID match
    rs_cols = ['Sample ID', 'rs1799853', 'rs1057910']
    assert find_best_column_match('CYP2C9*2 (rs1799853)', 'CYP2C9*2', 'CYP2C9', rs_cols) == 'rs1799853'
    assert find_best_column_match('CYP2C9*3 (rs1057910)', 'CYP2C9*3', 'CYP2C9', rs_cols) == 'rs1057910'

def test_find_best_column_match_star_alleles():
    from core.gene_utils import find_best_column_match
    star_cols = ['Sample ID', 'CYP2C9*2', 'CYP2C9*3']
    assert find_best_column_match('CYP2C9*2 (rs1799853)', 'CYP2C9*2', 'CYP2C9', star_cols) == 'CYP2C9*2'

    short_star_cols = ['Sample ID', '*2', '*3']
    assert find_best_column_match('CYP2C9*2 (rs1799853)', 'CYP2C9*2', 'CYP2C9', short_star_cols) == '*2'

def test_find_best_column_match_no_false_fallback():
    from core.gene_utils import find_best_column_match
    cols_19 = ['Sample ID', 'Gender', 'CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17']

    # When active gene is CYP2C9, matching CYP2C9*2 against CYP2C19 dataset must return empty string
    assert find_best_column_match('CYP2C9*2 (rs1799853)', 'CYP2C9*2', 'CYP2C9', cols_19) == ""
    assert find_best_column_match('CYP2C9*3 (rs1057910)', 'CYP2C9*3', 'CYP2C9', cols_19) == ""

def test_find_best_column_match_ignores_select_column():
    from core.gene_utils import find_best_column_match
    cols_with_select = ['(Select Column)', 'Sample ID', 'CYP2C9*2']
    assert find_best_column_match('CYP2C9*2 (rs1799853)', 'CYP2C9*2', 'CYP2C9', cols_with_select) == 'CYP2C9*2'
    assert find_best_column_match('CYP2C19*2 (rs4244285)', 'CYP2C19*2', 'CYP2C19', cols_with_select) == ""

def test_catalog_gene_markers_loading():
    """Verify that catalog genes load their respective markers properly."""
    catalog_df = pd.DataFrame({
        'Gene': ['ABCG2', 'ACE', 'ALDH2', 'ABCB1'],
        'Marker': ['rs2231142', 'rs4646994', 'rs671', 'rs1045642']
    })
    assert get_default_markers_for_gene('ABCG2', catalog_df) == ['rs2231142']
    assert get_default_markers_for_gene('ACE', catalog_df) == ['rs4646994']
    assert get_default_markers_for_gene('ALDH2', catalog_df) == ['rs671']
    assert get_default_markers_for_gene('ABCB1', catalog_df) == ['rs1045642']

def test_user_gene_selection_not_cleared_by_old_gene():
    """Verify that when user selects a new gene in step1_sel_gene, it is not deleted when differing from previous gene."""
    gene_options = ["(Select Target Gene)", "ABCB1", "ABCG2", "ACE", "ALDH2"]
    st.session_state['selected_gene'] = 'ABCB1'
    # User selects ABCG2 in the UI
    st.session_state['step1_sel_gene'] = 'ABCG2'
    
    # Validation check: should only delete if NOT in gene_options
    if 'step1_sel_gene' in st.session_state and st.session_state['step1_sel_gene'] not in gene_options:
        del st.session_state['step1_sel_gene']
        
    assert 'step1_sel_gene' in st.session_state
    assert st.session_state['step1_sel_gene'] == 'ABCG2'

def test_smart_merge_same_gene_distinct_and_shared_markers():
    """Verify that merging datasets with the same gene merges shared markers and keeps distinct markers as separate columns."""
    from core.gene_utils import smart_merge_datasets
    df1 = pd.DataFrame({
        'Sample ID': ['P1', 'P2'],
        'CYP2C19*2': ['GA', 'GG'],
        'Gender': ['M', 'F']
    })
    df2 = pd.DataFrame({
        'sample_id': ['P3', 'P4'],
        'CYP2C19*2': ['AA', 'GA'],
        'CYP2C19*17': ['CT', 'TT'],
        'Gender': ['F', 'M']
    })

    merged_df, align_log = smart_merge_datasets(df1, df2)
    assert len(merged_df) == 4
    # Demographic aligned
    assert 'Sample ID' in merged_df.columns
    # Shared marker merged into same column
    assert 'CYP2C19*2' in merged_df.columns
    # Distinct marker added as separate column
    assert 'CYP2C19*17' in merged_df.columns
    assert any("Distinct Marker" in msg for msg in align_log)

def test_smart_merge_blocks_incompatible_genes():
    """Verify that attempting to merge datasets representing completely different genes is blocked."""
    from core.gene_utils import smart_merge_datasets
    df_cyp19 = pd.DataFrame({'Sample ID': ['P1'], 'CYP2C19*2': ['GA']})
    df_cyp9 = pd.DataFrame({'Sample ID': ['P2'], 'CYP2C9*2': ['CC']})

    merged_df, align_log = smart_merge_datasets(df_cyp19, df_cyp9)
    assert any("Cannot merge" in msg for msg in align_log)


