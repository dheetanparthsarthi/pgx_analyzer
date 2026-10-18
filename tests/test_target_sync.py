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
