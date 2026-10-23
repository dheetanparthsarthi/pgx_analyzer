"""
Session State & Disk Cache Persistence Manager.
Handles saving and restoring session state, active datasets, and multi-dataset profiles across browser refreshes.
"""

import os
import pickle
import pandas as pd
import streamlit as st

DELETED_FLAG_FILE = ".data_deleted.flag"
ACTIVE_DATASET_CACHE = ".active_dataset.pkl"
ACTIVE_DATASETS_CACHE = ".active_datasets_dict.pkl"
ACTIVE_SESSION_CACHE = ".active_session_state.pkl"

def init_session_defaults():
    """Initializes default session state keys if not already present."""
    if os.path.exists(DELETED_FLAG_FILE):
        st.session_state['is_deleted'] = True

    defaults = {
        'is_deleted': False,
        'show_delete_confirm': False,
        'uploaded_df': None,
        'mapped_cols': {},
        'demographic_mapping': {},
        'datasets_dict': {},
        'selected_gene': None,
        'selected_markers': [],
        'workflow_step': '1 · Select Target',
        'custom_genes': [],
        'custom_markers': {},
        'extra_mapping_rows': [],
        'mapping_view_mode': '🗂️ Interactive Marker Cards (with Live QC)',
        'custom_snp_rules': {},
        'custom_translation_rules': {},
        'mapped_marker_details': {}
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

def save_active_session():
    """Saves workflow state, gene/marker selections, and dataset dictionary to disk so they survive page refreshes."""
    if st.session_state.get('is_deleted', False):
        return
    up_df = st.session_state.get('uploaded_df')
    has_df = up_df is not None and len(up_df) > 0
    state_to_save = {
        'workflow_step': st.session_state.get('workflow_step', '4 · Results' if has_df else '1 · Select Target'),
        'selected_gene': st.session_state.get('selected_gene'),
        'selected_markers': st.session_state.get('selected_markers', []),
        'mapped_cols': st.session_state.get('mapped_cols', {}),
        'demographic_mapping': st.session_state.get('demographic_mapping', {}),
        'custom_snp_rules': st.session_state.get('custom_snp_rules', {}),
        'custom_translation_rules': st.session_state.get('custom_translation_rules', {}),
        'mapped_marker_details': st.session_state.get('mapped_marker_details', {}),
        'custom_genes': st.session_state.get('custom_genes', []),
        'custom_markers': st.session_state.get('custom_markers', {}),
        'active_ds_name': st.session_state.get('active_ds_name'),
        'active_sheet_name': st.session_state.get('active_sheet_name'),
        'extra_mapping_rows': st.session_state.get('extra_mapping_rows', []),
        'mapping_view_mode': st.session_state.get('mapping_view_mode', '🗂️ Interactive Marker Cards (with Live QC)')
    }
    try:
        with open(ACTIVE_SESSION_CACHE, 'wb') as f:
            pickle.dump(state_to_save, f)
    except Exception:
        pass

    # Save all datasets in datasets_dict so multi-dataset profiles survive refresh
    d_dict = st.session_state.get('datasets_dict', {})
    if d_dict and any(isinstance(df, pd.DataFrame) and not df.empty for df in d_dict.values()):
        try:
            with open(ACTIVE_DATASETS_CACHE, 'wb') as f:
                pickle.dump(d_dict, f)
        except Exception:
            pass

def load_persisted_session():
    """Restores saved workflow step, gene, and mappings on page reload/refresh."""
    if os.path.exists(ACTIVE_SESSION_CACHE) and not st.session_state.get('is_deleted', False):
        try:
            with open(ACTIVE_SESSION_CACHE, 'rb') as f:
                saved = pickle.load(f)
            if isinstance(saved, dict):
                for k, v in saved.items():
                    if v is not None:
                        if k not in st.session_state or st.session_state[k] is None or st.session_state[k] == {} or st.session_state[k] == []:
                            st.session_state[k] = v
                        elif k == 'workflow_step' and st.session_state.get('workflow_step') == '1 · Select Target':
                            st.session_state['workflow_step'] = v
        except Exception:
            pass

def save_active_dataset(df: pd.DataFrame):
    """Saves active dataset and all datasets in dictionary to disk so they persist across page reloads and browser refreshes."""
    if df is not None and len(df) > 0:
        try:
            df.to_pickle(ACTIVE_DATASET_CACHE)
        except Exception:
            pass
    d_dict = st.session_state.setdefault('datasets_dict', {})
    active_name = st.session_state.get('active_ds_name', 'Active Dataset')
    if df is not None and len(df) > 0:
        d_dict[active_name] = df
    save_active_session()

def load_persisted_dataset() -> pd.DataFrame:
    """Loads dataset from persistent disk cache file if present."""
    if os.path.exists(ACTIVE_DATASET_CACHE):
        try:
            return pd.read_pickle(ACTIVE_DATASET_CACHE)
        except Exception:
            pass
    return None

def load_persisted_datasets_dict() -> dict:
    """Loads multiple datasets dictionary from disk cache file if present."""
    if os.path.exists(ACTIVE_DATASETS_CACHE):
        try:
            with open(ACTIVE_DATASETS_CACHE, 'rb') as f:
                d = pickle.load(f)
            if isinstance(d, dict):
                return d
        except Exception:
            pass
    return {}

def clear_active_dataset_cache():
    """Removes persistent dataset and session cache files."""
    for fn in [ACTIVE_DATASET_CACHE, ACTIVE_DATASETS_CACHE, ACTIVE_SESSION_CACHE]:
        if os.path.exists(fn):
            try:
                os.remove(fn)
            except Exception:
                pass

def set_active_target_gene(gene_symbol: str, markers: list = None, catalog_df: pd.DataFrame = None):
    """Sets active target gene and its markers across session state and synchronizes all widget keys."""
    from core.gene_utils import get_default_markers_for_gene, find_best_column_match
    gene_u = str(gene_symbol).strip().upper() if gene_symbol else None
    if markers is None and gene_u:
        def_m = get_default_markers_for_gene(gene_u, catalog_df)
        up_df = st.session_state.get('uploaded_df')
        if up_df is not None and not up_df.empty:
            matched_m = [m for m in def_m if find_best_column_match(m, m.split(' ')[0], gene_u, list(up_df.columns))]
            markers = matched_m if matched_m else def_m
        else:
            markers = def_m
    elif markers is None:
        markers = []

    st.session_state['selected_gene'] = gene_u
    st.session_state['selected_markers'] = list(markers)
    st.session_state['extra_mapping_rows'] = []

    for k in ['step1_sel_gene', 'step1_sel_markers', 'map_step_gene_sel', 'map_step_marker_sel']:
        if k in st.session_state:
            del st.session_state[k]
    for k in list(st.session_state.keys()):
        if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', 'ex_', '_mapped_', '_col_tracker_')):
            del st.session_state[k]
    save_active_session()


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
    st.session_state['uploaded_df'] = pd.DataFrame()
    st.session_state['selected_gene'] = None
    st.session_state['selected_markers'] = []
    st.session_state['mapped_cols'] = {}
    st.session_state['demographic_mapping'] = {}
    st.session_state['custom_translation_rules'] = {}
    st.session_state['datasets_dict'] = {}
    st.session_state['last_uploaded_file_id'] = None
    st.session_state['pending_new_df'] = None
    st.session_state['pending_file_id'] = None
    st.session_state['pending_file_name'] = None
    st.session_state['merge_alignment_log'] = None
    st.session_state['show_delete_confirm'] = False
    st.session_state['workflow_step'] = '1 · Select Target'
    for k in list(st.session_state.keys()):
        if k.startswith(('w_map_', 'sel_ref_', 'sel_var_', 'custom_ref_', 'custom_var_', 'ex_', 'step1_sel_', 'map_step_')):
            del st.session_state[k]

    clear_active_dataset_cache()
    try:
        with open(DELETED_FLAG_FILE, "w") as f:
            f.write("deleted")
    except Exception:
        pass

    try:
        from database.db_manager import DatabaseManager
        DatabaseManager().clear_database()
    except Exception:
        pass

    try:
        st.cache_data.clear()
    except Exception:
        pass


def init_session_state():
    """Initializes session defaults and restores persisted cache."""
    init_session_defaults()
    if not st.session_state.get('is_deleted', False):
        persisted_dict = load_persisted_datasets_dict()
        if persisted_dict:
            st.session_state['datasets_dict'] = persisted_dict

        load_persisted_session()

        cur_ds_name = st.session_state.get('active_ds_name')
        if cur_ds_name and cur_ds_name in st.session_state.get('datasets_dict', {}):
            st.session_state['uploaded_df'] = st.session_state['datasets_dict'][cur_ds_name]
        elif st.session_state.get('uploaded_df') is None or len(st.session_state['uploaded_df']) == 0:
            persisted_df = load_persisted_dataset()
            if persisted_df is not None and len(persisted_df) > 0:
                st.session_state['uploaded_df'] = persisted_df
                if not st.session_state.get('datasets_dict'):
                    st.session_state['datasets_dict'] = {'Active Dataset': persisted_df}
                    st.session_state['active_ds_name'] = 'Active Dataset'

    if st.session_state.get('uploaded_df') is None:
        st.session_state['uploaded_df'] = pd.DataFrame()


def has_active_data() -> bool:
    """Checks whether the active workspace has loaded dataset records."""
    df = st.session_state.get('uploaded_df')
    if df is not None and not df.empty and len(df) > 0:
        return True
    return False


