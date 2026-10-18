"""
Unit tests for Data Quality Control Engine (core.qc).
"""

import pytest
import pandas as pd
from core.qc import DataQCEngine


def test_data_qc_engine():
    data = {
        'sample ID': ['OD001', 'OD002', 'OD001', 'OD003'],
        'Gender': ['Female ', 'male', 'Female', 'Not Mentioned'],
        'State': ['Andhra Pradesh', 'Punjab', 'Andhra Pradesh', 'Maharashtra'],
        'CYP2C19*2 (rs4244285)': ['GA', 'GG', 'AA', 'INVALID']
    }
    df = pd.DataFrame(data)
    qc_engine = DataQCEngine(df)
    clean_df, report = qc_engine.run_qc()
    
    assert report['total_samples'] == 4
    assert 'OD001' in report['duplicate_sample_ids']
    assert clean_df.loc[clean_df['Sample_ID'] == 'OD001', 'Gender_Clean'].iloc[0] == 'Female'
    assert clean_df.loc[clean_df['Sample_ID'] == 'OD002', 'Gender_Clean'].iloc[0] == 'Male'
    assert clean_df.loc[clean_df['Sample_ID'] == 'OD003', 'Gender_Clean'].iloc[0] == 'Unknown'
    # Check invalid genotype converted to NaN
    assert pd.isna(clean_df.loc[clean_df['Sample_ID'] == 'OD003', 'CYP2C19*2'].iloc[0])


def test_ag_and_tc_are_invalid():
    data = {
        'Sample_ID': ['S1', 'S2', 'S3', 'S4'],
        'CYP2C19*2': ['GA', 'AG', 'GG', 'AA'],
        'CYP2C19*17': ['CT', 'TC', 'CC', 'TT']
    }
    df = pd.DataFrame(data)
    qc_engine = DataQCEngine(df)
    clean_df, report = qc_engine.run_qc()

    # AG is invalid for *2; only GA, GG, AA are valid
    assert clean_df.loc[clean_df['Sample_ID'] == 'S1', 'CYP2C19*2'].iloc[0] == 'GA'
    assert pd.isna(clean_df.loc[clean_df['Sample_ID'] == 'S2', 'CYP2C19*2'].iloc[0])
    assert clean_df.loc[clean_df['Sample_ID'] == 'S3', 'CYP2C19*2'].iloc[0] == 'GG'
    assert clean_df.loc[clean_df['Sample_ID'] == 'S4', 'CYP2C19*2'].iloc[0] == 'AA'

    # TC is invalid for *17; only CT, CC, TT are valid
    assert clean_df.loc[clean_df['Sample_ID'] == 'S1', 'CYP2C19*17'].iloc[0] == 'CT'
    assert pd.isna(clean_df.loc[clean_df['Sample_ID'] == 'S2', 'CYP2C19*17'].iloc[0])
    assert clean_df.loc[clean_df['Sample_ID'] == 'S3', 'CYP2C19*17'].iloc[0] == 'CC'
    assert clean_df.loc[clean_df['Sample_ID'] == 'S4', 'CYP2C19*17'].iloc[0] == 'TT'


def test_custom_permitted_genotypes_strict_flagging():
    """Verify that when a doctor permits only CC and CT, extra TT calls are strictly flagged as invalid."""
    data = {
        'Sample_ID': ['S1', 'S2', 'S3', 'S4'],
        'Raw_SNP_Col': ['CC', 'CT', 'TT', 'CC']
    }
    df = pd.DataFrame(data)
    custom_map = {'Raw_SNP_Col': 'Custom_Locus'}
    # Doctor allows only CC and CT; TT is unmapped/unexpected
    custom_rules = {
        'Custom_Locus': {
            'ref_allele': 'C',
            'var_allele': 'T',
            'valid_genotypes': ('CC', 'CT', 'TC')
        }
    }
    qc_engine = DataQCEngine(df, custom_col_mapping=custom_map, custom_snp_rules=custom_rules)
    clean_df, report = qc_engine.run_qc()

    # CC and CT must be retained
    assert clean_df.loc[clean_df['Sample_ID'] == 'S1', 'Custom_Locus'].iloc[0] == 'CC'
    assert clean_df.loc[clean_df['Sample_ID'] == 'S2', 'Custom_Locus'].iloc[0] == 'CT'
    assert clean_df.loc[clean_df['Sample_ID'] == 'S4', 'Custom_Locus'].iloc[0] == 'CC'

    # TT must be strictly converted to NaN and logged in invalid_genotypes_log
    assert pd.isna(clean_df.loc[clean_df['Sample_ID'] == 'S3', 'Custom_Locus'].iloc[0])
    assert 'Custom_Locus' in report['invalid_genotypes_log']
    assert report['invalid_genotypes_count'] == 1


def test_mixed_type_sample_ids_with_missing_and_duplicates():
    """Verify that sample IDs with mixed types (ints, floats, strings, NaNs) sort safely without TypeError."""
    data = {
        'Lab ID': [101, '102', 101, None, 'OD_005'],
        'CYP2C9*2': ['CC', 'CT', None, None, 'TT'],
        'CYP2C9*3': [None, 'AA', 'AA', 'AC', 'CC']
    }
    df = pd.DataFrame(data)
    qc_engine = DataQCEngine(df)
    clean_df, report = qc_engine.run_qc()

    assert report['total_samples'] == 5
    assert isinstance(report['samples_with_missing_ids'], list)
    assert isinstance(report['duplicate_sample_ids'], list)
    assert '101' in report['duplicate_sample_ids']
    # 101 has missing calls, None (auto-named SAMPLE_00004) has missing CYP2C9*2
    assert report['samples_with_missing_count'] == 2


