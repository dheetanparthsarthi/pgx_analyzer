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
