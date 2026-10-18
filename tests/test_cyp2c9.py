"""
Unit tests for CYP2C9 Translator (core.cyp2c9).
"""

import pytest
import pandas as pd
from core.cyp2c9 import CYP2C9Translator


def test_cyp2c9_diplotype_assignments():
    # Normal Metabolizer (*1/*1)
    row_nm = pd.Series({'CYP2C9*2': 'CC', 'CYP2C9*3': 'AA'})
    dip, pheno = CYP2C9Translator.assign_sample_diplotype(row_nm)
    assert dip == '*1/*1'
    assert pheno == 'Normal Metabolizer'

    # Intermediate Metabolizer (*1/*2)
    row_im1 = pd.Series({'CYP2C9*2': 'CT', 'CYP2C9*3': 'AA'})
    dip, pheno = CYP2C9Translator.assign_sample_diplotype(row_im1)
    assert dip == '*1/*2'
    assert pheno == 'Intermediate Metabolizer'

    # Intermediate Metabolizer (*1/*3)
    row_im2 = pd.Series({'CYP2C9*2': 'CC', 'CYP2C9*3': 'AC'})
    dip, pheno = CYP2C9Translator.assign_sample_diplotype(row_im2)
    assert dip == '*1/*3'
    assert pheno == 'Intermediate Metabolizer'

    # Intermediate Metabolizer (*2/*2)
    row_im3 = pd.Series({'CYP2C9*2': 'TT', 'CYP2C9*3': 'AA'})
    dip, pheno = CYP2C9Translator.assign_sample_diplotype(row_im3)
    assert dip == '*2/*2'
    assert pheno == 'Intermediate Metabolizer'

    # Poor Metabolizer (*2/*3)
    row_pm1 = pd.Series({'CYP2C9*2': 'CT', 'CYP2C9*3': 'AC'})
    dip, pheno = CYP2C9Translator.assign_sample_diplotype(row_pm1)
    assert dip == '*2/*3'
    assert pheno == 'Poor Metabolizer'

    # Poor Metabolizer (*3/*3)
    row_pm2 = pd.Series({'CYP2C9*2': 'CC', 'CYP2C9*3': 'CC'})
    dip, pheno = CYP2C9Translator.assign_sample_diplotype(row_pm2)
    assert dip == '*3/*3'
    assert pheno == 'Poor Metabolizer'
