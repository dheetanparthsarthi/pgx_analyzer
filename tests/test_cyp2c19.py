"""
Regression tests for CYP2C19 rules using the shared rule-based engine.
"""

import pytest
import pandas as pd
from config import GENE_TRANSLATION_RULES
from core.pgx_engine import RuleBasedPGxEngine


def test_diplotype_assignments():
    # Normal Metabolizer (*1/*1)
    row_nm = pd.Series({'CYP2C19*2': 'GG', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC'})
    dip, pheno = RuleBasedPGxEngine.assign_sample_diplotype(row_nm, GENE_TRANSLATION_RULES['CYP2C19'])
    assert dip == '*1/*1'
    assert pheno == 'Normal Metabolizer'

    # Poor Metabolizer (*2/*2)
    row_pm = pd.Series({'CYP2C19*2': 'AA', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC'})
    dip, pheno = RuleBasedPGxEngine.assign_sample_diplotype(row_pm, GENE_TRANSLATION_RULES['CYP2C19'])
    assert dip == '*2/*2'
    assert pheno == 'Poor Metabolizer'

    # Ultrarapid Metabolizer (*17/*17)
    row_um = pd.Series({'CYP2C19*2': 'GG', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'TT'})
    dip, pheno = RuleBasedPGxEngine.assign_sample_diplotype(row_um, GENE_TRANSLATION_RULES['CYP2C19'])
    assert dip == '*17/*17'
    assert pheno == 'Ultrarapid Metabolizer'

    # Indeterminate (*2/*17 or *17/*2)
    row_im = pd.Series({'CYP2C19*2': 'GA', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CT'})
    dip, pheno = RuleBasedPGxEngine.assign_sample_diplotype(row_im, GENE_TRANSLATION_RULES['CYP2C19'])
    assert dip == '*2/*17'
    assert pheno == 'Indeterminate'

    # A required marker call is missing, so the engine must not assume reference.
    missing_row = pd.Series({'CYP2C19*2': 'GG', 'CYP2C19*3': 'GG', 'CYP2C19*17': None})
    assert RuleBasedPGxEngine.assign_sample_diplotype(missing_row, GENE_TRANSLATION_RULES['CYP2C19']) == ('Indeterminate', 'Indeterminate')
