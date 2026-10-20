"""Tests for the gene-agnostic, rule-driven translation engine."""

import pandas as pd

from core.pgx_engine import RuleBasedPGxEngine


def test_rule_engine_supports_one_and_four_markers():
    one_marker_rules = {
        'markers': {
            'CUSTOM_rs1': {'genotype_dosage': {'AA': 0, 'AG': 1, 'GA': 1, 'GG': 2}},
        },
        'diplotype_rules': {'0': 'Ref/Ref', '1': 'Ref/Alt', '2': 'Alt/Alt'},
        'phenotype_map': {'Ref/Ref': 'Class 0', 'Ref/Alt': 'Class 1', 'Alt/Alt': 'Class 2'},
        'phenotype_order': ['Class 0', 'Class 1', 'Class 2'],
    }
    assert RuleBasedPGxEngine.assign_sample_diplotype(
        pd.Series({'CUSTOM_rs1': 'GA'}), one_marker_rules
    ) == ('Ref/Alt', 'Class 1')

    four_marker_rules = {
        'markers': {
            f'CUSTOM_rs{i}': {'genotype_dosage': {'CC': 0, 'CT': 1, 'TC': 1, 'TT': 2}}
            for i in range(1, 5)
        },
        'diplotype_rules': {'0,0,0,0': 'Reference Pattern', '1,0,0,0': 'Pattern A'},
        'phenotype_map': {'Reference Pattern': 'Reference Class', 'Pattern A': 'Alternate Class'},
    }
    sample = pd.Series({'CUSTOM_rs1': 'CT', 'CUSTOM_rs2': 'CC', 'CUSTOM_rs3': 'CC', 'CUSTOM_rs4': 'CC'})
    assert RuleBasedPGxEngine.assign_sample_diplotype(sample, four_marker_rules) == ('Pattern A', 'Alternate Class')


def test_unknown_gene_rules_are_not_guessed():
    called, summary = RuleBasedPGxEngine.process_dataframe(pd.DataFrame({'UNKNOWN_rs1': ['AA']}))
    assert called.loc[0, 'Diplotype'] == 'Not configured'
    assert called.loc[0, 'Phenotype'] == 'Not configured'
    assert summary['translation_configured'] is False
