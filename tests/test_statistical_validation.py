"""Unit tests for advanced statistical validation methods in core/stats.py."""

import numpy as np
import pandas as pd
import pytest
from core.stats import GeneticStatsEngine


def test_calculate_odds_ratio_haldane():
    # Regular 2x2 table
    # Group 1: 10 carriers, 90 non-carriers
    # Group 2: 5 carriers, 95 non-carriers
    res = GeneticStatsEngine.calculate_odds_ratio_haldane(10, 90, 5, 95)
    assert res['odds_ratio'] > 1.0
    assert res['ci_lower'] > 0
    assert res['ci_upper'] > res['ci_lower']
    assert 0.0 <= res['p_value'] <= 1.0
    assert not res['haldane_applied']

    # Zero-cell table (Haldane correction applied)
    res_zero = GeneticStatsEngine.calculate_odds_ratio_haldane(0, 10, 5, 50)
    assert res_zero['haldane_applied'] is True
    assert res_zero['odds_ratio'] > 0
    assert res_zero['ci_lower'] > 0
    assert res_zero['ci_upper'] > res_zero['ci_lower']
    assert 0.0 <= res_zero['p_value'] <= 1.0


def test_calculate_wilson_ci():
    # 50 variant alleles out of 1000 total alleles
    lower, upper = GeneticStatsEngine.calculate_wilson_ci(50, 1000, confidence=0.95)
    assert 0.0 < lower < 0.05
    assert 0.05 < upper < 0.10
    assert lower < upper

    # Edge cases
    assert GeneticStatsEngine.calculate_wilson_ci(0, 0) == (0.0, 0.0)


def test_calculate_fst():
    # 5 identical populations -> Fst must be 0
    identical_groups = {
        f"Pop_{i}": {"A": 90, "G": 10} for i in range(5)
    }
    res_ident = GeneticStatsEngine.calculate_fst(identical_groups)
    assert res_ident['fst'] == 0.0
    assert "Very low" in res_ident['interpretation'] or "Uniform" in res_ident['interpretation'] or "Low" in res_ident['interpretation']

    # Differentiated populations -> positive Fst
    diff_groups = {
        "Pop_1": {"A": 100, "G": 0},
        "Pop_2": {"A": 0, "G": 100},
    }
    res_diff = GeneticStatsEngine.calculate_fst(diff_groups)
    assert res_diff['fst'] == 1.0


def test_calculate_statistical_power():
    # Equal proportions -> power equals alpha
    pwr_equal = GeneticStatsEngine.calculate_statistical_power(50, 50, 0.20, 0.20, alpha=0.05)
    assert pwr_equal == 0.05

    # Large sample size and large difference -> high power
    pwr_high = GeneticStatsEngine.calculate_statistical_power(500, 500, 0.10, 0.30, alpha=0.05)
    assert pwr_high > 0.99

    # Small cohort (e.g. N=7 vs N=70 with modest diff) -> low power
    pwr_small = GeneticStatsEngine.calculate_statistical_power(7, 70, 0.286, 0.214, alpha=0.05)
    assert pwr_small < 0.20


def test_calculate_linkage_disequilibrium():
    # Create synthetic dataset with 2 markers
    np.random.seed(42)
    n = 200
    # Independent markers
    m1 = np.random.choice(["CC", "CT"], size=n, p=[0.9, 0.1])
    m2 = np.random.choice(["AA", "AC"], size=n, p=[0.85, 0.15])
    df = pd.DataFrame({"CYP2C9*2": m1, "CYP2C9*3": m2})

    cfg1 = {"wildtype_genotype": "CC", "het_genotype": "CT", "hom_var_genotype": "TT"}
    cfg2 = {"wildtype_genotype": "AA", "het_genotype": "AC", "hom_var_genotype": "CC"}

    ld_res = GeneticStatsEngine.calculate_linkage_disequilibrium(
        df, "CYP2C9*2", "CYP2C9*3", custom_info1=cfg1, custom_info2=cfg2
    )
    assert "r2" in ld_res
    assert 0.0 <= ld_res["r2"] <= 1.0
    assert ld_res["samples_analyzed"] == n
    assert "Independent segregation" in ld_res["interpretation"] or ld_res["r2"] < 0.1


def test_calculate_monte_carlo_chi2():
    # Contingency table with some sparse cells
    table = np.array([
        [7, 0],
        [55, 15],
        [36, 6],
        [8, 0],
        [293, 80]
    ])
    mc_res = GeneticStatsEngine.calculate_monte_carlo_chi2(table, n_sim=500, seed=42)
    assert mc_res["chi2_stat"] > 0
    assert mc_res["df"] == 4
    assert 0.0 <= mc_res["p_asymptotic"] <= 1.0
    assert 0.0 <= mc_res["p_monte_carlo"] <= 1.0
    assert isinstance(mc_res["cochran_violated"], (bool, np.bool_))
