"""PCA projections for genotype marker matrices."""

from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

import numpy as np
import pandas as pd


def _clean_genotype(value: Any) -> Optional[str]:
    if value is None or pd.isna(value):
        return None
    call = str(value).strip().upper()
    if call in {"", "NAN", "NONE", "NULL", "NA", "?", "MISSING", "UNKNOWN", "NOT PERFORMED"}:
        return None
    return call.replace("/", "").replace("|", "").replace(":", "").replace(" ", "")


def genotype_pca(
    dataframe: pd.DataFrame,
    markers: Iterable[str],
    marker_rules: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Return PCA scores from biallelic diploid genotype calls.

    Genotypes are converted to variant-allele dosages (0, 1, or 2). Alleles
    come from marker rules when available; otherwise they are inferred from
    the observed calls. Missing calls are mean-imputed after dosage encoding,
    then markers are centered and scaled by expected HWE standard deviation.
    """
    marker_rules = marker_rules or {}
    available = [marker for marker in markers if marker in dataframe.columns]
    if dataframe.empty or len(dataframe) < 2:
        return pd.DataFrame(index=dataframe.index), {"reason": "PCA needs at least two samples."}
    if len(available) < 2:
        return pd.DataFrame(index=dataframe.index), {"reason": "PCA needs at least two selected markers."}

    encoded: Dict[str, pd.Series] = {}
    skipped = []
    for marker in available:
        calls = dataframe[marker].map(_clean_genotype)
        rule = marker_rules.get(marker, {}) or {}
        ref = str(rule.get("ref_allele", "")).strip().upper()
        var = str(rule.get("var_allele", "")).strip().upper()

        if not ref or not var:
            observed_alleles = sorted({base for call in calls.dropna() if len(call) == 2 for base in call})
            if len(observed_alleles) != 2:
                skipped.append(marker)
                continue
            ref, var = observed_alleles
        if ref == var:
            skipped.append(marker)
            continue

        possible_dosages: Dict[str, set[int]] = {}
        for genotype, value in (
            (f"{ref}{ref}", 0),
            (f"{ref}{var}", 1),
            (f"{var}{ref}", 1),
            (f"{var}{var}", 2),
        ):
            possible_dosages.setdefault(genotype, set()).add(value)

        dosage = []
        for call in calls:
            choices = possible_dosages.get(call or "", set())
            if len(choices) != 1:
                dosage.append(np.nan)
            else:
                dosage.append(float(next(iter(choices))))
        series = pd.Series(dosage, index=dataframe.index, dtype=float)
        if series.notna().sum() < 2 or series.dropna().nunique() < 2:
            skipped.append(marker)
            continue
        encoded[marker] = series

    if len(encoded) < 2:
        return pd.DataFrame(index=dataframe.index), {
            "reason": "At least two variable biallelic markers with usable calls are required.",
            "skipped_markers": skipped,
        }

    matrix = pd.DataFrame(encoded, index=dataframe.index)
    matrix = matrix.apply(lambda col: col.fillna(col.mean()))
    matrix = matrix.loc[:, matrix.nunique(dropna=False) > 1]
    if matrix.shape[1] < 2:
        return pd.DataFrame(index=dataframe.index), {
            "reason": "Fewer than two variable markers remain after quality checks.",
            "skipped_markers": skipped,
        }

    centered = matrix.to_numpy(dtype=float)
    centered -= centered.mean(axis=0, keepdims=True)
    allele_freq = matrix.mean(axis=0).to_numpy(dtype=float) / 2.0
    expected_sd = np.sqrt(2.0 * allele_freq * (1.0 - allele_freq))
    valid_scale = expected_sd > np.finfo(float).eps
    centered = centered[:, valid_scale] / expected_sd[valid_scale]
    used_markers = matrix.columns[valid_scale].tolist()
    if centered.shape[1] < 2:
        return pd.DataFrame(index=dataframe.index), {
            "reason": "Fewer than two informative markers remain after scaling.",
            "skipped_markers": skipped,
        }

    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    eigenvalues = singular_values ** 2
    total_variance = float(eigenvalues.sum())
    if total_variance <= np.finfo(float).eps:
        return pd.DataFrame(index=dataframe.index), {
            "reason": "The selected genotypes do not contain variation for PCA.",
            "skipped_markers": skipped,
        }

    component_count = min(3, int(np.count_nonzero(eigenvalues > np.finfo(float).eps)))
    scores = u[:, :component_count] * singular_values[:component_count]
    columns = [f"PC{i + 1}" for i in range(component_count)]
    score_df = pd.DataFrame(scores, index=dataframe.index, columns=columns)
    variance_pct = (eigenvalues[:component_count] / total_variance * 100.0).tolist()
    return score_df, {
        "markers_used": used_markers,
        "skipped_markers": skipped,
        "explained_variance_percent": variance_pct,
    }
