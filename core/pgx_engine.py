"""Generic, data-driven pharmacogenomic diplotype and phenotype engine.

The engine contains no gene-specific calling logic. Marker definitions and
translation tables are supplied as data, so the same implementation can
process a gene with one marker or many markers.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Tuple

import pandas as pd


class RuleBasedPGxEngine:
    """Apply optional gene rules to sample rows and summarize the calls.

    A rule bundle has this shape::

        {
            "markers": {
                "GENE_marker1": {"genotype_dosage": {"AA": 0, "AG": 1, "GG": 2}},
                ...
            },
            "diplotype_rules": {"0,0": "*1/*1", "1,0": "*1/*2", ...},
            "phenotype_map": {"*1/*1": "Normal Metabolizer", ...},
            "phenotype_order": ["Normal Metabolizer", ...],
        }

    The marker order in ``markers`` defines the order of dosage values in
    ``diplotype_rules``. If any required call is absent or invalid, the sample
    is indeterminate. Without a rule bundle, marker statistics can still be
    calculated elsewhere in the pipeline, while clinical translation is
    explicitly marked as not configured.
    """

    @staticmethod
    def _normalize_call(value: Any) -> Optional[str]:
        if value is None or pd.isna(value):
            return None
        call = str(value).strip().upper()
        if call in {"", "NAN", "NONE", "NULL", "NA", "?", "MISSING", "UNKNOWN", "NOT PERFORMED"}:
            return None
        return call.replace("/", "").replace("|", "").replace(" ", "")

    @classmethod
    def assign_sample_diplotype(
        cls,
        row: pd.Series,
        gene_rules: Optional[Mapping[str, Any]] = None,
    ) -> Tuple[str, str]:
        """Return a rule-defined diplotype and phenotype for one sample."""
        if not gene_rules:
            return "Not configured", "Not configured"

        markers = gene_rules.get("markers", {})
        if not markers:
            return "Not configured", "Not configured"

        dosage_values = []
        for marker, marker_rules in markers.items():
            call = cls._normalize_call(row.get(marker))
            if call is None:
                return "Indeterminate", "Indeterminate"

            dosage_map = marker_rules.get("genotype_dosage", {})
            dosage = dosage_map.get(call)
            if dosage is None and len(call) == 2:
                dosage = dosage_map.get(call[::-1])
            if dosage is None:
                return "Indeterminate", "Indeterminate"
            dosage_values.append(int(dosage))

        pattern = ",".join(str(value) for value in dosage_values)
        diplotype = gene_rules.get("diplotype_rules", {}).get(pattern)
        if diplotype is None:
            return "Indeterminate", "Indeterminate"

        phenotype = gene_rules.get("phenotype_map", {}).get(diplotype, "Indeterminate")
        return str(diplotype), str(phenotype)

    @classmethod
    def process_dataframe(
        cls,
        df: pd.DataFrame,
        gene_rules: Optional[Mapping[str, Any]] = None,
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Add translation columns and return counts, frequencies, and percents."""
        result = df.copy()
        calls = [cls.assign_sample_diplotype(row, gene_rules) for _, row in result.iterrows()]
        result["Diplotype"] = [call[0] for call in calls]
        result["Phenotype"] = [call[1] for call in calls]

        total = len(result)
        configured = bool(gene_rules and gene_rules.get("markers") and gene_rules.get("diplotype_rules"))
        valid_mask = ~result["Phenotype"].isin(["Indeterminate", "Not configured"])
        valid_total = int(valid_mask.sum())

        diplotype_counts = result["Diplotype"].value_counts().to_dict()
        diplotype_summary = {
            name: {
                "count": int(count),
                "frequency": round(count / total, 4) if total else 0.0,
                "percentage": round(count / total * 100, 2) if total else 0.0,
            }
            for name, count in diplotype_counts.items()
        }

        phenotype_counts = result["Phenotype"].value_counts().to_dict()
        order = list(gene_rules.get("phenotype_order", [])) if gene_rules else []
        for name in phenotype_counts:
            if name not in order:
                order.append(name)
        phenotype_summary = {
            name: {
                "count": int(phenotype_counts.get(name, 0)),
                "frequency": round(phenotype_counts.get(name, 0) / total, 4) if total else 0.0,
                "percentage": round(phenotype_counts.get(name, 0) / total * 100, 2) if total else 0.0,
            }
            for name in order
        }

        summary = {
            "total_samples": total,
            "valid_pharmacogenomic_samples": valid_total,
            "diplotypes": diplotype_summary,
            "phenotypes": phenotype_summary,
            "translation_configured": configured,
        }
        return result, summary
