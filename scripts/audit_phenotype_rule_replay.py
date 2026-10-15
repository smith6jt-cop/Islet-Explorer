#!/usr/bin/env python3
"""Replay the checked-in scimap rules on a bounded frozen-cell sample.

This audit reads the canonical H5AD without loading its full expression matrix.
It samples up to ``--per-stratum`` cells from each phenotype-by-donor stratum,
replays ``data/phenotype_rules.csv`` on the saved rescaled ``X``, and writes only
aggregate results to an untracked JSON report. The aliases are deliberately
favorable comparisons, not evidence that those names were historically mapped.

Run with the scimap environment, for example::

    conda run -n scimap python scripts/audit_phenotype_rule_replay.py
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = REPO / "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
DEFAULT_RULES = REPO / "data/phenotype_rules.csv"
DEFAULT_REPORT = REPO / "tmp/phenotype_rule_replay.json"

# These generous name matches make the comparison an upper bound for this CSV.
# In particular, Dendritic cell -> APCs is a comparison convenience, not a
# recovered scientific rule.
COMPARISON_ALIASES = {
    "Smooth Muscle": "SMA+",
    "CD8 T cell": "CD8a Tcell",
    "CD4 T cell": "CD4 Tcell",
    "Blood vessel": "Blood Vessel",
    "Dendritic cell": "APCs",
}
TCELL_LABELS = {"CD8a Tcell", "CD4 Tcell", "T cell"}


def _decode(value: object) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def select_stratified_indices(
    phenotype_codes: np.ndarray,
    donor_codes: np.ndarray,
    phenotype_names: list[str],
    donor_names: list[str],
    per_stratum: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, object]]]:
    """Return sorted sampled rows, population weights, and all nonempty strata."""
    if per_stratum < 1:
        raise ValueError("per_stratum must be positive")
    if phenotype_codes.shape != donor_codes.shape:
        raise ValueError("phenotype and donor codes have different shapes")
    if np.any(phenotype_codes < 0) or np.any(donor_codes < 0):
        raise ValueError("missing phenotype or donor codes cannot be stratified")

    rng = np.random.default_rng(seed)
    chosen: list[int] = []
    weights_by_row: dict[int, float] = {}
    strata: list[dict[str, object]] = []
    for phenotype_index, phenotype in enumerate(phenotype_names):
        phenotype_rows = np.flatnonzero(phenotype_codes == phenotype_index)
        for donor_index, donor in enumerate(donor_names):
            pool = phenotype_rows[donor_codes[phenotype_rows] == donor_index]
            if pool.size == 0:
                continue
            take = rng.choice(pool, size=min(per_stratum, pool.size), replace=False)
            weight = float(pool.size / take.size)
            for row in take:
                chosen.append(int(row))
                weights_by_row[int(row)] = weight
            strata.append(
                {
                    "phenotype": phenotype,
                    "donor": donor,
                    "population_n": int(pool.size),
                    "sample_n": int(take.size),
                }
            )

    indices = np.array(sorted(chosen), dtype=np.int64)
    weights = np.array([weights_by_row[int(row)] for row in indices], dtype=np.float64)
    if len(indices) != len(set(chosen)):
        raise AssertionError("stratified sample contains duplicate rows")
    return indices, weights, strata


def apply_comparison_aliases(labels: np.ndarray) -> np.ndarray:
    """Apply generous display-name matches without changing the source data."""
    return np.array([COMPARISON_ALIASES.get(str(label), str(label)) for label in labels], dtype=object)


def apply_tcell_correction(labels: np.ndarray, expression: np.ndarray, markers: list[str]) -> np.ndarray:
    """Replay the canonical H5AD's later CD8a/CD4 subtype correction."""
    if "CD8a" not in markers or "CD4" not in markers:
        raise ValueError("CD8a and CD4 are required for T-cell correction")
    result = np.array(labels, dtype=object, copy=True)
    lineage = np.isin(result, list(TCELL_LABELS))
    cd8_positive = expression[:, markers.index("CD8a")] >= 0.5
    cd4_positive = expression[:, markers.index("CD4")] >= 0.5
    result[lineage & cd8_positive & ~cd4_positive] = "CD8a Tcell"
    result[lineage & cd4_positive & ~cd8_positive] = "CD4 Tcell"
    result[lineage & ~(cd8_positive ^ cd4_positive)] = "T cell"
    return result


def agreement(truth: np.ndarray, prediction: np.ndarray, weights: np.ndarray) -> dict[str, float | int]:
    """Score sampled agreement and the stratum-weighted cohort estimate."""
    if not (len(truth) == len(prediction) == len(weights)) or len(truth) == 0:
        raise ValueError("truth, prediction, and weights must have equal nonzero lengths")
    correct = np.asarray(truth) == np.asarray(prediction)
    return {
        "sample_agreement": float(np.mean(correct)),
        "population_weighted_agreement": float(np.average(correct, weights=weights)),
        "matches": int(np.count_nonzero(correct)),
        "mismatches": int(np.count_nonzero(~correct)),
    }


def _counts(values: np.ndarray) -> dict[str, int]:
    return dict(sorted((str(label), int(n)) for label, n in Counter(values).items()))


def run_audit(h5ad_path: Path, rules_path: Path, per_stratum: int, seed: int) -> dict[str, object]:
    # Import heavy scientific packages only for the real replay, so the small
    # fixture tests do not need to import scimap or open multi-GB artifacts.
    import anndata as ad
    import h5py
    import pandas as pd
    import scimap as sm
    from importlib.metadata import version

    h5ad_path = h5ad_path.resolve(strict=True)
    rules_path = rules_path.resolve(strict=True)
    rules = pd.read_csv(rules_path)
    rule_markers = [str(name) for name in rules.columns[2:]]
    rule_labels = [str(name) for name in rules.iloc[:, 1].dropna()]

    with h5py.File(h5ad_path, "r") as frozen:
        phenotype_names = [_decode(v) for v in frozen["obs/phenotype/categories"][:]]
        donor_names = [_decode(v) for v in frozen["obs/donor_id/categories"][:]]
        phenotype_codes = frozen["obs/phenotype/codes"][:]
        donor_codes = frozen["obs/donor_id/codes"][:]
        markers = [_decode(v) for v in frozen["var/_index"][:]]
        missing_markers = sorted(set(rule_markers) - set(markers))
        if missing_markers:
            raise ValueError(f"rule markers absent from H5AD: {missing_markers}")

        indices, weights, strata = select_stratified_indices(
            phenotype_codes, donor_codes, phenotype_names, donor_names, per_stratum, seed
        )
        expression = frozen["X"][indices, :]
        cell_ids = [_decode(v) for v in frozen["obs/_index"][indices]]
        image_ids = [_decode(v) for v in frozen["obs/imageid"][indices]]
        truth = np.asarray(phenotype_names, dtype=object)[phenotype_codes[indices]]
        donors = np.asarray(donor_names, dtype=object)[donor_codes[indices]]
        n_cells, n_markers = frozen["X"].shape
        frozen_counts = np.bincount(phenotype_codes, minlength=len(phenotype_names))

    sampled = ad.AnnData(
        X=expression,
        obs=pd.DataFrame({"imageid": image_ids}, index=cell_ids),
        var=pd.DataFrame(index=markers),
    )
    # scimap is verbose; the report records results without streaming every
    # phenotype step to the terminal.
    with contextlib.redirect_stdout(io.StringIO()):
        called = sm.tl.phenotype_cells(sampled, phenotype=rules, label="replayed_phenotype")
    literal = called.obs["replayed_phenotype"].astype(str).to_numpy()
    aliased = apply_comparison_aliases(literal)
    corrected = apply_tcell_correction(aliased, expression, markers)

    confusion = Counter(
        (str(saved), str(predicted))
        for saved, predicted in zip(truth, corrected)
        if saved != predicted
    )
    top_confusions = [
        {"saved": saved, "predicted": predicted, "count": count}
        for (saved, predicted), count in sorted(
            confusion.items(), key=lambda pair: (-pair[1], pair[0][0], pair[0][1])
        )[:12]
    ]
    report = {
        "audit": "frozen_phenotype_rule_replay",
        "inputs": {
            "h5ad": str(h5ad_path),
            "h5ad_size_bytes": h5ad_path.stat().st_size,
            "h5ad_mtime_ns": h5ad_path.stat().st_mtime_ns,
            "rules": str(rules_path),
            "rules_sha256": hashlib.sha256(rules_path.read_bytes()).hexdigest(),
            "scimap_version": version("scimap"),
        },
        "sampling": {
            "method": "up to N uniformly sampled without replacement per saved phenotype x donor",
            "seed": seed,
            "per_stratum": per_stratum,
            "sample_n": int(indices.size),
            "stratum_n": len(strata),
            "sample_indices_sha256": hashlib.sha256(indices.astype("<i8").tobytes()).hexdigest(),
            "sample_by_saved_phenotype": _counts(truth),
            "sample_by_donor": _counts(donors),
            "strata": strata,
        },
        "frozen": {
            "n_cells": n_cells,
            "n_markers": n_markers,
            "n_phenotypes": len(phenotype_names),
            "n_donors": len(donor_names),
            "phenotype_counts": dict(zip(phenotype_names, map(int, frozen_counts))),
        },
        "rules": {
            "n_rows": len(rule_labels),
            "row_labels": rule_labels,
            "saved_labels_absent_from_rules": sorted(set(phenotype_names) - set(rule_labels)),
        },
        "predictions": {
            "literal_counts": _counts(literal),
            "comparison_aliases": COMPARISON_ALIASES,
            "aliases_are_historical_evidence": False,
            "tcell_correction_predicted_lineage_n": int(np.isin(aliased, list(TCELL_LABELS)).sum()),
            "literal": agreement(truth, literal, weights),
            "after_aliases": agreement(truth, aliased, weights),
            "after_aliases_and_tcell_correction": agreement(truth, corrected, weights),
            "top_confusions_after_correction": top_confusions,
        },
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--per-stratum", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if not 1 <= args.per_stratum <= 100:
        parser.error("--per-stratum must be between 1 and 100 to keep this audit bounded")

    report = run_audit(args.h5ad, args.rules, args.per_stratum, args.seed)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Report: {args.report}")
    print(f"Sampled cells: {report['sampling']['sample_n']:,}")
    for name in ("literal", "after_aliases", "after_aliases_and_tcell_correction"):
        score = report["predictions"][name]
        print(f"{name}: {score['sample_agreement']:.2%} sampled; "
              f"{score['population_weighted_agreement']:.2%} population-weighted")


if __name__ == "__main__":
    main()
