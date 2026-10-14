#!/usr/bin/env python3
"""Owner-accepted protein rules, QC and historical gate-scale reconstruction.

These reconstruction rules were accepted on 2026-10-03; they are not recovered
historical labels. Unmatched
or conflicting top-level evidence yields Unknown. Acinar, ECAD+ and Structural
are not invented from negative-only evidence. Original labels remain untouched.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np

# Original notebook QC (SCVI_CODEX_v2.ipynb): keep cells with >=2 detected
# markers and total raw MFI >= the 1st percentile of per-cell totals over the
# UNFILTERED original cohort (/data/panc_CODEX/results/panc_CODEX.h5ad,
# 2,676,419 cells). With that threshold the rule reproduces the frozen
# phenotyped H5AD's 2,649,654 cells exactly. Re-derive with
# derive_original_qc_threshold(); the value is pinned so the release build does
# not depend on the out-of-repo unfiltered file.
QC_TOTAL_MFI = 5234.796272000001
QC_PERCENTILE = 0.01
ORIGINAL_UNFILTERED_CELLS = 2_676_419
ORIGINAL_RETAINED_CELLS = 2_649_654
RULES = Path(__file__).with_name('protein_phenotype_rules.csv')


def read_rules(path=RULES):
    with Path(path).open() as stream:
        rules = list(csv.DictReader(stream))
    seen = {'all'}
    for rule in rules:
        if rule['parent'] not in seen or rule['label'] in seen:
            raise ValueError('Rules must be a unique, parent-before-child tree')
        for key in ('all_positive', 'any_positive', 'all_negative'):
            rule[key] = [v for v in rule[key].split('|') if v]
        if not (rule['all_positive'] or rule['any_positive']):
            raise ValueError('Every proposed lineage requires positive evidence')
        seen.add(rule['label'])
    return rules


def derive_original_qc_threshold(path, expected_kept=ORIGINAL_RETAINED_CELLS, chunk=50_000):
    """Recompute QC_TOTAL_MFI from an unfiltered original H5AD (raw MFI in X).

    Returns (threshold, n_cells, n_retained) and raises if the >=2-markers plus
    total-MFI rule does not reproduce the frozen retained-cell count.
    """
    import h5py
    with h5py.File(path, 'r') as h5:
        matrix = h5['X']
        totals = np.empty(matrix.shape[0], dtype=np.float64)
        detected = np.empty(matrix.shape[0], dtype=np.int64)
        for start in range(0, len(totals), chunk):
            values = np.asarray(matrix[start:start + chunk], dtype=np.float64)
            totals[start:start + len(values)] = values.sum(axis=1)
            detected[start:start + len(values)] = (values > 0).sum(axis=1)
    threshold = float(np.quantile(totals, QC_PERCENTILE))
    retained = int(((detected >= 2) & (totals >= threshold)).sum())
    if expected_kept is not None and retained != expected_kept:
        raise ValueError('Original notebook QC does not reproduce the frozen cell count')
    return threshold, int(len(totals)), retained


def qc_mask(raw):
    raw = np.asarray(raw, dtype=np.float64)
    finite = np.isfinite(raw)
    complete, absent = finite.all(axis=1), ~finite.any(axis=1)
    if np.isinf(raw).any() or not np.all(complete | absent) or np.any(raw[finite] < 0):
        raise ValueError('Negative, infinite or partially missing protein measurements')
    return complete & ((raw > 0).sum(axis=1) >= 2) & (raw.sum(axis=1) >= QC_TOTAL_MFI)


def rescale_proteins(raw, gates):
    """Per-donor piecewise min-max scaling, with nearest observed gate anchors.

    Only QC-retained cells enter here. A constant feature maps to 0 or 1 by
    its actual gate; this avoids division by zero or false positive support.
    """
    logged = np.log1p(np.asarray(raw, dtype=np.float64))
    gates = np.asarray(gates, dtype=np.float64)
    if logged.ndim != 2 or not len(logged) or gates.shape != (logged.shape[1],):
        raise ValueError('Nonempty matrix and one gate per marker required')
    if not np.isfinite(logged).all() or not np.isfinite(gates).all():
        raise ValueError('Finite measured protein values and gates required')
    lower, upper = logged.min(axis=0), logged.max(axis=0)
    anchors = logged[np.argmin(np.abs(logged - gates), axis=0), np.arange(len(gates))]
    result = np.empty_like(logged)
    for j in range(len(gates)):
        if lower[j] == upper[j]:
            result[:, j] = float(lower[j] >= gates[j])
            continue
        # An out-of-range gate has no observed split. Preserve its actual
        # positive/negative meaning instead of declaring the extreme cell positive.
        if gates[j] > upper[j]:
            result[:, j] = 0.5 * (logged[:, j] - lower[j]) / (gates[j] - lower[j])
            continue
        if gates[j] < lower[j]:
            result[:, j] = 0.5 + 0.5 * (logged[:, j] - gates[j]) / (upper[j] - gates[j])
            continue
        lo = logged[:, j] <= anchors[j]
        result[lo, j] = (0.5 * (logged[lo, j] - lower[j]) / (anchors[j] - lower[j])
                         if anchors[j] > lower[j] else 0.5)
        result[~lo, j] = (0.5 + 0.5 * (logged[~lo, j] - anchors[j]) / (upper[j] - anchors[j])
                          if upper[j] > anchors[j] else 0.5)
    return result, logged, dict(minimum=lower.tolist(), maximum=upper.tolist(),
        anchors=anchors.tolist(), gates=gates.tolist(),
        gate_outside_observed_range=((gates < lower) | (gates > upper)).tolist())


def classify(expression, markers, rules=None):
    rules = read_rules() if rules is None else rules
    x = np.asarray(expression)
    if x.ndim != 2 or x.shape[1] != len(markers) or len(set(markers)) != len(markers):
        raise ValueError('Expression columns must match unique markers')
    needed = {m for r in rules for key in ('all_positive','any_positive','all_negative') for m in r[key]}
    if needed - set(markers):
        raise ValueError(f'Missing required markers: {sorted(needed - set(markers))}')
    idx = {m: i for i, m in enumerate(markers)}
    finite = np.isfinite(x).all(axis=1)
    positives = x >= 0.5
    labels = np.full(len(x), 'Unknown', dtype=object)
    reason = np.full(len(x), 'no_supported_lineage', dtype=object)
    reason[~finite] = 'missing_measurements'
    levels = ['all'] + list(dict.fromkeys(r['label'] for r in rules))
    for parent in levels:
        children = [r for r in rules if r['parent'] == parent]
        if not children:
            continue
        eligible = finite if parent == 'all' else finite & (labels == parent)
        matches = []
        for r in children:
            mask = eligible.copy()
            if r['all_positive']:
                mask &= positives[:, [idx[m] for m in r['all_positive']]].all(axis=1)
            if r['any_positive']:
                mask &= positives[:, [idx[m] for m in r['any_positive']]].any(axis=1)
            if r['all_negative']:
                mask &= ~positives[:, [idx[m] for m in r['all_negative']]].any(axis=1)
            matches.append(mask)
        support = np.sum(matches, axis=0)
        reason[eligible & (support > 1)] = 'conflicting_' + parent + '_evidence'
        for r, mask in zip(children, matches):
            use = mask & (support == 1)
            labels[use] = r['label']
            reason[use] = 'rule_supported'
    return labels, reason
