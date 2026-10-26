#!/usr/bin/env python3
"""Reproduce the exact H1 input CSV used for the 2026-10-25 paper-method models.

Read the corrected feature table, reuse its established eligibility function,
and add the six algebraic model variables. Existing CSVs are verified, never
overwritten. This script does not fit or modify any statistical model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import h1_analysis_20261022 as h1

REPO = Path(__file__).resolve().parents[2]
INPUT = REPO / 'output/h1_edge_corrected_20261024/h1_object_features.parquet'
OUTPUT = REPO / 'output/h1_h4_literature_methods_20261025/h1/analysis_objects.csv'
EXPECTED_CSV_SHA256 = 'befd01e84430e54bc8eaf970c5212134557766355aa3d58ca9ef6ca7f774dcf1'


def build_frame() -> pd.DataFrame:
    d = h1.eligible(pd.read_parquet(INPUT))
    d['allocation10'] = d.immune_localization_fraction / .1
    d['log_area'] = np.log(d.core_area_um2)
    d['log_area_ratio'] = np.log(d.peri_area_um2 / d.core_area_um2)
    d['log_nearby'] = np.log1p(d.n_immune_nearby)
    d['log_field'] = np.log1p(d.field_density)
    d['beta_density_log2'] = np.log2(1 + d.beta_density_mm2 / 1000)
    if len(d) != 2937 or d.case_id.nunique() != 20 or not d.islet_id.is_unique:
        raise ValueError('Corrected H1 population differs from the model input population')
    return d


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT,
                        help='Create only if absent; otherwise verify exact bytes')
    parser.add_argument('--check', action='store_true', help='Require an existing output; write no CSV')
    args = parser.parse_args()
    d = build_frame()
    payload = d.to_csv(index=False).encode('utf-8')
    generated_sha = hashlib.sha256(payload).hexdigest()
    if generated_sha != EXPECTED_CSV_SHA256:
        raise ValueError(f'Regenerated input differs from the recorded model input: {generated_sha}')
    if args.output.exists():
        if args.output.read_bytes() != payload:
            raise ValueError(f'Existing CSV differs from regenerated input; not overwriting {args.output}')
        action = 'verified_existing_byte_identical_csv'
    elif args.check:
        raise FileNotFoundError(args.output)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(payload)
        action = 'created_byte_identical_csv'
    report = dict(action=action, rows=len(d), columns=len(d.columns), donors=int(d.case_id.nunique()),
                  source=str(INPUT), source_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
                  eligibility_script=str(Path(h1.__file__).resolve()),
                  eligibility_script_sha256=hashlib.sha256(Path(h1.__file__).read_bytes()).hexdigest(),
                  csv=str(args.output), csv_sha256=generated_sha, byte_identical=True,
                  existing_models_and_results_modified=False)
    (args.output.parent / 'input_reproduction_check.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
