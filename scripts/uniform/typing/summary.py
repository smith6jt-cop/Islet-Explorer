#!/usr/bin/env python3
"""summary.py -- one row per image of the typing run, from the outputs of classify, refine and endocrine_checks.

Reports only; nothing here enters the pipeline. Written to `$OUT/typing/summary.csv`.

    python scripts/uniform/typing/summary.py --out data/development/uniform-20donor-20261015
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    typ = CT.typing_dir(a.out)
    checks = json.loads((typ / 'endocrine_checks.json').read_text())
    rows = []
    for r in CT.cohort_rows():
        t = r['image_name']
        audit = json.loads((typ / 'calls' / f'{t}.json').read_text())
        cls, imm = audit['cell_class'], audit['immune_definition']
        refined = pd.read_parquet(typ / 'refine' / f'{t}.parquet', columns=['refined_class'])['refined_class']
        ex = checks[t].get('Excluded', {})
        level = audit['hormone_level_in_objects']
        rows.append({
            'image': t, 'status': r['donor_status'],
            'cells': audit['n_cells'] + audit['n_not_analysis_eligible'], 'analysis_eligible': audit['n_cells'],
            'in_endocrine_objects': checks[t]['cells_in_objects'],
            'Beta': cls.get('Beta', 0), 'Alpha': cls.get('Alpha', 0), 'Delta': cls.get('Delta', 0),
            'Endocrine_unresolved': cls.get('Endocrine (unresolved subtype)', 0)
            + cls.get('Endocrine (ambiguous subtype)', 0),
            'Excluded': cls.get('Excluded', 0), 'Excluded_chga': ex.get('chga', 0),
            'Excluded_hormone_below_level': ex.get('hormone_below_level', 0),
            'Excluded_hormone_below_level_expected': ex.get('hormone_below_level_expected', 0),
            'Excluded_no_endocrine_evidence': ex.get('no_endocrine_evidence', 0),
            'Immune': cls.get('Immune', 0), 'T_cell': imm.get('Imm-T-cell', 0),
            'macrophage': imm.get('Imm-macrophage', 0), 'CD45_only': imm.get('Imm-CD45-leukocyte', 0),
            'B_cell': imm.get('Imm-B-cell', 0), 'Exocrine': cls.get('Exocrine', 0),
            'refined_Unresolved': int((refined == 'Unresolved (no lineage evidence)').sum()),
            'hormone_level_tissue_p99': round(level['tissue_p99'], 4), 'hormone_level': round(level['level'], 4)})
    pd.DataFrame(rows).to_csv(typ / 'summary.csv', index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
