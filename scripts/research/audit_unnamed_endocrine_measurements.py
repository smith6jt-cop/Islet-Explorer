#!/usr/bin/env python3
"""Compare unresolved object labels with source measurements without changing calls.

Raw QuPath means precede despeckling and REDSEA. Comparing them to the current
corrected-value level is a diagnostic of measurement changes, NOT an alternate
typing method or evidence that a cell should receive an endocrine name.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

MARKERS = ('INS', 'GCG', 'SST', 'CHGA')
UNRESOLVED = ('Islet', 'EndocrineCluster', 'EndocrineSingle')
PREVIOUS_NAMED = ('Alpha', 'Beta', 'Delta', 'CHGA+ hormone-negative')
CURRENT_NAMED = ('Alpha cell', 'Beta cell', 'Delta cell', 'CHGA+ hormone-negative')


def describe(frame):
    return {
        'cells': len(frame),
        'raw_above_current_corrected_level': int(frame.raw_any_above_level.sum()),
        'corrected_above_current_level': int(frame.corrected_any_above_level.sum()),
        'previous_named_endocrine': int(frame.previous_label.isin(PREVIOUS_NAMED).sum()),
        'cell_area_um2_median': float(frame.cell_area_um2.median()),
        'nucleus_cell_fraction_median': float(frame.nucleus_cell_fraction.median()),
        'nucleus_cell_fraction_p90': float(frame.nucleus_cell_fraction.quantile(.9)),
        'zero_cytoplasm_pixels': int((frame.cytoplasm_pixels == 0).sum()),
        'raw_positive_then_corrected_negative': {
            marker: int((frame['raw_above__' + marker] & ~frame['corrected_above__' + marker]).sum())
            for marker in MARKERS},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error('output already exists; choose a fresh directory')
    levels = pd.read_csv(args.run / 'typing/levels.csv').set_index(['image', 'marker'])
    frames, checks, reference = [], [], []
    for path in sorted((args.run / 'typing/identity').glob('*.parquet')):
        tag = path.stem
        ident = pd.read_parquet(path, columns=['object_id', 'label', 'status'] + ['pos__' + m for m in MARKERS])
        corrected = pd.read_parquet(args.run / 'typing/redsea' / path.name,
                                    columns=['object_id', 'Cell__n_px', 'Nucleus__n_px', 'Cytoplasm__n_px']
                                    + ['Cell__' + m for m in MARKERS])
        all_cells = ident.merge(corrected, on='object_id', validate='one_to_one', how='left', indicator=True)
        if not all_cells._merge.eq('both').all():
            raise ValueError(f'{tag}: identity missing from measurements')
        for m in MARKERS:
            level = levels.loc[(tag, m)]
            if level.measurement != 'whole_cell':
                raise ValueError(f'{tag}/{m}: unexpected measurement {level.measurement}')
            observed = np.log1p(all_cells['Cell__' + m].to_numpy(float)) > level.level_log
            if not np.array_equal(observed, all_cells['pos__' + m].to_numpy(bool)):
                raise ValueError(f'{tag}/{m}: stored positivity disagrees with measurements')
        checks.append({'image': tag, 'verified_cells': len(all_cells)})
        all_cells['cell_area_um2'] = all_cells.Cell__n_px * .3774 ** 2
        all_cells['nucleus_cell_fraction'] = all_cells.Nucleus__n_px / all_cells.Cell__n_px
        named = all_cells[all_cells.label.isin(CURRENT_NAMED)]
        reference.append({'image': tag, 'named_endocrine_cells': len(named),
                          'cell_area_um2_median': float(named.cell_area_um2.median()),
                          'nucleus_cell_fraction_median': float(named.nucleus_cell_fraction.median())})
        d = all_cells[all_cells.label.isin(UNRESOLVED)].drop(columns='_merge').copy()
        raw = pd.read_parquet(args.run / 'cells' / path.name,
                              columns=['Object ID'] + [f'Cell: {m}: Mean' for m in MARKERS])
        raw = raw.rename(columns={'Object ID': 'object_id', **{f'Cell: {m}: Mean': 'raw__' + m for m in MARKERS}})
        old = pd.read_parquet(args.run / 'superseded/typing-v17-20261028/calls' / path.name,
                              columns=['object_id', 'cell_class', 'excluded_reason'])
        old = old.rename(columns={'cell_class': 'previous_label', 'excluded_reason': 'previous_excluded_reason'})
        d = d.merge(raw, on='object_id', how='left', validate='one_to_one', indicator=True)
        if not d._merge.eq('both').all():
            raise ValueError(f'{tag}: unresolved cell missing from raw measurements')
        d = d.drop(columns='_merge').merge(old, on='object_id', how='left', validate='one_to_one', indicator=True)
        if not d._merge.eq('both').all():
            raise ValueError(f'{tag}: unresolved cell missing from previous calls')
        d = d.drop(columns='_merge').rename(columns={'Cytoplasm__n_px': 'cytoplasm_pixels'})
        for m in MARKERS:
            level = levels.loc[(tag, m)]
            d['level__' + m] = level.level_raw
            d['raw_above__' + m] = np.log1p(d['raw__' + m].to_numpy(float)) > level.level_log
            d['corrected_above__' + m] = d['pos__' + m].astype(bool)
        d['raw_any_above_level'] = d[['raw_above__' + m for m in MARKERS]].any(axis=1)
        d['corrected_any_above_level'] = d[['corrected_above__' + m for m in MARKERS]].any(axis=1)
        d['image'] = tag
        frames.append(d)
    if not frames:
        raise ValueError('No identity files found')
    combined = pd.concat(frames, ignore_index=True)
    report = {'run': str(args.run.resolve()), 'measurement_note': __doc__,
              'positivity_checks': checks, 'totals': describe(combined),
              'by_label': {label: describe(g) for label, g in combined.groupby('label')},
              'by_image': {tag: describe(g) for tag, g in combined.groupby('image')},
              'named_endocrine_reference': reference}
    args.out.mkdir(parents=True)
    combined.to_parquet(args.out / 'cells.parquet', index=False)
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['totals'], indent=2))


if __name__ == '__main__':
    main()
