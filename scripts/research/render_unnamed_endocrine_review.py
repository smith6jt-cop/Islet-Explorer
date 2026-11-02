#!/usr/bin/env python3
"""Render diagnostic evidence for cells retaining endocrine object names.

Read-only with respect to the run. For each image and unresolved label, select
the closest-to-level, median, and furthest-below-level cells by the maximum
INS/GCG/SST/CHGA fold. This is a boundary diagnostic, not a random accuracy
estimate. Reuse the established renderer's measured-scale grayscale panels,
cell outlines, DAPI, and coverage records. Never assign replacement identities.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import typing_accuracy_tiles_20261029 as renderer

LABELS = ('Islet', 'EndocrineCluster', 'EndocrineSingle')
MARKERS = ('INS', 'GCG', 'SST', 'CHGA')


def render_one(tag, rows, run):
    # Explicit channels avoid the original renderer's unrelated review strata.
    renderer.shown_markers = lambda *_: (list(MARKERS), [])
    return renderer.render_image(tag, rows, run)


def select_cells(frame):
    selected = []
    for label, group in frame[frame.label.isin(LABELS)].groupby('label'):
        group = group.copy()
        group['max_endocrine_fold'] = group[['fold__' + m for m in MARKERS]].max(axis=1)
        group = group.sort_values(['max_endocrine_fold', 'object_id'])
        for index, description in [(0, 'furthest below'), (len(group)//2, 'median'),
                                   (len(group)-1, 'closest to level')]:
            row = group.iloc[index].copy()
            row['group'] = label
            row['pair'] = description
            selected.append(row)
    return pd.DataFrame(selected).drop_duplicates('object_id') if selected else pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output already exists; choose a fresh directory')
    tasks = []
    for path in sorted((args.run / 'typing/identity').glob('*.parquet')):
        selected = select_cells(pd.read_parquet(path))
        if len(selected):
            selected['image'] = path.stem
            tasks.append((path.stem, selected))
    if not tasks:
        parser.error('no unresolved endocrine object labels found')
    args.out.mkdir(parents=True)
    pd.concat([rows for _, rows in tasks]).to_csv(args.out / 'selection.csv', index=False)
    done = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(render_one, tag, rows, args.run) for tag, rows in tasks]
        for (tag, _), future in zip(tasks, futures):
            done.extend(future.result())
            print(tag, 'rendered', flush=True)
    renderer.sheets(done, args.out)
    print('Rendered', len(done), 'cells from', len(tasks), 'images', flush=True)


if __name__ == '__main__':
    main()
