"""Append phenotype-specific raw overlap totals to an islet composition table."""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'app/shiny_app/python'))
import paga_pseudotime as engine

def augment(frame, cells_dir):
    names, counts = engine.phenotype_inputs(frame)
    labels = {name.split('_', 1)[1]: i for i, name in enumerate(names)}
    scopes = {label: names[i].split('_', 1)[0] for label, i in labels.items()}
    structures = engine.OVERLAP_STRUCTURES
    columns = ['crossed_frac_' + structure for structure in structures]
    mass = np.zeros((len(frame), len(names), len(columns)))
    if not frame.islet_id.is_unique:
        raise ValueError('Duplicate islet identifiers')
    for row, obj in enumerate(frame.itertuples(index=False)):
        cells = pd.read_csv(cells_dir / (obj.islet_id + '.csv'), usecols=['cell_region', 'phenotype', *columns])
        cls = cells.phenotype.map(engine.trajectory_class)
        scope = cls.map(scopes)
        retained = scope.eq('combined') | scope.eq(cells.cell_region)
        idx = cls[retained].map(labels).to_numpy(dtype=int)
        values = cells.loc[retained, columns].to_numpy(dtype=float)
        if not (np.isfinite(values) & (values >= 0) & (values <= 1)).all():
            raise ValueError(f'{obj.islet_id}: missing or invalid passthrough fraction')
        observed = np.bincount(idx, minlength=len(names))
        np.testing.assert_allclose(np.log1p(observed / max(obj.n_endocrine_core, 1)), counts[row], atol=1e-12)
        np.add.at(mass[row], idx, values)
    keys = ['overlap__' + name + '__' + structure for name in names for structure in structures]
    raw = pd.DataFrame(mass.reshape(len(frame), -1), columns=keys, index=frame.index)
    return pd.concat([frame.drop(columns=keys, errors='ignore'), raw], axis=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('source', 'cells', 'output'):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    augment(pd.read_csv(args.source, float_precision='round_trip'), args.cells).to_csv(args.output, index=False)
