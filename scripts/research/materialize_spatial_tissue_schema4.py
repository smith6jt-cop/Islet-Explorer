#!/usr/bin/env python3
"""Create schema-4 donor tissue tables from corrected cells and current typing."""
from pathlib import Path
import argparse
import pandas as pd

def main(root: Path):
    corrected = root / 'corrected'
    typed = root / 'typing' / 'cells'
    if not typed.is_dir():
        raise SystemExit('current schema-4 typing/cells is required; no legacy fallback')
    memberships = pd.read_parquet(root / 'analysis/memberships.parquet')
    memberships['case_id'] = memberships.case_id.astype(str)
    if 'band_um' in memberships:
        memberships = memberships[memberships.band_um.eq(20)].copy()
    islets = pd.read_parquet(root / 'analysis/islets.parquet', columns=['islet_id','islet_key','area_um2'])
    memberships = memberships.merge(islets, on='islet_id', how='left', validate='many_to_one')
    # One plotted point per raw cell. Prefer core membership, then the smallest
    # containing islet when regenerated candidate regions overlap.
    memberships['compartment_order'] = memberships.compartment.map({'core': 0, 'peri': 1}).fillna(2)
    memberships = memberships.sort_values(['case_id','object_id','compartment_order','area_um2','islet_id']).drop_duplicates(['case_id','object_id'])
    out = root / 'donors'; out.mkdir(exist_ok=True)
    files = sorted(corrected.glob('*-cells.parquet'))
    if not files: raise SystemExit('no corrected cell tables')
    for src in files:
        donor = src.stem.split('-')[0]
        c = pd.read_parquet(src, columns=['object_id','X_centroid','Y_centroid','archival_phenotype'])
        tpath = typed / f'{donor}.parquet'
        if tpath.exists():
            t = pd.read_parquet(tpath, columns=['object_id','base_type','phenotype_state'])
            t['object_key'] = t.object_id.astype(str).str.split(':').str[-1]
            c['object_key'] = c.object_id.astype(str).str.split(':').str[-1]
            c = c.merge(t[['object_key','base_type','phenotype_state']], on='object_key', how='left')
        else:
            c['base_type'] = None; c['phenotype_state'] = None
        base = c['base_type'].fillna('').astype(str)
        # Current typing is authoritative. Unassigned cells stay explicit rather
        # than silently inheriting the retired legacy phenotype labels.
        phenotype = base.where(base.str.len() > 0, 'Unresolved')
        member = memberships[memberships.case_id.eq(str(donor))][['object_id','compartment','islet_key']].copy()
        member['object_key'] = member.object_id.astype(str).str.split(':').str[-1]
        c = c.merge(member[['object_key','compartment','islet_key']], on='object_key', how='left')
        outdf = pd.DataFrame({
            'X_centroid': pd.to_numeric(c.X_centroid, errors='coerce'),
            'Y_centroid': pd.to_numeric(c.Y_centroid, errors='coerce'),
            'phenotype': phenotype,
            'cell_region': c.compartment.fillna('tissue'),
            'islet_name': c.islet_key.fillna('')
        })
        outdf = outdf.dropna(subset=['X_centroid','Y_centroid'])
        outdf.to_csv(out / f'{donor}.csv', index=False)
    print(f'wrote {len(files)} donor tissue tables to {out}')

if __name__ == '__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--input-root',type=Path,required=True)
    main(ap.parse_args().input_root)
