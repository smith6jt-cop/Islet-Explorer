#!/usr/bin/env python3
"""Build compact app summaries from canonical revised tables, never Excel."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import MARKERS, require_cohort, sha256, write_json


def read_table(root, relative):
    base = root / relative
    if base.with_suffix('.parquet').exists():
        return pd.read_parquet(base.with_suffix('.parquet'))
    return pd.read_csv(base.with_suffix('.csv'), dtype={'case_id': str})


def flags(values):
    return values.astype(str).str.lower().isin(['true', 't', '1'])


def build(root):
    root = Path(root).resolve(strict=True)
    if str(json.loads((root/'candidate.json').read_text())['schema_version']) != '4':
        raise ValueError('Canonical summary materialization requires schema 4')
    if (root/'preview/service.json').exists():
        raise ValueError('Refusing to mutate a candidate with a recorded preview')
    cohort = pd.read_csv(root/'cohort.csv', dtype={'case_id': str})
    meta = pd.read_csv(root/'donor_metadata.csv', dtype={'case_id': str})
    islets = pd.read_parquet(root/'analysis/islets.parquet')
    counts = read_table(root, 'trajectory/counts')
    endocrine = read_table(root, 'trajectory/endocrine_counts')
    frames = [meta, islets, counts, endocrine]
    for frame in frames:
        frame['case_id'] = frame.case_id.astype(str)
        require_cohort(frame, cohort.case_id)
    if islets.islet_id.duplicated().any() or meta.case_id.duplicated().any():
        raise ValueError('Duplicate canonical islet or donor identifier')
    if not set(counts.islet_id).issubset(islets.islet_id) or set(endocrine.islet_id) != set(islets.islet_id):
        raise ValueError('Canonical counts and geometry revisions do not align')
    if counts.duplicated(['case_id', 'islet_id', 'compartment', 'phenotype']).any():
        raise ValueError('Duplicate canonical phenotype counts')
    catalog = json.loads((root/'typing/catalog.json').read_text())
    labels = [x if isinstance(x, str) else x['name'] for x in catalog['base_types'] + catalog.get('states', [])]
    # Keep the complete delivered display vocabulary available to the app even
    # when the current scientific typing catalog contains only a subset.
    for entry in catalog.get('display_types', []):
        label = entry if isinstance(entry, str) else entry.get('name')
        if label and label not in labels:
            labels.append(label)
    if set(counts.phenotype) - set(labels):
        raise ValueError('Counts contain identities absent from the typing catalog')
    # Canonical metadata owns clinical fields; geometry owns size and identity.
    joined = islets.drop(columns=[x for x in meta if x != 'case_id' and x in islets], errors='ignore').merge(meta, on='case_id', validate='many_to_one')
    joined['Case ID'] = joined.case_id.astype(int)
    joined['Donor Status'] = joined.donor_status
    joined['islet_diam_um'] = 2 * np.sqrt(joined.area_um2 / np.pi)
    joined['core_region_um2'] = joined.area_um2
    joined['peri_region_um2'] = joined['peri_area_um2']
    if 'age' not in joined:
        joined['age'] = joined.get('age_years_npod', joined.get('age_years_key', np.nan))
    if 'gender' not in joined:
        joined['gender'] = joined.get('sex', joined.get('donor_gender', np.nan))
    lookup = joined.set_index('islet_id')
    counts = counts.copy()
    availability_column = 'validated' if 'validated' in counts else 'available'
    counts.loc[~flags(counts[availability_column]), 'count'] = np.nan
    summary = joined.copy()
    for region in ['core', 'peri']:
        subset = counts[counts.compartment == region]
        totals = subset.groupby('islet_id').n_measured.agg(['first', 'nunique'])
        if (totals['nunique'] > 1).any():
            raise ValueError('Phenotype rows disagree on measured-cell totals')
        summary['total_cells_'+region] = summary.islet_id.map(totals['first'])
        wide = subset.pivot(index='islet_id', columns='phenotype', values='count')
        for label in labels:
            values = summary.islet_id.map(wide[label]) if label in wide else np.nan
            total = summary['total_cells_'+region]
            prefix = '' if region == 'core' else 'peri_'
            summary[prefix+'count_'+label] = values
            summary[prefix+'prop_'+label] = np.where(total > 0, values / total, np.nan)
    summary['cells_total'] = summary.total_cells_core
    summary['endocrine_count'] = summary.islet_id.map(endocrine.set_index('islet_id').endocrine_count)
    # Load only the all-cell population, not the entire population × marker table.
    path = root/'trajectory/marker_summaries.parquet'
    marker = pd.read_parquet(path, filters=[('population', '=', 'All cells')]) if path.exists() else read_table(root, 'trajectory/marker_summaries').query('population == "All cells"')
    marker['case_id'] = marker.case_id.astype(str)
    if set(marker.marker) != set(MARKERS):
        raise ValueError('All 31 markers must be represented in canonical marker summaries')
    if marker.duplicated(['case_id','islet_id','compartment','marker']).any():
        raise ValueError('Duplicate canonical marker measurements')
    marker.loc[~flags(marker.positive_available), 'positive_count'] = np.nan
    for name, protein in [('Ins_any', 'INS'), ('Glu_any', 'GCG'), ('Stt_any', 'SST')]:
        subset = marker[(marker.compartment == 'core') & (marker.marker == protein)]
        summary[name] = summary.islet_id.map(subset.set_index('islet_id').positive_count)
    # No independent joint hormone counts can be inferred from marginals.
    for name in ['Ins_single','Glu_single','Stt_single','Multi_Pos','Triple_Neg']:
        if name not in summary:
            summary[name] = np.nan
    marker = marker.rename(columns={'positive_count': 'pos_count'})
    marker['region_type'] = marker.compartment.map({'core': 'islet_core', 'peri': 'islet_band'})
    marker['pos_frac'] = np.where(marker.n_cells > 0, marker.pos_count / marker.n_cells, np.nan)
    marker['pos_pct'] = marker.pos_frac * 100
    marker = marker.merge(joined.drop(columns=['case_id']), on='islet_id', validate='many_to_one', suffixes=('', '_islet'))
    targets_path = root/'analysis/spatial_measurements.parquet'
    if not targets_path.exists():
        raise ValueError('Canonical structural measurements are required; no old-cohort fallback')
    targets = pd.read_parquet(targets_path)
    if 'band_um' in targets:
        targets = targets[targets.band_um == 20]
    targets = targets[targets.compartment.isin(['core', 'peri'])].copy()
    if targets.duplicated(['case_id', 'islet_id', 'compartment', 'structure']).any():
        raise ValueError('Duplicate canonical structural areas')
    targets.loc[~flags(targets.available), 'positive_area_um2'] = np.nan
    keys = ['case_id', 'islet_id', 'structure']
    core = targets[targets.compartment == 'core'].set_index(keys)
    peri = targets[targets.compartment == 'peri'].set_index(keys)
    union = core[['positive_area_um2', 'valid_area_um2']].add(peri[['positive_area_um2', 'valid_area_um2']])
    union['compartment'] = 'both'; union['available'] = core.available & peri.available
    targets = pd.concat([targets, union.reset_index()], ignore_index=True)
    targets = targets.rename(columns={'structure':'class', 'positive_area_um2':'area_um2', 'valid_area_um2':'region_um2'})
    targets['type'] = targets.compartment.map({'core':'islet_core', 'peri':'islet_band', 'both':'islet_union'})
    targets['area_density'] = np.where(targets.region_um2 > 0, targets.area_um2 / targets.region_um2, np.nan)
    targets['count'] = np.nan  # structural area never masquerades as a cell count
    targets = targets.merge(joined.drop(columns=['case_id','area_um2']), on='islet_id', validate='many_to_one')
    pairs = pd.DataFrame([{'stem': label, 'label': label, 'core_col': 'prop_'+label, 'peri_col':'peri_prop_'+label} for label in labels])
    output = {'core_area': joined[['Case ID','Donor Status','islet_key','islet_diam_um','core_region_um2','peri_region_um2']],
              'islets':summary, 'markers':marker, 'targets':targets, 'combined_pairs':pairs}
    destination = root/'summary'; destination.mkdir(exist_ok=True)
    for name, table in output.items():
        path = destination/(name+'.parquet')
        temporary = path.with_suffix('.parquet.tmp')
        table.to_parquet(temporary, index=False, compression='zstd')
        temporary.replace(path)
    write_json(root/'summary/materialization.json', {'schema_version':'4', 'source':'canonical typed measurements; no Excel/H5AD fallback',
        'rows':{name:len(table) for name,table in output.items()},
        'sources':{str(path.relative_to(root)):sha256(path) for path in [root/'analysis/islets.parquet', root/'analysis/spatial_measurements.parquet',root/'typing/catalog.json',root/'donor_metadata.csv']}})
    return {name:len(table) for name,table in output.items()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, required=True)
    print(json.dumps(build(parser.parse_args().input_root), indent=2))
