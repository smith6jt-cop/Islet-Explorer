#!/usr/bin/env python3
"""Prepare an isolated project for the two source-reconciled acquisitions."""
import datetime
import json
from pathlib import Path
import shutil

import numpy as np
import tifffile

from convert_images import sha256, write_json


REPO = Path(__file__).resolve().parents[2]
BASE = REPO / 'data/development/expanded-20260911'
ROOT = BASE / 'recovery'


def main():
    if list((ROOT / 'logs').glob('*/*.done.json')):
        raise ValueError('Recovered processing has started; preserve its project/checkpoints')
    provenance_path = ROOT / 'recovered_project_provenance.json'
    if provenance_path.exists():
        raise ValueError('Recovered project already prepared; run its stages')
    matched = json.loads((BASE / 'reports/image_conversion/source_recovery_6090.json').read_text())
    rebuilt = json.loads((ROOT / 'source_reconstruction_6539.json').read_text())
    if not matched['all_planes_identical'] or len(matched['planes']) != 36 or rebuilt['status'] != 'validated':
        raise ValueError('Incomplete source recovery evidence')
    sources = {'6090': Path(matched['assembled_input']), '6539': Path(rebuilt['output'])}
    if sha256(sources['6539']) != rebuilt['output_sha256']:
        raise ValueError('Recovered 6539 changed after pixel validation')
    original_6090 = json.loads((BASE / 'reports/image_conversion/Aab_6090.json').read_text())
    if sha256(sources['6090']) != original_6090['source_sha256']:
        raise ValueError('6090 source changed after pixel validation')
    project_path = ROOT / 'qupath/project.qpproj'
    project = json.loads(project_path.read_text())
    entries = {39:'6090', 38:'6539'}
    project['images'] = [entry for entry in project['images'] if entry['entryID'] in entries]
    historical = ROOT / 'historical_annotations'
    historical.mkdir(exist_ok=True)
    cases = {}
    for entry in project['images']:
        case = entries[entry['entryID']]
        original = ROOT / 'qupath/data' / str(entry['entryID'])
        if original.exists():
            original.rename(historical / str(entry['entryID']))
        entry['serverBuilder'] = {'builderType':'uri', 'providerClassName':'qupath.imagej.images.servers.ImageJServerBuilder',
                                  'uri':sources[case].as_uri(), 'args':[]}
        cases[case] = {'source':str(sources[case]), 'source_sha256':sha256(sources[case]),
                       'historical_annotations':str(historical / str(entry['entryID'])),
                       'annotation_transfer':'None; retain old annotation data for review against its original pixels'}
    write_json(project_path, project)
    channels = []
    with tifffile.TiffFile(sources['6539']) as tif:
        labels = tif.imagej_metadata['Labels']
    stack = tifffile.memmap(sources['6539'], mode='r')
    for index, label in enumerate(labels):
        low, high = np.quantile(stack[index, ::8, ::8], [0.005, 0.998])
        if high <= low:
            raise ValueError(f'No usable display range for {label}')
        channels.append({'index':index, 'source_label':label, 'minDisplay':float(low), 'maxDisplay':float(high)})
    del stack
    write_json(ROOT / 'display_ranges_6539.json', {
        'source_sha256':rebuilt['output_sha256'], 'purpose':'SAM boundary refinement and diagnostic rendering only',
        'method':'0.5th and 99.8th percentiles of each original channel sampled every eighth source pixel in X and Y',
        'raw_measurements_changed':False, 'channels':channels})
    reports = ROOT / 'reports/image_conversion'
    reports.mkdir(parents=True, exist_ok=True)
    shutil.copy2(BASE / 'reports/image_conversion/Aab_6090.json', reports / 'Aab_6090.json')
    write_json(provenance_path, {'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'status':'source_reconciled', 'cases':cases,
        'source_match_report_sha256':sha256(BASE / 'reports/image_conversion/source_recovery_6090.json'),
        'reconstruction_report_sha256':sha256(ROOT / 'source_reconstruction_6539.json'),
        'display_ranges_sha256':sha256(ROOT / 'display_ranges_6539.json'),
        'limitation':'The old 6539 Ignore* annotation was made on duplicated 6090 pixels; it is preserved but cannot be transferred to the distinct recovered 6539 acquisition.'})
    print(project_path, flush=True)


if __name__ == '__main__':
    main()
