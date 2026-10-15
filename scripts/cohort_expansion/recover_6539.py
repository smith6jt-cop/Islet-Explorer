#!/usr/bin/env python3
"""Reassemble 6539's original marker planes without changing source pixels.

The existing PAN-6539 stack incorrectly contains 6090 pixels. The raw acquisition
has 60 single-plane files; retain the same 36 named channels as the selected
panel, excluding repeated DAPI and blank cycles. Preserve the misleading stack
and its annotations as historical evidence, never overwrite them.
"""
import datetime
import json
from pathlib import Path

import numpy as np
import tifffile

from convert_images import inspect_source, sha256, write_json


REPO = Path(__file__).resolve().parents[2]
BASE = Path('/data/QuPath_Projects/QuPath_project_IsletExplorer')
ROOT = REPO / 'data/development/expanded-20260911/recovery'
SOURCE = BASE / '2022-11-11  PAN-6539-ND test 34'
TEMPLATE = BASE / 'Images/New/6090AAbtest31_35p_led.tif'


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    folder = ROOT / 'sources'
    folder.mkdir(exist_ok=True)
    output = folder / 'PAN-6539-ND test 34.tif'
    report_path = ROOT / 'source_reconstruction_6539.json'
    if output.exists():
        report = json.loads(report_path.read_text())
        if report['status'] != 'validated' or sha256(output) != report['output_sha256']:
            raise ValueError('Recovered output differs from its checkpoint')
        for row in report['planes']:
            if sha256(Path(row['source'])) != row['source_sha256']:
                raise ValueError('Original channel plane changed after reconstruction')
        print('Recovered 6539 checkpoint verified', flush=True)
        return
    with tifffile.TiffFile(TEMPLATE) as template:
        labels = template.imagej_metadata['Labels']
    inputs = [SOURCE / label.removeprefix('reg001_') for label in labels]
    planes = []
    shape = None
    calibration = None
    for path in inputs:
        with tifffile.TiffFile(path) as tif:
            plane = tif.series[0]
            if len(tif.series) != 1 or plane.axes != 'YX' or plane.dtype != np.dtype('uint16'):
                raise ValueError(f'Unexpected source plane: {path}')
            scale = {axis: tif.pages[0].tags[f'{axis}Resolution'].value for axis in ['X', 'Y']}
            if (tif.imagej_metadata or {}).get('unit') != 'mm':
                raise ValueError('Missing original millimetre calibration')
            if shape is not None and (shape != plane.shape or calibration != scale):
                raise ValueError('Raw channel shapes/calibration do not agree')
            shape, calibration = plane.shape, scale
        planes.append({'source': str(path), 'bytes': path.stat().st_size, 'source_sha256': sha256(path)})
    partial = output.with_suffix('.partial.tif')
    if partial.exists():
        raise ValueError('Incomplete reconstruction exists; inspect it before restarting')
    # The source collection uses contiguous ImageJ hyperstacks, including its
    # >4-GiB extension. Use the same format and validate through an independent
    # read; the browser conversion below produces standard tiled OME-BigTIFF.
    stack = tifffile.memmap(partial, shape=(len(inputs), *shape), dtype='uint16', imagej=True,
                            resolution=(calibration['X'], calibration['Y']),
                            metadata={'axes': 'CYX', 'unit': 'mm', 'Labels': labels})
    for index, path in enumerate(inputs):
        stack[index] = tifffile.imread(path)
        stack.flush()
        print(f'Wrote {index + 1}/{len(inputs)}: {path.name}', flush=True)
    del stack
    verified = tifffile.memmap(partial, mode='r')
    for index, path in enumerate(inputs):
        raw = tifffile.imread(path)
        if not np.array_equal(verified[index], raw) or sha256(path) != planes[index]['source_sha256']:
            raise ValueError(f'Pixel/source mismatch in recovered {path}')
        planes[index]['pixels_verified'] = int(raw.size)
        planes[index]['identical_pixels'] = True
    del verified
    output_hash = sha256(partial)
    if output_hash == sha256(TEMPLATE):
        raise ValueError('Recovered 6539 unexpectedly duplicates 6090')
    partial.replace(output)
    info = inspect_source(output)
    report = {'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'status': 'validated',
              'case_id': '6539', 'raw_acquisition': str(SOURCE), 'output': str(output),
              'output_sha256': output_hash, 'output_bytes': output.stat().st_size, 'planes': planes,
              'shape_cyx': info['shape_cyx'], 'pixel_size_um': info['pixel_size_um'],
              'channel_selection': 'Same 36 cycle/channel labels; reg001_ prefix added only to metadata for existing adapters',
              'pixels_modified': False,
              'original_mislabeled_stack': str(BASE / 'Images/New/PAN-6539-ND test 34.tif'),
              'annotation_status': 'Original Ignore* annotation belongs to the duplicated stack; retained at its original path and not transplanted to different pixels'}
    write_json(report_path, report)
    print(f'Reconstructed and verified {output}', flush=True)


if __name__ == '__main__':
    main()
