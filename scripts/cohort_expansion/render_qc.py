#!/usr/bin/env python3
"""Render exported cell/nucleus outlines over a full-resolution review crop."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np


def rings(geometry, origin):
    if not geometry:
        return []
    kind = geometry['type']
    coordinates = geometry.get('coordinates', [])
    if kind == 'Polygon':
        polygons = [coordinates]
    elif kind == 'MultiPolygon':
        polygons = coordinates
    elif kind == 'GeometryCollection':
        return [ring for child in geometry['geometries'] for ring in rings(child, origin)]
    else:
        raise ValueError(f'Unexpected review geometry: {kind}')
    return [np.asarray(ring)[:, :2] - origin for polygon in polygons for ring in polygon]


def render(root, case):
    folder = root / 'exports' / str(int(case)).zfill(4) / 'qc'
    metadata = json.loads((folder / 'crop.json').read_text())
    features = json.loads((folder / 'objects.geojson').read_text())['features']
    source = plt.imread(folder / 'DAPI_INS_GCG.png')
    if source.shape[:2] != (metadata['height_px'], metadata['width_px']):
        raise ValueError('Review image dimensions differ from crop metadata')
    origin = np.array([metadata['x_px'], metadata['y_px']])
    outlines = {'Cells': [], 'Nuclei': [], 'Islet core': [], 'Peri-islet': []}
    for feature in features:
        properties = feature.get('properties', {})
        if properties.get('objectType') == 'cell':
            outlines['Cells'].extend(rings(feature['geometry'], origin))
            outlines['Nuclei'].extend(rings(feature.get('nucleusGeometry'), origin))
        else:
            classification = properties.get('classification', {})
            name = classification.get('name', '') if isinstance(classification, dict) else classification
            label = {'Islet': 'Islet core', 'IsletExpanded': 'Peri-islet'}.get(name)
            if label:
                outlines[label].extend(rings(feature['geometry'], origin))
    if not outlines['Cells'] or not outlines['Nuclei']:
        raise ValueError('Actual cell and nucleus geometries are required for review')
    colors = {'Cells': '#ffff66', 'Nuclei': '#00e5ff', 'Islet core': '#ff5686', 'Peri-islet': '#ffffff'}
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    for ax in axes:
        ax.imshow(source, interpolation='nearest')
        ax.set_xlim(0, metadata['width_px'])
        ax.set_ylim(metadata['height_px'], 0)
        ax.set_axis_off()
        length = 100 / metadata['pixel_size_um']
        y = metadata['height_px'] - 40
        ax.plot([30, 30 + length], [y, y], color='white', linewidth=3)
        ax.text(30, y - 12, '100 µm', color='white', fontsize=10)
    axes[0].set_title('Source channels: DAPI / INS / GCG')
    axes[1].set_title('InstanSeg 0.1.1 cell and nucleus outlines')
    for label, lines in outlines.items():
        if lines:
            axes[1].add_collection(LineCollection(lines, colors=colors[label],
                linewidths=0.45 if label in {'Cells', 'Nuclei'} else 1.0, alpha=0.8))
    axes[1].legend(handles=[Line2D([0], [0], color=color, label=label) for label, color in colors.items()],
                   loc='upper right', fontsize=8)
    fig.suptitle(f"Case {metadata['case_id']} · {metadata['representative_islet']} · "
                 f"{metadata['cell_count']:,} cells intersect crop\nVisual review aid; biological acceptance remains pending", fontsize=12)
    destination = root / 'reports/segmentation' / f'{str(int(case)).zfill(4)}_cell_review.png'
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=160)
    plt.close(fig)
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('data/development/expanded-20260911'))
    parser.add_argument('--cases', nargs='+', required=True)
    args = parser.parse_args()
    for case in args.cases:
        print(render(args.root.resolve(), case), flush=True)
