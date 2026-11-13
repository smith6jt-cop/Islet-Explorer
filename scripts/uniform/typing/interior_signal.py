"""Compare competing marker sets at the same native cell-interior pixels."""
from __future__ import annotations

import json

import numpy as np
from scipy.ndimage import maximum_filter, minimum_filter

import area_coverage as AC
import codex_typing as CT
import redsea_cells as RS


def compare_signals(left, right):
    """Return exclusive pixel counts and median brightness for two calibrated groups.

    Input matrices are markers by interior pixels, on the uncapped scale
    (native signal - marker level) / marker span. A group is detectable at a
    pixel if any of its active, group-specific markers clears its level.
    """
    left, right = np.asarray(left, float), np.asarray(right, float)
    if (left.ndim != 2 or right.ndim != 2 or left.shape[1] != right.shape[1]
            or not left.size or not right.size or not np.isfinite(left).all()
            or not np.isfinite(right).all()):
        return None
    a, b = left.mean(axis=0), right.mean(axis=0)
    a_on, b_on = (left >= 0).any(axis=0), (right >= 0).any(axis=0)
    a_wins = a_on & (~b_on | (a > b))
    b_wins = b_on & (~a_on | (b > a))
    return a.size, int(a_wins.sum()), int(b_wins.sum()), float(np.median(a)), float(np.median(b))


class InteriorSignal:
    """Read the native channels for eligible two-way conflicts on demand per image."""

    def __init__(self, run, tag, reg, object_ids, counts, markers, levels, tops):
        self.run, self.tag, self.reg = run, tag, reg
        self.object_ids = np.asarray(object_ids)
        self.counts = np.asarray(counts)
        self.markers, self.levels, self.spans = markers, levels, tops - levels
        self.owners = None
        self.image_path = None
        self.positions = {}
        self.samples = {}
        self.evidence = {}

    def _index_pixels(self, rows):
        missing = sorted(set(rows) - self.positions.keys())
        if not missing:
            return
        if self.owners is None:
            entry = CT.manifest(self.run)[self.tag]
            self.image_path = entry['image_path']
            self.channels = CT.channels(self.run, entry['entry_id'])
            img = RS.Channels(self.image_path)
            try:
                height, width = img.shape[1:]
            finally:
                img.close()
            features = json.loads((self.run / 'cell_geojson' / f'{self.tag}.geojson').read_text())
            if isinstance(features, dict):
                features = features['features']
            mask, nucleus, ids, _ = RS.rasterize(features, (height, width))
            del nucleus, features
            self.labels = {oid: i for i, oid in enumerate(ids, 1)}
            if len(ids) != len(self.labels):
                raise ValueError(f'{self.tag}: duplicate native cell geometry IDs')
            kernel = 2 * AC.INTERIOR_PX + 1
            interior = (minimum_filter(mask, size=kernel) == mask) & (maximum_filter(mask, size=kernel) == mask)
            mask[~interior] = 0
            self.owners = mask.ravel()
        labels = [self.labels[self.object_ids[i]] for i in missing]
        index = np.flatnonzero(np.isin(self.owners, labels, kind='table'))
        index = index[np.argsort(self.owners[index], kind='stable')]
        owner = self.owners[index]
        for i, label in zip(missing, labels):
            start, stop = np.searchsorted(owner, [label, label + 1])
            if stop - start != self.counts[i]:
                raise ValueError(f'{self.tag}/{self.object_ids[i]}: native interior count differs from coverage')
            self.positions[i] = index[start:stop]

    def measure(self, requests, stage):
        """Requests are (source row, left marker indices, right marker indices)."""
        active = []
        needed = {}
        for row, left, right in requests:
            shared = set(left) & set(right)
            left = tuple(m for m in left if m not in shared)
            right = tuple(m for m in right if m not in shared)
            cols = list(left + right)
            if not left or not right or self.counts[row] == 0:
                continue
            if not (np.isfinite(self.spans[cols]) & (self.spans[cols] > 0)
                    & np.isfinite(self.levels[cols])).all():
                continue
            active.append((row, left, right))
            for m in cols:
                if (row, m) not in self.samples:
                    needed.setdefault(m, set()).add(row)
        if needed:
            self._index_pixels(row for row, _, _ in active)
            img = RS.Channels(self.image_path)
            try:
                for m, rows in needed.items():
                    name = self.markers[m]
                    channel = self.reg.table.loc[name, 'channel']
                    plane = RS.despeckle(img.native(self.channels.index(channel if channel in self.channels else name)))
                    for row in rows:
                        raw = plane.ravel()[self.positions[row]].astype(float)
                        self.samples[row, m] = (raw - self.levels[m]) / self.spans[m]
            finally:
                img.close()
        result = {}
        for row, left, right in active:
            measured = compare_signals([self.samples[row, m] for m in left],
                                       [self.samples[row, m] for m in right])
            if measured is None:
                continue
            result[row] = measured
            n, a, b, ma, mb = measured
            self.evidence.setdefault(row, {})[stage] = {
                'markers': [[self.markers[m] for m in left], [self.markers[m] for m in right]],
                'interior_pixels': n, 'winner_pixels': [a, b], 'median_scaled': [ma, mb]}
        return result
