#!/usr/bin/env python3
"""Cross-site equivalence of the HPG copy (plan step E3); run on HPG inside the site container. Reads only besides --out.

pdpn       the 20 despeckled PDPN rasters written on HPG by an unmodified copy of despeckle_pdpn.py, against the sha256
           the workstation recorded (output/typing_20261030/lymphatic_despeckle/pdpn/<tag>.json)
detection  the Lymphatic objects HPG's copy of lymphatic_despeckle_run.py drew, against the workstation's, as normalized
           geometry (QuPath writes new object IDs and order; the driver compares the same way)
phases     every output recorded by HPG's run_v20 journals of the structures and typing phases, against the sha256 the
           10-30 workstation journals recorded for the same paths

    python scripts/hpg/equivalence.py --copy-dir output/hpg_equivalence_20261031 --reference-journals A B \
        --hpg-journals C D --out output/hpg_equivalence_20261031/equivalence.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
T30 = ROOT / 'output/typing_20261030'
sys.path.insert(0, str(ROOT / 'scripts/uniform'))
import structure_generation as SG  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 23), b''):
            h.update(block)
    return h.hexdigest()


def pdpn(copy_dir: Path) -> dict:
    rows = {}
    for ref in sorted((T30 / 'lymphatic_despeckle/pdpn').glob('*.json')):
        tag = ref.stem
        expected = json.loads(ref.read_text())['sha256']
        produced = copy_dir / 'lymphatic_despeckle/pdpn' / f'{tag}.u16'
        rows[tag] = {'equal': produced.is_file() and sha256(produced) == expected}
    return rows


def detection(copy_dir: Path) -> dict:
    rows = {}
    for mode in ('despeckled', 'original'):
        for produced in sorted((copy_dir / 'lymphatic_despeckle/runs' / mode).glob('*/Lymphatic.geojson')):
            tag = produced.parent.name
            reference = (T30 / 'lymphatic_despeckle/runs' / mode / tag / 'Lymphatic.geojson')
            a, b = SG._normalized_polygons(produced), SG._normalized_polygons(reference)
            rows[f'{mode}/{tag}'] = {'objects_hpg': len(a), 'objects_workstation': len(b), 'equal_geometry': a == b}
    return rows


def phases(reference: list[Path], hpg: list[Path]) -> dict:
    rows = {}
    for ref, new in zip(reference, hpg):
        r, n = json.loads(ref.read_text()), json.loads(new.read_text())
        if r['phase'] != n['phase'] or n.get('status') != 'complete':
            raise ValueError(f'{new}: not a complete journal of phase {r["phase"]}')
        ro, no = r['outputs'], n['outputs']
        rows[r['phase']] = {'reference_journal': str(ref), 'hpg_journal': str(new), 'outputs': len(no),
                            'same_paths': set(ro) == set(no),
                            'identical': sorted(p for p in no if ro.get(p, {}).get('sha256') == no[p]['sha256']),
                            'different': sorted(p for p in no if p in ro and ro[p]['sha256'] != no[p]['sha256']),
                            'only_hpg': sorted(set(no) - set(ro)), 'only_workstation': sorted(set(ro) - set(no))}
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--copy-dir', required=True)
    ap.add_argument('--reference-journals', nargs=2, required=True)
    ap.add_argument('--hpg-journals', nargs=2, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        raise FileExistsError(out)
    copy_dir = Path(a.copy_dir).resolve()
    result = {'pdpn': pdpn(copy_dir), 'detection': detection(copy_dir),
              'phases': phases([Path(p) for p in a.reference_journals], [Path(p) for p in a.hpg_journals])}
    out.write_text(json.dumps(result, indent=1) + '\n')
    p = result['pdpn']; d = result['detection']
    print(f"pdpn equal {sum(r['equal'] for r in p.values())}/{len(p)}; "
          f"detection equal {sum(r['equal_geometry'] for r in d.values())}/{len(d)}")
    for phase, r in result['phases'].items():
        print(f"{phase}: {len(r['identical'])}/{r['outputs']} outputs identical, {len(r['different'])} different, "
              f"same paths {r['same_paths']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
