#!/usr/bin/env python3
"""Read-only dump of the typing launch gates, to compare the workstation with HPG and before/after a code change.

Writes one JSON: run_remaining_qc_pipeline.code_inventory(), redsea_ready(), structure_generation.validate(RUN) and
run_v20.py's printed stage plans, with the site identity hook's state. Takes no lock and writes nothing else.

    python scripts/hpg/gate_dump.py --out /path/new.json
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
T30 = ROOT / 'output/typing_20261030'
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(T30))
sys.path.insert(0, str(ROOT / 'scripts/uniform'))
import site_identity as SI  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True)
    out = Path(ap.parse_args().out)
    if out.exists():
        raise FileExistsError(out)
    import run_v20 as D
    import structure_generation as SG
    P = D.P
    a = D.options()
    plans = {phase: subprocess.run([sys.executable, str(T30 / 'run_v20.py'), '--phase', phase], cwd=ROOT,
                                   capture_output=True, text=True, check=True).stdout.splitlines()
             for phase in D.PHASES}
    dump = {'site_identity': SI.active(), 'code_inventory': P.code_inventory(), 'redsea_ready': P.redsea_ready(a),
            'structure_generation': SG.validate(a.run), 'run_v20_plans': plans}
    out.write_text(json.dumps(dump, indent=1, sort_keys=True, default=str) + '\n')
    print(f"ready {dump['redsea_ready']['all20_current']}; wrote {out}")
    return 0 if dump['redsea_ready']['all20_current'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
