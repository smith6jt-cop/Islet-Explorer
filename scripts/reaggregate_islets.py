#!/usr/bin/env python3
"""Offline trajectory entry point for the explicit all-donor candidate.

The former scVI/INS-root aggregation pipeline is preserved in
archive/trajectory-pre-20260916/reaggregate_islets.py. For new preparation use
scripts/research/build_candidate.py; this command consumes its islet features.
--input-root, --output-root and --islet-review are mandatory. See scripts/research/README.md.
"""
import runpy
import sys
from pathlib import Path

if __name__ == '__main__':
    runner=Path(__file__).parent/'research'/'trajectory.py'
    sys.path.insert(0,str(runner.parent))
    runpy.run_path(str(runner),run_name='__main__')
