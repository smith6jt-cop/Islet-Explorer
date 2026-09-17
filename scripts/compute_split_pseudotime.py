#!/usr/bin/env python3
"""Compatibility CLI: size is now a display filter, not a separate fitted axis.

Pass --input-root, --output-root and --islet-review explicitly; see research/trajectory.py.
Historical implementation is preserved in archive/trajectory-pre-20260916.
"""
import runpy
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'research'))
runpy.run_path(str(Path(__file__).parent / 'research/trajectory.py'), run_name='__main__')
