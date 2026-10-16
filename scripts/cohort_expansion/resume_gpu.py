#!/usr/bin/env python3
"""Validate local GPU access, stop only the exact unfinished CPU cell job, and resume."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
ROOT = REPO / 'data/development/expanded-20260911'
QUPATH = '/home/smith6jt/QuPath/bin/QuPath'

parser = argparse.ArgumentParser()
parser.add_argument('--cases', nargs='+', choices=['6232', '6558', '6520'], default=['6232', '6558', '6520'])
args = parser.parse_args()
os.chdir(REPO)
env = dict(os.environ)
env.update(PYTORCH_LIBRARY_PATH='/home/smith6jt/.djl.ai/pytorch/2.7.1-cu128-linux-x86_64',
           PYTORCH_VERSION='2.7.1', PYTORCH_FLAVOR='cu128',
           JAVA_TOOL_OPTIONS='-Xmx4g -Djava.awt.headless=true -Dai.djl.pytorch.num_threads=1 -Dai.djl.pytorch.num_interop_threads=1')
probe = 'def m=ai.djl.ndarray.NDManager.newBaseManager(ai.djl.Device.gpu(0));assert m.ones(new ai.djl.ndarray.types.Shape(2,2)).sum().getFloat()==4.0f;m.close();println "GPU_PROBE_PASS"'
result = subprocess.run([QUPATH, 'script', '--cmd', probe], env=env, capture_output=True, text=True)
print(result.stdout, flush=True)
if result.returncode or 'GPU_PROBE_PASS' not in result.stdout:
    print(result.stderr, file=sys.stderr)
    raise SystemExit('GPU probe failed; any existing CPU job remains running')
expected = [str(ROOT / 'qupath/project.qpproj'), 'CX_PAN-6232-ND_test33_35p_led.tif',
            str(ROOT / 'qupath/scripts/3_Islet_Counter_v2.groovy')]
stopped = []
if not (ROOT / 'logs/6232/cells.done.json').exists():
    for process in Path('/proc').glob('[0-9]*'):
        try:
            if process.stat().st_uid != os.getuid():
                continue
            command = [s.decode() for s in (process / 'cmdline').read_bytes().split(b'\0') if s]
            if command and command[0] == QUPATH and all(value in command for value in expected):
                pid = int(process.name)
                print(f'Stopping unfinished owned CPU QuPath cell job PID {pid}', flush=True)
                os.kill(pid, signal.SIGTERM)
                stopped.append(pid)
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    deadline = time.monotonic() + 30
    while any(Path(f'/proc/{pid}').exists() for pid in stopped) and time.monotonic() < deadline:
        time.sleep(0.2)
    if any(Path(f'/proc/{pid}').exists() for pid in stopped):
        raise SystemExit('Owned CPU cell job did not terminate; GPU processing was not started')
    time.sleep(1)
subprocess.run([sys.executable, str(HERE / 'prepare_qupath.py')], check=True)
for case in args.cases:
    run_env = dict(os.environ)
    if case == '6520':
        evidence = ROOT / 'reports/channel_mapping/krt8_workflow_alias_evidence.json'
        if not evidence.is_file():
            raise SystemExit('Missing accepted KRT8 workflow-alias evidence')
        json.loads(evidence.read_text())
        run_env['ISLET_KRT8_PANEL_VERIFIED'] = '1'
    command = [sys.executable, str(HERE / 'run_qupath.py'), '--case', case, '--device', 'gpu']
    if case == '6232':
        command += ['--stages', 'cells', 'markers', 'targets', 'export', 'qc']
    print('Running case ' + case, flush=True)
    subprocess.run(command, env=run_env, check=True)
