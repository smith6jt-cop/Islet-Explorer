#!/usr/bin/env python3
"""Run isolated QuPath stages with durable logs and explicit completion markers."""
import argparse
from contextlib import contextmanager
import datetime
import fcntl
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys

CASES = {'6232': 'CX_PAN-6232-ND_test33_35p_led.tif', '6558': 'PAN-6558-AAb_test35_35p-LED.tif',
         '6520': 'CX_PANT6520-01A_T1D_Test22_34p-16LED50.tif', '0078': 'CX_PANT-HPAP-078-T1D_test28_34p_16LED50.tif',
         '6090': '6090AAbtest31_35p_led.tif', '6539': 'PAN-6539-ND test 34.tif'}


@contextmanager
def project_lock(root):
    """QuPath saves the shared project record even for a single-image run."""
    with (root / '.processing.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f'Another processing runner holds the project lock: {root}') from error
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def main():
    parser = argparse.ArgumentParser()
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--case', choices=CASES)
    selection.add_argument('--cases', choices=CASES, nargs='+', help='Run cases sequentially under one project lock')
    parser.add_argument('--root', type=Path, default=Path('data/development/expanded-20260911'))
    parser.add_argument('--device', choices=['cpu', 'gpu'], default='cpu')
    parser.add_argument('--stages', nargs='+', default=['normalize', 'tissue', 'sam_prepare', 'sam', 'sam_import', 'cells', 'markers', 'targets', 'export', 'qc'])
    args = parser.parse_args()
    with project_lock(args.root.resolve()):
        for case in args.cases or [args.case]:
            args.case = case
            run_stages(args)


def run_stages(args):
    root = args.root.resolve()
    if args.case in {'6090', '6539'}:
        recovered = json.loads((root / 'recovered_project_provenance.json').read_text())
        if recovered.get('status') != 'source_reconciled' or args.case not in recovered.get('cases', {}):
            raise SystemExit('Recovered cases require a reconciled source project')
    scripts = Path(__file__).resolve().parent
    project = root / 'qupath/project.qpproj'
    exports = root / 'exports' / args.case
    logs = root / 'logs' / args.case
    exports.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env['ISLET_EXPORT_DIR'] = str(exports)
    env['ISLET_CASE_ID'] = args.case
    preferences = root / 'runtime' / 'preferences'
    preferences.mkdir(parents=True, exist_ok=True)
    copied_preferences = preferences / '.java' / '.userPrefs'
    if not (preferences / 'qupath_preferences_initialized').exists():
        original_preferences = Path('/home/smith6jt/.java/.userPrefs')
        for original in original_preferences.glob('*/prefs.xml'):
            if original.parent.name.startswith('_!'):
                target = copied_preferences / original.relative_to(original_preferences)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(original, target)
        for original in original_preferences.glob('*/*/prefs.xml'):
            if 'QuPath/v0.7' in original.read_text():
                target = copied_preferences / original.relative_to(original_preferences)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(original, target)
        (preferences / 'qupath_preferences_initialized').write_text('Copied existing QuPath 0.7 configuration into isolated preferences.\n')
    env['JAVA_TOOL_OPTIONS'] = f'-Xmx64g -Djava.awt.headless=true -Djava.util.concurrent.ForkJoinPool.common.parallelism=32 -Djava.util.prefs.userRoot={preferences} -Djava.util.prefs.systemRoot={preferences}'
    env['OMP_NUM_THREADS'] = '32'
    env['MKL_NUM_THREADS'] = '32'
    env.setdefault('PYTORCH_VERSION', '2.5.1')
    env.setdefault('PYTORCH_FLAVOR', 'cpu')
    env.setdefault('PYTORCH_LIBRARY_PATH', '/home/smith6jt/.djl.ai/pytorch/2.5.1-cpu-linux-x86_64')
    env['ISLET_INSTANSEG_DEVICE'] = args.device
    env['ISLET_INSTANSEG_THREADS'] = '2' if args.device == 'gpu' else '32'
    if args.device == 'gpu':
        env['PYTORCH_VERSION'] = '2.7.1'
        env['PYTORCH_FLAVOR'] = 'cu128'
        env['PYTORCH_LIBRARY_PATH'] = '/home/smith6jt/.djl.ai/pytorch/2.7.1-cu128-linux-x86_64'
    stages = {'normalize': scripts / 'normalize_channels.groovy',
              'tissue': root / 'qupath/scripts/1_Panc_Objects.groovy',
              'sam_prepare': scripts / 'prepare_sam.groovy',
              'sam_import': scripts / 'import_sam.groovy',
              'cells': root / 'qupath/scripts/3_Islet_Counter_v2.groovy',
              'markers': root / 'qupath/scripts/4_Islet_Cells.groovy',
              'targets': root / 'qupath/scripts/5_Islet_Objects.groovy',
              'export': scripts / 'export_measurements.groovy',
              'qc': scripts / 'export_qc_crop.groovy'}
    if args.case in {'6090', '6539'}:
        stages['normalize'] = scripts / 'normalize_recovered.groovy'
    for stage in args.stages:
        done = logs / f'{stage}.done.json'
        if stage == 'sam':
            command = ['/home/smith6jt/miniconda3/envs/samapi/bin/python', str(scripts / 'predict_sam.py'), str(exports)]
        else:
            command = ['/home/smith6jt/QuPath/bin/QuPath', 'script', '--project', str(project),
                       '--image', CASES[args.case], '--save', str(stages[stage])]
        script = scripts / 'predict_sam.py' if stage == 'sam' else stages[stage]
        script_hash = hashlib.sha256(script.read_bytes()).hexdigest()
        if done.exists():
            previous = json.loads(done.read_text())
            if previous.get('script_sha256') != script_hash:
                raise SystemExit(f'{stage} script changed after completion; explicitly review checkpoint before rerunning')
            print(f'{args.case} {stage}: already complete', flush=True)
            continue
        stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
        print(f'{stamp} {args.case} {stage}: starting', flush=True)
        logfile = logs / f'{stage}.log'
        if logfile.exists():
            logfile.rename(logs / f'{stage}.{datetime.datetime.now().strftime("%Y%m%dT%H%M%S")}.log')
        with logfile.open('w') as log:
            log.write(f'\nSTART {stamp}\n')
            log.flush()
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        text = logfile.read_text(errors='replace')
        if result.returncode or 'Script failed' in text or '[ERROR]' in text or 'Exception in' in text:
            print(f'{args.case} {stage}: failed; inspect {logfile}', file=sys.stderr, flush=True)
            raise SystemExit(result.returncode or 1)
        done.write_text(json.dumps({'case': args.case, 'stage': stage, 'started_utc': stamp,
                                    'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                    'command': command, 'script_sha256': script_hash,
                                    'device': args.device, 'tile_workers': int(env['ISLET_INSTANSEG_THREADS'])}) + '\n')
        print(f'{args.case} {stage}: complete', flush=True)


if __name__ == '__main__':
    main()
