#!/usr/bin/env python3
"""Start an isolated loopback Shiny Server, preserving existing services."""
import argparse
import getpass
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
from audit_original_parity import audit as audit_original_parity


def require_original_parity(baseline_root, candidate_root):
    """Recheck source values before creating files or stopping any service."""
    report = audit_original_parity(baseline_root, candidate_root)
    if not report['passed']:
        raise ValueError(
            'Original-donor parity failed; candidate launch is blocked. '
            'Run audit_original_parity.py for the full read-only report. '
            'No preview files or services have been changed.')
    return report


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')

def choose_port(start=3838):
    for port in range(start,65536):
        with socket.socket() as s:
            s.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
            try: s.bind(('127.0.0.1',port))
            except OSError: continue
            return port
    raise RuntimeError('No free loopback port')

def stop_recorded_preview(data_root):
    """Stop only a process matching this candidate's exact recorded config."""
    runtime=Path(data_root).resolve(strict=True)/'preview'
    record=runtime/'service.json'
    if not record.exists():raise ValueError('No recorded preview to replace')
    old=json.loads(record.read_text());procpath=Path(f'/proc/{old["pid"]}/cmdline')
    if not procpath.exists():return old
    expected=str(runtime/'shiny-server.conf').encode()
    if expected not in procpath.read_bytes().split(b'\0'):
        raise ValueError('Recorded PID no longer belongs to this preview; leaving it untouched')
    os.kill(old['pid'],15)
    import time
    for _ in range(100):
        try:
            if not procpath.read_bytes():break
        except FileNotFoundError:
            break
        time.sleep(.05)
    return old

def main():
    p=argparse.ArgumentParser();p.add_argument('--data-root',required=True);p.add_argument('--app-root',required=True)
    p.add_argument('--baseline-root',required=True,type=Path,
                   help='Restored original data root; exact preservation is checked before launch')
    p.add_argument('--python',required=True);p.add_argument('--restart',action='store_true')
    p.add_argument('--port',type=int,help='Fixed loopback port for the preview; omit to choose the first free port.')
    p.add_argument('--replace-preview',type=Path,help='Replace only the recorded service for this previous candidate');a=p.parse_args()
    data=Path(a.data_root).resolve(strict=True);source=Path(a.app_root).resolve(strict=True)
    metadata=data/'manifest.json' if (data/'manifest.json').exists() else data/'candidate.json'
    schema=str(json.loads(metadata.read_text())['schema_version'])
    if schema not in {'2','3','4'}:raise ValueError('Unsupported candidate runtime schema')
    require_original_parity(a.baseline_root, data)
    runtime=data/'preview';runtime.mkdir(exist_ok=True)
    app=runtime/'app'
    if source==app or source.is_relative_to(app):
        raise ValueError('The app source must be outside the preview destination; copy a rollback snapshot to a separate source directory first')
    previous=runtime/'service.json'
    previous_record={}
    if previous.exists():
        old=json.loads(previous.read_text());previous_record=old;procpath=Path(f'/proc/{old["pid"]}/cmdline')
        if procpath.exists() and str(runtime/'shiny-server.conf').encode() in procpath.read_bytes():
            if not a.restart: raise RuntimeError('Candidate is already running; use --restart to restart only its recorded process')
            stop_recorded_preview(data)
    if app.exists(): shutil.rmtree(app)
    shutil.copytree(source,app,ignore=shutil.ignore_patterns('local_images','node_modules','.git','tests'))
    (app/'www/local_images').symlink_to(data/'images',target_is_directory=True)
    replaced=(stop_recorded_preview(a.replace_preview) if a.replace_preview
              else previous_record.get('replaced_preview'))
    port=a.port if a.port is not None else choose_port()
    if not (1 <= port <= 65535): raise ValueError('Port must be between 1 and 65535')
    # Fail before replacing the service record if an explicitly requested port
    # is occupied by an unrelated process.
    with socket.socket() as probe:
        probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try: probe.bind(('127.0.0.1', port))
        except OSError as exc: raise RuntimeError(f'Fixed preview port {port} is not free') from exc
    log=runtime/'logs';log.mkdir(exist_ok=True)
    for sub in ['sockets','bookmarks']:(runtime/sub).mkdir(exist_ok=True)
    config=f'run_as {getpass.getuser()};\npreserve_logs true;\nsanitize_errors true;\nsocket_dir {runtime}/sockets;\nbookmark_state_dir {runtime}/bookmarks;\nserver {{\n listen {port} 127.0.0.1;\n'
    for alias in ['local_images','images','img_data']:
        config+=f' location /islet-explorer/{alias}/ {{ site_dir {data}/images; log_dir {log}; directory_index off; }}\n'
    config+=f' location /islet-explorer/ {{ app_dir {app}; log_dir {log}; app_init_timeout 240; app_idle_timeout 0; directory_index off; }}\n}}\n'
    conf=runtime/'shiny-server.conf';conf.write_text(config)
    env=os.environ.copy()
    for key in ['KEY','OPENAI_API_KEY']:env.pop(key,None)
    env.update(ISLET_DATA_ROOT=str(data),LOCAL_IMAGE_ROOT=str(data/'images'),ISLET_RUNTIME_MODE='development',
               ISLET_DATA_SCHEMA=schema,RETICULATE_PYTHON=a.python,R_ENVIRON_USER='/dev/null',
               ISLET_ENABLE_AI='FALSE',SHOW_STATISTICS_TAB='0',SHOW_INFLAMMATION_TAB='0',
               DEBUG_CREDENTIALS='0',VIEWER_DEBUG='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    # Shiny Server resets worker environments. The isolated application snapshot
    # receives an explicit, non-secret profile before its ordinary R sources.
    allowed=['ISLET_DATA_ROOT','LOCAL_IMAGE_ROOT','ISLET_RUNTIME_MODE','ISLET_DATA_SCHEMA','RETICULATE_PYTHON',
             'ISLET_ENABLE_AI','SHOW_STATISTICS_TAB','SHOW_INFLAMMATION_TAB','DEBUG_CREDENTIALS','VIEWER_DEBUG',
             'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']
    profile='Sys.setenv('+',\n'.join(k+' = '+json.dumps(env[k]) for k in allowed)+')\n'
    (app/'R/00_00_candidate_profile.R').write_text(profile)
    with (runtime/'server.log').open('a') as stream:
        proc=subprocess.Popen(['/opt/shiny-server/bin/shiny-server',str(conf)],env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
    import time
    for _ in range(100):
        try:
            with socket.create_connection(('127.0.0.1',port),timeout=.2):break
        except OSError:
            if proc.poll() is not None: raise RuntimeError('Preview server exited; inspect server.log')
            time.sleep(.1)
    else: raise RuntimeError('Preview did not open its loopback listener; inspect server.log')
    restart_argv=[str(Path(sys.executable).resolve()),str(Path(__file__).resolve()),
                  '--data-root',str(data),'--app-root',str(source),'--python',a.python,
                  '--baseline-root',str(a.baseline_root.resolve(strict=True)),
                  '--port',str(port),'--restart']
    stop_argv=['kill','-TERM',str(proc.pid)]
    record=dict(url=f'http://127.0.0.1:{port}/islet-explorer/',pid=proc.pid,logs=str(log),server_log=str(runtime/'server.log'),
                app_snapshot=str(app),data_root=str(data),schema_version=schema,ufrc_accepted=False,
                replaced_preview=replaced,
                fixed_port=port,restart_argv=restart_argv,stop_argv=stop_argv,
                restart_command=shlex.join(restart_argv),stop_command=shlex.join(stop_argv))
    # Network exposure is configured separately from this loopback service.
    # Preserve its operator record across an app-only restart.
    network_record=previous_record or (replaced if isinstance(replaced,dict) else {})
    for key in ['public_url','network_proxy','network_config','network_restart_argv']:
        if key in network_record:record[key]=network_record[key]
    write_json(runtime/'service.json',record);print(json.dumps(record,indent=2))

if __name__=='__main__':main()
