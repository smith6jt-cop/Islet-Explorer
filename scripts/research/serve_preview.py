#!/usr/bin/env python3
"""Start an isolated loopback Shiny Server, preserving existing services."""
import argparse
import getpass
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
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
    p.add_argument('--python',required=True);p.add_argument('--restart',action='store_true')
    p.add_argument('--replace-preview',type=Path,help='Replace only the recorded service for this previous candidate');a=p.parse_args()
    data=Path(a.data_root).resolve(strict=True);source=Path(a.app_root).resolve(strict=True)
    metadata=data/'manifest.json' if (data/'manifest.json').exists() else data/'candidate.json'
    schema=str(json.loads(metadata.read_text())['schema_version'])
    if schema not in {'2','3'}:raise ValueError('Unsupported candidate runtime schema')
    runtime=data/'preview';runtime.mkdir(exist_ok=True)
    previous=runtime/'service.json'
    if previous.exists():
        old=json.loads(previous.read_text());procpath=Path(f'/proc/{old["pid"]}/cmdline')
        if procpath.exists() and str(runtime/'shiny-server.conf').encode() in procpath.read_bytes():
            if not a.restart: raise RuntimeError('Candidate is already running; use --restart to restart only its recorded process')
            stop_recorded_preview(data)
    app=runtime/'app'
    if app.exists(): shutil.rmtree(app)
    shutil.copytree(source,app,ignore=shutil.ignore_patterns('local_images','node_modules','.git','tests'))
    (app/'www/local_images').symlink_to(data/'images',target_is_directory=True)
    replaced=stop_recorded_preview(a.replace_preview) if a.replace_preview else None
    port=choose_port();log=runtime/'logs';log.mkdir(exist_ok=True)
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
    record=dict(url=f'http://127.0.0.1:{port}/islet-explorer/',pid=proc.pid,logs=str(log),server_log=str(runtime/'server.log'),
                app_snapshot=str(app),data_root=str(data),schema_version=schema,ufrc_accepted=False,
                replaced_preview=replaced,
                restart_argv=['python3',str(Path(__file__).resolve()),'--data-root',str(data),'--app-root',str(source),'--python',a.python,'--restart'],
                stop_argv=['kill','-TERM',str(proc.pid)])
    write_json(runtime/'service.json',record);print(json.dumps(record,indent=2))

if __name__=='__main__':main()
