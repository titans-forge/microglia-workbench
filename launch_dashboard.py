#!/usr/bin/env python3
"""Start or reuse the matching local dashboard, wait for readiness, then open it."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import os
import signal
import shlex
from urllib.request import urlopen
from dashboard_server import APP, VERSION, build_fingerprint

ROOT=Path(__file__).resolve().parent
URL='http://127.0.0.1:8879'


def probe():
    try:
        with urlopen(URL+'/health',timeout=.6) as response:return json.load(response)
    except Exception:return None


def restart_owned_dashboard():
    """Only terminate the listener whose command identifies this exact project server."""
    result=subprocess.run(['/usr/sbin/lsof','-nP','-iTCP:8879','-sTCP:LISTEN','-t'],capture_output=True,text=True)
    pids=set(result.stdout.split())
    if len(pids)!=1:raise RuntimeError('Cannot identify the dashboard listener safely. Close its window and retry.')
    pid=int(pids.pop())
    command=subprocess.run(['/bin/ps','-p',str(pid),'-o','command='],capture_output=True,text=True,check=True).stdout.strip()
    owner=subprocess.run(['/bin/ps','-p',str(pid),'-o','uid='],capture_output=True,text=True,check=True).stdout.strip()
    args=shlex.split(command)
    if owner!=str(os.getuid()) or len(args)<2 or Path(args[1]).resolve()!=ROOT/'dashboard_server.py':
        raise RuntimeError('Port 8879 process identity does not match this project; close the other instance manually.')
    os.kill(pid,signal.SIGTERM)
    for _ in range(50):
        if probe() is None:return
        time.sleep(.1)
    raise RuntimeError('The old dashboard has not stopped yet. Retry after it exits.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--no-open',action='store_true');args=parser.parse_args()
    expected=str((ROOT/'data/review').resolve())
    health=probe();process=None
    if health and (health.get('app')!=APP or health.get('workspace')!=expected):
        raise RuntimeError('Port 8879 belongs to another app or workspace. Close that instance before retrying.')
    if health and (health.get('version')!=VERSION or health.get('build')!=build_fingerprint()):
        restart_owned_dashboard();health=None
    if not health:
        logs=ROOT/'data/dashboard.log';logs.parent.mkdir(exist_ok=True)
        with logs.open('ab') as log:
            process=subprocess.Popen([sys.executable,str(ROOT/'dashboard_server.py')],cwd=ROOT,stdout=log,stderr=log,start_new_session=True)
        for _ in range(100):
            if process.poll() is not None:raise RuntimeError(f'Dashboard could not start. See {logs}')
            health=probe()
            if health and health.get('app')==APP and health.get('workspace')==expected and health.get('version')==VERSION and health.get('build')==build_fingerprint():break
            time.sleep(.1)
        else:
            process.terminate();raise RuntimeError(f'Dashboard readiness timed out. See {logs}')
    print(f'Titans Drug Discovery ready: {URL}',flush=True)
    if not args.no_open:subprocess.run(['/usr/bin/open',URL],check=True)


if __name__=='__main__':
    try:main()
    except Exception as error:print(f'Unable to open dashboard: {error}',file=sys.stderr);sys.exit(1)
