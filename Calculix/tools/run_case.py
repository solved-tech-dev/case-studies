# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Run ONE fixed input with an existing PanguLU CCX binary; no scheduler submit."""
import argparse
import gzip
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('case',choices=['dam','ax6ht','capacitor','oneel20fi','pan'])
parser.add_argument('--binary',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--profile',choices=['fixed','small'],default='fixed',help='Fixed tested mesh, or original upstream-sized integration check')
parser.add_argument('--reference-frd',type=Path,help='Optional independently supplied reference for this exact mesh')
args=parser.parse_args()
binary=args.binary.resolve(strict=True)
if not binary.is_file() or not os.access(binary,os.X_OK):parser.error('binary must be an executable file')
if not os.environ.get('PANGULU_LICENSE_PATHS'):parser.error('set your authorized PANGULU_LICENSE_PATHS first; no key is bundled')
if shutil.disk_usage(args.output.parent.resolve(strict=True)).free<20*1024**3:parser.error('at least20GiB free required for full application/evidence outputs')
out=args.output.resolve();out.mkdir()  # New run directory only; never overwrite.
case=root/'cases'/args.case
params=json.loads((case/'parameters.json').read_text())
if args.profile=='small':
    from prepare_small import prepare
    prepare(args.case,out/'model.inp')
    params['measured_system'].update(n=None,stored_full_nnz=None)
else:
    with gzip.open(case/'input/model.inp.gz','rb') as source,(out/'model.inp').open('xb') as dest:
        shutil.copyfileobj(source,dest)
params['run_profile']=args.profile
(out/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
env=os.environ.copy()
env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='4',
    CCX_NPROC_EQUATION_SOLVER='1',CCX_NPROC_STIFFNESS='4',CCX_NPROC_RESULTS='4',
    HULL_CANONICAL=str(out/'canonical.csr.gz'),HULL_CREATE_CANONICAL='1',
    HULL_ACTUAL_PREFIX=str(out/'actual'))
env.pop('HULL_CONDITION_DIAG',None)
started=time.monotonic();timed_out=False
with (out/'application.stdout').open('x') as stdout,(out/'application.stderr').open('x') as stderr:
    command=[str(binary),'-i','model']
    if Path('/usr/bin/time').is_file():command=['/usr/bin/time','-v','-o',str(out/'resources.txt')]+command
    process=subprocess.Popen(command,cwd=out,env=env,
        stdout=stdout,stderr=stderr,start_new_session=True)
    try:rc=process.wait(timeout=240 if args.profile=='small' else 600)
    except subprocess.TimeoutExpired:
        timed_out=True;os.killpg(process.pid,signal.SIGTERM)
        try:process.wait(timeout=10)
        except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
        rc=124
with (out/'application-exit.json').open('x') as f:
    json.dump(dict(exit_code=rc,timeout=timed_out,whole_application_wall_s=time.monotonic()-started,
        binary=str(binary),case=args.case,profile=args.profile),f,indent=2);f.write('\n')
check=subprocess.run([sys.executable,str(root/'tools/check_log.py'),str(out),str(out/'parameters.json')])
physical=None
if rc==0:
    command=[sys.executable,str(root/'tools/validate_fields.py'),args.case,str(out)]
    if args.reference_frd:command+=['--reference-frd',str(args.reference_frd.resolve())]
    physical=subprocess.run(command)
qualified=rc==0 and check.returncode==0 and physical is not None and physical.returncode==0
(out/'validation.json').write_text(json.dumps(dict(status='passed' if qualified else 'failed',
    profile=args.profile,log_checks_passed=check.returncode==0,physical_checks_passed=physical is not None and physical.returncode==0,
    qualification='integration_and_physical_correctness',
    reference_field_comparison=args.reference_frd is not None),indent=2)+'\n')
print('Run preserved:',out)
raise SystemExit(0 if qualified else 1)
