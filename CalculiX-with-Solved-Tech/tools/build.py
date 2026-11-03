# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Build CalculiX2.23 + PanguLU from source; never download PanguLU or a license."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pangulu-root',type=Path,required=True,help='User-provided external SDK directory')
    p.add_argument('--cuda-root',type=Path,required=True,help='User-provided CUDA include/ and lib64/ or lib/')
    p.add_argument('--spooles-root',type=Path,required=True,help='External SPOOLES2.2 headers and spooles.a for original MPC helpers')
    p.add_argument('--cuda-include',type=Path,action='append',default=[],help='Extra CUDA include directory, e.g. split nvcc headers')
    p.add_argument('--cuda-lib',type=Path,help='Override CUDA library directory')
    p.add_argument('--blas-libs',default='-lopenblas',help='User-provided BLAS/LAPACK linker flags; LP64 ABI')
    p.add_argument('--arpack-libs',default='-larpack',help='User-provided sequential ARPACK linker flags')
    p.add_argument('--workspace',type=Path,default=ROOT/'_work/build')
    p.add_argument('--jobs',type=int,default=min(8,os.cpu_count() or 1))
    p.add_argument('--cc',default=os.environ.get('CC','gcc'))
    p.add_argument('--cxx',default=os.environ.get('CXX','g++'))
    p.add_argument('--fc',default=os.environ.get('FC','gfortran'))
    p.add_argument('--timing',action='store_true',help='Apply the supplied optional Newton timing patch')
    p.add_argument('--resume',action='store_true',help='Resume only an identical source/flag workspace')
    a=p.parse_args();assert 1<=a.jobs<=128
    pg=a.pangulu_root.resolve(strict=True);cuda=a.cuda_root.resolve(strict=True);spooles=a.spooles_root.resolve(strict=True)
    for file in (pg/'include/pangulu.h',pg/'lib/libpangulu.so',pg/'lib/libpangulu_platform.so',cuda/'include/cuda_runtime_api.h'):
        if not file.is_file():p.error('Missing external dependency: '+str(file))
    for file in (spooles/'misc.h',spooles/'spooles.a'):
        if not file.is_file():p.error('Missing external SPOOLES dependency: '+str(file))
    cuda_lib=a.cuda_lib.resolve(strict=True) if a.cuda_lib else next((cuda/x for x in ('lib64','lib') if (cuda/x).is_dir()),None)
    if cuda_lib is None:p.error('CUDA library directory not found')
    # Do not copy or inspect license contents. The library validates runtime authorization.
    workspace=a.workspace.resolve();config=dict(pangulu_root=str(pg),cuda_root=str(cuda),cuda_include=list(map(str,a.cuda_include)),
        cuda_lib=str(cuda_lib),spooles_root=str(spooles),blas=a.blas_libs,arpack=a.arpack_libs,cc=a.cc,cxx=a.cxx,fc=a.fc,timing=a.timing)
    if workspace.exists():
        if not a.resume:p.error('Workspace exists; choose a new directory, or --resume for identical settings')
        if json.loads((workspace/'build-config.json').read_text())!=config:p.error('Cannot resume with different dependencies/flags')
    else:
        workspace.mkdir(parents=True)
        (workspace/'build-config.json').write_text(json.dumps(config,indent=2)+'\n')
    src=workspace/'CalculiX/ccx_2.23/src'
    if not src.exists():
        with tarfile.open(ROOT/'provenance/upstream/ccx_2.23.src.tar.bz2','r:bz2') as t:
            for m in t.getmembers():
                name=Path(m.name)
                if name.is_absolute() or '..' in name.parts or m.issym() or m.islnk():raise ValueError('Unsafe source archive member')
            t.extractall(workspace,filter='data')
        shutil.copyfile(ROOT/'adapter/spooles_entry.h',src/'spooles.h')
        if a.timing:subprocess.run(['patch','-p1','-i',str(ROOT/'patches/nonlingeo-timing.patch')],cwd=src,check=True)
    inc='\n'.join(line.split('#',1)[0] for line in (src/'Makefile.inc').read_text().splitlines())
    def sources(name):
        match=re.search(r'^'+name+r'\s*=\s*(.*?)(?=^\w+\s*=|\Z)',inc,re.M|re.S)
        return re.findall(r'[\w.]+\.(?:f|c|cpp)\b',match[1]) if match else []
    fs=[x for x in sources('SCCXF') if x!='dgesv.f']
    # The fixed 2.23 archive indexes mafillmm.c, but contains only mafillmm.f.
    # Keep the real Fortran implementation; do not create a placeholder C file.
    cs=[x for x in sources('SCCXC') if x not in ('spooles.c','mafillmm.c')]+['ccx_2.23.c']
    # Like the archive's Makefile, do not use SCCXCXX: its umat_dl.cpp is not
    # included in this source distribution. external.c remains compiled.
    cxx=[]
    assert fs and cs and len(set(fs+cs+cxx))==len(fs+cs+cxx)
    obj=workspace/'obj';logs=workspace/'logs';obj.mkdir(exist_ok=True);logs.mkdir(exist_ok=True)
    defines=['-DSPOOLES','-DARPACK','-DMATRIXSTORAGE','-DNETWORKOUT','-DARCH=Linux']
    tasks=[]
    for filename in fs+cs+cxx:
        source=src/filename;assert source.is_file(),source
        out=obj/(filename+'.o')
        flags=['-O2','-fopenmp','-I'+str(src)]
        if filename.endswith('.f'):compiler=a.fc;flags+=['-cpp','-std=legacy','-fallow-argument-mismatch']
        elif filename.endswith('.cpp'):compiler=a.cxx;flags+=['-std=c++17']+defines
        else:compiler=a.cc;flags+=['-std=gnu99','-I'+str(spooles)]+defines
        tasks.append((source,out,[compiler]+flags+['-c',str(source),'-o',str(out)]))
    adapter=ROOT/'adapter/ccx_pangulu.cpp';ao=obj/'ccx_pangulu.cpp.o'
    tasks.append((adapter,ao,[a.cxx,'-O2','-std=c++17','-fopenmp','-I'+str(pg/'include'),'-I'+str(cuda/'include')]+
        ['-I'+str(i.resolve()) for i in a.cuda_include]+['-c',str(adapter),'-o',str(ao)]))
    header_mtime=max((src/'CalculiX.h').stat().st_mtime,(src/'spooles.h').stat().st_mtime,Path(__file__).stat().st_mtime)
    def compile_one(task):
        source,out,cmd=task
        if a.resume and out.exists() and out.stat().st_mtime>max(source.stat().st_mtime,header_mtime):return
        log=logs/(source.name+'.log')
        with log.open('w') as f:
            f.write(shlex.join(cmd)+'\n');f.flush()
            cp=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
        if cp.returncode:raise RuntimeError('Compilation failed: '+source.name+'; see '+str(log))
    print('Compiling',len(tasks),'source files; no prebuilt application objects',flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for i,future in enumerate(as_completed([pool.submit(compile_one,t) for t in tasks]),1):
            future.result()
            if i%100==0:print('Compiled',i,'/',len(tasks),flush=True)
    binary=workspace/'ccx_pangulu'
    main_object=obj/'ccx_2.23.c.o'
    archive=workspace/'ccx-source.a'
    # Match upstream linkage: main + newly compiled source archive. Unused
    # alternative driver objects must not be forced into the executable.
    subprocess.run(['ar','rcs',str(archive)]+[str(t[1]) for t in tasks if t[1] not in (main_object,ao)],check=True)
    libs=['-L'+str(pg/'lib'),'-lpangulu','-lpangulu_platform','-L'+str(cuda_lib),'-l:libcudart.so.12']
    libs+=[str(spooles/'spooles.a')]+shlex.split(a.arpack_libs)+shlex.split(a.blas_libs)+['-lz','-lstdc++','-lpthread','-ldl','-lrt','-lm']
    libs+=['-Wl,-rpath,'+str(pg/'lib'),'-Wl,-rpath,'+str(cuda_lib)]
    cmd=[a.fc,'-O2','-fopenmp',str(main_object),str(ao),str(archive)]+libs+['-o',str(binary)]
    with (workspace/'link.log').open('w') as f:
        f.write(shlex.join(cmd)+'\n');f.flush();cp=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    if cp.returncode:raise RuntimeError('Link failed; see '+str(workspace/'link.log'))
    report=dict(status='compiled',compiled_source_files=len(tasks),prebuilt_application_objects_used=False,
        old_test_stubs_used=False,original_spooles_implementation_linked=False,external_ARPACK_enabled=True,
        original_MPC_helpers='external open-source SPOOLES2.2; CCX spooles.c linear-solver wrapper is not compiled',
        source_changes=['build-tree declaration-only spooles.h']+(['optional nonlingeo timing patch'] if a.timing else []),
        internal_dgesv_omitted='use user-provided external BLAS/LAPACK to avoid symbol interposition',
        absent_upstream_index_entries_omitted=['mafillmm.c (mafillmm.f is compiled)','unused SCCXCXX umat_dl.cpp'],
        binary=str(binary),runtime_validation='not implied by successful compilation')
    (workspace/'build-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Built',binary,flush=True)

if __name__=='__main__':main()
