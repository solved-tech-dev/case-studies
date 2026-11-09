#!/usr/bin/env python3
"""Build all case programs afresh. Production bridge requires Linux/CUDA/PanguLU."""
import argparse, hashlib, json, os, platform, shlex, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def run(cmd,cwd,log):
    print('BUILD',log.name,flush=True)
    with log.open('w') as out:r=subprocess.run(cmd,cwd=cwd,stdout=out,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError(f'Build failed: {log.name} (exit {r.returncode})')
def main():
    p=argparse.ArgumentParser()
    p.add_argument('--workspace',default='_work/clean')
    p.add_argument('--reference-only',action='store_true')
    p.add_argument('--jobs',default='4')
    p.add_argument('--blas')
    p.add_argument('--lapack')
    p.add_argument('--pangulu',default=os.environ.get('PANGULU_ROOT'))
    p.add_argument('--cuda',default=os.environ.get('CUDA_HOME','/usr/local/cuda-12.9'))
    p.add_argument('--cuda-include')
    p.add_argument('--cuda-nvcc-include')
    p.add_argument('--cuda-lib')
    args=p.parse_args()
    if not args.reference_only and (platform.system()!='Linux' or platform.machine()!='x86_64' or not args.pangulu):
        p.error('PanguLU build requires Linux x86-64 and --pangulu or PANGULU_ROOT')
    w=Path(args.workspace).resolve();phg=w/'phg-0.9.7';solid=w/'solid';logs=w/'logs';logs.mkdir(exist_ok=True);out=w/'bin';out.mkdir(exist_ok=True)
    cc=os.environ.get('CC','cc');cxx=os.environ.get('CXX','c++')
    if not (phg/'Makefile.inc').exists():
        flags=['--disable-'+s for s in ['mpi','omp','shared','hypre','petsc','mumps','superlu','trilinos','zoltan','metis','parmetis','eigen-solver','arpack','jdbsym','blopex','slepc','spooles','laspack','minres','hsparse','papi','p4est','matheval','gmp','mpfr','quadmath']]
        cmd=['./configure',*flags,'CC='+cc,'CXX='+cxx,'CFLAGS=-O2 -std=gnu99 -Wno-implicit-int -Wno-implicit-function-declaration -Wno-incompatible-pointer-types','CXXFLAGS=-O2 -std=gnu++11']
        if args.blas:cmd+=['--with-blas='+args.blas]
        if args.lapack:cmd+=['--with-lapack='+args.lapack]
        if os.environ.get('FC'):cmd+=['FC='+os.environ['FC']]
        if os.environ.get('F77'):cmd+=['F77='+os.environ['F77']]
        run(cmd,phg,logs/'configure.log')
    run(['make','-j'+args.jobs,'-C','src'],phg,logs/'phg-build.log')
    def val(k):return shlex.split(subprocess.check_output(['make','-s','-f','Makefile.inc',k],cwd=phg,text=True).strip())
    cpp=val('CPPFLAGS');cf=val('CFLAGS');cxf=val('CXXFLAGS');ld=val('LDFLAGS');libs=val('LIBS')
    for name in ['interface']:
        obj=out/(name+'.o');run([cc,*cf,*cpp,'-c',str(phg/'examples'/(name+'.c')),'-o',str(obj)],w,logs/(name+'-compile.log'))
        run([cxx,*ld,'-o',str(out/name),str(obj),*libs],w,logs/(name+'-link.log'))
    objects=[]
    for src in sorted((solid/'src').glob('*.cxx')):
        obj=out/('solid-'+src.stem+'.o');objects.append(str(obj))
        run([cxx,*cxf,*cpp,'-std=gnu++11','-include',str(ROOT/'integration/phg/solid-compat.hpp'),'-I'+str(solid/'include'),'-c',str(src),'-o',str(obj)],w,logs/('solid-'+src.stem+'.log'))
    run([cxx,*ld,'-o',str(out/'solid'),*objects,*libs,'-lexpat'],w,logs/'solid-link.log')
    shared=['-dynamiclib'] if platform.system()=='Darwin' else ['-shared','-fPIC']
    run([cxx,'-std=c++11','-O2',*shared,'-I'+str(ROOT/'adapter'),str(ROOT/'tests/reference_bridge.cpp'),'-o',str(out/'libtest_reference.so')],w,logs/'reference-bridge.log')
    run([cxx,'-std=c++11','-O2','-I'+str(ROOT/'adapter'),str(ROOT/'tests/test_csc.cpp'),'-o',str(out/'test_csc')],w,logs/'unit-build.log')
    run([str(out/'test_csc')],w,logs/'unit-run.log')
    if not args.reference_only:
        pang=Path(args.pangulu).resolve();cuda=Path(args.cuda).resolve()
        inc=Path(args.cuda_include).resolve() if args.cuda_include else cuda/'include'
        lib=Path(args.cuda_lib).resolve() if args.cuda_lib else cuda/'lib64'
        extra_inc=['-I'+str(Path(args.cuda_nvcc_include).resolve())] if args.cuda_nvcc_include else []
        run([cxx,'-std=c++11','-O3','-fopenmp',*shared,'-I'+str(ROOT/'adapter'),'-I'+str(pang/'include'),'-I'+str(inc),
             *extra_inc,str(ROOT/'adapter/pangulu_bridge.cpp'),'-L'+str(pang/'lib'),'-Wl,-rpath,'+str(pang/'lib'),
             '-Wl,--no-as-needed','-lpangulu','-lpangulu_platform','-Wl,--as-needed',
             '-L'+str(lib),'-Wl,-rpath,'+str(lib),'-l:libcudart.so.12','-o',str(out/'libphg_pangulu.so')],w,logs/'pangulu-bridge.log')
        # Load vendor dependencies at process startup, as in its official example.
        # Recent glibc rejects dlopen of DSOs requiring an executable stack.
        startup=['-L'+str(out),'-Wl,-rpath,'+str(out),'-Wl,--no-as-needed','-lphg_pangulu','-Wl,--as-needed']
        for name in ['interface']:
            run([cxx,*ld,'-o',str(out/name),str(out/(name+'.o')),*libs,*startup],w,logs/(name+'-production-link.log'))
        run([cxx,*ld,'-o',str(out/'solid'),*objects,*libs,'-lexpat',*startup],w,logs/'solid-production-link.log')
    expected=json.loads((w/'case-source-manifest.json').read_text())
    for name,base in [('phg',phg),('solid',solid)]:
        for f,sha in expected[name].items():assert hashlib.sha256((base/f).read_bytes()).hexdigest()==sha,f
    (w/'build-report.json').write_text(json.dumps({'all_case_programs_compiled':True,'case_sources_unchanged':True,'reference_only':args.reference_only,'real_pangulu_bridge_compiled':not args.reference_only,'source_origin':'pinned official archives','prebuilt_objects_used':False},indent=2)+'\n')
    print('Build complete:',out)
if __name__=='__main__':main()
