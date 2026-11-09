#!/usr/bin/env python3
"""Fetch pinned official archives, create a fresh build tree, apply library-only integration."""
import argparse, difflib, hashlib, json, shutil, sys, tarfile, tempfile, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(base, patterns):
    return {str(p.relative_to(base)):digest(p) for pattern in patterns for p in sorted(base.glob(pattern)) if p.is_file()}
def main():
    if sys.version_info < (3,12):raise SystemExit('Python 3.12+ is required; select a newer interpreter before preparing sources.')
    a=argparse.ArgumentParser();a.add_argument('--workspace',default='_work/clean');args=a.parse_args()
    work=Path(args.workspace).resolve();work.mkdir(parents=True,exist_ok=False)
    downloads=ROOT/'_work/downloads';downloads.mkdir(parents=True,exist_ok=True)
    lock=json.loads((ROOT/'upstream.lock.json').read_text())
    for meta in lock.values():
        archive=downloads/meta['archive']
        if not archive.exists():
            with tempfile.TemporaryDirectory(dir=downloads) as scratch:
                pending=Path(scratch)/archive.name
                with urllib.request.urlopen(meta['url'],timeout=45) as r,pending.open('xb') as out:
                    shutil.copyfileobj(r,out)
                if digest(pending)!=meta['sha256']:
                    raise RuntimeError('Official archive differs from pinned input: '+meta['archive'])
                pending.replace(archive)
        if digest(archive)!=meta['sha256']:raise RuntimeError('Official archive differs from pinned input: '+meta['archive'])
        with tarfile.open(archive) as t:
            members=[]
            for member in t.getmembers():
                # Build from source only; never copy upstream platform object files.
                target=(work/member.name).resolve()
                if not target.is_relative_to(work):raise RuntimeError('Unsafe archive path')
                if member.issym() or member.islnk() or not (member.isfile() or member.isdir()):continue
                if 'obj' in Path(member.name).parts or Path(member.name).suffix in ['.o','.so','.a','.dll']:continue
                members.append(member)
            t.extractall(work,members=members,filter='data')
    phg=work/'phg-0.9.7';solid=work/'solid'
    protected={'phg':manifest(phg,['examples/*.c','examples/*.h','test/*']),
               'solid':manifest(solid,['src/*.cxx','include/*','examples/*'])}
    if protected!=json.loads((ROOT/'provenance/case-files.json').read_text()):
        raise RuntimeError('Official case files do not match the repository source-protection manifest')
    (work/'case-source-manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    patches=[]
    def edit(relative,old,new):
        p=phg/relative;before=p.read_text();assert before.count(old)==1,(relative,old)
        after=before.replace(old,new);p.write_text(after)
        patches.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='a/'+relative,tofile='b/'+relative))
    edit('src/solver.c','/* list of solvers */','extern OEM_SOLVER phgSolverPanguLU_;\n\n/* list of solvers */')
    edit('src/solver.c','    "diagonal",','    "pangulu",\n    "diagonal",')
    edit('src/solver.c','    SOLVER_DIAGONAL\n','    &phgSolverPanguLU_,\n    SOLVER_DIAGONAL\n')
    edit('src/Makefile','OBJS =\t','OBJS = solver-pangulu.o repro-capture.o\t')
    edit('src/dof.c','    MagicCheck(DOF, dof)\n','    MagicCheck(DOF, dof)\n\n    if (dof != NULL) { extern void phgReproCapture(DOF *); phgReproCapture(dof); }\n')
    for source in ['solver-pangulu.c','repro-capture.c']:shutil.copyfile(ROOT/'integration/phg'/source,phg/'src'/source)
    shutil.copyfile(ROOT/'adapter/bridge.h',phg/'src/bridge.h')
    (work/'library-integration.patch').write_text(''.join(patches))
    current={'phg':manifest(phg,['examples/*.c','examples/*.h','test/*']),
             'solid':manifest(solid,['src/*.cxx','include/*','examples/*'])}
    assert current==protected,'Case source changed during integration'
    (work/'prepare-report.json').write_text(json.dumps({'official_archives_verified':True,'case_sources_unchanged':True,'prebuilt_objects_used':False,'changed_library_files':['src/solver.c','src/Makefile','src/dof.c'],'added_library_files':['src/solver-pangulu.c','src/repro-capture.c','src/bridge.h']},indent=2)+'\n')
    print('Prepared clean official source tree:',work)
if __name__=='__main__':main()
