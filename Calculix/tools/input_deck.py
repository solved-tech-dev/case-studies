# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Read only the node/element/set/boundary cards needed by the five examples."""
from pathlib import Path

WIDTH={'CPS8':8,'CAX6':6,'C3D20':20,'C3D20R':20}

def flat_lines(path,root=None,chain=()):
    path=Path(path).resolve();root=path.parent if root is None else Path(root).resolve()
    if not path.is_relative_to(root) or path in chain:raise ValueError('Unsafe/cyclic include')
    with path.open(encoding='latin1') as f:
        for line in f:
            if line.strip().upper().startswith('*INCLUDE'):
                _,a=header(line)
                yield from flat_lines(path.parent/a['INPUT'],root,chain+(path,))
            else:yield line.rstrip('\r\n')

def header(line):
    fields=[x.strip() for x in line.lstrip()[1:].split(',')]
    a={}
    for field in fields[1:]:
        k,sep,v=field.partition('=');a[k.upper()]=v.strip() if sep else True
    return fields[0].upper(),a

def parse(path):
    nodes={};elements={};ns={};es={};bc=[];film=[];flux=[]
    kind='';a={};pending=[];width=0
    for line in flat_lines(path):
        line=line.strip()
        if not line or line.startswith('**'):continue
        if line.startswith('*'):
            if pending:raise ValueError('Incomplete element')
            kind,a=header(line)
            if kind=='ELEMENT':width=WIDTH[a['TYPE'].upper()]
            continue
        v=[x.strip() for x in line.split(',') if x.strip()]
        if kind=='NODE':
            i=int(v[0]);nodes[i]=tuple(float(x.replace('D','E')) for x in v[1:4])
            if 'NSET' in a:ns.setdefault(a['NSET'].upper(),[]).append(i)
        elif kind=='ELEMENT':
            pending+=list(map(int,v))
            if len(pending)==width+1:
                i=pending[0];elements[i]=pending[1:];pending=[]
                if 'ELSET' in a:es.setdefault(a['ELSET'].upper(),[]).append(i)
            elif len(pending)>width+1:raise ValueError('Unexpected element width')
        elif kind in ('NSET','ELSET'):
            target=(ns if kind=='NSET' else es).setdefault(a[kind].upper(),[])
            if 'GENERATE' in a:
                lo,hi=int(v[0]),int(v[1]);step=int(v[2]) if len(v)>2 else 1
                target.extend(range(lo,hi+1,step))
            else:target.extend(int(x) if x.isdigit() else x.upper() for x in v)
        elif kind=='BOUNDARY':bc.append(v)
        elif kind=='FILM':film.append(v)
        elif kind=='DFLUX':flux.append(v)
    def resolve(sets,key,stack=()):
        if key in stack:raise ValueError('Cyclic set')
        result=set()
        for item in sets[key]:
            if isinstance(item,int):result.add(item)
            else:result.update(resolve(sets,item,stack+(key,)))
        return result
    return dict(nodes=nodes,elements=elements,nsets={k:resolve(ns,k) for k in ns},
        esets={k:resolve(es,k) for k in es},boundary=bc,film=film,flux=flux)
