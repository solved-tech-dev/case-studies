# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Geometry-preserving mesh preparation for the five supported thermal examples.

Extracted from the original model-preparation routines. No solver is executed.
Dependencies: NumPy. Fixed tested meshes are supplied separately in cases/.
"""
import argparse
import collections
import itertools
import json
from pathlib import Path
import numpy as np
from input_deck import flat_lines
ROOT=Path(__file__).resolve().parents[1]
OUT=None
Q8=np.array([[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[1,0],[0,1],[-1,0]])
H20=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1],
 [0,-1,-1],[1,0,-1],[0,1,-1],[-1,0,-1],[0,-1,1],[1,0,1],[0,1,1],[-1,0,1],[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]])


def shape(points, ref):
    points=np.atleast_2d(points); d=ref.shape[1]; out=[]
    if not np.any(ref==0):
        return np.array([np.prod(1+points*v,axis=1)/2**d for v in ref]).T
    for v in ref:
        zero=np.flatnonzero(v==0)
        if len(zero)==0:
            out.append(np.prod(1+points*v,axis=1)*(points@v-(d-1))/2**d)
        else:
            j=zero[0]; ids=np.flatnonzero(v!=0)
            out.append((1-points[:,j]**2)*np.prod(1+points[:,ids]*v[ids],axis=1)/2**(d-1))
    return np.array(out).T

def blocks(path):
    result=[]
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith('**'): continue
        if line.startswith('*'):
            p=[s.strip().upper() for s in line.split(',')]
            result.append((p[0][1:],dict(s.split('=',1) if '=' in s else (s,True) for s in p[1:]),[]))
        else: result[-1][2].append(line)
    return result

def vals(lines): return [s.strip() for line in lines for s in line.split(',') if s.strip()]

def nums(lines): return [int(v) for v in vals(lines)]

def indices(lines, gen=False):
    if not gen: return nums(lines)
    out=[]
    for line in lines:
        x=[int(v) for v in line.split(',') if v.strip()]
        out.extend(range(x[0],x[1]+1,x[2] if len(x)>2 else 1))
    return out

def output(case,level,nodes,elements,kind,nsets,esets,cards,meta):
    directory=OUT; directory.mkdir(parents=True,exist_ok=True)
    assert not (directory/'model.inp').exists(), 'do not overwrite a generated model'
    lines=['*HEADING',case+' mesh-refinement comparison; provenance in model.json','*NODE,NSET=ALLN']
    # CCX nodes.f reads each coordinate from textpart(field)(1:20).
    # Fixed scientific width preserves tiny negative coordinates without
    # truncating an exponent (e.g. -9.99...e-16 to -9.99...e-1).
    coordinate_lines=[str(n)+','+','.join(f'{v:.12e}' for v in p) for n,p in sorted(nodes.items())]
    assert all(len(v)<=20 for line in coordinate_lines for v in line.split(','))
    lines += coordinate_lines
    lines += ['*ELEMENT,TYPE='+kind+',ELSET=ALLE']
    for eid,conn in elements.items():
        if len(conn)==20:
            lines += [str(eid)+','+','.join(map(str,conn[:10])), ','.join(map(str,conn[10:]))]
        else: lines.append(str(eid)+','+','.join(map(str,conn)))
    for tag, ids in nsets.items():
        if tag=='ALLN': continue
        lines.append('*NSET,NSET='+tag)
        ids=sorted(ids); lines += [','.join(map(str,ids[i:i+12])) for i in range(0,len(ids),12)]
    for tag,ids in esets.items():
        if tag=='ALLE': continue
        lines.append('*ELSET,ELSET='+tag)
        ids=sorted(ids); lines += [','.join(map(str,ids[i:i+12])) for i in range(0,len(ids),12)]
    lines += cards
    (directory/'model.inp').write_text('\n'.join(lines)+'\n')
    meta.update(case=case,level=level,kind=kind,input_nodes=len(nodes),elements=len(elements),
        used_nodes=len({n for conn in elements.values() for n in conn}),
        source_files_retained=True,actual_n=None,actual_full_nnz=None,
        output='complete native NT/RFL DAT and NDTEMP/FLUX FRD; original scalar physics/control cards retained',
        coordinate_serialization='12-decimal scientific notation, every CCX field <=20 characters')
    (directory/'model.json').write_text(json.dumps(meta,indent=2)+'\n')
    (directory/'mesh.json').write_text(json.dumps(dict(nodes=nodes,elements=elements,nsets={k:sorted(v) for k,v in nsets.items()})))
    print(json.dumps(meta),flush=True)

def refine(case,f,level):
    source=ROOT/'cases'/case/'source'/f'{case}.inp'; bb=blocks(source)
    nodes={}; elems={}; kind=None; nsets={}; eset0={}
    for k,a,data in bb:
        if k=='NODE':
            for line in data:
                p=[float(v) for v in line.split(',') if v.strip()]; nodes[int(p[0])]=p[1:4]
            if 'NSET' in a: nsets.setdefault(a['NSET'],set()).update(int(line.split(',')[0]) for line in data)
        elif k=='ELEMENT':
            assert kind in (None,a['TYPE']); kind=a['TYPE']; count=20 if kind.startswith('C3D20') else 8
            raw=nums(data); assert len(raw)%(count+1)==0
            for i in range(0,len(raw),count+1): elems[raw[i]]=raw[i+1:i+count+1]
            if 'ELSET' in a: eset0.setdefault(a['ELSET'],set()).update(raw[::count+1])
        elif k=='NSET': nsets.setdefault(a['NSET'],set()).update(indices(data,'GENERATE' in a))
        elif k=='ELSET':
            if 'GENERATE' in a: ids=indices(data,True)
            else:
                ids=[]
                for x in vals(data): ids.extend(eset0[x.upper()] if not x.isdigit() else [int(x)])
            eset0.setdefault(a['ELSET'],set()).update(ids)
    assert kind in ('C3D20','C3D20R','S8R','CAX8R','C3D8','CPS8')
    d=3 if kind.startswith('C3D') else 2; ref=(H20[:8] if kind=='C3D8' else H20) if d==3 else Q8
    fs=np.array([f]*d if isinstance(f,int) else f,dtype=int);assert len(fs)==d and np.min(fs)>0
    assert np.max(abs(shape(ref,ref)-np.eye(len(ref))))<1e-12
    loadsets={line.split(',')[0].strip().upper() for k,a,data in bb if k=='CLOAD' for line in data}
    oldnodes=nodes.copy(); lookup={tuple(np.round(p,10)):n for n,p in nodes.items()}
    assert len(lookup)==len(nodes), 'coincident original nodes must not be silently merged'
    linear=kind=='C3D8'
    grid=[p for p in itertools.product(*(range(2*q+1) for q in fs)) if sum(x%2 for x in p)<=(0 if linear else 1)]
    param=np.array(grid)/fs-1; N=shape(param,ref); where={p:i for i,p in enumerate(grid)}
    active=np.abs(N)>1e-12
    nextnode=max(nodes)+1; newelements={}; children={}; childposition={}; newsets={k:set(v) for k,v in nsets.items()}
    for parent,conn in elems.items():
        xyz=np.array([oldnodes[n] for n in conn]); coords=N@xyz; ids=[]
        localsets={tag:~np.any(active[:,[n not in members for n in conn]],axis=1)
                   for tag,members in nsets.items() if tag not in loadsets}
        for index,p in enumerate(coords):
            key=tuple(np.round(p,10)); node=lookup.get(key)
            if node is None:
                node=nextnode; nextnode+=1; lookup[key]=node; nodes[node]=p.tolist()
            ids.append(node)
            for tag,mask in localsets.items():
                if mask[index]: newsets[tag].add(node)
        children[parent]=[]
        for pos in itertools.product(*(range(q) for q in fs)):
            center=2*np.array(pos)+1
            child=[ids[where[tuple(center+q)]] for q in ref]
            eid=len(newelements)+1; newelements[eid]=child; children[parent].append(eid); childposition[eid]=pos
    # Confirm each child mapping is non-degenerate/orientation-preserving at2^d Gauss points.
    # Finite differences differentiate only known polynomial shape functions, not the physics.
    xyz=np.array([[nodes[n] for n in c] for c in newelements.values()]); minjac=float('inf')
    for gp in itertools.product((-1/np.sqrt(3),1/np.sqrt(3)),repeat=d):
        D=[]
        for j in range(d):
            e=np.eye(d)[j]*1e-6; D.append(((shape(np.array(gp)+e,ref)-shape(np.array(gp)-e,ref))/(2e-6))[0])
        J=np.einsum('ai,eib->eab',np.array(D),xyz)
        det=np.linalg.det(J) if d==3 else np.linalg.norm(np.cross(J[:,0],J[:,1]),axis=1)
        assert np.all(np.isfinite(det)) and np.min(det)>1e-18
        minjac=min(minjac,float(np.min(det)))
    esets={tag:{c for eid in members for c in children[eid]} for tag,members in eset0.items()}
    cards=[]
    face_axes={1:(2,0),2:(2,fs[-1]-1),3:(1,0),4:(0,fs[0]-1),5:(1,fs[1]-1),6:(0,0)}
    skip={'HEADING','NODE','ELEMENT','NSET','ELSET','NODE FILE','EL FILE','NODE PRINT','EL PRINT','END STEP'}
    for k,a,data in bb:
        if k in skip: continue
        assert k not in ('EQUATION','TRANSFORM','NODAL THICKNESS','INCLUDE'), 'unsupported original card must not be silently dropped'
        assert not (k=='STEP' and 'NLGEOM' in a), 'no silent linearization of official verification inputs'
        header='*'+k+''.join(','+key+('' if value is True else '='+value) for key,value in a.items())
        cards.append(header)
        if k=='SURFACE':
            assert d==3
            for line in data:
                oldid,face=[x.strip() for x in line.split(',')]; fn=int(face[1:]); axis,target=face_axes[fn]
                cards += [f'{c},S{fn}' for c in children[int(oldid)] if childposition[c][axis]==target]
        elif k=='FILM':
            assert case=='oneel20fi' and d==3
            for line in data:
                oldid,face,ambient,coefficient=[x.strip() for x in line.split(',')]
                assert face.upper()=='F1' and float(ambient)==-12 and float(coefficient)==10
                axis,target=face_axes[1]
                selected=[c for c in children[int(oldid)] if childposition[c][axis]==target]
                assert len(selected)==fs[0]*fs[1]
                cards += [f'{c},{face},{ambient},{coefficient}' for c in selected]
        else: cards += data
    cards += ['*NODE PRINT,NSET=ALLN','NT,RFL','*EL PRINT,ELSET=ALLE,TOTALS=ONLY','EVOL',
              '*NODE FILE,OUTPUT=2D','NT','*EL FILE,OUTPUT=2D','HFL','*END STEP']
    used={n for conn in elems.values() for n in conn}
    output(case,level,nodes,newelements,kind,newsets,esets,cards,dict(source=str(source.relative_to(ROOT)),
        subdivision_factor=f,mesh_dimension=d,original_nodes=len(oldnodes),original_elements=len(elems),
        original_reference_nodes=sorted(set(oldnodes)-used),minimum_child_Gauss_Jacobian=minjac,
        geometry='exact subdivision of original isoparametric parent elements, coordinate merge1e-10 source units',
        boundary_policy='retain original scalar boundary values, sets and descendant boundary faces',
        original_physics_cards_unchanged=True))
    meta_path=OUT/'model.json'
    saved=json.loads(meta_path.read_text())
    saved['output']='complete native NT/RFL DAT and NDTEMP/FLUX FRD'
    meta_path.write_text(json.dumps(saved,indent=2)+'\n')

def tri_shape(p):
    p=np.atleast_2d(p);a=1-p[:,0]-p[:,1];b=p[:,0];c=p[:,1]
    return np.array([a*(2*a-1),b*(2*b-1),c*(2*c-1),4*a*b,4*b*c,4*c*a]).T

def tri_refine(f,level):
    bb=blocks(ROOT/'cases/ax6ht/source/ax6ht.inp');oldnodes={};elems={}
    for k,a,data in bb:
        if k=='NODE':
            for line in data:
                v=[float(x) for x in line.split(',') if x.strip()];oldnodes[int(v[0])]=v[1:4]
        elif k=='ELEMENT':
            assert a['TYPE']=='CAX6';v=nums(data);assert len(v)%7==0
            for i in range(0,len(v),7):elems[v[i]]=v[i+1:i+7]
    ref=np.array([[0,0],[1,0],[0,1],[.5,0],[.5,.5],[0,.5]])
    assert np.max(abs(tri_shape(ref)-np.eye(6)))<1e-12
    nodes=oldnodes.copy();lookup={tuple(np.round(x,10)):n for n,x in nodes.items()};assert len(lookup)==len(nodes)
    grid=[(i,j) for i in range(2*f+1) for j in range(2*f+1-i)];N=tri_shape(np.array(grid)/(2*f));where={p:i for i,p in enumerate(grid)}
    elements={}
    def child(ids,a,b,c):
        mids=[tuple((np.array(x)+y)//2) for x,y in ((a,b),(b,c),(c,a))]
        elements[len(elements)+1]=[ids[where[p]] for p in (a,b,c,*mids)]
    for conn in elems.values():
        xyz=np.array([oldnodes[n] for n in conn]);coords=N@xyz;ids=[]
        for p in coords:
            key=tuple(np.round(p,10));node=lookup.get(key)
            if node is None:node=max(oldnodes)+len(nodes)-len(oldnodes)+1;lookup[key]=node;nodes[node]=p.tolist()
            ids.append(node)
        for i in range(f):
            for j in range(f-i):
                child(ids,(2*i,2*j),(2*i+2,2*j),(2*i,2*j+2))
                if i+j<f-1:child(ids,(2*i+2,2*j),(2*i+2,2*j+2),(2*i,2*j+2))
    assert len(elements)==len(elems)*f*f
    xyz=np.array([[nodes[n][:2] for n in conn] for conn in elements.values()]);minjac=float('inf')
    for gp in ((1/6,1/6),(2/3,1/6),(1/6,2/3)):
        D=[]
        for j in range(2):
            e=np.eye(2)[j]*1e-6;D.append(((tri_shape(np.array(gp)+e)-tri_shape(np.array(gp)-e))/(2e-6))[0])
        det=np.linalg.det(np.einsum('ai,eib->eab',np.array(D),xyz));assert np.all(np.isfinite(det)) and np.min(det)>1e-18
        minjac=min(minjac,float(np.min(det)))
    cards=[]
    for k,a,data in bb:
        if k in ('NODE','ELEMENT','EL PRINT','END STEP'):continue
        assert k in ('MATERIAL','DENSITY','CONDUCTIVITY','SPECIFIC HEAT','SOLID SECTION','INITIAL CONDITIONS','STEP','HEAT TRANSFER','BOUNDARY')
        cards.append('*'+k+''.join(','+key+('' if value is True else '='+value) for key,value in a.items()));cards+=data
    cards+=['*NODE PRINT,NSET=ALLN','NT,RFL','*EL PRINT,ELSET=ALLE,TOTALS=ONLY','EVOL','*NODE FILE,OUTPUT=2D','NT','*EL FILE,OUTPUT=2D','HFL','*END STEP']
    output('ax6ht',level,nodes,elements,'CAX6',{'NALL':set(nodes)}, {'EALL':set(elements)},cards,
        dict(source='cases/ax6ht/source/ax6ht.inp',subdivision_factor=f,mesh_dimension=2,original_nodes=len(oldnodes),original_elements=len(elems),
             original_reference_nodes=[],minimum_child_Gauss_Jacobian=minjac,original_physics_cards_unchanged=True,
             geometry='exact quadratic T6 parent subdivision; r,z axisymmetric geometry retained',
             boundary='original point temperatures ONLY at node1=0 and46=100; not expanded to edges; mesh-dependent localization is not a converged engineering prediction'))
    assert len(nodes)==(16*f+1)*(4*f+1), 'conforming T6 interfaces must match exact structured source topology'

def blocks_lines(lines):
    result=[]
    for line in lines:
        line=line.strip()
        if not line or line.startswith('**'):continue
        if line.startswith('*'):
            items=[x.strip().upper() for x in line.split(',')]
            result.append((items[0][1:],dict(x.split('=',1) if '=' in x else (x,True) for x in items[1:]),[]))
        else:result[-1][2].append(line)
    return result

EDGES=((0,1),(1,2),(2,0));OMIT=(2,0,1)
def pan_shape(p):
    p=np.atleast_2d(p);a=1-p.sum(axis=1);b=p[:,0];c=p[:,1]
    return np.column_stack((a*(2*a-1),b*(2*b-1),c*(2*c-1),4*a*b,4*b*c,4*a*c))

def pan_template(f):
    tri=[]
    for i in range(f):
        for j in range(f-i):
            tri.append([(i,j),(i+1,j),(i,j+1)])
            if i+j<f-1:tri.append([(i+1,j),(i+1,j+1),(i,j+1)])
    points={};con=[]
    for t in np.array(tri)/f:
        ids=[]
        for q in np.concatenate([t,np.array([(t[i]+t[j])/2 for i,j in EDGES])]):
            key=tuple(np.rint(q*2*f).astype(int));points.setdefault(key,len(points));ids.append(points[key])
        con.append(ids)
    P=np.array(list(points))/(2*f);C=np.array(con,dtype=int);B={i:[] for i in (1,2,3)};seen=collections.Counter()
    assert len(C)==f*f
    for child,c in enumerate(C):
        for face,loc in enumerate(EDGES,1):
            ids=c[list(loc)];seen[tuple(sorted(ids))]+=1
            L=np.column_stack((1-P[ids].sum(axis=1),P[ids]))
            for original,omit in enumerate(OMIT,1):
                if abs(L[:,omit]).max()<1e-12:B[original].append((child,face))
    assert all(n in (1,2) for n in seen.values()) and sum(n==1 for n in seen.values())==3*f
    assert all(len(x)==f for x in B.values())
    return P,C,B

def pan_refine(f,level):
    folder=ROOT/'cases/pan/source';bb=blocks_lines(flat_lines(folder/'_steadystate.inp',folder));nodes={};elems={};ns={};es={}
    for k,a,d in bb:
        if k=='NODE':
            for l in d:
                v=[x.strip() for x in l.split(',') if x.strip()];nodes[int(v[0])]=list(map(float,v[1:4]))
            if 'NSET' in a:ns.setdefault(a['NSET'],set()).update(int(l.split(',')[0]) for l in d)
        elif k=='ELEMENT':
            assert a['TYPE']=='CAX6';v=[int(x) for l in d for x in l.split(',') if x.strip()];assert len(v)%7==0
            for i in range(0,len(v),7):elems[v[i]]=v[i+1:i+7]
            if 'ELSET' in a:es.setdefault(a['ELSET'],set()).update(v[::7])
        elif k in ('NSET','ELSET'):
            sets=ns if k=='NSET' else es
            if 'GENERATE' in a:ids=indices(d,True)
            else:
                ids=[]
                for x in (x.strip().upper() for l in d for x in l.split(',') if x.strip()):ids+=([int(x)] if x.isdigit() else list(sets[x]))
            sets.setdefault(a[k],set()).update(ids)
    old=nodes.copy();lookup={tuple(np.round(p,10)):i for i,p in nodes.items()};assert len(lookup)==len(nodes)
    P,C,B=pan_template(f);N=pan_shape(P);active=abs(N)>1e-12;new={};children={};newns={k:set(v) for k,v in ns.items()};nextnode=max(nodes)+1
    for parent,conn in elems.items():
        xyz=N@np.array([old[i] for i in conn]);ids=[];masks={k:~np.any(active[:,[i not in v for i in conn]],axis=1) for k,v in ns.items()}
        for j,p in enumerate(xyz):
            key=tuple(np.round(p,10));node=lookup.get(key)
            if node is None:node=nextnode;nextnode+=1;lookup[key]=node;nodes[node]=p.tolist()
            ids.append(node)
            for tag,mask in masks.items():
                if mask[j]:newns[tag].add(node)
        childids=[]
        for row in C:
            eid=len(new)+1;new[eid]=[ids[j] for j in row];childids.append(eid)
        children[parent]=childids
    newes={k:{c for p in v for c in children[p]} for k,v in es.items()};cards=[];film=[];flux=[]
    skip={'HEADING','NODE','ELEMENT','NSET','ELSET','NODE PRINT','NODE FILE','EL PRINT','EL FILE','END STEP'}
    for k,a,d in bb:
        if k in skip:continue
        assert k in ('SURFACE','MATERIAL','CONDUCTIVITY','SOLID SECTION','STEP','HEAT TRANSFER','FILM','DFLUX')
        cards.append('*'+k+''.join(','+x+('' if y is True else '='+y) for x,y in a.items()))
        if k in ('FILM','DFLUX','SURFACE'):
            for line in d:
                v=[x.strip() for x in line.split(',')];of=int(v[1][1:]);assert 1<=of<=3
                targets=[int(v[0])] if v[0].isdigit() else es[v[0].upper()]
                for parent in targets:
                    for child,face in B[of]:
                        eid=children[parent][child]
                        if k=='FILM':cards.append(f'{eid},F{face},{v[2]},{v[3]}');film.append((eid,face,float(v[2]),float(v[3])))
                        elif k=='DFLUX':cards.append(f'{eid},S{face},{v[2]}');flux.append((eid,face,float(v[2])))
                        else:cards.append(f'{eid},S{face}')
        else:cards+=d
    cards+=['*NODE PRINT,NSET=ALLN','NT,RFL','*EL PRINT,ELSET=ALLE,TOTALS=ONLY','EVOL','*NODE FILE,OUTPUT=2D','NT','*EL FILE,OUTPUT=2D','HFL','*END STEP']
    minjac=float('inf');connlist=list(new.values())
    for start in range(0,len(new),30000):
        xyz=np.array([[nodes[i][:2] for i in row] for row in connlist[start:start+30000]])
        for gp in ((1/6,1/6),(2/3,1/6),(1/6,2/3)):
            D=np.array([((pan_shape(np.array(gp)+np.eye(2)[j]*1e-6)-pan_shape(np.array(gp)-np.eye(2)[j]*1e-6))/(2e-6))[0] for j in range(2)])
            det=np.linalg.det(np.einsum('ai,eib->eab',D,xyz));assert np.isfinite(det).all() and det.min()>1e-18
            minjac=min(minjac,float(det.min()))
    output('pan',level,nodes,new,'CAX6',newns,newes,cards,dict(source='cases/pan/source/_steadystate.inp and original includes',
       subdivision_factor=f,mesh_dimension=2,original_nodes=len(old),original_elements=len(elems),original_reference_nodes=[],
       minimum_child_Gauss_Jacobian=minjac,geometry='exactquadratictriangleparentsubdivision,originalaxisymmetricironpanprofile',
       original_physics_cards_unchanged=True,source_revision='a079253398c05f96b019e8c07ce076e370ddeb75',
       material='originalironconductivity50',boundary='originalairh10/foodh150sink293.15,originalbottomflux31500',
       dimension_disclosure='2Daxisymmetric,not3Dmesh;nativesectorfields/volume retained;independentbalanceusesfull360axisweight'))
    (OUT/'film.json').write_text(json.dumps(film)+'\n');(OUT/'flux.json').write_text(json.dumps(flux)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',choices=['dam','ax6ht','capacitor','oneel20fi','pan'])
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--factor',type=int);group.add_argument('--factors',type=int,nargs=3)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();OUT=args.output.resolve()
    if OUT.exists():parser.error('Choose a new output directory; existing models are never overwritten')
    factor=args.factors if args.factors else args.factor
    if any(x<1 for x in (factor if isinstance(factor,list) else [factor])):parser.error('Positive subdivisions required')
    if args.case in ('ax6ht','pan') and not args.factor:parser.error('This triangular mesh takes one --factor')
    if args.case=='pan':pan_refine(factor,'generated')
    elif args.case=='ax6ht':tri_refine(factor,'generated')
    else:refine(args.case,factor,'generated')
