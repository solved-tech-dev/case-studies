# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Validate complete native temperature/flux fields and case-specific physics.

No NumPy, private mesh metadata, SDK or license is needed for this validator.
An optional independently supplied FRD enables full-field reference comparison.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import re
from input_deck import parse

CASES=['dam','ax6ht','capacitor','oneel20fi','pan']

def text(path):
    path=Path(path)
    if not path.exists():path=Path(str(path)+'.gz')
    return gzip.open(path,'rt') if path.suffix=='.gz' else path.open()

def fields(path):
    result={};active=None
    with text(path) as f:
        for line in f:
            if line.startswith(' -4'):
                active=next((k for k in ('NDTEMP','FLUX') if k in line),None)
                if active:result[active]={}
            elif active and line.startswith(' -1'):
                result[active][int(line[3:13])]=tuple(float(v.replace('D','E')) for v in
                    re.findall(r'[-+]?\d+\.\d+(?:[EeDd][-+]?\d+)?',line[13:]))
            elif line.startswith(' -3'):active=None
    return result

def scalar_table(dat,label):
    result={};active=False
    for line in dat.splitlines():
        if label.lower() in line.lower():active=True;result={};continue
        if active:
            v=line.split()
            if len(v)==2 and v[0].isdigit():result[int(v[0])]=float(v[1])
            elif result:active=False
    return result

def check(case,run,reference=None):
    model=parse(run/'model.inp');used={n for row in model['elements'].values() for n in row}
    f=fields(run/'model.frd');log=(run/'application.stdout').read_text(errors='replace')
    status=json.loads((run/'application-exit.json').read_text());assert status['exit_code']==0
    assert 'Job finished' in log and '*ERROR' not in log
    warnings=[x for x in log.splitlines() if '*WARNING' in x and
        '*WARNING reading *HEAT TRANSFER: a nonlinear geometric analysis is requested' not in x]
    assert not warnings,warnings
    report=dict(status='passed',case=case,used_nodes=len(used),full_fields={},reference_comparison_performed=reference is not None)
    for name,width in (('NDTEMP',1),('FLUX',3)):
        assert used<=set(f[name]),(name,'missing original nodes')
        assert all(len(v)==width and all(math.isfinite(x) for x in v) for v in f[name].values())
        report['full_fields'][name]=dict(nodes=len(f[name]),finite=True)
    if reference:
        ref=fields(reference)
        for name,tolerance in (('NDTEMP',1e-5),('FLUX',1e-4)):
            assert set(ref[name])==set(f[name]),'Reference must use exactly the same mesh'
            ds=rs=0.
            for node,row in f[name].items():
                for x,y in zip(row,ref[name][node]):ds+=(x-y)**2;rs+=y*y
            rel=math.sqrt(ds/(rs or 1));assert rel<tolerance,(name,rel)
            report['full_fields'][name]['reference_relative_L2']=rel
    with text(run/'model.dat') as stream:dat=stream.read()
    temp=scalar_table(dat,'temperatures for set ALLN');heat=scalar_table(dat,'heat generation for set ALLN')
    assert used<=set(temp) and used<=set(heat)
    assert all(math.isfinite(x) for x in list(temp.values())+list(heat.values()))
    report['temperature_range']=[min(temp.values()),max(temp.values())]
    # All prescribed scalar degrees of freedom, including set expansion.
    for row in model['boundary']:
        assert len(row)>=4 and int(row[1])==int(row[2])==11
        nodes=[int(row[0])] if row[0].isdigit() else model['nsets'][row[0].upper()]
        assert max(abs(temp[n]-float(row[3])) for n in nodes)<1e-5
    if case!='pan':
        balance=abs(sum(heat.values()))/max(sum(abs(x) for x in heat.values()),1e-30)
        assert balance<1e-5,balance;report['boundary_heat_balance_relative']=balance
    if case=='dam':assert 3.9<min(temp.values())<=max(temp.values())<8.6
    elif case=='ax6ht':
        assert abs(temp[1])<1e-5 and abs(temp[46]-100)<1e-5
        assert -.01<min(temp.values()) and max(temp.values())<100.01
        report['interpretation']='Original point temperatures; mesh-dependent localization, not an edge-temperature engineering model'
    elif case=='capacitor':
        assert .99<min(temp.values()) and max(temp.values())<2.01
        rms=math.sqrt(sum((temp[i]-(2-math.log(math.hypot(*model['nodes'][i][1:3]))/math.log(2)))**2 for i in used)/len(used))
        assert rms<.002;report['potential_RMS_vs_radial_analytic']=rms
        actual=float(re.search(r'total surface flux \(q\).*?\n\s*([-+0-9.Ee]+)',dat,re.I)[1])
        expected=-8.8541878176e-12*(math.pi/180)*10/math.log(2)
        error=abs(actual/expected-1);assert error<.02
        report['surface_flux_relative_error_vs_analytic']=error
    elif case=='oneel20fi':
        rms=math.sqrt(sum((temp[i]+2*model['nodes'][i][0])**2 for i in used)/len(used))
        qrms=math.sqrt(sum(sum((v[j]-(100 if j==0 else 0))**2 for j in range(3)) for i,v in f['FLUX'].items() if i in used)/len(used))
        assert rms<1e-4 and qrms<.01,(rms,qrms)
        report.update(temperature_RMS_vs_T_minus2x=rms,flux_RMS_vs_100_0_0=qrms)
    else:
        assert min(temp.values())>=293.15-.01
        quadrature=[]
        for x,w in ((-math.sqrt(3/5),5/9),(0.,8/9),(math.sqrt(3/5),5/9)):
            quadrature.append(([x*(x-1)/2,x*(x+1)/2,1-x*x],[x-.5,x+.5,-2*x],w))
        def integrate(eid,face,value,h=None):
            conn=model['elements'][eid];loc=list(((0,1),(1,2),(2,0))[face-1])+[face+2]
            ids=[conn[j] for j in loc];xyz=[model['nodes'][i] for i in ids];total=0.
            for N,D,w in quadrature:
                r=sum(N[j]*xyz[j][0] for j in range(3));assert r>=-1e-10
                ds=math.sqrt(sum(sum(D[j]*xyz[j][k] for j in range(3))**2 for k in (0,1)))
                q=value if h is None else h*(value-sum(N[j]*temp[ids[j]] for j in range(3)))
                total+=q*2*math.pi*max(r,0)*ds*w
            return total
        def targets(row):return [int(row[0])] if row[0].isdigit() else model['esets'][row[0].upper()]
        qin=sum(integrate(e,int(row[1][1:]),float(row[2])) for row in model['flux'] for e in targets(row))
        qout=sum(integrate(e,int(row[1][1:]),float(row[2]),float(row[3])) for row in model['film'] for e in targets(row))
        assert qin>0 and qout<0
        balance=abs(qin+qout)/qin;assert balance<1e-3,balance
        report.update(axisymmetric_full360_input=qin,Robin_output=qout,Robin_heat_balance_relative=balance)
    volume=float(re.search(r'total volume.*?\n\s*([-+0-9.Ee]+)',dat,re.I)[1]);assert math.isfinite(volume) and volume>0
    report['native_volume']=volume
    report['Newton_iterations']=list(map(int,re.findall(r'^\s*iteration\s+(\d+)\s*$',log,re.M)))
    assert len(re.findall(r'increment 1 attempt 1',log))==1
    assert report['Newton_iterations']==list(range(1,len(report['Newton_iterations'])+1)) and report['Newton_iterations']
    report['physical_time']='original single steady increment, end time 1'
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case',choices=CASES);p.add_argument('run',type=Path)
    p.add_argument('--reference-frd',type=Path);a=p.parse_args()
    try:r=check(a.case,a.run,a.reference_frd)
    except (AssertionError,ValueError,KeyError,FileNotFoundError,IndexError) as e:
        r=dict(status='failed',case=a.case,error=repr(e))
    with (a.run/'field-validation.json').open('x') as stream:json.dump(r,stream,indent=2);stream.write('\n')
    print(json.dumps(r,indent=2));raise SystemExit(0 if r['status']=='passed' else 1)
