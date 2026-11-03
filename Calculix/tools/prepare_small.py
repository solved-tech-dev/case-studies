# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Prepare an upstream-sized integration check; changes output requests only."""
import argparse
from pathlib import Path
from input_deck import flat_lines,header,parse

ROOT=Path(__file__).resolve().parents[1]

def prepare(case,output):
    source=ROOT/'cases'/case/'source'/('_steadystate.inp' if case=='pan' else case+'.inp')
    model=parse(source);lines=[];skip=False;injected=False
    for line in flat_lines(source):
        if line.strip().startswith('*') and not line.strip().startswith('**'):
            kind,_=header(line);skip=kind in ('NODE PRINT','NODE FILE','EL PRINT','EL FILE')
            if kind=='STEP':
                lines+=['*NSET,NSET=ALLN',','.join(map(str,sorted(model['nodes']))),
                    '*ELSET,ELSET=ALLE',','.join(map(str,sorted(model['elements'])))]
            if kind=='END STEP':
                lines+=['*NODE PRINT,NSET=ALLN','NT,RFL','*EL PRINT,ELSET=ALLE,TOTALS=ONLY','EVOL',
                    '*NODE FILE,OUTPUT=2D','NT','*EL FILE,OUTPUT=2D','HFL']
                injected=True
        if not skip:lines.append(line)
    assert injected
    # Limit set rows to CalculiX's normal input line length.
    result=[]
    for line in lines:
        if (len(line)>120 or len(line.split(','))>16) and all(c in '0123456789,' for c in line):
            values=line.split(',');result+= [','.join(values[i:i+12]) for i in range(0,len(values),12)]
        else:result.append(line)
    with output.open('x') as f:f.write('\n'.join(result)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case',choices=['dam','ax6ht','capacitor','oneel20fi','pan'])
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();prepare(a.case,a.output)
