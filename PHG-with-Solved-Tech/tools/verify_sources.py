#!/usr/bin/env python3
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--workspace',default='_work/clean');a=p.parse_args();w=Path(a.workspace).resolve()
expected=json.loads((w/'case-source-manifest.json').read_text());count=0
root=Path(__file__).resolve().parents[1]
assert expected==json.loads((root/'provenance/case-files.json').read_text()),'Workspace manifest differs from pinned official case manifest'
for group,base in [('phg',w/'phg-0.9.7'),('solid',w/'solid')]:
    for file,digest in expected[group].items():
        assert hashlib.sha256((base/file).read_bytes()).hexdigest()==digest,'Modified case input/source: '+file
        count+=1
print('PASS: original official case files unchanged:',count)
