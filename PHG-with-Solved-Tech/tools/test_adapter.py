#!/usr/bin/env python3
"""Host-only API contract check with the external PanguLU 26.44 public header."""
import argparse,json,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--api-include',required=True);p.add_argument('--output',default='_work/api-contract');a=p.parse_args();out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
command=[os.environ.get('CXX','c++'),'-std=c++11','-O2','-DPHG_PANGULU_TEST_API','-I'+str(ROOT/'adapter'),'-I'+str(Path(a.api_include).resolve()),'-I'+str(ROOT/'tests/stubs'),str(ROOT/'adapter/pangulu_bridge.cpp'),str(ROOT/'tests/test_api_contract.cpp'),'-o',str(out/'test_api_contract')]
r=subprocess.run(command,capture_output=True,text=True);(out/'compile.log').write_text(r.stdout+r.stderr);r.check_returncode()
r=subprocess.run([str(out/'test_api_contract')],capture_output=True,text=True);(out/'test.log').write_text(r.stdout+r.stderr);r.check_returncode()
(out/'report.json').write_text(json.dumps({'production_bridge_api_contract_pass':True,'real_cuda_runtime_tested':False,'real_pangulu_runtime_tested':False,'external_header_version':'26.44 / ABI 2'},indent=2)+'\n');print(r.stdout)
