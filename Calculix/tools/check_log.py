# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Check actual adapter logs, not a replacement for physical-field validation."""
import argparse
import json
import math
import re
from pathlib import Path

def check(run,params):
    log=(run/'application.stdout').read_text(errors='replace')
    entries=re.findall(r'HULL_ENTRY call=(\d+) n=(\d+) lower_nnz=(\d+)',log)
    residuals=[float(x) for x in re.findall(r'HULL_RESIDUAL relative_l2=([^\s]+)',log)]
    phases=[dict(phase=p,wall_s=float(w),CPU_s=float(c)) for p,w,c in re.findall(r'HULL_TIME native (\S+) wall=([^\s]+) cpu=([^\s]+)',log)]
    expected=params['measured_system'];errors=[]
    status=json.loads((run/'application-exit.json').read_text())
    if status['exit_code']!=0:errors.append(f'application_exit={status["exit_code"]}')
    if len(entries)!=expected['linear_calls_per_complete_application']:errors.append('unexpected_number_of_linear_calls')
    for ordinal,(call,n,nz) in enumerate(entries,1):
        if int(call)!=ordinal:errors.append('entry_call_mismatch')
        if expected.get('n') is not None and (int(n)!=expected['n'] or 2*int(nz)+int(n)!=expected['stored_full_nnz']):errors.append('entry_dimension_mismatch')
    if len(residuals)!=len(entries) or not residuals:errors.append('missing_independent_original_residual')
    if any(not math.isfinite(x) or x>1e-8 or x<0 for x in residuals):errors.append('invalid_original_residual')
    if not phases or any(not math.isfinite(p[k]) or p[k]<0 for p in phases for k in ('wall_s','CPU_s')):errors.append('invalid_or_missing_native_times')
    if 'Job finished' not in log:errors.append('missing_native_job_finished')
    if 'ordering=L0 reuse=L2' not in log or 'HULL_PROVIDER pangulu_execute ' not in log:errors.append('missing_PanguLU_L0_L2_or_provider_evidence')
    if 'HULL_FAILURE' in log:errors.append('adapter_failure')
    if len(re.findall('HULL_ACTUAL_INPUT_RETAINED ',log))!=len(entries):errors.append('missing_complete_original_CSR_RHS')
    if len([x for x in phases if x['phase']=='analysis_order_symbolic_preprocess'])!=1:errors.append('expected_one_real_analysis')
    if len([x for x in phases if x['phase'] in ('numeric','refactor_numeric')])!=len(entries):errors.append('numeric_factorization_count_mismatch')
    for call in range(2,len(entries)+1):
        if f'HULL_REUSE call={call} same_entire_graph=1 native_analysis_skipped=1' not in log:errors.append('missing_actual_graph_reuse')
    fields={name:(run/name).is_file() and (run/name).stat().st_size>0 for name in ('model.frd','model.dat')}
    if not all(fields.values()):errors.append('missing_native_output')
    return dict(status='log_checks_passed' if not errors else 'failed',errors=errors,
        actual_linear_calls=len(entries),residuals=residuals,native_phases=phases,
        measured_dimensions=[dict(n=int(n),full_nnz=2*int(nz)+int(n)) for _,n,nz in entries],
        native_wall_sum_s=sum(p['wall_s'] for p in phases),nonempty_native_output=fields,
        full_physical_field_validation='not_performed_by_this_tool',
        qualified_performance_sample=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path);parser.add_argument('parameters',type=Path)
    args=parser.parse_args()
    result=check(args.run,json.loads(args.parameters.read_text()))
    with (args.run/'log-validation.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['status']=='log_checks_passed' else 1)
