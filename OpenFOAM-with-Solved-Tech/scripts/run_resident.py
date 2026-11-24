#!/usr/bin/env python3
"""Run one resident process using case-specific GS ordering and sealed configurations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from preload_resident_configs import PreparedConfigurations, document


def validate_gs(path):
    _, _, config = document(path)
    solver = config.get('solver', {})
    expected = {'method': 'smooth', 'smoother': 'gauss_seidel',
                'stopping': 'fixed_sweeps', 'preconditioner': 'none',
                'max_iterations': 1, 'initial_guess': 'provided'}
    if config.get('schema_version') != 2 or any(solver.get(k) != v for k, v in expected.items()):
        raise ValueError('GS configuration must use the resident-r1 provided-initial-guess one-sweep route')
    offsets = config.get('parameters', {}).get('partition_offsets')
    if offsets is not None:
        if (not isinstance(offsets, list) or len(offsets) < 2
                or any(type(v) is not int or not 0 <= v <= 2147483647 for v in offsets)
                or offsets[0] != 0 or any(a >= b for a, b in zip(offsets, offsets[1:]))):
            raise ValueError('partition_offsets must be strictly increasing int32 offsets starting at zero')
        return offsets[-1]
    return None


def digest(path):
    checksum = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    parser.add_argument('gs_config', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--check-case', action='store_true', help='validate complete input without creating GPU state')
    parser.add_argument('--dump-fields', action='store_true', help='write raw FP64 diagnostic field arrays')
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    case = args.case.resolve(strict=True)
    gs = args.gs_config.resolve(strict=True)
    cells = validate_gs(gs)
    executable = Path(os.environ['RESIDENT_BUILD']) / 'bin/rawResident'
    if not executable.is_file():
        parser.error('build rawResident with scripts/build-resident.sh first')
    provider = Path(os.environ['LINGLONGAMG_SDK']) / 'lib/liblinglongamg.so'
    output = args.output.absolute()
    if output.resolve().is_relative_to(case):
        parser.error('output must be outside the input case')
    output.mkdir(parents=True, exist_ok=False)
    output = output.resolve()
    provider_config = output / 'linglongamg.json'
    provider_config.write_bytes((package / 'examples/resident-linglong-vcycle-lowfill.json').read_bytes())
    pressure = json.loads((package / 'examples/resident-pressure-amg.json').read_text())
    pressure['extended_solver'].update(library_path=str(provider.resolve()), config_path=str(provider_config))
    pressure_config = output / 'pressure.json'
    pressure_config.write_text(json.dumps(pressure, indent=2) + '\n')
    command = [str(executable), '-case', str(case)]
    if cells is not None:
        command += ['-expected-cells', str(cells)]
    if args.check_case:
        command += ['-check-case']
    libraries = {
        'preparation': Path(os.environ['PREPARATION_SDK']) / 'lib/libLinearSystemPreparation.so',
        'core': Path(os.environ['XIANZHIKRYLOV_SDK']) / 'lib/libXianzhiKrylov.so',
        'linglongamg': provider,
    }
    identity = {name: {'path': str(path.resolve()), 'sha256': digest(path)} for name, path in libraries.items()}
    identity['executable'] = {'path': str(executable), 'sha256': digest(executable)}
    (output / 'runtime-identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    with PreparedConfigurations(gs, pressure_config) as prepared:
        environment = prepared.environment()
        if args.dump_fields:
            environment['RESIDENT_STATE_DUMPS'] = str(output / 'fields')
        (output / 'configuration-snapshots.json').write_text(json.dumps(prepared.metadata, indent=2) + '\n')
        with (output / 'application.log').open('w') as log, (output / 'application.err').open('w') as errors:
            status = subprocess.run(command, env=environment, pass_fds=tuple(prepared.fds),
                                    stdout=log, stderr=errors, cwd=output).returncode
    if status:
        print((output / 'application.err').read_text(), end='')
        raise SystemExit(status if status > 0 else 1)
    text = (output / 'application.log').read_text()
    marker = r'^RAW_RESIDENT_INPUT_VALIDATED gpu_initialized=0$' if args.check_case else r'^RAW_RESIDENT_COMPLETE steps=(\d+) requested_d2h_bytes=\d+$'
    match = re.search(marker, text, re.MULTILINE)
    if not match:
        raise RuntimeError('application exited without its expected completion marker')
    if not args.check_case:
        records = re.findall(r'^RAW_RESIDENT_STEP outer=(\d+) ', text, re.MULTILINE)
        if [int(value) for value in records] != list(range(int(match.group(1)))):
            raise RuntimeError('incomplete resident outer-step sequence')
    print('Input validated.' if args.check_case else 'Resident run completed.')
    print('Logs: ' + str(output))


if __name__ == '__main__':
    main()
