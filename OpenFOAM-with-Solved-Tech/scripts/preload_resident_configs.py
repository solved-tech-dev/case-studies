#!/usr/bin/env python3
"""Preload local26.47 experiment configurations; never changes solver parameters."""
import argparse
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess


def document(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    def nonfinite(value):
        raise ValueError('nonfinite JSON configuration: ' + value)
    def object_pairs(pairs):
        result = {}
        for name, value in pairs:
            if name in result:
                raise ValueError('duplicate JSON option: ' + name)
            result[name] = value
        return result
    value = json.loads(data, parse_constant=nonfinite, object_pairs_hook=object_pairs)
    if not isinstance(value, dict):
        raise ValueError('configuration must be a JSON object')
    return path, data, value


class PreparedConfigurations:
    def __init__(self, gs, pressure):
        self.fds = []
        self.metadata = {'input_data_preloaded': True,
            'configuration_validation_in_compute': True,
            'authorization_and_library_initialization_in_compute': True,
            'public_API_ABI_changed': False, 'performance_claim': False}
        try:
            gs_path, gs_bytes, gs_doc = document(gs)
            pressure_path, pressure_bytes, pressure_doc = document(pressure)
            solver = gs_doc.get('solver', {})
            if gs_doc.get('schema_version') != 2 or solver.get('method') != 'smooth' or solver.get('smoother') != 'gauss_seidel':
                raise ValueError('this experiment requires the selected GS route')
            ext = pressure_doc.get('extended_solver', {})
            if ext.get('method') != 'bicgstab_amg' or ext.get('provider') != 'linglongamg':
                raise ValueError('this experiment requires BiCGStab/LinglongAMG')
            for name in ['library_path', 'config_path']:
                if not isinstance(ext.get(name), str) or not Path(ext[name]).is_absolute():
                    raise ValueError('explicit absolute dependency paths required')
            provider_path, provider_bytes, provider_doc = document(ext['config_path'])
            self.provider = self.snapshot('provider', provider_doc)
            resolved = copy.deepcopy(pressure_doc)
            resolved['extended_solver']['config_path'] = self.provider
            self.gs = self.snapshot('GS', gs_doc)
            self.pressure = self.snapshot('pressure', resolved)
            self.metadata['inputs'] = {name: {'path': str(path), 'bytes': len(data),
                'sha256': hashlib.sha256(data).hexdigest()} for name, path, data in [
                    ('GS', gs_path, gs_bytes), ('pressure', pressure_path, pressure_bytes),
                    ('provider', provider_path, provider_bytes)]}
            self.metadata['rewritten_options'] = ['extended_solver.config_path only']
        except BaseException:
            self.close()
            raise

    def snapshot(self, name, value):
        encoded = (json.dumps(value, allow_nan=False, separators=(',', ':')) + '\n').encode()
        fd = os.memfd_create('xk2647-' + name, os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING)
        self.fds.append(fd)
        remaining = memoryview(encoded)
        while remaining:
            count = os.write(fd, remaining)
            if not count:
                raise OSError('configuration snapshot short write')
            remaining = remaining[count:]
        os.lseek(fd, 0, os.SEEK_SET)
        seals = fcntl.F_SEAL_WRITE | fcntl.F_SEAL_GROW | fcntl.F_SEAL_SHRINK | fcntl.F_SEAL_SEAL
        fcntl.fcntl(fd, fcntl.F_ADD_SEALS, seals)
        assert fcntl.fcntl(fd, fcntl.F_GET_SEALS) == seals
        path = '/proc/self/fd/' + str(fd)
        assert json.loads(Path(path).read_bytes()) == value
        self.metadata.setdefault('snapshots', {})[name] = {'descriptor': fd,
            'memory_path': path, 'bytes': len(encoded), 'seals': seals,
            'sha256': hashlib.sha256(encoded).hexdigest()}
        return path

    def environment(self):
        return {**os.environ, 'RESIDENT_GS_CONFIG': self.gs,
            'RESIDENT_PRESSURE_CONFIG': self.pressure,
            'RESIDENT_PROVIDER_CONFIG': self.provider,
            'RESIDENT_CONFIG_PRELOADED': '1'}

    def close(self):
        for fd in self.fds:
            os.close(fd)
        self.fds.clear()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gs', type=Path, required=True)
    parser.add_argument('--pressure', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('child command required after --')
    with PreparedConfigurations(args.gs, args.pressure) as prepared:
        with args.metadata.open('x') as output:
            json.dump(prepared.metadata, output, indent=2)
            output.write('\n')
        return subprocess.run(command, env=prepared.environment(), pass_fds=tuple(prepared.fds)).returncode


if __name__ == '__main__':
    raise SystemExit(main())
