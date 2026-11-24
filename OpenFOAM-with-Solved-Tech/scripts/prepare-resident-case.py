#!/usr/bin/env python3
"""Stage an already prepared serial case while preserving its cell/field ordering."""
import argparse
import json
from pathlib import Path
import shutil

from run_resident import validate_gs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--gs-config', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=5)
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    destination = args.destination.absolute()
    if not 1 <= args.steps <= 2147483647:
        parser.error('--steps must be a positive int32 value')
    if destination.exists() or destination.is_symlink():
        parser.error('destination already exists')
    if destination.resolve().is_relative_to(source):
        parser.error('destination must be outside the input case')
    for name in ('U', 'p', 'k', 'omega', 'nut'):
        if not any((source / '0' / (name + suffix)).is_file() for suffix in ('', '.gz')):
            parser.error('missing serial initial field 0/' + name)
    for name in ('points', 'faces', 'owner', 'neighbour', 'boundary'):
        if not any((source / 'constant/polyMesh' / (name + suffix)).is_file() for suffix in ('', '.gz')):
            parser.error('missing serial mesh file ' + name)
    for name in ('fvSchemes', 'fvSolution'):
        if not (source / 'system' / name).is_file():
            parser.error('missing system/' + name)
    for name in ('transportProperties', 'turbulenceProperties'):
        if not (source / 'constant' / name).is_file():
            parser.error('missing constant/' + name)
    gs_path = args.gs_config.resolve(strict=True)
    validate_gs(gs_path)
    destination.mkdir(parents=True)
    shutil.copytree(source / 'system', destination / 'system')
    (destination / '0').symlink_to((source / '0').resolve(), target_is_directory=True)
    (destination / 'constant').mkdir()
    for item in (source / 'constant').iterdir():
        if item.is_file():
            shutil.copy2(item, destination / 'constant' / item.name)
    (destination / 'constant/polyMesh').symlink_to((source / 'constant/polyMesh').resolve(), target_is_directory=True)
    # Input fields stay at time 0; runtime controls are local to this staged case.
    (destination / 'system/controlDict').write_text('''FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object controlDict;
}
application rawResident;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime %d;
deltaT 1;
adjustTimeStep no;
writeControl timeStep;
writeInterval 1000000000;
runTimeModifiable no;
libs ();
functions {};
''' % args.steps)
    shutil.copy2(gs_path, destination / 'gs-config.json')
    (destination / 'input-source.json').write_text(json.dumps({
        'serial_input': str(source), 'gs_configuration': str(gs_path),
        'outer_steps': args.steps, 'input_mesh_and_fields': 'read-only references; preserve source',
    }, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
