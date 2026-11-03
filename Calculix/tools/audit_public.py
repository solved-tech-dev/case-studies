# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 解得科技
"""Check the public source tree, not private build/run workspaces."""
import gzip
import json
from pathlib import Path
import re
import tarfile

ROOT=Path(__file__).resolve().parents[1]
SKIP={'.git','_work','__pycache__'}
DENY={'.key','.pem','.p12','.pfx','.lic','.license','.o','.a','.so','.dll','.dylib','.exe','.sbatch','.slurm'}
# Construct these strings so that this checker does not flag its own source.
private_paths=['/'+x+'/' for x in ('Users','home')]
private_markers=['ssslab'+'-login','/mnt/'+'NAS03', 'BEGIN '+'PRIVATE KEY','BEGIN '+'OPENSSH PRIVATE KEY']

def main():
    errors=[];count=0;largest=0;total=0
    for p in ROOT.rglob('*'):
        rel=p.relative_to(ROOT)
        if any(x in SKIP for x in rel.parts):continue
        if p.is_symlink():errors.append(dict(file=str(rel),issue='symlink_not_allowed'));continue
        if not p.is_file():continue
        count+=1;size=p.stat().st_size;largest=max(largest,size);total+=size
        if p.suffix.lower() in DENY or p.name.startswith(('libpangulu','pangulu-','.env')) or p.name=='pangulu.h':
            errors.append(dict(file=str(rel),issue='private_or_generated_binary_file'))
        if size>=100*1024**2:errors.append(dict(file=str(rel),issue='exceeds_normal_GitHub_file_limit'))
        if p.suffix=='.png':continue
        if p.name=='ccx_2.23.src.tar.bz2':
            with tarfile.open(p,'r:bz2') as archive:
                for m in archive:
                    if m.name.startswith('/') or '..' in Path(m.name).parts or m.issym() or m.islnk():
                        errors.append(dict(file=str(rel),issue='unsafe_source_archive_member'))
            continue
        if p.suffix=='.gz':
            if not re.fullmatch(r'cases/(dam|ax6ht|capacitor|oneel20fi|pan)/input/model\.inp\.gz',rel.as_posix()):
                errors.append(dict(file=str(rel),issue='unexpected_compressed_payload'))
                continue
            stream=gzip.open(p,'rt',encoding='latin1')
        else:stream=p.open(encoding='utf8',errors='replace')
        with stream:
            for number,line in enumerate(stream,1):
                if any(x in line for x in private_paths+private_markers):
                    errors.append(dict(file=str(rel),line=number,issue='private_path_or_credential_marker'))
                    break
    report=dict(status='passed' if not errors else 'failed',files=count,total_bytes=total,
        largest_file_bytes=largest,errors=errors,
        limitations='Pattern/content check only; not a legal review or complete secret detector')
    print(json.dumps(report,indent=2));return 0 if not errors else 1

if __name__=='__main__':raise SystemExit(main())
