#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: verify-public-tree.sh [OPTIONS]

Scan OPM-with-Solved-Tech/ for forbidden binaries, closed SDKs, license files,
secrets, private absolute paths, oversized files, and broken Markdown links.

Options:
  -h, --help            Show this help message and exit
EOF
}

for arg in "$@"; do
    case "$arg" in
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

python3 - "$package_dir" <<'PY'
import os
import re
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
errors = []

forbidden_exts = {
    ".so", ".a", ".o", ".obj", ".cubin", ".ptx", ".exe",
    ".lic", ".key", ".receipt",
    ".unsmry", ".smspec", ".esmry", ".dbg", ".prt", ".unrst", ".egrid", ".init", ".rft",
}
forbidden_dirs = {"sdk", "third_party", "lib", "bin", "license", "work", "build", "install", "__pycache__"}
max_bytes = 512 * 1024

private_prefixes = (
    "/ho" + "me/",
    "/mnt/" + "NAS",
    "/wo" + "rk2/",
    "/opt/" + "spack/",
    "/var/" + "tmp/",
)
private_path_re = re.compile(
    r"(?:^|[\s\"'=:,;])(" + "|".join(re.escape(p) for p in private_prefixes) + r")"
)
secret_prefixes = ("ghp_", "gho_", "github_pat_", "AKIA")
secret_re = re.compile(
    r"(-----BEGIN [A-Z ]*PRIVATE KEY-----|"
    r"(?:ghp_|gho_|github_pat_)[A-Za-z0-9_]{20,}|AKIA[0-9A-Z]{16})"
)
md_link_re = re.compile(r"\[[^\]]+\]\(([^)]+)\)")

file_count = 0
for dirpath, dirnames, filenames in os.walk(root):
    dirnames[:] = [d for d in dirnames if d != ".git"]
    rel_dir = Path(dirpath).relative_to(root)
    for d in list(dirnames):
        if d in forbidden_dirs:
            errors.append(f"Forbidden directory present: {rel_dir / d}")
    for fname in filenames:
        if fname == ".gitkeep":
            errors.append(f"Placeholder .gitkeep must be removed: {rel_dir / fname}")
        fpath = Path(dirpath) / fname
        rel_path = fpath.relative_to(root)
        file_count += 1
        lower = fname.lower()
        if ".so." in lower or any(lower.endswith(ext) for ext in forbidden_exts) or lower.endswith(".license.key"):
            errors.append(f"Forbidden file extension: {rel_path}")
        size = fpath.stat().st_size
        if size > max_bytes:
            errors.append(f"File exceeds {max_bytes} bytes ({size}): {rel_path}")
        raw = fpath.read_bytes()
        if raw.startswith(b"\x7fELF"):
            errors.append(f"ELF binary detected: {rel_path}")
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            errors.append(f"Non-UTF8 binary file detected: {rel_path}")
            continue
        if secret_re.search(text):
            errors.append(f"Potential secret/key material detected in: {rel_path}")
        for idx, line in enumerate(text.splitlines(), 1):
            m = private_path_re.search(line)
            if m:
                errors.append(f"Private absolute path '{m.group(1)}' in {rel_path}:{idx}")
        if fpath.suffix == ".md":
            for target in md_link_re.findall(text):
                if target.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                target_clean = target.split("#", 1)[0]
                if not target_clean:
                    continue
                resolved = (fpath.parent / target_clean).resolve()
                if not resolved.exists():
                    errors.append(f"Broken relative markdown link '{target}' in {rel_path}")

if errors:
    for err in errors:
        print(f"FAIL: {err}", file=sys.stderr)
    sys.exit(1)

print(f"verify-public-tree: OK ({file_count} files checked, 0 violations)")
PY
