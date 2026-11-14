#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Build against the exact PETSc used by the MOOSE executable."""

import argparse, os, shlex, subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument(
    "--petsc-dir",
    default=os.environ.get("PETSC_DIR"),
    required=not os.environ.get("PETSC_DIR"),
)
p.add_argument("--petsc-arch", default=os.environ.get("PETSC_ARCH", ""))
p.add_argument(
    "--pangulu-root",
    default=os.environ.get("PANGULU_ROOT"),
    required=not os.environ.get("PANGULU_ROOT"),
)
p.add_argument(
    "--cuda-root",
    default=os.environ.get("CUDA_HOME"),
    required=not os.environ.get("CUDA_HOME"),
)
p.add_argument("--cxx", default=os.environ.get("CXX", "mpicxx"))
p.add_argument("--output", default="build/libmoose_pangulu.so")
p.add_argument(
    "--native-test",
    action="store_true",
    help="Build the GPU reuse/block-RHS test executable instead",
)
p.add_argument(
    "--petsc-test",
    action="store_true",
    help="Build the PETSc callback/block-RHS test executable instead",
)
a = p.parse_args()
if a.native_test and a.petsc_test:
    p.error("Choose only one test target")
if a.output == "build/libmoose_pangulu.so":
    if a.native_test:
        a.output = "build/test-native"
    if a.petsc_test:
        a.output = "build/test-petsc"
root = Path(__file__).resolve().parents[1]
petsc = Path(a.petsc_dir).resolve()
arch = petsc / a.petsc_arch
pg = Path(a.pangulu_root).resolve()
cuda = Path(a.cuda_root).resolve()
out = Path(a.output).resolve()
out.parent.mkdir(parents=True, exist_ok=True)
includes = [
    petsc / "include",
    arch / "include",
    pg / "include",
    cuda / "include",
    root / "adapter",
]
clib = next(
    (
        x
        for x in [cuda / "lib64", cuda / "lib", cuda / "targets/x86_64-linux/lib"]
        if (x / "libcudart.so").exists() or (x / "libcudart.so.12").exists()
    ),
    None,
)
if clib is None:
    p.error("CUDA runtime library not found under --cuda-root")
libs = [arch / "lib", pg / "lib", clib]
cmd = (
    shlex.split(a.cxx)
    + ["-std=c++17", "-O2", "-Wall", "-Wextra"]
    + ([] if (a.native_test or a.petsc_test) else ["-fPIC", "-shared"])
)
cmd += [f"-I{x}" for x in includes] + [
    str(
        root
        / (
            "tests/test_native.cpp"
            if a.native_test
            else (
                "tests/test_petsc.cpp"
                if a.petsc_test
                else "integration/petsc/solver_pangulu.cpp"
            )
        )
    )
]
cmd += [f"-L{x}" for x in libs] + [f"-Wl,-rpath,{x}" for x in libs]
cmd += [
    "-lpetsc",
    "-Wl,--no-as-needed",
    "-lpangulu",
    "-lpangulu_platform",
    "-Wl,--as-needed",
    "-l:libcudart.so.12",
    "-ldl",
    "-o",
    str(out),
]
subprocess.run(cmd, check=True)
print(out)
