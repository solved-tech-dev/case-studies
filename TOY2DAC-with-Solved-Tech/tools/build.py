#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Build the Marmousi 3D inversion cases with PanguLU 26.46."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pangulu", type=Path, required=True, help="PanguLU 26.46 SDK root")
    parser.add_argument("--cuda", type=Path, required=True, help="CUDA 12.9 toolkit or runtime development root")
    parser.add_argument("--cuda-nvcc-include", type=Path, help="additional CUDA headers for split installations")
    parser.add_argument("--output", type=Path, default=ROOT / "_work/build", help="new build directory")
    args = parser.parse_args()
    if sys.platform != "linux":
        parser.error("use Linux with GCC and G++")
    sdk, cuda, output = (p.resolve() for p in (args.pangulu, args.cuda, args.output))
    sdk_lib = sdk / "lib"
    cuda_lib = next((cuda / d for d in ("lib64", "lib", "targets/x86_64-linux/lib")
                     if (cuda / d / "libcudart.so.12").is_file()), None)
    header = sdk / "include/pangulu.h"
    for path in (header, sdk_lib / "libpangulu.so", sdk_lib / "libpangulu_platform.so",
                 cuda / "include/cuda_runtime_api.h"):
        if not path.is_file():
            parser.error(f"missing {path}")
    if cuda_lib is None:
        parser.error(f"libcudart.so.12 is missing under {cuda}")
    if not re.search(r'^\s*#\s*define\s+PANGULU_VERSION_STRING\s+"26\.46"\s*$', header.read_text(), re.M):
        parser.error("use the PanguLU 26.46 SDK")
    output.mkdir(parents=True)
    (output / "bin").mkdir()
    adapter = output / "pangulu_fwi_adapter.o"
    app = ROOT / "apps/marmousi3d"
    includes = [f"-I{sdk / 'include'}", f"-I{cuda / 'include'}", f"-I{ROOT / 'adapter'}"]
    if args.cuda_nvcc_include:
        includes.append(f"-I{args.cuda_nvcc_include.resolve()}")
    subprocess.run(["g++", "-O3", "-std=c++17", *includes, "-c",
                    str(ROOT / "adapter/pangulu_fwi_adapter.cpp"), "-o", str(adapter)], check=True)
    for grid in json.loads((ROOT / "cases/fwi-grids.json").read_text()).values():
        h = grid["spacing_m"]
        model = output / f"model-h{h}.o"
        binary = output / "bin" / f"marmousi3d_fwi_h{h}"
        grid_flag = f"-DFWI_SPACING={h}"
        subprocess.run(["gcc", "-O3", "-std=c11", grid_flag, "-c",
                        str(app / "model.c"), "-o", str(model)], check=True)
        subprocess.run(["g++", "-O3", "-std=c++17", "-fopenmp", grid_flag, f"-I{ROOT / 'adapter'}",
                        str(app / "fwi.cpp"), str(model), str(adapter), f"-L{sdk_lib}", "-lpangulu",
                        "-Wl,--no-as-needed", "-lpangulu_platform", "-Wl,--as-needed",
                        f"-L{cuda_lib}", "-l:libcudart.so.12", "-ldl", "-lpthread", "-lm",
                        "-o", str(binary)], check=True)
    library_paths = [str(sdk_lib), str(cuda_lib)]
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = ":".join(filter(None, [*library_paths, env.get("LD_LIBRARY_PATH")]))
    for binary in sorted((output / "bin").iterdir()):
        result = subprocess.run(["ldd", "-r", str(binary)], env=env, text=True, capture_output=True)
        report = result.stdout + result.stderr
        log = output / f"loader-{binary.name}.txt"
        log.write_text(report)
        if result.returncode or "not found" in report or "undefined symbol" in report:
            raise SystemExit(f"unresolved runtime dependencies; see {log}")
        print(binary)
    (output / "runtime.json").write_text(json.dumps({"library_paths": library_paths}, indent=2) + "\n")


if __name__ == "__main__":
    main()
