#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Build the prepared SU2 8.5.0 tree and link XianzhiKrylov 26.45."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent.parent
SOURCE_CONTROL_FILES = ("meson.build", "Common/include/CConfig.hpp", "Common/src/CConfig.cpp",
                        "Common/include/parallelization/omp_structure.hpp")
INTERNAL_SYMBOLS = (
    "xk_internal_su2_htd_alloc_f64", "xk_internal_su2_htd_upload_f64",
    "xk_internal_su2_htd_alloc_u64", "xk_internal_su2_htd_upload_u64",
    "xk_internal_su2_htd_alloc_i32", "xk_internal_su2_htd_upload_i32",
    "xk_internal_su2_htd_free", "xk_internal_get_last_stats",
)


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run(command, cwd, env, log):
    log.with_suffix(".command.json").write_text(json.dumps(command, indent=2) + "\n")
    print(shlex.join(command), flush=True)
    with log.open("w") as output:
        subprocess.run(command, cwd=cwd, env=env, stdout=output,
                       stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True,
                        help="Prepared SU2 8.5.0 tree produced by prepare-source.sh.")
    parser.add_argument("--sdk-root", type=Path, default=ROOT / "sdk",
                        help="SDK with include/ and lib/, default: package sdk/.")
    parser.add_argument("--sdk-library", type=Path,
                        help="Defaults to SDK_ROOT/lib/libXianzhiKrylov.so.2.")
    parser.add_argument("--cuda-root", type=Path, required=True)
    parser.add_argument("--cc", default=os.environ.get("CC", "gcc"))
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--ninja", default="ninja")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--workspace", type=Path, default=ROOT / "build/xianzhikrylov")
    args = parser.parse_args()
    if sys.version_info < (3, 12):
        parser.error("Python 3.12 or newer is required.")
    if args.jobs < 1:
        parser.error("--jobs must be positive.")
    if not sys.platform.startswith("linux"):
        parser.error("This integration builds on Linux with the GNU linker.")

    sdk = args.sdk_root.resolve()
    cuda = args.cuda_root.resolve()
    source = args.source_dir.resolve()
    for relative in (*SOURCE_CONTROL_FILES, "adapter/SU2SolveWrapper.cpp",
                     "adapter/Xianzhi45Adapter.cpp", "adapter/Xianzhi45Adapter.hpp",
                     "su2preconfig.timestamp", "externals/meson/mesonbuild/__init__.py",
                     "externals/eigen/Eigen/Core", "externals/mel/mel.hpp"):
        if not (source / relative).is_file():
            parser.error(f"Missing prepared SU2 source file: {source / relative}")
    if not re.search(r"version:\s*['\"]8\.5\.0(?:\s|['\"])", (source / "meson.build").read_text()):
        parser.error("The prepared source must be SU2 8.5.0.")
    if 'GetXianzhiPreconditionerLevel' not in (source / "Common/include/CConfig.hpp").read_text():
        parser.error("Apply the integration overlay with scripts/prepare-source.sh first.")
    library = (args.sdk_library or sdk / "lib/libXianzhiKrylov.so.2").resolve()
    for path in (library, sdk / "include/XianzhiKrylovAPI.h",
                 sdk / "include/XianzhiKrylovData.h", cuda / "bin/nvcc"):
        if not path.is_file():
            parser.error(f"Missing dependency: {path}")
    package_library = sdk / "lib/libXianzhiKrylov.so.2"
    if not package_library.is_file() or package_library.resolve() != library:
        parser.error("Retain lib/libXianzhiKrylov.so.2 in the matching external SDK.")
    compiler = shutil.which(args.cc)
    cxx = shutil.which(args.cxx)
    ninja = shutil.which(args.ninja)
    if not all((compiler, cxx, ninja, shutil.which("nm"))):
        parser.error("Provide working gcc, g++, ninja and GNU nm executables.")
    exports = subprocess.check_output(["nm", "-D", "--defined-only", str(library)], text=True)
    exported_names = {line.split()[-1].split("@")[0] for line in exports.splitlines() if line.split()}
    missing = sorted(set(INTERNAL_SYMBOLS) - exported_names)
    if missing:
        parser.error("This adapter requires the matching 26.45 SDK build. Missing symbols: "
                     + ", ".join(missing))

    workspace = args.workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    with (workspace / ".build.lock").open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("Another build is using this workspace.")
        settings = {
            "source_dir": str(source),
            "source_control_sha256": {name: sha256(source / name) for name in SOURCE_CONTROL_FILES},
            "adapter_sha256": {name: sha256(source / "adapter" / name) for name in
                               ("SU2SolveWrapper.cpp", "Xianzhi45Adapter.cpp", "Xianzhi45Adapter.hpp")},
            "sdk_library": str(library), "sdk_root": str(sdk), "cuda_root": str(cuda),
            "cc": compiler, "cxx": cxx, "ninja": ninja,
        }
        settings_path = workspace / "build-settings.json"
        if settings_path.exists() and json.loads(settings_path.read_text()) != settings:
            parser.error("Build settings changed. Select a new --workspace.")
        settings_path.write_text(json.dumps(settings, indent=2) + "\n")
        runtime_sdk = workspace / "sdk"
        if runtime_sdk.exists() or runtime_sdk.is_symlink():
            if runtime_sdk.resolve() != sdk:
                parser.error("The build bundle references another SDK. Select a new workspace.")
        else:
            runtime_sdk.symlink_to(sdk, target_is_directory=True)

        compile_dir = workspace / "compile"
        install_dir = workspace / "install"
        logs = workspace / "logs"
        for directory in (logs, install_dir / "bin", workspace / "lib",
                          workspace / "objects", workspace / "generated"):
            directory.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env.update(CC=compiler, CXX=cxx, NINJA=ninja, CUDACXX=str(cuda / "bin/nvcc"),
                   CUDA_HOME=str(cuda), CUDA_ROOT=str(cuda), SU2_SOURCE_DIR=str(source),
                   PYTHONDONTWRITEBYTECODE="1")
        env["PATH"] = os.pathsep.join((str(Path(compiler).parent), str(Path(ninja).parent),
                                       str(cuda / "bin"), env.get("PATH", "")))
        runtime_dirs = [cuda / "lib64", library.parent]
        compiler_runtime = subprocess.check_output([cxx, "-print-file-name=libstdc++.so"], text=True).strip()
        if Path(compiler_runtime).is_file():
            runtime_dirs.append(Path(compiler_runtime).resolve().parent)
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            [str(path) for path in runtime_dirs] + [env.get("LD_LIBRARY_PATH", "")])
        if not (compile_dir / "build.ninja").exists():
            run([sys.executable, str(ROOT / "scripts/meson_driver.py"), "setup",
                 str(compile_dir), str(source), "--prefix=" + str(install_dir),
                 "--buildtype=release", "-Dcpu-arch=x86-64-v3", "-Dwith-mpi=disabled",
                 "-Dwith-omp=true", "-Denable-cuda=true", "-Denable-mixedprec=false",
                 "-Denable-autodiff=false", "-Denable-directdiff=false",
                 "-Denable-pywrapper=false"], source, env, logs / "configure.log")
        run([ninja, "-C", str(compile_dir), "-j" + str(args.jobs), "SU2_CFD/src/SU2_CFD"],
            source, env, logs / "compile.log")
        commands = json.loads(subprocess.check_output(
            [ninja, "-C", str(compile_dir), "-t", "compdb"], env=env, text=True))
        entry = next(item for item in commands
                     if item.get("file", "").endswith("linear_algebra/CSysSolve.cpp"))
        native_args = shlex.split(entry["command"])
        native_object = Path(entry["directory"]) / native_args[native_args.index("-o") + 1]
        symbols = subprocess.check_output(["nm", "-g", str(native_object)], text=True)
        candidates = [line.split()[-1] for line in symbols.splitlines()
                      if line.split() and line.split()[-1].startswith("_ZN9CSysSolveIdE5SolveE")]
        if len(candidates) != 1:
            raise RuntimeError(f"Expected one SU2 Solve symbol, found {candidates}")
        symbol = candidates[0]
        (workspace / "generated/wrapped_symbol.hpp").write_text(
            '#pragma once\n#define SU2_WRAP_SOLVE_SYMBOL "__wrap_' + symbol + '"\n'
            '#define SU2_REAL_SOLVE_SYMBOL "__real_' + symbol + '"\n')
        adapter = workspace / "lib/libSU2Xianzhi45.so"
        run([cxx, "-std=c++17", "-O2", "-march=x86-64-v3", "-fPIC", "-fno-fast-math",
             "-shared", str(source / "adapter/Xianzhi45Adapter.cpp"),
             "-I" + str(source / "adapter"), "-I" + str(sdk / "include"),
             str(library), "-ldl", "-pthread",
             "-Wl,-soname,libSU2Xianzhi45.so", "-Wl,-rpath,$ORIGIN/../sdk/lib",
             "-o", str(adapter)], workspace, env, logs / "adapter.log")
        wrapper = workspace / "objects/SU2SolveWrapper.o"
        wrapper_args = []
        skip = False
        for arg in native_args:
            if skip:
                skip = False
                continue
            if arg in ("-MQ", "-MF", "-o"):
                if arg == "-o":
                    wrapper_args += ["-o", str(wrapper)]
                skip = True
            elif arg not in ("-MD", "-ffast-math", "-fno-finite-math-only"):
                wrapper_args.append(str(source / "adapter/SU2SolveWrapper.cpp")
                                    if arg == entry["file"] else arg)
        wrapper_args += ["-fno-fast-math", "-I" + str(source / "adapter"),
                         "-I" + str(workspace / "generated"),
                         "-I" + str(source / "Common/include/linear_algebra")]
        run(wrapper_args, Path(entry["directory"]), env, logs / "wrapper.log")
        link_commands = subprocess.check_output(
            [ninja, "-C", str(compile_dir), "-t", "commands", "SU2_CFD/src/SU2_CFD"],
            env=env, text=True)
        link = shlex.split(link_commands.splitlines()[-1])
        binary = install_dir / "bin/SU2_CFD"
        link[link.index("-o") + 1] = str(binary)
        link.insert(1, str(wrapper))
        link += ["-Xlinker=--wrap=" + symbol, str(adapter), "-Xlinker=--enable-new-dtags",
                 "-Xlinker=-rpath,$ORIGIN/../../lib", "-Xlinker=-rpath,$ORIGIN/../../sdk/lib"]
        for directory in runtime_dirs:
            if directory == library.parent:
                continue
            link.append("-Xlinker=-rpath," + str(directory))
        run(link, compile_dir, env, logs / "link.log")
        loader = subprocess.check_output(["ldd", str(binary)], env=env, text=True)
        (logs / "ldd.log").write_text(loader)
        if "not found" in loader:
            raise RuntimeError("Unresolved runtime libraries. See logs/ldd.log.")
        settings.update(binary=str(binary), wrapped_symbol=symbol,
                        binary_sha256=sha256(binary))
        (workspace / "build-report.json").write_text(json.dumps(settings, indent=2) + "\n")
        print(f"Built {binary}")


if __name__ == "__main__":
    main()
