# Build and run

[English](BUILD-AND-RUN.md) | [简体中文](BUILD-AND-RUN.zh-CN.md) · [Home](../README.md)

Run the commands below from `Calculix/`. On a cluster, obtain a GPU allocation before running the examples.

## 1. Requirements

| Component | Requirement |
|---|---|
| Platform | Linux x86_64 |
| Build tools | GCC, G++, GFortran, GNU Make, ar, patch; Python 3.12 or later |
| PanguLU | 26.42, ABI2, CUDA 12.9; SDK headers and shared libraries |
| GPU runtime | Compatible NVIDIA driver and CUDA runtime development files |
| Numerical libraries | LP64 BLAS/LAPACK, sequential ARPACK, SPOOLES 2.2 headers and sequential library |
| Other library | zlib development files |
| Mesh generation | NumPy, only for `tools/generate_mesh.py` |

The reference runtime configuration is one A100 40 GB GPU, 40 host CPU cores, and 100 GiB host memory. The runner requires at least 20 GiB of free disk space.

CalculiX uses SPOOLES helper routines for constraint processing, and ARPACK is a build dependency for its associated modules. The main linear-system solver entry is provided by the PanguLU adapter.

To purchase PanguLU and obtain a license, visit the [company website](https://solved-tech.com/). Install the SDK and license separately, using a version compatible with the interface above.

## 2. Build CalculiX with the adapter

Replace the example paths with your installed dependency locations:

```bash
python3 tools/build.py \
  --pangulu-root /path/to/external/pangulu \
  --cuda-root /path/to/external/cuda \
  --spooles-root /path/to/external/SPOOLES.2.2 \
  --arpack-libs "/path/to/external/libarpack.a" \
  --blas-libs "-L/path/to/external/openblas/lib -lopenblas" \
  --workspace _work/build --jobs 16
```

The script extracts `provenance/upstream/ccx_2.23.src.tar.bz2` into the workspace, compiles CalculiX and the adapter, and creates `_work/build/ccx_pangulu`. It uses the default 32-bit CalculiX integer interface and LP64 numerical libraries.

| Option | Purpose |
|---|---|
| `--pangulu-root` | Directory containing `include/pangulu.h` and PanguLU libraries under `lib/` |
| `--cuda-root` | CUDA directory containing `include/` and `lib64/` or `lib/` |
| `--spooles-root` | Directory containing SPOOLES headers and `spooles.a` |
| `--arpack-libs` | Linker flags or archive path for sequential ARPACK |
| `--blas-libs` | Linker flags for LP64 BLAS/LAPACK |
| `--jobs` | Number of concurrent compiler processes |
| `--timing` | Enable Newton-stage timing output |
| `--resume` | Resume a workspace with the same dependency paths and options |

For split CUDA installations, repeat `--cuda-include` to add header directories and use `--cuda-lib` for the runtime library directory. For oneMKL, provide GNU Fortran-compatible LP64/OpenMP linker options through `--blas-libs`.

Compiler logs, `link.log`, and `build-report.json` are written to the workspace. Use a separate workspace when changing build settings.

## 3. Configure the runtime

Configure the license location and the search path for external runtime libraries:

```bash
export PANGULU_LICENSE_PATHS=/path/to/your/authorized-license.key
export LD_LIBRARY_PATH=/path/to/external/runtime-libraries:$LD_LIBRARY_PATH
```

The executable includes local runtime search paths for the PanguLU and CUDA directories selected during compilation. BLAS, Fortran, and other runtime libraries may require additional paths.

## 4. Run a small example

```bash
python3 -m unittest discover -s tests -v
mkdir -p _work/runs
python3 tools/run_case.py dam --profile small \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam-small
```

The `small` profile uses the upstream mesh with temperature, heat-flux, and balance output enabled. Replace `dam` with another case name to check its small model.

## 5. Run the supplied refined models

```bash
python3 tools/run_case.py dam \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam-full
```

Available case names are `dam`, `ax6ht`, `capacitor`, `oneel20fi`, and `pan`. The runner decompresses the input model and launches CalculiX. Each run needs a new output directory.

Use an exclusive GPU and run the examples sequentially. The runner uses 4 assembly/result threads, and the adapter uses 40 host threads. Application limits are 240 seconds for `small` and 600 seconds for refined models, followed by a 10-second termination grace period.

### Output files

| File | Contents |
|---|---|
| `model.inp` | Input model used by CalculiX |
| `model.frd`, `model.dat`, `model.sta` | Native fields, tabular results, and increment information |
| `application.stdout`, `application.stderr` | Application messages and solver timings |
| `application-exit.json` | Exit status and application elapsed time |
| `resources.txt` | GNU time CPU/memory statistics, when available |
| `actual.call*.csr.gz` | Complete matrix and right-hand side for each call |
| `log-validation.json` | Dimensions, solver calls, residuals, and analysis/refactorization sequence |
| `field-validation.json` | Temperature, heat-flux, boundary, and balance checks |
| `validation.json` | Combined run-check status |

To compare with a reference field on the same mesh, add `--reference-frd /path/to/reference.frd.gz`. See [Checking a calculation](VALIDATION.md) for the residual and physical-field checks.

## 6. Optional: generate a different mesh

The supplied refined models are ready to run. To generate meshes from the upstream examples:

```bash
python3 tools/generate_mesh.py dam --factor 32 --output _work/mesh-dam
python3 tools/generate_mesh.py ax6ht --factor 125 --output _work/mesh-ax6ht
python3 tools/generate_mesh.py capacitor --factors 3 16 400 --output _work/mesh-capacitor
python3 tools/generate_mesh.py oneel20fi --factor 44 --output _work/mesh-oneel20fi
python3 tools/generate_mesh.py pan --factor 16 --output _work/mesh-pan
```

The tool writes `model.inp` and mesh metadata to the selected directory. Refinement preserves the parent-element geometry and maps sets and boundary faces to their child elements. A modified mesh is a separate test configuration.

## 7. Optional: check coefficient updates

`tests/changed_values.cpp` solves two 3×3 systems with the same graph and different coefficients and right-hand sides. Both have the known solution (1,2,3).

```bash
c++ -O2 -std=c++17 -fopenmp tests/changed_values.cpp \
  _work/build/obj/ccx_pangulu.cpp.o \
  -L/path/to/external/pangulu/lib -lpangulu -lpangulu_platform \
  -L/path/to/external/cuda/lib64 -l:libcudart.so.12 -lz -ldl \
  -o _work/changed_values
mkdir _work/reuse-check
HULL_CANONICAL="$PWD/_work/reuse-check/canonical.csr.gz" \
HULL_ACTUAL_PREFIX="$PWD/_work/reuse-check/actual" \
HULL_CREATE_CANONICAL=1 _work/changed_values
```

The output reports one analysis, two numerical factorizations, and `CHANGED_VALUES_PASSED` when the solutions satisfy the check.
