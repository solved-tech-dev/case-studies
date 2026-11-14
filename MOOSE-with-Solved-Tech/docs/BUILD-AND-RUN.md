# Build and run

## Prerequisites

Python 3.10 or newer is required for the preparation and validation scripts. Install their dependencies from `requirements.txt`.

Use Linux x86_64 with a supported NVIDIA GPU, a CUDA 12.9-compatible driver and complete CUDA toolkit, and **PanguLU 26.46, ABI 2**. This integration follows the updated 26.46 package distributed on 2026-11-12 (archive SHA-256 `6aeb44706e01827e56924e2ea4d97068276ddbf50418fdf03103f69582dfff29`). The package is a separately obtained dependency; no authenticated download link is embedded here.

Use a built MOOSE environment with real, double-precision PETSc and compatible MPI/libMesh. Development uses PETSc 3.25.4. The plugin must match the PETSc build linked by your application. CUDA include files, runtime, cuBLAS and other SDK dependencies must be visible to the compiler/loader. Obtain the device license separately and install it at `~/.pangulu/license.key` with access restricted to your account. Do not put the key inside this repository.

MOOSE dependency installation is documented by [MOOSE](https://mooseframework.inl.gov/getting_started/installation/index.html). This package does not alter an existing installation or automatically upgrade PETSc, MPI or CUDA.

## 1. Source and application

Run these commands from the package root. Replace `/path/to/...` with your installation paths.

```bash
python3 -m pip install -r requirements.txt
python3 tools/fetch_upstream.py _upstream/moose
export MOOSE_DIR="$PWD/_upstream/moose"
# Build/install compatible PETSc and libMesh following MOOSE's source instructions.
# Point to those exact libraries before building the application:
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-your-build  # leave empty for a prefix install
export LIBMESH_DIR=/path/to/libmesh/installed
make -C app -j 8 MOOSE_DIR="$MOOSE_DIR" METHOD=opt
```

`fetch_upstream.py` selects commit `07b49de938ed357c136b487078065406d838aa26` in a new checkout. It does not overwrite existing source. MOOSE dependency/submodule setup is a prerequisite, as in a normal MOOSE build. The small application in `app/` links the original Solid Mechanics and Heat Transfer modules. You can also use an existing real-FP64 MOOSE application that registers these modules and matches the input syntax.

## 2. Plugin

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/path/to/cuda-12.9
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$PANGULU_ROOT/lib:$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
python3 tools/build_plugin.py
```

The default compiler is `mpicxx`; use `--cxx /path/to/the/matching/mpicxx` if needed. `--petsc-dir`, `--petsc-arch`, `--pangulu-root` and `--cuda-root` can replace their environment variables. The output is `build/libmoose_pangulu.so`. When using a system that splits CUDA across packages, supply the additional include/library paths through your compiler wrapper and loader environment.

Do not mix PETSc headers from a different MOOSE installation. A `not found` dependency reported by `ldd build/libmoose_pangulu.so` must be resolved before running.

## 3. A complete case

Run on an allocated GPU. The commands below use one MPI rank; they do not launch `mpirun`.

```bash
export OMP_NUM_THREADS=8
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
python3 tools/run_case.py thermomech \
  --exe "$PWD/app/moosesolved-opt" \
  --profile small --output runs/thermomech-small
```

The runner creates a **new** output directory, prepares the mesh where needed, launches the complete original-physics application, checks convergence/physical times and fields, and parses PanguLU's independent full-CSR residuals. Outputs include `run.json`, mesh and application command records, `application.log`, an Exodus file, CSV output and `validation.json`.

Run the supplied enlarged grid by selecting `--profile million`. The million profiles require substantially more memory than the small checks; select resources appropriate to your installation. No scheduler, cluster name, private path, GPU model, performance target or automatic time limit is embedded in the runner.

All six cases use the same entry point:

```bash
for case in thermomech annular-torsion mug pinched-cylinder gap-heat hex20-heat; do
  python3 tools/run_case.py "$case" --exe "$PWD/app/moosesolved-opt" \
    --profile million --output "runs/$case-million"
done
```

This runs each application once, sequentially. `--prepare-only` creates the mesh and `run.json` without solving the physical problem. Mesh-only conversion still invokes MOOSE.

## 4. Reference result and tests

For a small grid, PETSc's built-in sequential LU requires no additional commercial solver:

```bash
python3 tools/run_case.py thermomech --exe "$PWD/app/moosesolved-opt" \
  --profile small --solver reference --reference-solver petsc \
  --output runs/thermomech-reference
python3 tools/validate.py runs/thermomech-small --reference runs/thermomech-reference
c++ -std=c++17 -O2 -Iadapter tests/test_sparse.cpp -o /tmp/moose-test-sparse
/tmp/moose-test-sparse
python3 tools/build_plugin.py --native-test --output build/test-native
./build/test-native
python3 tools/build_plugin.py --petsc-test --output build/test-petsc
./build/test-petsc -dll_prepend "$PWD/build/libmoose_pangulu.so"
```

The native test needs a licensed GPU. It checks explicit zeros becoming nonzero and vice versa, unchanged-matrix factor reuse, graph-change reanalysis and a true two-column RHS block. It is a correctness test, not a performance benchmark. Reference `petsc` is intended for small grids; other available PETSc LU providers can be selected explicitly.

See [validation](VALIDATION.md) for what each check establishes and [integration](INTEGRATION.md) for supported matrix types and timing definitions.
