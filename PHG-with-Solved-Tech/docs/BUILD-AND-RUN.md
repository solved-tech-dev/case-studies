# Build and run

[简体中文](BUILD-AND-RUN.zh-CN.md)

Use Linux x86-64 with GNU C, C++, and Fortran compilers, Make, Expat development files, BLAS/LAPACK, and Python 3.12 or later. Obtain PanguLU 26.44 with ABI 2 for CUDA 12.9 and a license for the target NVIDIA GPU. Each case uses one PHG process and one GPU.

```bash
python3 tools/prepare.py
export PANGULU_ROOT=/path/to/pangulu-26.44
export CUDA_HOME=/usr/local/cuda-12.9
python3 tools/build.py --pangulu "$PANGULU_ROOT" --cuda "$CUDA_HOME"
```

The preparation tool downloads the two official archives pinned in `upstream.lock.json`. For offline setup, place the archives in `_work/downloads`. The default build tree is `_work/clean`. To build elsewhere, pass the same `--workspace` to the prepare, build, and run commands. Binaries and build logs are under `bin/` and `logs/` in that tree.

Use `--blas` and `--lapack` for custom link flags, such as `--blas="-L/path/to/openblas/lib -lopenblas"`. Use `--cuda-include`, `--cuda-nvcc-include`, and `--cuda-lib` for split CUDA installations. Configure other shared-library dependencies as described by the SDK.

Run from the project root on a GPU node. Use a new output directory for each run.

```bash
export PANGULU_LICENSE_PATHS=/path/to/your-gpu.license.key
python3 tools/run_case.py solid-adaptive-stress --output _work/runs/adaptive --vtk
python3 tools/run_case.py solid-static-p2 --output _work/runs/static-p2 --vtk
python3 tools/run_case.py interface-sphere-dg1 --output _work/runs/interface
```

The default profile uses the full case. Use `--profile connection-test` for a small grid and `--solver-threads` to set the PanguLU host thread count, which defaults to 40. `--vtk` enables solid displacement visualization. The interface case uses its native VTK output. `run-report.json` records matrix sizes, original-system residuals, native phase times, and field files.
