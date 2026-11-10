# Build and run

[English](BUILD-AND-RUN.md) · [简体中文](BUILD-AND-RUN.zh-CN.md)

Use a Linux x86-64 system with GCC, G++ with C++17 support, Python 3.9 or later, a compatible NVIDIA GPU driver, CUDA 12.9 development files and the PanguLU 26.46 SDK. The 3D application uses a single process. Obtain the SDK and a license for your GPU from [Solved Tech](https://solved-tech.com/).

Run from the project root with your installation paths. Install the license following the SDK instructions. The default license location used for the benchmark was `~/.pangulu/license.key`.

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9

python3 tools/build.py --pangulu "$PANGULU_ROOT" --cuda "$CUDA_HOME"

for h in 85 75 65; do
    python3 tools/run_case.py "marmousi3d-h${h}" --output "_work/runs/marmousi3d-h${h}"
done
```

The default build directory is `_work/build`. To use another directory, pass `--output` to the builder and `--build` to the runner. For a split CUDA installation, pass `--cuda` with the runtime development root and `--cuda-nvcc-include` with the additional header directory. Make the CUDA libraries required by the SDK available on `LD_LIBRARY_PATH` at build and run time. The builder records its loader check in `loader-marmousi3d_pangulu.txt`. On a cluster, run the application inside a GPU allocation.

Use a new output directory for each run. To stage the input model and case metadata without a build or GPU, add `--prepare-only`.

| Output | Contents |
| --- | --- |
| `case.json`, `metadata.txt` | Case parameters, grid dimensions and application residual |
| `model.vtk` | 3D P-wave and S-wave speeds in m/s |
| `field.vtk` | Real and imaginary displacement vectors and their combined magnitude |
| `field_complex.bin` | Full displacement field. Three little-endian 32-bit node counts precede six float64 values per node, ordered Re/Im for x, y and z. The x index varies fastest. |
| `midplane_y.csv` | Velocity and complex displacement on the middle crossline plane |
| `stdout.log`, `stderr.log`, `validation.json` | Application and solver logs, residual, field checks and native stage timings |

The runner checks the expected matrix order and COO entry count, a finite original-system relative residual at most `1e-8`, and the dimensions, size and finite values of the full binary displacement field. Core solver time sums native analysis, numeric factorization and solve. It excludes model generation, assembly, explicit transfers, output and validation. For a failed run, the original phase and status are available in the logs.
