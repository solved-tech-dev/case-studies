# Build and run

[English](BUILD-AND-RUN.md) · [简体中文](BUILD-AND-RUN.zh-CN.md)

Use Linux x86-64, GCC, G++ with C++17 and OpenMP, Python 3.9 or newer, a compatible NVIDIA driver, CUDA 12.9 development files, and the PanguLU 26.46 SDK. The application runs in one process on one GPU. Obtain the SDK and GPU license from [Solved Tech](https://solved-tech.com/).

Set the installation paths from the project directory and configure the license as described by the SDK.

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9

python3 tools/build.py --pangulu "$PANGULU_ROOT" --cuda "$CUDA_HOME"

for h in 85 75 65; do
    python3 tools/run_case.py "marmousi3d-h${h}" --output "_work/runs/marmousi3d-h${h}"
done
```

The default build directory is `_work/build`. To change it, pass `--output` when building and `--build` when running. For split CUDA installations, pass the runtime development directory to `--cuda` and the additional header directory to `--cuda-nvcc-include`. The builder writes a dynamic-library check for each grid. On a cluster, run inside a GPU allocation and supply the required SDK runtime libraries through `LD_LIBRARY_PATH`.

Use a new output directory for each run. The runner generates observations at 3 and 5 Hz using the truth model and then inverts from the smoothed starting model. Use `--prepare-only` to copy the case configuration and input without running the solver.

| File | Contents |
| --- | --- |
| `case.json` | Grid, frequencies, matrix order, and COO storage count |
| `observations/vp-true.bin`, `observations/vp-initial.bin`, `vp-final.bin` | True, initial, and recovered P-wave velocity in m/s. Headerless little-endian float64 node values with x varying fastest. |
| `observations/data-{0,1}.bin`, `predicted-{0,1}.bin` | Three-component receiver displacement at 3 and 5 Hz, stored receiver by receiver. Each component contains its real and imaginary parts as little-endian float64 values. |
| `generate.*.log`, `invert.*.log`, `validation.json` | Generation and inversion logs, validation, iteration history, and native phase timings |

The runner checks every original COO residual (finite and ≤ 1e-8), the forward receiver traces, a finite-difference gradient check, 12 velocity updates, and the waveform and velocity errors. Core time sums native analysis, factorization, and solve during optimization. Observation generation, gradient and forward checks, and application wall time are tracked separately.
