# Prepare, build and run

[简体中文](BUILD-AND-RUN.zh-CN.md) · [Solver interface](INTEGRATION.md)

Run the commands from the package root. Input preparation needs Python 3 and NumPy. The tested target is Linux x86-64, GCC with C++17/OpenMP, FEBio 4.13.0 with matching FECore headers and shared library, and a separately licensed PanguLU 26.46 CUDA 12.9 installation (ABI 2). The application runs used one A100 40 GB. The plugin sets 40 OpenMP threads during solver calls and restores the caller's setting afterwards.

## Generate the three input decks

The script verifies the pinned upstream file, subdivides the parent mesh, maps its boundary conditions, selects the plugin, and checks the output against the exact FEBio input used in the validated run. Select just one case by giving its name.

```bash
git clone --filter=blob:none --no-checkout https://github.com/febiosoftware/TestSuite.git TestSuite
git -C TestSuite sparse-checkout set --no-cone /Verify/bp10.feb /Verify/tu01.feb /Verify/mi06.feb
git -C TestSuite checkout --detach 986005fc1865f5de7f2c789b8153d6f8fe2c0acd
python3 -m pip install numpy
python3 scripts/prepare.py TestSuite bp10 tu01 mi06
```

## Build the plugin

Replace the example paths with the local installations. Point FEBIO_SRC at a directory containing FECore/ and use an executable matching those headers and libraries. Configure PanguLU's license and runtime environment as instructed by its supplier.

```bash
export FEBIO_SRC=/opt/FEBio/include
export FEBIO_LIB=/opt/FEBio/lib
export FEBIO_BIN=/opt/FEBio/bin/febio4
export PANGULU_ROOT=/opt/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9
export LD_LIBRARY_PATH="$FEBIO_LIB:$PANGULU_ROOT/lib:$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
bash scripts/build.sh
```

The output is `FEBio-with-Solved-Tech/build/pangulu.so`. For split CUDA installations, override CUDA_INCLUDE, CUDA_NVCC_INCLUDE and CUDA_LIB with the runtime headers, compiler headers and library directories respectively.

## Run and inspect

```bash
mkdir -p runs
bash scripts/run.sh bp10 runs/bp10
bash scripts/run.sh tu01 runs/tu01
bash scripts/run.sh mi06 runs/mi06
```

Each call creates a new output directory containing console.log, application.log and physical.xplt. The wrapper checks FEBio's normal termination, presence of the per-solve residual record, and a nonempty plot. Open the XPLT in FEBio Studio to inspect fields over the converged time steps. Each FEBIO_RESIDUAL line reports ‖Ax−b‖₂/‖b‖₂ on the original system, with an absolute residual for zero RHS. The plugin enforces a threshold of 10⁻⁸.

To check against an independently generated reference, use the same mesh, material, loading and field requests with the alternate solver. For example:

```bash
python3 scripts/compare_fields.py cases/tu01/tu01-8.feb \
    /path/to/reference/physical.xplt runs/tu01/physical.xplt
```

The reader supports the uncompressed FEBio 3.5 plot format used here. It loads both plots into memory, checks the mesh identities, requested variables and converged times, then checks every field component using |candidate−reference| / max(1, |reference|) ≤ 10⁻⁶. A failed check returns a nonzero exit code.
