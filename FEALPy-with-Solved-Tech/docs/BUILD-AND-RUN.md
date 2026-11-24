# Build and run

[简体中文](BUILD-AND-RUN.zh-CN.md) · [Overview](../README.md) · [Solver interface](INTEGRATION.md)

## Environment

Use Linux x86-64, a C++17 compiler, a CUDA GPU and a valid PanguLU 26.46 license for that device. Install the CUDA toolkit matching your SDK (the case runs used the CUDA 12.9 edition on an A100 PCIe 40 GB). The full-size runs used jobs with 100 GiB of host memory allocated. This is a reference allocation, not a measured minimum.

Use Python 3.11 with a NumPy backend. Start in the extracted `FEALPy-with-Solved-Tech` directory. The following installs upstream in a sibling directory and places the Python environment outside the case package. On a headless Linux system, install the system OpenGL/GLU libraries required by the Gmsh Python distribution.

```bash
python3.11 -m venv ../fealpy-venv
source ../fealpy-venv/bin/activate
git clone https://github.com/weihuayi/fealpy.git ../fealpy-upstream
git -C ../fealpy-upstream checkout --detach e033533bd8e24f796f9c5284685b20a98b4d11cb
python -m pip install numpy==1.26.4 scipy==1.16.2
env -u WITH_PANGULU python -m pip install -e ../fealpy-upstream
```

## Build the adapter

Set the installation paths below to your own directories. Configure the PanguLU license as directed by the SDK vendor. For a license file at a custom location, export `PANGULU_LICENSE_PATHS` with its absolute path. Use `CUDA_VISIBLE_DEVICES` to select the device when necessary; the adapter opens visible device 0.

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/usr/local/cuda-12.9
export LD_LIBRARY_PATH="$PANGULU_ROOT/lib:$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
bash scripts/build.sh
```

The script writes `build/libfealpy_pangulu_real.so` and `build/libfealpy_pangulu_complex.so`. For a split CUDA installation, set `CUDA_INCLUDE`, `CUDA_NVCC_INCLUDE` and `CUDA_LIB` to the runtime headers, compiler headers and runtime library directory before building. Add the matching runtime directory and any SDK-specific dependencies to the loader path when running.

## Run and check the solution

The following commands start with the small meshes. Each program calls the actual GPU backend once. Omit `--small` to select the mesh documented in the case README.

```bash
python cases/elastic3d/run.py --small
python cases/maxwell2d/run.py --small
python cases/advdiff2d/run.py --small
```

To save coefficients from the three full-sized models:

```bash
python cases/elastic3d/run.py --output output/elastic3d.npy
python cases/maxwell2d/run.py --output output/maxwell2d.npy
python cases/advdiff2d/run.py --output output/advdiff2d.npy
```

The terminal reports the matrix order, CSR entry count, full relative residual and case-specific field error. A successful run ends with `validation=passed`. Output is optional; specify a separate `.npy` file per run to retain it. The file contains one NumPy vector in the original finite-element degree-of-freedom order. See the individual [elasticity](../cases/elastic3d/README.md), [Maxwell](../cases/maxwell2d/README.md) and [transport](../cases/advdiff2d/README.md) guides for the interpretation.
