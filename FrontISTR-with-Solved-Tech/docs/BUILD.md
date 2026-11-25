# 构建、运行与校验 · Build and run

在 Linux 上准备 GNU C/C++/Fortran、CMake、Git、Python 3、OpenMP、MPI，以及与 MPI 匹配的 MUMPS（双精度实数）、ScaLAPACK、BLAS 和 LAPACK。GPU 端使用有授权的 PanguLU 26.48 ABI 2、兼容的 NVIDIA 驱动和 CUDA 12 运行库与头文件。MUMPS is a build dependency because the optional integration uses its FrontISTR wrapper. The three example runs use one process. Original GPU measurements used an A100 PCIe 40 GB.

## 1. Prepare the upstream source

Run the following from the root of this source package. Choose a fresh sibling directory for the upstream checkout.

```bash
repo=$PWD
git clone https://github.com/FrontISTR/FrontISTR.git ../FrontISTR-5.9
export FRONTISTR_SRC="$(cd ../FrontISTR-5.9 && pwd)"
git -C "$FRONTISTR_SRC" checkout 8ad8993860151e7070f0ca6078b8a0157966561e
git -C "$FRONTISTR_SRC" apply --check "$repo/patches/frontistr-pangulu.patch"
git -C "$FRONTISTR_SRC" apply "$repo/patches/frontistr-pangulu.patch"
```

## 2. Build the adapter and application

将以下示例路径改成本机安装路径，并按 SDK 的说明配置授权和依赖。Keep the SDK's runtime and license configuration active. Use matching MPI compiler wrappers and MPI libraries throughout the build.

```bash
export PANGULU_ROOT=/path/to/pangulu-26.48
export CUDA_HOME=/usr/local/cuda
export MUMPS_ROOT=/path/to/mumps
export LD_LIBRARY_PATH="$PANGULU_ROOT/lib:$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
scripts/build-adapter.sh
export LD_LIBRARY_PATH="$repo/build:$LD_LIBRARY_PATH"
scripts/build-frontistr.sh
```

Split CUDA installations can set `CUDA_INCLUDE`, `CUDA_NVCC_INCLUDE` and `CUDA_LIB`. For non-default linear algebra installations, pass CMake cache entries to the second script, e.g. `-DSCALAPACK_LIBRARIES=/path/to/libscalapack.so`, `-DBLAS_LIBRARIES=/path/to/libopenblas.so` and `-DLAPACK_LIBRARIES=/path/to/libopenblas.so`. Include the corresponding directories on the runtime library path. The build uses `BUILD_JOBS` (default 4), `CXX` and `CMAKE` if supplied. Check `ldd build/frontistr/fistr1/fistr1` for missing dependencies before running.

## 3. Run the examples

在有 GPU 的计算节点上运行。复现推送配置时，为一个 MPI 进程分配一张 A100 和 40 个 CPU 核，再从包根目录执行。The example runner defaults to `OMP_NUM_THREADS=4` and `FRONTISTR_PANGULU_THREADS=40`. Existing environment settings take precedence. Set the latter to a positive integer within your CPU allocation. When calling the patched executable directly, the adapter inherits the caller’s OpenMP setting unless this solver-specific variable is set. Each call logs the entry, solver and restored settings. The configuration log reports the SDK header version and ABI, and `FRONTISTR_PROVIDER` identifies the loaded library.

Each command creates a new directory, unpacks and hashes its mesh, copies native control files, enables PanguLU, starts FrontISTR and checks results. A failed run leaves its logs for inspection. Use a new directory when retrying.

```bash
python3 scripts/run.py conrod runs/conrod --fistr1 build/frontistr/fistr1/fistr1
python3 scripts/run.py hinge  runs/hinge  --fistr1 build/frontistr/fistr1/fistr1
python3 scripts/run.py heat   runs/heat   --fistr1 build/frontistr/fistr1/fistr1
python3 scripts/validate.py heat runs/heat
```

## 4. Interpret the outputs

`application.log` contains one `FRONTISTR_SOLVE` record per linear solve. Its residual is computed against the assembled matrix and the original right-hand side, with threshold `1e-8`. For a zero right-hand side the denominator is one. `FRONTISTR_TIME` records labelled `scope=native` measure solver analysis, numerical (re)factorization and solve. Each model's `model.json` stores its decompressed mesh SHA-256 and final field reference norms. The validator checks the field names, widths, node IDs and finite values, then compares each full field's L2 norm to the reference (relative tolerance `1e-6`, absolute tolerance `1e-12`). This scalar regression check complements the per-system residual.

| Case | Final native field relative to the run directory | Field names |
|---|---|---|
| conrod | `fstrRES/STEP1/conrod.res.0.1` | DISPLACEMENT, NodalSTRESS, NodalMISES |
| hinge | `hinge.res.0.1` | DISPLACEMENT, NodalSTRESS, NodalMISES |
| heat | `block.res.0.1` | TEMPERATURE |

The tested PanguLU 26.48 installation used CUDA 12.9. Follow the SDK instructions to configure a license valid for the selected device. For a custom license location, set `PANGULU_LICENSE_PATHS` to your own authorized file. Status 4 from `pangulu_create` means the SDK did not accept a license. For builds on a login node without the NVIDIA driver libraries, the final executable link may require `-DCMAKE_EXE_LINKER_FLAGS=-Wl,--allow-shlib-undefined`. At run time the matching libraries must resolve on the compute node.
