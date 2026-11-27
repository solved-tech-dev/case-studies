# 构建指南 (BUILDING)

本文说明如何基于上游 **OPM Flow 2026.04** 应用本仓库的 overlay 与补丁，并按两条路线分别编译。

---

## 1. 通用前置依赖

在任一 GPU 路线编译前，先准备标准 OPM Flow 2026.04 构建环境：

- Linux x86_64，GCC（支持 C++20），CMake ≥ 3.23，Ninja
- 上游依赖：Dune 2.11（`dune-common`, `dune-geometry`, `dune-grid`, `dune-istl`）、`opm-common` 2026.04、`opm-grid` 2026.04、Boost、SuiteSparse/UMFPACK、Zoltan、METIS、Open MPI
- 上游 `opm-simulators` 2026.04 源码目录（`OPM_SOURCE_ROOT`）
- CUDA Toolkit（`CUDA_HOME`，RTX 5090 使用 `--cuda-arch=120`；XianzhiKrylov `26.48` 需 CUDA 13 runtime `libcudart.so.13`）

---

## 2. 按路线准备 SDK 与环境变量

### 2.1 `linglongamg`（压力子系统路线）

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_ROOT=/path/to/LinglongAMG
# 可选覆盖：
# export LINGLONGAMG_INCLUDE_DIR=/path/to/LinglongAMG/include
# export LINGLONGAMG_LIBRARY=/path/to/LinglongAMG/lib/liblinglongamg.so
# export LINGLONGAMG_LICENSE=/path/to/linglongamg.license.key
```

> 若使用早于 `26.46` 异步同步约定版本的底层库，可按需编译 `overlay/linglongamg/compat/linglong2646_compat.c` 作为 `dlopen` 包装层，并在运行时通过 `LINGLONGAMG_2646_LIBRARY` 指定底层动态库路径。

### 2.2 `xianzhikrylov`（全耦合路线）

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda-13.x
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov-26.48
# 需包含 include/XianzhiKrylovAPI.h、include/XianzhiKrylovData.h 与 lib/libXianzhiKrylov.so
# 可选覆盖：
# export XIANZHIKRYLOV_INCLUDE_DIR=/path/to/XianzhiKrylov-26.48/include
# export XIANZHIKRYLOV_LIBRARY=/path/to/XianzhiKrylov-26.48/lib/libXianzhiKrylov.so
# export XIANZHIKRYLOV_LICENSE=/path/to/xianzhikrylov.license.key
```

---

## 3. 执行配置与构建

```bash
# 1) 依赖与 SDK 预检
bash scripts/check-env.sh --route=<linglongamg|xianzhikrylov>

# 2) 准备独立源码副本并应用补丁与路线文件（默认位于 work/source-<route>/）
bash scripts/prepare-source.sh --route=<linglongamg|xianzhikrylov>

# 3) 配置 CMake（如有额外依赖目录，通过 --prefix-path 传入）
bash scripts/configure.sh \
  --route=<linglongamg|xianzhikrylov> \
  --prefix-path="/path/to/opm-deps" \
  --cuda-arch=120

# 4) 编译 flow 并执行动态库隔离校验
bash scripts/build.sh --jobs=16
```

`configure.sh` 默认关闭无关上游可选项（`-DUSE_GPU_BRIDGE=OFF -DUSE_HYPRE=OFF -DUSE_AMGX=OFF -DUSE_OPENCL=OFF -DBUILD_TESTING=OFF -DBUILD_FLOW_VARIANTS=OFF`），并按路线向 `CMakeLists.txt` 传入 `-DOPM_SOLVED_TECH_ROUTE=<route>`：

- `--route=linglongamg`：定义 `OPM_SOLVED_TECH_ENABLE_LINGLONGAMG=1` 并链接 `dl`。
- `--route=xianzhikrylov`：定义 `OPM_SOLVED_TECH_ENABLE_XIANZHIKRYLOV=1`，包含 `${XIANZHIKRYLOV_INCLUDE_DIR}` 并链接 `libXianzhiKrylov.so`、`CUDA::cudart` 与 `dl`。

> 编译提示：若 `--prefix-path` 依赖路径中包含不同版本的 CUDA 头文件，补丁已通过 `target_include_directories(opmsimulators BEFORE SYSTEM PRIVATE ${CUDAToolkit_INCLUDE_DIRS})` 优先使用 `CUDA_HOME` 头文件；若宿主机使用较新版本的 GCC 编译 OPM 自带 `.cu` 文件，可在配置前设置 `export CUDAFLAGS="-allow-unsupported-compiler -I$CUDA_HOME/targets/x86_64-linux/include"`。
