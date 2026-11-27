# 构建指南 (BUILDING)

本文说明如何基于上游 **OPM Flow 2026.04** 应用本仓库的 overlay 与补丁，并按两条公开 SDK 路线分别编译。

---

## 1. 通用前置依赖

在任意 GPU 路线编译前，先准备标准 OPM Flow 2026.04 构建环境：

- Linux x86_64，GCC（支持 C++20），CMake ≥ 3.23，Ninja
- 上游依赖：Dune 2.11（`dune-common`, `dune-geometry`, `dune-grid`, `dune-istl`）、`opm-common` 2026.04、`opm-grid` 2026.04、Boost、SuiteSparse/UMFPACK、Zoltan、METIS、Open MPI
- 未修改的上游 `opm-simulators` 2026.04 源码目录（`OPM_SOURCE_ROOT`）
- CUDA Toolkit（`CUDA_HOME`，RTX 5090 使用 `--cuda-arch=120`；XianzhiKrylov 26.48 发行包需要 CUDA 13 runtime `libcudart.so.13`）

---

## 2. 按路线准备 SDK 与环境变量

### 2.1 `linglongamg`（无需 XianzhiKrylov）

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_ROOT=/path/to/LinglongAMG
# 可选覆盖：
# export LINGLONGAMG_INCLUDE_DIR=/path/to/LinglongAMG/include
# export LINGLONGAMG_LIBRARY=/path/to/LinglongAMG/lib/liblinglongamg.so
# export LINGLONGAMG_LICENSE=/path/to/linglongamg.license.key
```

> 若使用了不支持 `linglongamg_solve` 后隐式同步消费结果的早期 26.46 SDK，可按需编译 `overlay/linglongamg/compat/linglong2646_compat.c` 作为 `dlopen` 包装层，并在运行时设置 `LINGLONGAMG_2646_LIBRARY` 指向底层 SDK 动态库。

### 2.2 `xianzhikrylov`（无需 LinglongAMG；使用公开 `XianzhiKrylovAPI.h`）

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda-13.x
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov-26.48
# 确认包含 include/XianzhiKrylovAPI.h、include/XianzhiKrylovData.h 与 lib/libXianzhiKrylov.so
# 可选覆盖：
# export XIANZHIKRYLOV_INCLUDE_DIR=/path/to/XianzhiKrylov-26.48/include
# export XIANZHIKRYLOV_LIBRARY=/path/to/XianzhiKrylov-26.48/lib/libXianzhiKrylov.so
# export XIANZHIKRYLOV_LICENSE=/path/to/xianzhikrylov.license.key
```

---

## 3. 执行配置与构建

```bash
# 如需在无现成 OPM 依赖环境时先做预检：
bash scripts/check-env.sh --route=<linglongamg|xianzhikrylov>

# 1) 准备带补丁的独立源码副本（默认在 work/source-<route>/）
bash scripts/prepare-source.sh --route=<linglongamg|xianzhikrylov>

# 2) 配置 CMake（如有额外 Spack/系统依赖目录，通过 --prefix-path 传入）
bash scripts/configure.sh \
  --route=<linglongamg|xianzhikrylov> \
  --prefix-path="/path/to/opm-deps" \
  --cuda-arch=120

# 3) 编译 flow 并自动运行 ldd 隔离校验
bash scripts/build.sh --jobs=16
```

`configure.sh` 默认关闭无关上游可选项（`-DUSE_GPU_BRIDGE=OFF -DUSE_HYPRE=OFF -DUSE_AMGX=OFF -DUSE_OPENCL=OFF -DBUILD_TESTING=OFF -DBUILD_FLOW_VARIANTS=OFF`），并按路线向 `CMakeLists.txt` 传入 `-DOPM_SOLVED_TECH_ROUTE=<route>`：

- `--route=linglongamg`：定义 `OPM_SOLVED_TECH_ENABLE_LINGLONGAMG=1` 并链接 `dl`。
- `--route=xianzhikrylov`：定义 `OPM_SOLVED_TECH_ENABLE_XIANZHIKRYLOV=1`，包含 `${XIANZHIKRYLOV_INCLUDE_DIR}` 并链接 `libXianzhiKrylov.so`、`CUDA::cudart` 与 `dl`。

> 环境提示：若 `--prefix-path` 中混有旧版 CUDA 头文件（例如较早版本的 Spack 视图），补丁已通过 `target_include_directories(opmsimulators BEFORE SYSTEM PRIVATE ${CUDAToolkit_INCLUDE_DIRS})` 优先使用 `CUDA_HOME` 的头文件；当宿主机使用 GCC 15 编译 `opm-simulators` 自带 `gpuistl` `.cu` 文件时，如遇 Spack 视图旧头文件拦截，可在配置前额外设置 `export CUDAFLAGS="-allow-unsupported-compiler -I$CUDA_HOME/targets/x86_64-linux/include"`。
