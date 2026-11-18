# Code_Saturne 9.2.0 接入 LinglongAMG

在 Code_Saturne 的压力修正方程中，使用 LinglongAMG 26.46 的库内 GPU GMRES（广义最小残量法）和 AMG（代数多重网格）预条件。Code_Saturne 在 CPU 上组装方程，适配层把标量压力系统转换为 CSR（压缩稀疏行格式），上传到一张 GPU 求解。支持单 MPI 进程、FP64 和 32 位 CSR 索引。其他方程沿用算例的求解器设置。

| 位置 | 内容 |
| --- | --- |
| `overlay/adapter/` | 压力求解适配层和可选的原生压力计时组件 |
| `patches/` | 宿主补丁，让 Code_Saturne 以 host HYPRE ParCSR 组装压力矩阵 |
| `examples/` | BENCH_C016_03、mixing_tee 的 SDK 配置和 `cs_user_linear_solvers` 接入示例 |
| `scripts/` | 环境加载、编译、宿主接口检查和迁移安装的启动工具 |

## 准备 Code_Saturne

本包对应 [Code_Saturne v9.2.0，提交 5315e9f](https://github.com/code-saturne/code_saturne/commit/5315e9f08d37982213d7ee5efb603f89655d8b05)。使用 Linux、Bash、Git、GNU 构建工具和 Python 3，并按 [Code_Saturne 构建说明](https://github.com/code-saturne/code_saturne/blob/5315e9f08d37982213d7ee5efb603f89655d8b05/INSTALL.md)准备算例所需的 MPI、网格读取器及其他依赖。为宿主启用 CPU 版 HYPRE，参数为 `--with-hypre="$HYPRE_ROOT"`。HYPRE 与宿主使用同一套 MPI 和双精度构建。

在本包根目录执行以下命令，准备锁定版本的源码。应用补丁后再编译和安装 Code_Saturne。

```bash
bundle="$PWD"
git clone --branch v9.2.0 --depth 1 https://github.com/code-saturne/code_saturne.git /path/to/code_saturne
export CS_SOURCE_DIR=/path/to/code_saturne
git -C "$CS_SOURCE_DIR" apply --check "$bundle/patches/0001-linglongamg-host-parcsr.patch"
git -C "$CS_SOURCE_DIR" apply "$bundle/patches/0001-linglongamg-host-parcsr.patch"
```

补丁中的 `cs_sles_linglong_t` 是现有宿主识别矩阵格式用的内部名称。构建脚本从宿主源码读取该名称，并在安装的 `libsaturne.so` 中核对。对外的注册入口为 `cs_linglongamg_define`。HYPRE 在这条路线中负责提供已组装的矩阵，线性求解由 LinglongAMG 完成。

## 编译适配层

通过 SDK 交付渠道获取 LinglongAMG 26.46 / CUDA 13.1 版及当前 GPU 的授权。SDK 的 `lib/` 下需要 `linglongamg.h`、`liblinglongamg.so` 和 `delivery.receipt`。服务器环境使用 NVIDIA A100、CUDA Toolkit 13.2、GCC 15.2、Open MPI 5.0 和 HYPRE 3.1。路径替换为本机安装位置。

```bash
export CS_PREFIX=/path/to/patched-code-saturne-install
export HYPRE_ROOT=/path/to/cpu-hypre
export MPI_HOME=/path/to/matching-openmpi
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_ROOT=/path/to/LinglongAMG
export LINGLONGAMG_LICENSE_FILE=/path/to/device-authorization.key

bash scripts/build.sh
source scripts/environment.sh
```

编译器 `CXX` 默认是 `g++`，请使用与 Code_Saturne ABI 兼容的编译器。如果使用已经配置好头文件路径的 MPI 编译器包装器，可以将 `CXX` 设为该命令并省略 `MPI_HOME`。若动态库装在额外的目录，将目录以冒号连接后写入 `CS_LINGLONGAMG_DEP_LIBPATH`。输出默认为 `build/`，或由 `CS_LINGLONGAMG_BUILD` 指定。脚本会把 SDK 运行文件和授权复制到该本地构建目录，授权副本权限为 `0600`。

## 启用压力求解器

准备已可运行的 Code_Saturne 算例及网格。有现成的 `cs_user_linear_solvers` 时，把示例中的动态加载与注册逻辑合入该函数。新算例可直接使用 `examples/SRC/cs_user_parameters.cpp`，放在算例的 `SRC/` 下。函数会对 `CS_F_(p)->id` 注册压力求解器。

```bash
case_dir=/path/to/case
cp -n examples/SRC/cs_user_parameters.cpp "$case_dir/SRC/"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
"$CS_PREFIX/bin/code_saturne" run --case "$case_dir" --id linglongamg -n 1 --nt "$OMP_NUM_THREADS"
```

| 配置 | SDK 内部容差 | `CS_LINGLONGAMG_ROW_SCALE` |
| --- | --- | --- |
| `examples/linglongamg.json`（默认，BENCH_C016_03） | `5e-9` | `1` |
| `examples/linglongamg_mixing_tee.json` | `5e-6` | `0` |

表中为已有实验配置。运行 mixing_tee 时，把 `CS_LINGLONGAMG_CONFIG` 设为该 JSON 的绝对路径。行缩放默认关闭，BENCH_C016_03 需显式设 `CS_LINGLONGAMG_ROW_SCALE=1`。缩放同时作用于矩阵行和右端项，使当前矩阵的行范数与建层时一致。其他算例须根据压力容差和数值变化选择参数。

对于迁移过的宿主安装，可用 `python3 scripts/run_saturne.py run ...`。它通过 `CS_ORIGINAL_BUILD_ROOT` 和 `CS_BUILD_ROOT` 替换已迁移的路径，并可从 `CS_LINGLONGAMG_GCC_RUNTIME` 读取存放版本化 Fortran 运行库的目录。原有安装配置保留在安装前缀下，启动工具使用临时副本。

## 运行检查

稀疏结构不变时，适配层原位上传新矩阵数值，复用 GPU 缓冲区和 AMG 层级。结构变化时重新建层。每次求解下载解向量后，用当前未缩放矩阵在 CPU 上复算 `||b - A x||₂`，阈值为 Code_Saturne 传入的 `precision * r_norm`。超过阈值会终止计算。日志中的 `CS_LINGLONGAMG_DEVICE`、`CS_GPU_LINEAR` 和 `CS_LINGLONGAMG_SUMMARY` 分别记录设备、逐次残差和建层次数。

服务器原始 A100 实验在 BENCH_C016_03 和 mixing_tee 上各做了 20 个 CFD 外迭代步，两例各 40 次 LinglongAMG 压力求解均通过真实残差验收。源码包中的 JSON 来自这两条实验路线。计时组件 `libcs_pressure_timing.so` 用于单独测量原生 CPU 压力求解，通过 `LD_PRELOAD` 加载到原生路线的求解进程。

## 来源与许可

宿主补丁基于 Code_Saturne 的 `cs_sles_default.cpp`，上游版权归 EDF S.A.。上游代码采用 GPL-2.0-or-later。 [上游许可正文](overlay/COPYING)

LinglongAMG SDK、运行库、设备授权和算例网格由各自提供渠道交付。
