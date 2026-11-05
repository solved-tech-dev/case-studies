# OpenCAEPoroX的SPE10、HYPRE和LinglongAMG接入包

本包只保留以下内容：

- OpenCAEPoroX内部生成SPE10r1-r5压力矩阵所需代码。
- HYPRE的CPU/GPU求解器接入代码。
- LinglongAMGCUDA求解器接入代码。
- HYPRE和LinglongAMG的CMake链接配置。
- SPE10数据卡和Slurm测试脚本。

本包不包含HYPRE、LinglongAMG的源码、头文件、动态库、静态库或授权文件。运行时必须提供已经安装好的第三方库。

## 代码结构

```text
overlay/                        需要覆盖到原版源码树的文件
examples/linglongamg.json       LinglongAMG参数模板
scripts/prepare-source.sh       从干净源码树生成工作副本
scripts/build-hypre-cpu.sh      HYPRE CPU构建脚本
scripts/build-hypre-gpu.sh      HYPRE GPU构建脚本
scripts/build-linglongamg.sh    LinglongAMG构建脚本
scripts/test-hypre-cpu.sbatch   HYPRE CPU的r1-r5测试
scripts/test-hypre-gpu.sbatch   HYPRE GPU的r1-r5测试
scripts/test-linglongamg.sbatch LinglongAMG的r1-r5测试
scripts/test-all.sh             一次提交三组r1-r5测试
scripts/summarize-results.sh    汇总日志为CSV
```

## 准备源码

进入解压后的包目录，通过`OCP_UPSTREAM_DIR`指定未修改的OpenCAEPoroX源码树：

```bash
export OCP_UPSTREAM_DIR=../OpenCAEPoroX
bash scripts/prepare-source.sh
```

默认生成：

```text
work/OpenCAEPoroX
```

如果目标目录已经存在，脚本会停止，不会覆盖现有源码。

## 配置第三方库

包内不写死任何第三方库路径。编译和运行前，用户需要根据实际安装位置设置以下环境变量：

```text
HYPRE_CPU_CONFIG       与CPU节点匹配的HYPREConfig.cmake所在目录
HYPRE_GPU_CONFIG       GPU版HYPREConfig.cmake所在目录
LINGLONGAMG_ROOT       LinglongAMG发布包根目录
CPU_MPI_ROOT           CPU节点MPI安装目录
GPU_MPI_ROOT           GPU节点MPI安装目录
CUDA_ROOT              CUDA安装目录
OPENBLAS_LIB_DIR       运行时OpenBLAS库目录
BLAS_LIB               GPU节点BLAS库路径
LAPACK_LIB             GPU节点LAPACK库路径
CPU_BLAS_LIB           CPU节点BLAS库路径
CPU_LAPACK_LIB         CPU节点LAPACK库路径
METIS_ROOT             METIS头文件根目录
PARMETIS_ROOT          ParMETIS头文件根目录
```

`CMAKE_BIN`可选，默认使用`PATH`中的`cmake`。用户也可以直接提供`OPENBLAS_ROOT`，脚本会据此推导`OPENBLAS_LIB_DIR`及默认BLAS/LAPACK库路径。

环境变量示例：

```bash
export CPU_MPI_ROOT=...
export GPU_MPI_ROOT=...
export CUDA_ROOT=...
export METIS_ROOT=...
export PARMETIS_ROOT=...
export HYPRE_CPU_CONFIG=...
export HYPRE_GPU_CONFIG=...
export LINGLONGAMG_ROOT=...
export OPENBLAS_ROOT=...
export CPU_BLAS_LIB=...
export CPU_LAPACK_LIB=...
```

所有实际计算和编译都必须通过Slurm提交。不同集群的分区和GPU资源名可能不同，使用前应检查脚本中的`--partition`、`--nodelist`、`--gres`与`--qos`设置。

如果LinglongAMG头文件位于`${LINGLONGAMG_ROOT}/lib`，CMake会自动搜索该目录。

r5含953,048,225个非零元。若HYPRE版本的ParCSR转置仍使用32位局部计数，需要在外部库侧修复该整数溢出；本包只负责链接修复后的HYPRE。

## 构建

提交三个构建作业：

```bash
bash scripts/submit-builds.sh
```

也可以分别提交：

```bash
sbatch scripts/build-hypre-cpu.sbatch
sbatch scripts/build-hypre-gpu.sbatch
sbatch scripts/build-linglongamg.sbatch
```

构建结果分别位于：

```text
build/hypre-cpu/main/testOpenCAEPoro
build/hypre-gpu/main/testOpenCAEPoro
build/linglongamg/main/testOpenCAEPoro
```

## 测试

构建完成后提交三组r1-r5测试：

```bash
bash scripts/test-all.sh
```

测试脚本使用`normal`QoS和`00:20:00`上限。HYPRECPU的r1-r4使用`cpu9462`的64个物理核；r5单独申请116G内存并使用8个OpenMP线程，以避开外部HYPRE高线程临时内存膨胀。HYPREGPU和LinglongAMG各使用一块A100。

如需调整r5线程数：

```bash
OCP_CPU_R5_THREADS=8 bash scripts/test-all.sh
```

测试完成后汇总：

```bash
bash scripts/summarize-results.sh | tee results/summary.csv
```

主要计时字段：

- `wall_s`：Slurm任务墙钟时间。
- `common_total_seconds`：从程序启动到基础一阶压力矩阵完成的公共时间。
- `matrix_generation_seconds`：在OpenCAEPoroX内生成目标r矩阵的时间。
- `solver_preprocess_seconds`：求解器CSR整理、数据传输和插入矩阵的时间。
- `setup_seconds`：AMG层次结构建立时间。
- `solve_seconds`：GMRES求解时间。

## 矩阵参数

核心环境变量：

```text
OCP_SPE10_REFINE_FACTOR=1..5
OCP_SINGLE_SOLVE_BENCHMARK=1
OCP_LINEAR_SOLVER=hypre
OCP_LINEAR_SOLVER=linglongamg
```

矩阵由OpenCAEPoroX在内存中直接生成，不读取外部矩阵文件。

本实现生成的矩阵规模为：

| 算例 | n | nnz |
| --- | ---: | ---: |
| r1 | 1,094,421 | 7,478,141 |
| r2 | 8,755,368 | 60,556,352 |
| r3 | 29,549,367 | 205,200,315 |
| r4 | 70,042,944 | 487,375,712 |
| r5 | 136,802,625 | 953,048,225 |

## R5说明

当前原生R5时间步比例尚未由软件直接给出。R5沿用R4的固定时间步比例扩展，测试脚本会显式设置：

```text
OCP_SPE10_R5_USE_R4_TIME_RATIO=1
```

日志中的`native_r5_time_step=0`表示R5不是由原生R5时间步比例生成。

## HYPRE参数

默认配置：

```text
OCP_HYPRE_GMRES_RESTART=10
OCP_HYPRE_INTERP_MAX_ELEMENTS=4
OCP_HYPRE_STRONG_THRESHOLD=0.25
OCP_HYPRE_AGGRESSIVE_LEVELS=0
OCP_HYPRE_RAP2=1
```

HYPRECPU/GPU版本由链接的HYPRE库决定，不通过OpenCAEPoroX运行参数切换。`r4`和`r5`要求HYPRE库能够处理对应的两亿级和一亿级矩阵。本包不修改也不附带HYPRE源码。

## LinglongAMG参数

默认读取：

```text
examples/linglongamg.json
```

可以通过以下环境变量替换：

```bash
export OCP_LINGLONGAMG_CONFIG=linglongamg.json
```

LinglongAMG通过外部安装包中的`linglongamg.h`、`liblinglongamg.so`和CUDA runtime链接，这些文件均不随包分发。
