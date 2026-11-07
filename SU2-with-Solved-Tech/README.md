# SU2 五算例与 XianzhiKrylov 接入包

本目录提供 SU2 8.5.0 的 XianzhiKrylov 接入覆盖层、五个算例的配置与网格，以及构建、测试和日志汇总脚本。

五个算例为 Cavity、Kenics、IncTurbulentFlatPlate、IncNACA 和 Urban。SU2 负责网格、方程组装和流场更新，适配器在隐式线性求解入口调用 XianzhiKrylov。

## 代码结构

```text
overlay/                        覆盖到 SU2 8.5.0 源码树的接入文件
overlay/Common/                 预条件等级配置项与 NVCC/OpenMP 兼容修改
overlay/adapter/                XianzhiKrylov 适配器及链接包装器
examples/<case>/run_100.cfg      五例原配置
examples/<case>/*.su2.gz         五例网格，无损 gzip 压缩
scripts/prepare-source.sh       从外部源码树生成工作副本并应用覆盖层
scripts/build-xianzhikrylov.sh  编译 SU2 与适配器
scripts/build-xianzhikrylov.sbatch  Slurm 构建作业
scripts/run-case.sh             单算例运行入口
scripts/test-cpu.sbatch         9462 CPU 的五例测试
scripts/test-xianzhikrylov.sbatch  A100 的五例测试
scripts/test-all.sh             提交两组五例测试
scripts/summarize-results.sh    汇总运行记录为 CSV
```

## 准备源码

使用外部的 [SU2 v8.5.0 源码](https://github.com/su2code/SU2/tree/v8.5.0)，先按上游脚本初始化构建依赖。下面命令在外部源码目录执行，依赖下载需要网络连接。

```bash
git clone --depth 1 --branch v8.5.0 \
  https://github.com/su2code/SU2.git /path/to/SU2-8.5.0
cd /path/to/SU2-8.5.0
python3 preconfigure.py --with-own-meson \
  --no-codi --no-medi --no-opdi --no-mpp --no-coolprop --no-fado --no-mlpcpp
```

回到 `SU2-with-Solved-Tech/` 目录，用 `SU2_UPSTREAM_DIR` 指定已准备的源码树。

```bash
export SU2_UPSTREAM_DIR=/path/to/SU2-8.5.0
bash scripts/prepare-source.sh
```

默认生成 `work/SU2/`。脚本复制外部源码和依赖，再覆盖 `overlay/` 中的文件。目标目录已存在时会停止。可用 `SU2_SOURCE_DIR` 指定其他工作副本位置。

## 配置第三方库

XianzhiKrylov SDK、头文件与授权由用户单独安装。获取匹配的 26.45 版本和许可证，请访问[解得科技商城](https://store.solved-tech.com/)。SDK 根目录应包含

```text
include/XianzhiKrylovAPI.h
include/XianzhiKrylovData.h
lib/libXianzhiKrylov.so.2
license.key
```

保留 SDK 自带的库文件和 soname 链接。适配器除 ABI 3 公开接口外，还使用匹配库中的 `xk_internal_su2_htd_*` 和 `xk_internal_get_last_stats`，构建工具会检查所需的八个符号。

设置本机依赖位置。

```bash
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov
export CUDA_ROOT=/path/to/cuda
export CC=/path/to/gcc
export CXX=/path/to/g++
export NINJA_BIN=/path/to/SU2-8.5.0/ninja
export PYTHON_BIN=python3
```

环境要求为 Linux x86_64、Python 3.12 或更新版本、支持 C++17/OpenMP 的 GCC/G++、GNU 链接器及 `nm`、Ninja 和 CUDA 开发工具。交接包的构建环境为 GCC 12.4、CUDA 13.2.1，CUDA 编译目标为 A100 SM80。Meson、Eigen、MEL 由外部 SU2 源码的依赖初始化提供。

所有实际编译和计算通过 Slurm 提交。使用前按集群检查脚本中的分区、节点、GPU 资源名、QoS、内存和时限。

## 构建

```bash
bash scripts/submit-builds.sh
```

也可直接提交构建作业。

```bash
sbatch scripts/build-xianzhikrylov.sbatch
```

默认输出为

```text
build/xianzhikrylov/install/bin/SU2_CFD
build/xianzhikrylov/lib/libSU2Xianzhi45.so
build/xianzhikrylov/sdk -> 外部 SDK 根目录
build/xianzhikrylov/logs/
build/xianzhikrylov/build-report.json
```

构建先生成 SU2 CUDA 版本，再通过 GNU 链接器 `--wrap` 接入 `CSysSolve<double>::Solve`。同一二进制通过 `SU2_XIANZHI_BACKEND=cpu|gpu` 选择后端，算例中保留 `LINEAR_SOLVER=FGMRES`。

SDK 加载按二进制所在包解析 `sdk/lib/libXianzhiKrylov.so.2`。需要整体迁移运行包时，保留 `install/`、`lib/`、`sdk/` 的相对布局，并将 SDK 链接替换为实际 SDK 文件。目标机器还需提供兼容的 CUDA、驱动及编译器运行库。

`SU2_BUILD_ROOT` 可改变构建根目录。相同设置可在原工作目录继续构建，改变源文件、编译器或依赖设置时使用新的构建目录。

## 测试

构建完成后提交 CPU、GPU 两组五例作业。

```bash
bash scripts/test-all.sh
```

也可分别提交。

```bash
sbatch scripts/test-cpu.sbatch
sbatch scripts/test-xianzhikrylov.sbatch
```

两组均为单进程，五例按作业数组串行运行。CPU 默认使用 `cpu9462` 的 64 个物理核，GPU 默认使用一块 A100 和 8 个主机线程。运行目录保存实际配置、展开网格、应用日志和逐次线性调用记录。

默认开启独立残差复核。计时运行可显式关闭额外复核，并单独保存新一轮结果。

```bash
SU2_TIMING=1 bash scripts/test-all.sh
```

输出分别写入 `logs/` 和 `results/`。汇总命令为

```bash
bash scripts/summarize-results.sh > results/summary.csv
```

主要计时字段

- `application_wall_seconds` 是应用从启动到退出的墙钟时间。
- `whole_linear_seconds` 是完整线性求解调用的累计墙钟时间，CPU、GPU 均记录。
- `gpu_execute_seconds` 是 GPU 适配器内部分析、数值准备与求解区间的累计墙钟时间。
- `gpu_prepare_seconds` 是 GPU 结构分析与数值准备的累计时间。
- `gpu_solve_seconds` 是 GPU 迭代求解区间的累计时间。

GPU 内部时间属于完整调用的一部分，两者不能相加。比较完整线性调用时使用双方的 `whole_linear_seconds`。这些字段由主机墙钟测量，解释时应保留其计时边界。

## 算例与迭代参数

| 算例 | 物理场景 | 维度 | 单元 | 节点 | 复测内迭代上限 |
|---|---|---:|---:|---:|---:|
| Cavity | 方腔自然对流 | 2 | 65,536 | 66,049 | 30 |
| Kenics | 静态混合器、组分输运 | 3 | 128,790 | 138,324 | 20 |
| IncTurbulentFlatPlate | 湍流平板边界层 | 2 | 208,896 | 209,825 | 20 |
| IncNACA | NACA0012 湍流绕流 | 2 | 229,376 | 230,336 | 10 |
| Urban | 城市街区风场 | 2 | 786,472 | 488,212 | 10 |

运行器默认设置 `ITER=20`、`CONV_STARTITER=21`、`RESTART_SOL=NO`。内迭代上限采用上表，CPU 与 GPU 一致。Kenics 原输入写 5，最新复测脚本实际使用 20，本包运行器保留该覆盖。

`XIANZHI_PRECONDITIONER_LEVEL` 直接在各例配置中设置，允许 0–3，缺省 2，五例当前均为 1。适配器默认每次线性调用刷新数值预条件器。

## 当前复测记录

以下数值来自维护者提供的更新交接包，SHA256 为 `f6e56c9ca2a78769c2916902cc883fb0a0347c501f2fe6f25bcbe041d3533b30`。平台为 A100-PCIE-40GB、单进程、8 个主机线程、20 步外迭代，预条件等级均为 1。

| 算例 | GPU 内部区间累计（s） | 完整调用累计（s） | 内迭代累计 | 平均残差 |
|---|---:|---:|---:|---:|
| Cavity | 0.642 | 1.252 | 600 | 8.403e-4 |
| Kenics | 1.218 | 3.885 | 601 | 1.144e-4 |
| IncTurbulentFlatPlate | 0.946 | 2.950 | 800 | 1.748e-3 |
| IncNACA | 0.738 | 2.934 | 400 | 4.554e-2 |
| Urban | 1.687 | 7.828 | 329 | 3.853e-3 |

这些是交接包已有的复测记录。源码整理后尚未完成独立的完整 CUDA/SDK 重编译。当前源树提供 CPU 完整调用计时，已有成品中的 CPU 分阶段插桩没有包含在交接源码中。

20 步用于固定性能工作负载。检查结果时结合退出码、内迭代上限、线性残差及物理场残差判断，达到内迭代上限也可能伴随应用退出码 0。修改预条件等级后需重新检查收敛表现。

## 许可

SU2 覆盖文件和适配器保留 LGPL-2.1-or-later 声明，许可全文见 [overlay/COPYING](overlay/COPYING)。本目录新增工具和文档采用相同许可。算例配置保留原作者和上游说明。SDK 及其授权按供应方条款单独安装使用。
