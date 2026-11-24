# OpenFOAM v2412 接入 LinglongAMG 和 XianzhiKrylov

本目录保留原有的两条 26.46 独立求解路线，并追加 26.47-resident-r1 GPU 常驻方案。三条路线的入口和依赖分别如下。

| 路线 | 求解范围 | 运行入口 |
| --- | --- | --- |
| 26.46 LinglongAMG 独立路线 | LinglongAMG 求解压力方程 `p`，其余方程沿用 OpenFOAM | `scripts/run-case.sh linglongamg ...` |
| 26.46 XianzhiKrylov 独立路线 | XianzhiKrylov 求解 `Ux`、`Uy`、`Uz`、`k`、`omega`，压力仍由 OpenFOAM GAMG 求解 | `scripts/run-case.sh xianzhikrylov ...` |
| 26.47 GPU 常驻路线 | 准备库负责 GPU 几何、装配和 SIMPLE/SST 更新；XianzhiKrylov 求解线性系统，LinglongAMG 为压力提供 AMG 预条件 | `scripts/run-resident.sh ...` |

原有适配源码、OCC 配置和运行脚本保持原样。新增源码位于 `overlay/resident/application/`，新增脚本和配置通过 `resident` 命名与旧入口区分。

| 目录 | 内容 |
| --- | --- |
| `overlay/linglongamg/`、`overlay/xianzhikrylov/` | 原有 26.46 适配层、调用入口及 `wmake` 配置 |
| `overlay/resident/application/` | 新增 26.47 原始输入转换、常驻上下文调用和 `wmake` 配置 |
| `examples/` | 原有 OCC 配置及 SDK 配置，另有新增的 `resident-*.json` 模板 |
| `scripts/` | 原有脚本，以及新增的常驻方案编译、准备和运行入口 |
| `resident-dependencies.json` | 常驻方案所需的配套 SDK 版本和身份摘要 |

## 原有 26.46 独立路线

### 环境与编译

- Linux、OpenFOAM v2412 的双精度 32 位索引构建（`DPInt32`）、`wmake`、Open MPI 和 Python 3。编译适配层时使用与 OpenFOAM 一致的编译器和 MPI。
- CUDA 工具包、支持所用 SDK 的 NVIDIA 驱动和一块 GPU。两条路线各自用 20 个 MPI 进程，均由第 0 个进程调用 GPU。为 65M 网格准备足够的主机内存、显存和磁盘空间。
- 通过 SDK 交付渠道获取匹配的安装包及当前 GPU 的授权。LinglongAMG 目录下需要 `lib/liblinglongamg.so` 和 `linglongamg.h`。XianzhiKrylov 需要 26.46 头文件和匹配的动态库，动态库须导出 `xk_bench_set_tolerances` 与 `xk_internal_get_last_stats`，构建脚本会检查这两个符号。

在包根目录执行以下命令。将路径替换为实际安装位置。

```bash
source /path/to/OpenFOAM-v2412/etc/bashrc
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_SDK=/path/to/LinglongAMG
export XIANZHIKRYLOV_SDK=/path/to/XianzhiKrylov

bash scripts/prepare-source.sh
bash scripts/build.sh linglongamg
bash scripts/build.sh xianzhikrylov
```

默认构建目录为 `work/source/`。当 CUDA 或 SDK 布局不同时，可设置 `CUDA_INCLUDE_DIR`、`CUDA_LIB_DIR`、`LINGLONGAMG_INCLUDE_DIR` 或 `XIANZHIKRYLOV_LIB_DIR`。如果授权文件在 SDK 根目录之外，用 `LINGLONGAMG_LICENSE` 或 `XIANZHIKRYLOV_LICENSE` 指向相应文件。

### 下载网格并准备算例

从 OpenFOAM HPC Challenge 在 Zenodo 发布的网格中下载 `polyMesh_65M.tar.gz`。下载包约 3 GiB。以下命令使用 GNU/Linux 工具校验并解包。 [Zenodo 网格记录](https://zenodo.org/records/15012221)

```bash
mkdir -p mesh-download/unpacked
curl --fail --location --continue-at - \
  'https://zenodo.org/records/15012221/files/polyMesh_65M.tar.gz?download=1' \
  --output mesh-download/polyMesh_65M.tar.gz
printf '%s\n' '50393a15d84f53255fda85ac702390e9  mesh-download/polyMesh_65M.tar.gz' \
  | md5sum --check
tar -xzf mesh-download/polyMesh_65M.tar.gz -C mesh-download/unpacked
find mesh-download/unpacked -type f \( -name points -o -name points.gz \)
```

将下面的 `/path/to/polyMesh` 替换为上一步找到的 `points` 或 `points.gz` 所在目录。在分配了 20 个 CPU 核的计算节点上准备一次算例。

```bash
bash scripts/prepare-case.sh /path/to/polyMesh work/occ-prepared
```

脚本按 `5 × 2 × 2` 分区并重排网格，用 `potentialFoam` 初始化流场，再用 `applyBoundaryLayer` 初始化边界层。两次运行都会读取这个目录下的同一份网格和初场。

如果已有按相同网格和分区准备好的输入，可以向算例提供方取得 `constant/` 及 `processor0/` 至 `processor19/` 下的 `0/`、`constant/polyMesh/`，并直接指定该输入目录。可选的 LinglongAMG 行排列文件 `processor0.bin` 至 `processor19.bin` 也由算例提供方另行提供。已有对应文件时可设置 `LINGLONGAMG_ORDERING_DIR`，文件必须与所用的 20 个分区一一匹配。

### 分别运行

在单节点、20 个 CPU 核和一块 GPU 的计算资源上，从同一入口提交两个作业。将 `GPU_PARTITION` 设为集群的 GPU 分区。最后一个参数为 SIMPLE 外迭代步数，默认 5。

```bash
export GPU_PARTITION=a100
sbatch --partition="$GPU_PARTITION" --gres=gpu:1 --exclusive --mem=0 --time=00:30:00 \
  scripts/run-case.sbatch linglongamg work/occ-prepared work/occ-linglongamg 5
sbatch --partition="$GPU_PARTITION" --gres=gpu:1 --exclusive --mem=0 --time=00:30:00 \
  scripts/run-case.sbatch xianzhikrylov work/occ-prepared work/occ-xianzhikrylov 5
```

已有作业分配时，也可以在作业内用 `bash scripts/run-case.sh` 传入同样的参数。若集群需要专用 MPI 启动器，可用 `MPI_RUN` 指向兼容 `mpirun -np` 的启动脚本。脚本会创建各自的运行目录和求解器字典。运行完成时脚本检查实际步数，运行日志位于各自目录下的 `log.openfoam`。

## 新增 26.47 GPU 常驻路线

常驻入口读取原始网格、初场和字典，将它们转换为准备库的数据结构，构造一次 `ResidentContext`，再调用 `advance()` 推进 SIMPLE/SST 外迭代。准备库、求解核心和压力预条件器配合运行，使用**一个进程、一块 GPU**；CPU 仍负责输入解析、控制和资源管理。

### 依赖与编译

需要 Linux、Python 3.9 或更新版本、OpenFOAM v2412 `DPInt32`，以及与 OpenFOAM 一致的编译器和 MPI 环境。本适配对应 CUDA 13.1、A100/SM80 的 26.47-resident-r1 交付物。

| 变量 | 配套目录必须包含 |
| --- | --- |
| `PREPARATION_SDK` | `include/XianzhiPreparation2647.hpp`、`lib/libLinearSystemPreparation.so` |
| `XIANZHIKRYLOV_SDK` | `include/XianzhiKrylovAPI.h`、匹配的头文件和 `lib/libXianzhiKrylov.so` |
| `LINGLONGAMG_SDK` | `lib/liblinglongamg.so` |
| `RESIDENT_CUDA` | CUDA 13.1 的 `include/` 和 `lib64/` |

通过 SDK 交付渠道取得匹配的准备库、求解核心、LinglongAMG 及当前设备的有效授权。准备库使用配套的实验性 C++ 接口，头文件、动态库和 C++ 工具链必须匹配；原有 26.46 SDK 不能直接用于此入口。SDK 源码、动态库和授权文件由交付渠道另行提供。

以下命令均从 `OpenFOAM-with-Solved-Tech/` 目录执行：

```bash
source /path/to/OpenFOAM-v2412/etc/bashrc
export PREPARATION_SDK=/path/to/XianzhiKrylov-26.47-resident-r1/preparation
export XIANZHIKRYLOV_SDK=/path/to/XianzhiKrylov-26.47-resident-r1/core
export LINGLONGAMG_SDK=/path/to/LinglongAMG
export RESIDENT_CUDA=/path/to/cuda-13.1

bash scripts/build-resident.sh
```

默认生成 `work/resident/build/bin/rawResident`，可通过 `RESIDENT_BUILD` 指定其他构建目录。常驻适配入口不依赖交付示例中的 CUPTI 追踪组件。

### 输入要求与准备

输入必须是已经初始化、单元编号与场顺序一致的**串行 OpenFOAM 算例**，包含以下内容：

```text
prepared-serial-case/
├── 0/                          # U、p、k、omega、nut
├── constant/
│   ├── polyMesh/               # points、faces、owner、neighbour、boundary
│   ├── transportProperties
│   └── turbulenceProperties
└── system/
    ├── fvSchemes
    ├── fvSolution
    └── ...                    # 字典引用的 include 文件等
```

OCC 原始网格的获取方式见前文。常驻路线还需要配套的完整串行初场、边界和离散设置，可向算例提供方取得 resident-r1 所用输入及对应的 `gs-config.json`。

若输入来自原有的 20 个分区，必须先重建完整网格与初场，并核对重建后的单元编号、场顺序和 GS 分区偏移。GS 配置中的 `parameters.partition_offsets` 应保留原分区在重建单元序中的连续区间，末项应等于总单元数。不能直接将原有的 `processor*` 输入传给常驻入口。

`examples/resident-gs-one-sweep-provided.json` 是单分区模板；需要保留原分区扫描顺序时，使用对应算例的 GS 配置。下面的脚本复制控制字典，将原始网格和初场作为只读引用保留，因此源目录需要一直存在。

```bash
python3 scripts/prepare-resident-case.py \
  /path/to/prepared-serial-case work/occ-resident \
  --gs-config /path/to/case-specific-gs-config.json --steps 5

bash scripts/run-resident.sh \
  work/occ-resident work/occ-resident/gs-config.json work/resident-input-check \
  --check-case
```

`--steps` 设置固定外迭代次数。`--check-case` 会解析网格、初场、边界和离散设置，完成数据转换，并在提供分区偏移时核对总单元数。该检查不创建 GPU 上下文，设备运行与数值验证仍需单独完成。

### 运行与输出

在分配了匹配 GPU 的计算节点上运行：

```bash
export RESIDENT_CPU_THREADS=20
bash scripts/run-resident.sh \
  work/occ-resident work/occ-resident/gs-config.json work/resident-output
```

Slurm 示例；按集群实际情况设置 GPU 分区、内存和时限：

```bash
sbatch --partition=a100 --gres=gpu:1 --mem=100G --time=00:20:00 \
  scripts/run-resident.sbatch \
  work/occ-resident work/occ-resident/gs-config.json work/resident-output
```

输入与输出目录应分开，输出目录必须是新目录。运行器填写本机依赖路径，将三份求解配置封存到只读内存后交给应用；运行日志、实际配置和依赖身份记录在输出目录。脚本检查返回值和外迭代完成标记。

需要导出诊断场时添加 `--dump-fields`，原始 FP64 数组写到 `fields/initial/` 和 `fields/stepN/`，数组顺序由配套 SDK 定义。这些 `.f64` 文件需要配套网格与数组解释才能可视化，此入口不写出 OpenFOAM 原生时间目录。

### 支持范围

当前适配限于固定拓扑、不可压 Newtonian 流体、稳态 `kOmegaSST`、consistent SIMPLE、零次非正交附加修正及源码登记的离散格式和边界组合。速度使用 15 扫 Gauss–Seidel，`k` 与 `omega` 使用 8 扫；模板 `max_iterations=1` 为每扫入口，完整扫描次数由准备库控制。

压力使用 XianzhiKrylov BiCGStab 与 LinglongAMG 预条件器。压力模板保留配套交付示例的 `1e-12` 相对容差和 20000 次上限，这属于所选算例配置。当前按固定外迭代步数运行，尚未接入 SIMPLE `residualControl` 收敛终止，完成标记只表示请求的步骤执行完成。动态网格、耦合 patch、未登记的模型、边界或 `fvOptions` 需要另行适配。

## 来源与许可

`overlay/` 中基于 OpenFOAM 的代码保留了 GPL-3.0-or-later 标记及原有许可声明。见 [GPL 许可正文](overlay/COPYING)。

新增常驻适配来自 26.47-resident-r1 的输入读取与生命周期调用示例。准备库对应源码标识 `7c9a616b7b3800765a48ace0171de19aabc99cdc`，求解核心对应 `b39ea72a3c84bc94b93ad3f3d06023be65119750`，配套交付物身份见 [resident-dependencies.json](resident-dependencies.json)。SDK 配套 C++ 接口不作为稳定公共 ABI 承诺，外部 SDK 按其交付许可使用。

OCC 算例原始设置由 Charles Mockett、Hendrik Hetmann 和 Felix Kramer（Upstream CFD GmbH）编写，由 Mark Wasserman（Huawei）和 Sergey Lesnik（Wikki GmbH）为 OpenFOAM HPC Challenge 修改。版权归 Upstream CFD GmbH（2022–2023）、Huawei Technologies Co., Ltd.（2024–2025）及 Wikki GmbH（2024–2025）。本包将算例改为 20 个分区，并提供单独的求解器设置和工作流。算例配置沿用上游的 CC BY-SA 4.0 许可，网格数据使用 Zenodo 记录标注的 CC BY 4.0 许可。 [上游算例](https://develop.openfoam.com/committees/hpc/-/tree/develop/incompressible/simpleFoam/occDrivAerStaticMesh) · [算例配置许可](https://creativecommons.org/licenses/by-sa/4.0/) · [网格作者及许可](https://zenodo.org/records/15012221)
