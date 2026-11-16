# OpenFOAM v2412 接入 LinglongAMG 和 XianzhiKrylov

本包使用 `occDrivAerStaticMesh` 的 65M 单元网格，在同一份初始流场上分别运行两套配置。LinglongAMG 求解压力方程 `p`。XianzhiKrylov 求解 `Ux`、`Uy`、`Uz`、`k` 和 `omega`，压力方程仍由 OpenFOAM 的 GAMG 求解。

| 目录 | 内容 |
| --- | --- |
| `overlay/` | 两套 OpenFOAM 适配层、调用入口及 `wmake` 配置 |
| `examples/` | OCC 初始和边界条件、离散参数、分区设置及两套 SDK 的 JSON 配置 |
| `scripts/` | 准备源码、编译、初始化算例和运行作业的脚本 |

## 环境

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

## 下载网格并准备算例

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

## 分别运行

在单节点、20 个 CPU 核和一块 GPU 的计算资源上，从同一入口提交两个作业。将 `GPU_PARTITION` 设为集群的 GPU 分区。最后一个参数为 SIMPLE 外迭代步数，默认 5。

```bash
export GPU_PARTITION=a100
sbatch --partition="$GPU_PARTITION" --gres=gpu:1 --exclusive --mem=0 --time=00:30:00 \
  scripts/run-case.sbatch linglongamg work/occ-prepared work/occ-linglongamg 5
sbatch --partition="$GPU_PARTITION" --gres=gpu:1 --exclusive --mem=0 --time=00:30:00 \
  scripts/run-case.sbatch xianzhikrylov work/occ-prepared work/occ-xianzhikrylov 5
```

已有作业分配时，也可以在作业内用 `bash scripts/run-case.sh` 传入同样的参数。若集群需要专用 MPI 启动器，可用 `MPI_RUN` 指向兼容 `mpirun -np` 的启动脚本。脚本会创建各自的运行目录和求解器字典。运行完成时脚本检查实际步数，运行日志位于各自目录下的 `log.openfoam`。

## 来源与许可

`overlay/` 中基于 OpenFOAM 的代码保留了 GPL-3.0-or-later 标记及原有许可声明。见 [GPL 许可正文](overlay/COPYING)。

OCC 算例原始设置由 Charles Mockett、Hendrik Hetmann 和 Felix Kramer（Upstream CFD GmbH）编写，由 Mark Wasserman（Huawei）和 Sergey Lesnik（Wikki GmbH）为 OpenFOAM HPC Challenge 修改。版权归 Upstream CFD GmbH（2022–2023）、Huawei Technologies Co., Ltd.（2024–2025）及 Wikki GmbH（2024–2025）。本包将算例改为 20 个分区，并提供单独的求解器设置和工作流。算例配置沿用上游的 CC BY-SA 4.0 许可，网格数据使用 Zenodo 记录标注的 CC BY 4.0 许可。 [上游算例](https://develop.openfoam.com/committees/hpc/-/tree/develop/incompressible/simpleFoam/occDrivAerStaticMesh) · [算例配置许可](https://creativecommons.org/licenses/by-sa/4.0/) · [网格作者及许可](https://zenodo.org/records/15012221)
