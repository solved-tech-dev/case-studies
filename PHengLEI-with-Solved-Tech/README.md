# 风雷 PHengLEI 接入 LinglongAMG

本地发布候选包。包含新版适配层、三例单区计算输入及构建、运行和验证脚本。拟公开文件及再分发核对项见[发布清单](docs/发布清单.md)和[许可与来源](LICENSE-NOTICE.md)。

风雷负责网格、有限体积离散、组装方程及更新流场。不可压求解器的压力修正方程可切换到 LinglongAMG，在一张 GPU 上使用 GMRES（广义最小残量法）和 AMG（代数多重网格）预条件。动量、湍流、能量和组分方程按算例设置调用 CPU HYPRE。

```text
overlay/             新增适配器、CSR 组装、计时和会话测试源码
patches/             风雷求解入口、CMake、压力参考及 HYPRE 兼容补丁
examples/            H11、H12、H13 运行参数、压缩网格及校验值
scripts/             准备源码、构建、运行、会话测试和 CSV 汇总
upstream.lock.json   基线提交及补丁涉及文件的 SHA-256
```

| 算例 | 工况 | 三维单元数 |
| --- | --- | ---: |
| H11 | 三维不可压气体泄漏与混合 | 640,178 |
| H12 | 三维不可压六组分输运 | 225,610 |
| H13 | SUBOFF 标模三维不可压绕流 | 695,490 |

三个算例都已转换成单网格区。各例 `case.json` 记录原始教程名、网格及壁面距离文件的校验值。脚本在新的运行目录解压文件并核对哈希。风雷以 `./grid/mesh.fts` 为入口名，单区读入实际文件 `mesh_0.fts`。

## 准备环境

- 在 Linux x86-64 上构建和运行。准备 Bash、Git、Python 3.9+、CMake 3.17+（建议 3.25+）、GCC/G++ 与 MPI C++ 编译器。风雷上游所需的 HDF5、网格工具及其他依赖按其构建说明准备。
- 使用有权取得的风雷源码。基线为 NAS03 上的风雷 2606 派生提交 `9da7a4e5c48e05daf26b332d298c99004291be2e`。准备脚本优先从 Git 提取该提交，也接受补丁所需文件与锁定哈希相符的源码快照。[风雷官方文档及获取说明](https://docs.osredm.com/phenglei/)
- 准备 Eigen 3.2.5。若授权源码工作树已包含 `PHengLEI/3rdparty/eigen/Eigen/Eigen`，准备脚本会沿用其中的头文件。否则把 `EIGEN_ROOT` 指向 Eigen 解包根目录。[Eigen 3.2.5](https://gitlab.com/libeigen/eigen/-/tags/3.2.5)
- 安装 CPU 版 HYPRE 3.2.1，采用双精度及默认整数模式。HYPRE 与风雷使用同一 MPI 工具链，启用 OpenMP 时使用相容的 GNU 运行时。[HYPRE 上游仓库](https://github.com/hypre-space/hypre)
- 玲珑构建额外依赖 CUDA 13.1 或更高的 13.x 工具包、匹配的 NVIDIA 驱动，以及 LinglongAMG 26.44 / CUDA 13.1 持久会话版 SDK。SDK 与对应设备的授权按供应方流程单独获取。SDK 根目录下需有 `lib/linglongamg.h` 与 `lib/liblinglongamg.so`。[解得科技商城](https://store.solved-tech.com/)

```bash
export PHENGLEI_UPSTREAM_DIR=/path/to/authorized/PHengLEI
export HYPRE_ROOT=/path/to/hypre-cpu-install
export MPICXX=/path/to/mpicxx
# 源码树缺少 Eigen 时设置下一行。
# export EIGEN_ROOT=/path/to/eigen-3.2.5

bash scripts/prepare-source.sh
bash scripts/build.sh hypre

export LINGLONGAMG_ROOT=/path/to/LinglongAMG
export CUDA_ROOT=/path/to/cuda-13.x
export PHENGLEI_LINGLONG_TESTS=ON
bash scripts/build.sh linglong
```

`PHENGLEI_SOURCE_DIR` 可指定一个尚不存在的工作副本目录。构建默认写入 `build/hypre/` 和 `build/linglong/`，`BUILD_JOBS` 控制编译并行度。运行前按工具链和 SDK 说明加载动态库搜索路径。如果使用集群，在已分配的计算节点执行下列运行命令。

## 运行算例

```bash
export OMP_NUM_THREADS=1

for case in H11 H12 H13; do
  python3 scripts/run-case.py "$case" hypre --steps 20 --run-dir "runs/$case-hypre"
  python3 scripts/run-case.py "$case" linglong --steps 20 --run-dir "runs/$case-linglong"
done

python3 scripts/summarize.py runs/H{11,12,13}-{hypre,linglong} > summary.csv
```

每次启动一个 MPI 进程。玲珑路线使用单区、单 GPU、FP64 及 32 位 CSR 索引，`--device` 取作业内可见的 CUDA 编号，默认 0。HYPRE 路线沿用算例的原生压力求解参数。可通过 `PHENGLEI_EXE` 指定已构建的二进制，通过 `--prepare-only` 仅解压并生成运行目录。运行脚本拒绝覆盖已有算例目录。

## 会话和验证

压力适配器初始化时创建 SDK 求解句柄。每次压力方程装配后，转换成 CSR（压缩稀疏行格式），上传矩阵和向量，再执行 `setup → solve → get_solve_info`。更换线性系统时对原句柄执行 `reset`，后续仍重新构建 AMG 层级。每次结果查询由 SDK 按其契约执行授权检查。退出时释放句柄。授权或运行错误会终止求解，失效会话需要重新启动进程。

```bash
bash scripts/test-session.sh "$PWD/runs/session-test"
```

先用 `PHENGLEI_LINGLONG_TESTS=ON` 构建，在有有效授权的 GPU 作业中运行会话测试。测试转发层调用真实 SDK，记录调用计数与句柄，并验证 SDK 查询被拒时不会把结果写回风雷。

`run.json` 给出进程退出状态与墙钟。玲珑的 `results/linglong_pressure.csv` 逐次记录独立重算的残差、收敛标记、上传、setup、solve、结果查询、下载和校验时间。主机侧 setup/solve 含设备同步，`sdk_setup_ms` 和 `sdk_solve_ms` 记录 SDK 自报值。压力收敛标记和外层 CFD 收敛分开核对。历史 CUDA 12.4 子进程实验的耗时不用于本直接适配版本。详见[验证记录](docs/验证记录.md)。
