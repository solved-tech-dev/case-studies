# OPM Flow 2026.04 接入 LinglongAMG 与 XianzhiKrylov

本目录提供基于 **OPM Flow 2026.04**（`release/2026.04/final`）的公开应用集成层（Overlay + 最小统一补丁），允许用户通过统一入口按需选择以下**三条互不耦合的构建与求解路线**：

```text
                        OPM Flow 2026.04 (Black-Oil BSR3 / CPRW)
                                          │
          ┌───────────────────────────────┼───────────────────────────────┐
          ▼                               ▼                               ▼
   route=linglongamg              route=xianzhikrylov                 route=cpr
   (LinglongAMG 独立路线)          (XianzhiKrylov 独立路线)          (双后端 CPR 组合路线)
          │                               │                               │
   CPU 外层 BiCGStab               GPU 全耦合 BSR3 BiCGStab        GPU 外层 BiCGStab (XK)
   + CPU paroverilu0               + GPU BSR3 块 ILU(0) (XK)       + True-IMPES 压力 CSR
   + CPRW True-IMPES               + GPU 独立标准井 Schur 算子     + GPU GMRES(10) (XK)
   + GPU GMRES(10)/AMG (LL)        (无 CPR / 无 AMG)               + 单次 AMG V-cycle (LL)
   [不链接 XianzhiKrylov]          [不链接 LinglongAMG]            + GPU 细层块 ILU(0) (XK)
                                                                   + GPU 独立标准井 Schur 算子
```

> [!IMPORTANT]
> 三条路线对应不同的数学算法路径、预条件器边界和对比基线，**严禁**将 `linglongamg` 的压力 API 耗时、`xianzhikrylov` 的全耦合 ILU(0) 耗时与 `cpr` 的双后端分阶段耗时放进同一组直接加速比中对比。

---

## 1. 三条路线与线性系统分工

| 路线 (`--route`) | 外层迭代与矩阵 | 细层预条件器 | 压力子系统 / AMG | 井贡献处理 | 编译/链接依赖 |
| --- | --- | --- | --- | --- | --- |
| `linglongamg` | 原生 OPM CPU BiCGStab（全系统，默认 `tol=0.005, maxiter=20`） | 原生 OPM CPU `paroverilu0` (`relaxation=1.0`) | CPRW True-IMPES 压力 CSR；GPU 原生 GMRES(10) + AMG V-cycle (`tol=1e-6, maxiter=100`) 由 `liblinglongamg.so` 承担 | 原生 OPM `WellModelOperator` | 仅 CUDA + LinglongAMG（**不链接** XianzhiKrylov，不需要其授权） |
| `xianzhikrylov` | XianzhiKrylov GPU 全耦合 `BSR3` BiCGStab (`tol=0.005, maxiter=200`) | XianzhiKrylov GPU 原序 `BSR3` 块 ILU(0)（仅储层矩阵） | 不启用 CPR / 不调用 AMG | GPU 储层 `BSR3` SpMV + 修复版 `WellContributionsCuda` 独立标准井 Schur 补 | 仅 CUDA + XianzhiKrylov（**不链接** LinglongAMG，不需要其授权） |
| `cpr` | XianzhiKrylov GPU 外层 BiCGStab (`tol=0.005, maxiter=20`) | XianzhiKrylov GPU 细层原序 `BSR3` 块 ILU(0)（`pre_smooth=0, post_smooth=1`） | CPRW True-IMPES + 36 个井压力未知量；XianzhiKrylov GPU GMRES(10) (`tol=1e-6, maxiter=100`) + LinglongAMG 单次 V-cycle | GPU 储层 `BSR3` SpMV + 修复版 `WellContributionsCuda` 独立标准井 Schur 补 | CUDA + XianzhiKrylov + LinglongAMG（缺任一后端均在配置阶段报错，不静默降级） |

各路线的源码边界与调度契约详见 [docs/ROUTES.md](docs/ROUTES.md)。

---

## 2. 测试平台与已核对组件版本

下表所有版本均从项目构建清单、运行日志和头文件直接核对：

| 组件 | 已核对版本 / 标识 | 核对来源 |
| --- | --- | --- |
| OPM (`opm-simulators`) | `2026.04`（标签 `release/2026.04/final`，提交 `b82f21db`） | Git 提交与 Flow 运行日志 |
| Dune-ISTL | `2.11.0` | `dune-istl-config.hh` |
| CUDA 工具链 | `12.8.1` / `13.1` / `13.2.1` | 归档 `environment.log` 与编译日志 |
| NVIDIA 驱动（RTX 5090 节点） | `610.57.04`（CUDA Driver `13.3`） | 节点 `nvidia-smi` 日志 |
| LinglongAMG | `26.46` / 支持 SM120 单 V-cycle 的 NVIDIA SDK | 运行日志与头文件契约 |
| XianzhiKrylov | `26.48` OPM CPR/ILU 组件库 | 构建清单与 `XianzhiKrylov_CprComponents.h` |
| HYPRE（对照） | `3.2.0` (`HYPRE_RELEASE_NUMBER=30200`) | `HYPRE_config.h` |
| AMGX（对照） | `2.5.0` | `amgx/run.log` |
| Ginkgo（对照） | `1.11.0` (`libginkgo.so.1.11.0`) | `OPM-NORNE-ILU-ONLY` 归档校验日志 |
| CPU 测试平台 | 双路 Intel Xeon CPU Max 9462（64 核） | `lscpu` 运行环境日志 |
| GPU 测试平台 | 单张 NVIDIA GeForce RTX 5090（约 32 GiB 显存） | `nvidia-smi` 与 Slurm 日志 |

---

## 3. 目录结构

本仓库仅保存轻量级 overlay 源码、最小统一补丁、公开配置、构建/验证脚本与脱敏文档，**不复制整套上游 OPM 源码，也不包含任何闭源 `.so`、`.a`、可执行文件或授权文件**：

```text
OPM-with-Solved-Tech/
  README.md
  .gitignore
  docs/
    BUILDING.md                         # 各路线依赖与详细编译说明
    ROUTES.md                           # 三条路线数学路径、源码边界与接口契约
    TESTING.md                          # 静态检查、组件烟雾测试、首日验算与计时口径
    RESULTS.md                          # 分路线历史可追溯实测结果与作业号索引
  overlay/
    COPYING                             # 上游 OPM GPL-3.0 许可正文
    common/                             # 公共 GPU 井算子适配 (NorneWellDevice.hpp)
    linglongamg/                        # LinglongAMG 压力预条件器与 26.46 兼容层
    xianzhikrylov/                      # 全耦合 BSR3 BiCGStab + 块 ILU(0) 适配器
    cpr/                                # 双后端 GPU CPRW 适配器
    patches/opm-2026.04-solved-tech.patch
  examples/
    configs/                            # 公开示例 JSON（仅用相对路径与环境变量）
    norne/                              # Norne N64 输入下载、模型说明与首日参考摘要
  manifests/
    source.json                         # 上游 OPM 与各对比组件版本清单
    sdk-contract.json                   # 三条路线所需的 SDK 头文件、库符号与隔离契约
  scripts/
    common.sh                           # 公共校验与路线辅助函数
    check-env.sh                        # 按路线检查工具链、SDK、符号与授权
    prepare-source.sh                   # 复制上游源码、打统一补丁并按路线覆盖文件
    configure.sh                        # 统一配置入口 (--route=linglongamg|xianzhikrylov|cpr)
    build.sh                            # 统一构建入口并验证 ldd 无无关后端链接
    run-case.sh                         # 统一运行入口
    run-tests.sh                        # 宿主静态校验与烟雾测试入口
    verify-public-tree.sh               # 公开仓库禁传文件、密钥与私有路径扫描
  tests/
    config/validate_configs.py          # JSON 配置与清单校验
    smoke/                              # CUDA 单组件与隔离验证测试
    norne-day1/                         # Norne 首日摘要校验脚本
```

---

## 4. 环境准备与分模式依赖

**SDK 动态库与授权获取说明**：LinglongAMG 与 XianzhiKrylov 的 SDK 头文件、动态库（`.so`）和机器绑定授权（`.license.key` / `.lic`）均由用户通过正式交付渠道单独取得，本公开仓库不分发任何闭源库或授权文件。

所有配置与脚本均通过标准环境变量引用外部路径（禁止在公开配置中硬编码私有绝对路径）：

- `OPM_SOURCE_ROOT`：未修改的上游 `opm-simulators` 2026.04 源码根目录
- `CUDA_HOME`：CUDA Toolkit 安装根目录
- `LINGLONGAMG_ROOT`：LinglongAMG SDK 根目录（含 `include/` 或 `lib/linglongamg.h`、`lib/liblinglongamg.so` 及授权）
- `XIANZHIKRYLOV_ROOT`：XianzhiKrylov SDK 根目录（含 `include/XianzhiKrylov_CprComponents.h`、`lib/libXianzhiKrylov.so` 及授权）

### 各模式依赖矩阵

| 使用模式 | 必需依赖 | 不需要 / 不链接 |
| --- | --- | --- |
| **无 GPU（纯宿主静态检查）** | Bash、Python 3、Git/`patch` | 无需 CUDA、无需 LinglongAMG、无需 XianzhiKrylov、无需授权 |
| **LinglongAMG-only** | `OPM_SOURCE_ROOT`、OPM 依赖（Dune 2.11、MPI 等）、`CUDA_HOME`、`LINGLONGAMG_ROOT` + LinglongAMG 授权 | 不需要 `XIANZHIKRYLOV_ROOT`，不链接 `libXianzhiKrylov`，不要求 XianzhiKrylov 授权 |
| **XianzhiKrylov-only** | `OPM_SOURCE_ROOT`、OPM 依赖、`CUDA_HOME`、`XIANZHIKRYLOV_ROOT` + XianzhiKrylov 授权 | 不需要 `LINGLONGAMG_ROOT`，不链接 `liblinglongamg`，不要求 LinglongAMG 授权 |
| **CPR 双后端** | `OPM_SOURCE_ROOT`、OPM 依赖、`CUDA_HOME`、`LINGLONGAMG_ROOT`、`XIANZHIKRYLOV_ROOT` 及两套授权 | 缺任一后端在 `configure.sh` / CMake 阶段直接失败，不静默降级 |

---

## 5. 源码预览与 Overlay 应用方式

本仓库不复制完整的 `opm-simulators` 代码树。`scripts/prepare-source.sh` 执行以下步骤：

1. 将 `$OPM_SOURCE_ROOT` 复制到 `work/source-<route>/`（保持用户原始 OPM 源码只读）。
2. 应用唯一最小补丁 `overlay/patches/opm-2026.04-solved-tech.patch`，仅改动 5 个上游文件：
   - `CMakeLists.txt`：增加 `OPM_SOLVED_TECH_ROUTE`（`none|linglongamg|xianzhikrylov|cpr`）选项，按路线隔离编译宏、`.cu` 源码与链接库。
   - `opm/simulators/linalg/StandardPreconditioners_serial.hpp`：仅在 `linglongamg` 路线下向 `PreconditionerFactory` 注册 `"linglongamg"` 压力预条件器。
   - `opm/simulators/linalg/ISTLSolver.hpp` 与 `ISTLSolver.cpp`：按编译宏和运行时环境变量显式路由至 `NorneIluSolver`（`xianzhikrylov`）或 `NorneCprwSolver`（`cpr`），未启用时完全保持 OPM 原有默认路径与公共 ABI。
   - `opm/simulators/linalg/gpubridge/cuda/cuWellContributions.cu`：修复每口井射孔条目超过 warp 大小时的循环覆盖及同单元原子累加（`atomicAdd`）。
3. 仅复制所选路线对应的 `overlay/<route>/`（以及 `xianzhikrylov`/`cpr` 共用的 `overlay/common/`）头文件与设备代码，确保源码树内也不存在无关后端文件。

---

## 6. 配置、构建与运行命令

运行前先[下载并解压 Norne N64 输入包](examples/norne/README.md#simulation-input)。将 `NORNE_N64.DATA` 和 `GRID_PROPS.INC` 保存在同一目录。默认运行完整历史调度，短时测试可按算例说明设置结束日期。

### 纯宿主静态检查（无需 GPU 与 SDK）

```bash
bash scripts/verify-public-tree.sh
bash scripts/run-tests.sh --static-only
```

如本机同时有干净的 `OPM_SOURCE_ROOT`，可附加 `--opm-source="$OPM_SOURCE_ROOT"` 验证三条路线的补丁干净应用与目录隔离。

### 路线一：LinglongAMG-only

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_ROOT=/path/to/LinglongAMG
# 注意：无需设置 XIANZHIKRYLOV_ROOT

bash scripts/configure.sh --route=linglongamg
bash scripts/build.sh
bash scripts/run-case.sh --route=linglongamg \
  --deck=/path/to/Norne_N64/NORNE_N64.DATA \
  --config=examples/configs/linglongamg.json
```

### 路线二：XianzhiKrylov-only

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov
# 注意：无需设置 LINGLONGAMG_ROOT

bash scripts/configure.sh --route=xianzhikrylov
bash scripts/build.sh
bash scripts/run-case.sh --route=xianzhikrylov \
  --deck=/path/to/Norne_N64/NORNE_N64.DATA \
  --config=examples/configs/xianzhikrylov.json \
  --profile --true-residual
```

### 路线三：CPR 双后端

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_ROOT=/path/to/LinglongAMG
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov

bash scripts/configure.sh --route=cpr
bash scripts/build.sh
bash scripts/run-case.sh --route=cpr \
  --deck=/path/to/Norne_N64/NORNE_N64.DATA \
  --config=examples/configs/cpr.json \
  --no-output --profile --true-residual
```

详细构建参数（`CMAKE_PREFIX_PATH`、CUDA 架构 `--cuda-arch=120`、许可证路径覆盖等）见 [docs/BUILDING.md](docs/BUILDING.md)。

---

## 7. 测试方法与计时口径

1. **宿主静态检查**：`scripts/run-tests.sh --static-only` 检查 Shell 语法、JSON 配置有效性、`--help` 入口、非法路线报错、无 SDK 预检失败信息以及 `scripts/verify-public-tree.sh` 合规扫描。
2. **单组件与隔离检查**：`tests/smoke/check-linglong-cycle.cu` 与 `tests/smoke/check-xianzhikrylov-ilu.cu` 在 GPU 上分别验证单 V-cycle 与 BSR3 BiCGStab/GMRES/ILU(0)，并通过 `/proc/self/maps` 确认未加载对方动态库。
3. **计时口径严格分离**：
   - `linglongamg` 路线记录 OPM 累计 `Linear setup`/`Linear solve`、端到端 `WALL_SECONDS`（含磁盘输出）以及单独抽出的纯压力 API 耗时（`OPM_AMG_TIMING`，已同步 GPU 且不含宿主 CPR/ILU 与数据传输）。
   - `xianzhikrylov` 路线（无输出计时）记录禁用 Eclipse 文件写入下的端到端时间、OPM 线性总计（含准备与传输）、外层迭代及 `BSR3` 块 ILU(0) 分解/应用耗时。
   - `cpr` 路线（配对无输出计时）在同一 RTX 5090 节点单核固定掩码下对比 1 MPI × 1 线程 CPU CPRW，记录无文件端到端、线性总计及嵌套分阶段耗时（注意：AMG V-cycle 计时嵌套于压力 GMRES，压力 GMRES 嵌套于外层纯求解，各分项不得直接相加）。

完整测试规程、残差验证标准与故障处理见 [docs/TESTING.md](docs/TESTING.md)。

---

## 8. 已核对历史测试结果摘要（按路线分开记录）

算例均为 **Norne N64（4×4×4 网格加细）**：2,843,196 个活跃单元，8,529,588 阶黑油主系统（`BSR3`, FP64）。以下所有数据均源自各自任务归档日志，完整明细与作业号见 [docs/RESULTS.md](docs/RESULTS.md)。

原始 Norne 角点网格的上游位置，以及 N64 `.EGRID` 的规格、校验和与发布位置，见 [Norne N64 算例说明](examples/norne/README.md#corner-point-grid-files)。

### A. LinglongAMG 路线（混合路径：CPU 外层 + CPU ILU(0) + GPU 压力 GMRES/AMG）

- **首日（Jobs `352448` / `353062`）**：OPM 线性阶段 `setup 2.83 s + solve 1.55 s = 4.38 s`（32 次 GMRES 迭代，9 层 AMG，4 次压力调用），同口径 CPU64/HYPRE CPU 参考为 `3.03 + 3.99 = 7.02 s`。
- **176 天重测（压力 API 独立计时，Jobs `354839` / `354933` / `354826`）**：
  - 纯压力 API（Setup + GMRES + AMG，排除 CPU 装配/细层 ILU/传输）：CPU64 HYPRE CPU `225.33 s`（6,385 次迭代） vs RTX 5090 LinglongAMG `33.73 s`（5,611 次迭代） vs RTX 5090 AMGX `63.46 s`（5,043 次迭代） vs RTX 5090 HYPRE GPU `43.16 s`（4,761 次迭代）。
  - 同次运行 OPM 整套应用墙钟（含 1 线程 CPU 细层 ILU(0) 与读/写开销）：CPU64 `547.26 s`，5090 LinglongAMG `1029.62 s`，5090 AMGX `1042.70 s`，5090 HYPRE `948.74 s`。

### B. XianzhiKrylov 路线（全耦合 `BSR3` BiCGStab + 单级 ILU(0)，无 CPR/AMG）

- **首日无输出配对计时（Jobs `364339` / `364340` / `364338`，3 时间步，11 Newton）**：
  - OPM 线性总计：CPU64（64 MPI Xeon Max 9462）`26.56 s`（450 次迭代） vs RTX 5090 Ginkgo 标量 CSR ILU(0) `76.03 s`（363 次迭代） vs RTX 5090 XianzhiKrylov `BSR3` 块 ILU(0) `20.71 s`（366 次迭代）。
  - 端到端无输出耗时：`194.95 s`（CPU64） vs `193.31 s`（Ginkgo） vs `138.81 s`（XianzhiKrylov）。
  - ILU 分解 / 应用：Ginkgo `10.82 s / 61.20 s` vs XianzhiKrylov `0.34 s / 16.16 s`。
  - 独立真相对残差（审核作业 `364337`/`364329`/`364328`）均收敛且物性汇总一致。
- **176 天状态**：CPU64（Job `364343`）完整跑完 176 天（1347.65 s，线性 996.87 s）；两条 GPU 路线受 RTX 5090 分区 20 分钟墙钟上限（1170 s 进程超时，退出码 124）阻塞，分别推进至 56.00 天（XianzhiKrylov，Job `364341`）与 14.93 天（Ginkgo，Job `364342`），未取得完整 176 天 GPU 对比样本。

### C. CPR 双后端路线（GPU 外层 BiCGStab + GPU 压力 GMRES + LinglongAMG V-cycle + GPU 细层块 ILU(0)）

- **首日同节点无输出配对计时（Job `356802`，RTX 5090 节点 1 MPI × 1 线程，1 时间步，3 Newton，3 个线性系统）**：
  - 无模拟输出端到端：`99.09 s`（CPU CPRW）→ `83.39 s`（GPU CPR，`1.19×`）。
  - OPM 线性阶段（含准备）：`19.10 s` → `3.67 s`（`5.20×`）。
  - 纯外层求解（含 CPRW 内层应用）：`7.656 s` → `0.910 s`（`8.42×`）。
  - GPU 嵌套阶段拆解：压力构建 `1.256 s`，压力 GMRES `0.715 s`（含 LinglongAMG V-cycle API `0.657 s`），AMG setup `0.287 s`，细层块 ILU(0) 分解 `0.099 s` / 应用 `0.152 s`，储层 `BSR3` SpMV `0.032 s`，GPU 井算子 `0.0027 s`。
  - 数值验算：独立含井算子相对误差 `8.96e-17`，CPRW 完整修正相对误差 `2.99e-15`，3 次外层真相对残差依次为 `3.7489e-4`、`8.3250e-5`、`1.5460e-3`（均 `< 0.005`），首日末物性 6 项汇总在摘要精度下完全一致。
- **未验证项**：`cpr` 路线尚未完成 176 天长窗实测；本仓库整合补丁在本次开源整理后仅完成了宿主静态检查与目录隔离检查，未在带真实 SDK 的 GPU 节点重新执行编译与首日复测。

---

## 9. 支持范围与已知限制

- **适用问题范围**：`xianzhikrylov` 与 `cpr` 的 GPU 路径当前面向单 MPI 进程、FP64 `BSR3`（3 相黑油）及标准井（`WellContributionsCuda`，遇到多段井时显式报错而不静默回退）。
- **LinglongAMG-only 限制**：`LinglongAMGPreconditioner` 仅处理单节点串行 `1x1` 标量 FP64 压力系统，32 位 CSR 索引。
- **CPR 单 V-cycle 契约**：`cpr` 路线优先调用 `linglongamg_apply_cycle`；若所用 SDK 未导出该符号，可设置 `NORNE_AMG_ALLOW_SINGLE_STEP=1` 配合 `examples/configs/linglong-vcycle.json`（`method=amg, cycle=v, max_iters=1`）并通过 `linglongamg_get_solve_info` 校验每次调用恰好执行 1 步。

---

## 10. 上游组件来源与许可说明

- **OPM (`opm-simulators` 2026.04)**：基于 GNU General Public License v3.0 or later（GPL-3.0-or-later）发布；本仓库 `overlay/` 中基于 OPM 的修改保留 GPL-3.0 声明，许可正文见 [overlay/COPYING](overlay/COPYING)。[OPM 项目仓库](https://github.com/OPM/opm-simulators)
- **Dune-ISTL (`2.11.0`)**：GPL-2.0 with runtime exception。[Dune 项目](https://dune-project.org/)
- **HYPRE (`3.2.0`)**：Apache-2.0 OR MIT（用于历史 CPU/GPU 压力对照）。[HYPRE](https://github.com/hypre-space/hypre)
- **NVIDIA AMGX (`2.5.0`)**：BSD-3-Clause（用于历史 GPU 压力对照）。[AMGX](https://github.com/NVIDIA/AMGX)
- **Ginkgo (`1.11.0`)**：BSD-3-Clause（用于历史单级 ILU(0) GPU 对照）。[Ginkgo](https://github.com/ginkgo-project/ginkgo)
- **LinglongAMG / XianzhiKrylov**：闭源商业/专有加速库及其授权协议，不包含在本仓库中，需由用户从正式交付渠道获取并按对应许可使用。
