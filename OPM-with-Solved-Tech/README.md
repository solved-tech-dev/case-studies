# OPM Flow 2026.04 接入 LinglongAMG 与 XianzhiKrylov

本目录提供基于 **OPM Flow 2026.04**（`release/2026.04/final`）的公开应用集成层（Overlay + 最小统一补丁），允许用户通过统一入口按需选择以下**两条互不耦合的构建与求解路线**：

```text
                        OPM Flow 2026.04 (Black-Oil BSR3 / CPRW)
                                          │
                  ┌───────────────────────┴───────────────────────┐
                  ▼                                               ▼
         route=linglongamg                               route=xianzhikrylov
       (LinglongAMG 混合路线)                        (XianzhiKrylov 公开 API 全耦合路线)
                  │                                               │
        CPU 外层 BiCGStab                            GPU 显式全耦合 BSR3 BiCGStab
        + CPU paroverilu0                            + GPU 原生 BSR3 块 ILU(0) L0
        + CPRW True-IMPES                            + OPM 显式井 Schur 矩阵注入
        + GPU GMRES(10)/AMG (LL)                     (仅调用 XianzhiKrylovAPI.h)
        [不链接 XianzhiKrylov]                       [不链接 LinglongAMG]
```

> [!IMPORTANT]
> **版本与路线说明**：此前公开版本中曾出现的预条件器单独 `apply` 接口与 `cpr` 组合路线依赖于未包含在 XianzhiKrylov `26.48` 发行包中的内部实验接口。当前仓库已移除 `cpr` 路线，**OPM + XianzhiKrylov 26.48 仅支持通过公开 API（`XianzhiKrylovAPI.h`）直接求解显式全耦合矩阵**；预条件器单独 `apply` 功能请等待 `26.49` 版本更新后再行使用。此外，两条路线对应不同的数学路径与对比基线，禁止跨路线直接对比耗时。

---

## 1. 两条路线与线性系统分工

| 路线 (`--route`) | 外层迭代与矩阵 | 细层预条件器 | 压力子系统 / AMG | 井贡献处理 | 编译/链接依赖 |
| --- | --- | --- | --- | --- | --- |
| `linglongamg` | 原生 OPM CPU BiCGStab（全系统，默认 `tol=0.005, maxiter=20`） | 原生 OPM CPU `paroverilu0` (`relaxation=1.0`) | CPRW True-IMPES 压力 CSR；GPU 原生 GMRES(10) + AMG V-cycle (`tol=1e-6, maxiter=100`) 由 `liblinglongamg.so` 承担 | 原生 OPM `WellModelOperator` (`--matrix-add-well-contributions=false`) | 仅 CUDA + LinglongAMG（**不链接** XianzhiKrylov，不需要其授权） |
| `xianzhikrylov` | XianzhiKrylov 公开 C API v3 GPU 全耦合 `BSR3` BiCGStab (`tol=0.005, max_iterations=400`) | XianzhiKrylov GPU 原生 `BSR3` 块 ILU(0) `level=0` | 不启用 CPR / 不调用 AMG | OPM 显式将标准井 Schur 补注入主矩阵（`--matrix-add-well-contributions=true`） | 仅 CUDA + XianzhiKrylov `26.48`（**不链接** LinglongAMG，不需要其授权） |

各路线的源码边界与调度契约详见 [docs/ROUTES.md](docs/ROUTES.md)。

---

## 2. 测试平台与已核对组件版本

| 组件 | 已核对版本 / 标识 | 核对来源 |
| --- | --- | --- |
| OPM (`opm-simulators`) | `2026.04`（标签 `release/2026.04/final`，提交 `b82f21db`） | Git 提交与 Flow 运行日志 |
| Dune-ISTL | `2.11.0` | `dune-istl-config.hh` |
| CUDA 工具链 | `12.8.1` / `13.1` / `13.2.1` | 编译与节点运行日志 |
| NVIDIA 驱动（RTX 5090 节点） | `610.57.04`（CUDA Driver `13.3`） | 节点 `nvidia-smi` 日志 |
| LinglongAMG | `26.46` / 支持 SM120 单 V-cycle 的 NVIDIA SDK | 运行日志与头文件契约 |
| XianzhiKrylov | `26.48` 公开 C API v3（`XianzhiKrylovAPI.h` + `XianzhiKrylovData.h`） | `26.48` 发行包与头文件 |
| HYPRE（对照） | `3.2.0` (`HYPRE_RELEASE_NUMBER=30200`) | `HYPRE_config.h` |
| AMGX（对照） | `2.5.0` | `amgx/run.log` |
| CPU 测试平台 | 双路 Intel Xeon CPU Max 9462（64 核） | `lscpu` 运行环境日志 |
| GPU 测试平台 | 单张 NVIDIA GeForce RTX 5090（约 32 GiB 显存） | `nvidia-smi` 与 Slurm 日志 |

---

## 3. 目录结构

本目录包含 overlay 源码、统一补丁、公开配置、构建和验证脚本、文档以及 Norne N64 公开算例数据，**不复制整套上游 OPM 源码，也不包含任何闭源 `.so`、`.a`、可执行文件或授权文件**：

```text
OPM-with-Solved-Tech/
  README.md
  .gitignore
  docs/
    BUILDING.md                         # 各路线依赖与详细编译说明
    ROUTES.md                           # 两条路线数学路径、源码边界与接口契约
    TESTING.md                          # 静态检查、组件烟雾测试、首日验算与计时口径
    RESULTS.md                          # 分路线实测结果与作业号索引
  overlay/
    COPYING                             # 上游 OPM GPL-3.0 许可正文
    linglongamg/                        # LinglongAMG 压力预条件器与 26.46 兼容层
    xianzhikrylov/                      # 公开 API 显式全耦合 BSR3 适配器 (NorneXianzhi.hpp)
    patches/opm-2026.04-solved-tech.patch
  examples/
    configs/                            # 公开示例 JSON（仅用相对路径与环境变量）
    norne/                              # Norne N64 模型说明与参考摘要
  norne-n64-grid-v1/                    # Norne N64 Eclipse 输入与角点网格下载
  manifests/
    source.json                         # 上游 OPM 与各对比组件版本清单
    sdk-contract.json                   # 两条路线所需的公开 SDK 头文件、库符号与隔离契约
  scripts/
    common.sh                           # 公共校验与路线辅助函数
    check-env.sh                        # 按路线检查工具链、公开 SDK、符号与授权
    prepare-source.sh                   # 复制上游源码、打统一补丁并按路线覆盖文件
    configure.sh                        # 统一配置入口 (--route=linglongamg|xianzhikrylov)
    build.sh                            # 统一构建入口并验证 ldd 无无关后端链接
    run-case.sh                         # 统一运行入口
    run-tests.sh                        # 宿主静态校验与烟雾测试入口
    verify-public-tree.sh               # 公开仓库禁传文件、密钥与私有路径扫描
  tests/
    config/validate_configs.py          # JSON 配置与清单校验
    smoke/                              # 单组件与隔离验证测试
    norne-day1/                         # Norne 摘要校验脚本
```

---

## 4. 环境准备与分模式依赖

**SDK 动态库与授权获取说明**：LinglongAMG 与 XianzhiKrylov 的 SDK 头文件、动态库（`.so`）和机器绑定授权（`.license.key` / `.key` / `.lic`）均由用户通过正式交付渠道单独取得，本公开仓库不分发任何闭源库或授权文件。

所有配置与脚本均通过标准环境变量引用外部路径（禁止在公开配置中硬编码私有绝对路径）：

- `OPM_SOURCE_ROOT`：未修改的上游 `opm-simulators` 2026.04 源码根目录
- `CUDA_HOME`：CUDA Toolkit 安装根目录（`26.48` 需 CUDA 13 runtime）
- `LINGLONGAMG_ROOT`：LinglongAMG SDK 根目录（含 `include/` 或 `lib/linglongamg.h`、`lib/liblinglongamg.so` 及授权）
- `XIANZHIKRYLOV_ROOT`：XianzhiKrylov `26.48` SDK 根目录（含 `include/XianzhiKrylovAPI.h`、`include/XianzhiKrylovData.h`、`lib/libXianzhiKrylov.so` 及授权）

### 各模式依赖矩阵

| 使用模式 | 必需依赖 | 不需要 / 不链接 |
| --- | --- | --- |
| **无 GPU（纯宿主静态检查）** | Bash、Python 3、Git/`patch` | 无需 CUDA、无需 LinglongAMG、无需 XianzhiKrylov、无需授权 |
| **LinglongAMG-only** | `OPM_SOURCE_ROOT`、OPM 依赖（Dune 2.11、MPI 等）、`CUDA_HOME`、`LINGLONGAMG_ROOT` + LinglongAMG 授权 | 不需要 `XIANZHIKRYLOV_ROOT`，不链接 `libXianzhiKrylov`，不要求 XianzhiKrylov 授权 |
| **XianzhiKrylov-only** | `OPM_SOURCE_ROOT`、OPM 依赖、`CUDA_HOME`、`XIANZHIKRYLOV_ROOT` (`26.48`) + XianzhiKrylov 授权 | 不需要 `LINGLONGAMG_ROOT`，不链接 `liblinglongamg`，不要求 LinglongAMG 授权 |

---

## 5. 源码预览与 Overlay 应用方式

本仓库不复制完整的 `opm-simulators` 代码树。`scripts/prepare-source.sh` 执行以下步骤：

1. 将 `$OPM_SOURCE_ROOT` 复制到 `work/source-<route>/`（保持用户原始 OPM 源码只读）。
2. 应用唯一最小补丁 `overlay/patches/opm-2026.04-solved-tech.patch`，仅改动 4 个上游文件：
   - `CMakeLists.txt`：增加 `OPM_SOLVED_TECH_ROUTE`（`none|linglongamg|xianzhikrylov`）选项，按路线隔离编译宏、CUDA 优先头文件路径与链接库（无需开启 `USE_GPU_BRIDGE`）。
   - `opm/simulators/linalg/StandardPreconditioners_serial.hpp`：仅在 `linglongamg` 路线下向 `PreconditionerFactory` 注册 `"linglongamg"` 压力预条件器。
   - `opm/simulators/linalg/ISTLSolver.hpp` 与 `ISTLSolver.cpp`：提供可选的性能与独立残差校验统计，并在 `xianzhikrylov` 路线下显式路由至 `NorneXianzhi`（要求单 MPI 与 `--matrix-add-well-contributions=true`）；未启用时完全保持 OPM 原有默认路径与公共 ABI。
3. 仅复制所选路线对应的 `overlay/<route>/` 头文件，确保源码树内也不存在无关后端文件。

---

## 6. 配置、构建与运行命令

运行前先[下载并解压 Norne N64 输入包](norne-n64-grid-v1/README.md#downloads)。将 `NORNE_N64.DATA` 和 `GRID_PROPS.INC` 保存在同一目录。默认运行完整历史调度，短时测试可按算例说明设置结束日期。

### 纯宿主静态检查（无需 GPU 与 SDK）

```bash
bash scripts/verify-public-tree.sh
bash scripts/run-tests.sh --static-only
```

如本机同时有干净的 `OPM_SOURCE_ROOT`，可附加 `--opm-source="$OPM_SOURCE_ROOT"` 验证两条路线的补丁干净应用与目录隔离。

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

### 路线二：XianzhiKrylov-only（公开 API 显式全耦合）

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda-13.x
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov-26.48
# 注意：无需设置 LINGLONGAMG_ROOT

bash scripts/configure.sh --route=xianzhikrylov
bash scripts/build.sh
bash scripts/run-case.sh --route=xianzhikrylov \
  --deck=/path/to/Norne_N64/NORNE_N64.DATA \
  --config=examples/configs/xianzhikrylov.json \
  --xk-config=examples/configs/xk-bsr3-bicgstab-ilu0.json \
  --no-output --profile --true-residual
```

详细构建参数（`CMAKE_PREFIX_PATH`、CUDA 架构 `--cuda-arch=120`、许可证路径覆盖等）见 [docs/BUILDING.md](docs/BUILDING.md)。

---

## 7. 测试方法与计时口径

1. **宿主静态检查**：`scripts/run-tests.sh --static-only` 检查 Shell 语法、JSON 配置有效性、`--help` 入口、非法路线（含已移除的 `cpr`）报错、无 SDK 预检失败信息以及 `scripts/verify-public-tree.sh` 合规扫描。
2. **单组件与隔离检查**：`tests/smoke/check-linglong-cycle.cu` 与纯 C++ 的 `tests/smoke/check-xianzhikrylov-bsr3.cpp` 在 GPU 上分别验证单 V-cycle 与公开 API `BSR3` BiCGStab/ILU(0)，并通过 `/proc/self/maps` 确认未加载对方动态库。
3. **计时口径严格分离**：
   - `linglongamg` 路线记录 OPM 累计 `Linear setup`/`Linear solve`、端到端 `WALL_SECONDS`（含磁盘输出）以及单独抽出的纯压力 API 耗时（`OPM_AMG_TIMING`，已同步 GPU 且不含宿主 CPR/ILU 与数据传输）。
   - `xianzhikrylov` 路线通过 `OPM_XK_ANALYSIS`、`OPM_XK_PREPARE`、`OPM_XK_LINEAR` 与 `OPM_SOLVED_TECH_TRUE_RESIDUAL` 分项输出结构分析、打包、传输、ILU 分解、纯求解及独立真相对残差。

完整测试规程、残差验证标准与故障处理见 [docs/TESTING.md](docs/TESTING.md)。

---

## 8. 已核对实测结果摘要（按路线分开记录）

算例均为 **Norne N64（4×4×4 网格加细）**：2,843,196 个活跃单元，8,529,588 阶黑油主系统（`BSR3`, FP64）。以下数据均源自实机 Slurm 日志，完整明细见 [docs/RESULTS.md](docs/RESULTS.md)。

原始 Norne 角点网格的上游位置，以及 N64 `.EGRID` 的规格、校验和与发布位置，见 [Norne N64 算例说明](examples/norne/README.md#corner-point-grid-files)。

### A. LinglongAMG 路线（混合路径：CPU 外层 + CPU ILU(0) + GPU 压力 GMRES/AMG）

- **首日（Jobs `352448` / `353062`）**：OPM 线性阶段 `setup 2.83 s + solve 1.55 s = 4.38 s`（32 次 GMRES 迭代，9 层 AMG，4 次压力调用），同口径 CPU64/HYPRE CPU 参考为 `3.03 + 3.99 = 7.02 s`。
- **176 天重测（压力 API 独立计时，Jobs `354839` / `354933` / `354826`）**：
  - 纯压力 API（Setup + GMRES + AMG，排除 CPU 装配/细层 ILU/传输）：CPU64 HYPRE CPU `225.33 s`（6,385 次迭代） vs RTX 5090 LinglongAMG `33.73 s`（5,611 次迭代） vs RTX 5090 AMGX `63.46 s`（5,043 次迭代） vs RTX 5090 HYPRE GPU `43.16 s`（4,761 次迭代）。
  - 同次运行 OPM 整套应用墙钟（含 1 线程 CPU 细层 ILU(0) 与读/写开销）：CPU64 `547.26 s`，5090 LinglongAMG `1029.62 s`，5090 AMGX `1042.70 s`，5090 HYPRE `948.74 s`。

### B. XianzhiKrylov 路线（26.48 公开 API 显式全耦合 `BSR3` BiCGStab + 块 ILU(0) L0）

- **纯 C++ 烟雾测试 `check-xianzhikrylov-bsr3`（Job `367086` / `367095`，RTX 5090）**：
  - 使用 `26.48` 发行包 `libXianzhiKrylov.so.26.48.0`（SHA-256 `8661b5e5...`）编译并运行通过，两次数值制备+求解独立真实相对残差 `1.40e-16` / `1.62e-16`。
- **OPM Flow Norne N64 首日实测（Jobs `367087` / `367095`，RTX 5090，`--matrix-add-well-contributions=true`，`--no-output`）**：
  - 含显式井 Schur 填充后主矩阵为 `2,843,196` 块行、`55,065,494` 非零块。
  - 作业 `367095`（`max_iterations=400`）零报错跑通 OPM→`XianzhiKrylovAPI.h` 全流程，完整收敛时间步 0（`0.33` 天，3 次 Newton，37 次线性迭代）及时间步 1 前两次 Newton，独立真相对残差（`0.002393`–`0.004943`）与库报告残差吻合至 12 位以上。
  - 因 5,506 万块显式含井矩阵的 `BSR3` ILU(0) `level=0` 单次数值分解耗时约 `70.82 s`（结构分析 `0.71 s`、打包 `0.19 s`、上传 `0.17 s`），作业于时间步 1 第 11 次线性系统准备阶段达到 `1050 s` 超时上限（退出码 `124`）。

---

## 9. 支持范围与已知限制

- **适用问题范围**：`xianzhikrylov` 的 GPU 路径当前面向单 MPI 进程、FP64 `BSR3`（3 相黑油），要求 `--matrix-add-well-contributions=true`（由 OPM 将标准井 Schur 补直接注入主矩阵）。
- **LinglongAMG-only 限制**：`LinglongAMGPreconditioner` 仅处理单节点串行 `1x1` 标量 FP64 压力系统，32 位 CSR 索引。
- **CPR / 单独 `apply` 说明**：已移除依赖内部符号的 `cpr` 路线，待 XianzhiKrylov `26.49` 正式提供预条件器单独 `apply` 接口后再行评估接入。

---

## 10. 上游组件来源与许可说明

- **OPM (`opm-simulators` 2026.04)**：基于 GNU General Public License v3.0 or later（GPL-3.0-or-later）发布；本仓库 `overlay/` 中基于 OPM 的修改保留 GPL-3.0 声明，许可正文见 [overlay/COPYING](overlay/COPYING)。[OPM 项目仓库](https://github.com/OPM/opm-simulators)
- **Dune-ISTL (`2.11.0`)**：GPL-2.0 with runtime exception。[Dune 项目](https://dune-project.org/)
- **HYPRE (`3.2.0`)**：Apache-2.0 OR MIT（用于历史 CPU/GPU 压力对照）。[HYPRE](https://github.com/hypre-space/hypre)
- **NVIDIA AMGX (`2.5.0`)**：BSD-3-Clause（用于历史 GPU 压力对照）。[AMGX](https://github.com/NVIDIA/AMGX)
- **LinglongAMG / XianzhiKrylov**：闭源商业/专有加速库及其授权协议，不包含在本仓库中，需由用户从正式交付渠道获取并按对应许可使用。
