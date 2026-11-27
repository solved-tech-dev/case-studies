# OPM Flow 2026.04 接入 LinglongAMG 与 XianzhiKrylov

本目录提供基于 **OPM Flow 2026.04**（`release/2026.04/final`）的公开应用集成层（Overlay + 最小统一补丁），通过统一入口支持两条相互独立的 GPU 线性求解加速路线：

```text
                        OPM Flow 2026.04 (Black-Oil BSR3 / CPRW)
                                          │
                  ┌───────────────────────┴───────────────────────┐
                  ▼                                               ▼
         route=linglongamg                               route=xianzhikrylov
     (LinglongAMG 压力子系统路线)                      (XianzhiKrylov 全耦合求解路线)
                  │                                               │
        CPU 外层 BiCGStab                            GPU 显式全耦合 BSR3 BiCGStab
        + CPU paroverilu0                            + GPU 原生 BSR3 块 ILU(0) L0
        + CPRW True-IMPES                            + 显式含井全耦合矩阵组装
        + GPU GMRES + AMG (LL)                       (基于 XianzhiKrylovAPI.h)
        [不链接 XianzhiKrylov]                       [不链接 LinglongAMG]
```

> [!NOTE]
> `linglongamg` 与 `xianzhikrylov` 两条路线面向不同的数学求解路径与评估阶段：前者针对 CPRW 压力子系统（GMRES + AMG），后者针对含井全耦合线性系统（BSR3 块 ILU(0) + BiCGStab），请按对应口径分别构建与对比。

---

## 1. 两条路线与线性系统分工

| 路线 (`--route`) | 外层迭代与矩阵 | 细层预条件器 | 压力子系统 / AMG | 井贡献处理 | 编译/链接依赖 |
| --- | --- | --- | --- | --- | --- |
| `linglongamg` | OPM CPU BiCGStab（默认 `tol=0.005, maxiter=20`） | OPM CPU `paroverilu0` (`relaxation=1.0`) | CPRW True-IMPES 压力 CSR；GPU GMRES + AMG V-cycle (`tol=1e-6, maxiter=100`) 由 `liblinglongamg.so` 求解 | OPM `WellModelOperator` (`--matrix-add-well-contributions=false`) | CUDA + LinglongAMG（**不链接** XianzhiKrylov） |
| `xianzhikrylov` | XianzhiKrylov GPU 全耦合 `BSR3` BiCGStab (`tol=0.005, max_iterations=400`) | XianzhiKrylov GPU 原生 `BSR3` 块 ILU(0) `level=0` | 不启用 CPR / AMG | OPM 显式组装含井全耦合矩阵（`--matrix-add-well-contributions=true`） | CUDA + XianzhiKrylov `26.48`（**不链接** LinglongAMG） |

各路线的源码边界与接口说明详见 [docs/ROUTES.md](docs/ROUTES.md)。

---

## 2. 环境与组件版本

| 组件 | 版本 | 说明 |
| --- | --- | --- |
| OPM (`opm-simulators`) | `2026.04`（`release/2026.04/final`） | 上游开源黑油油藏模拟器 |
| Dune-ISTL | `2.11.0` | OPM 线性求解接口基础库 |
| CUDA Toolkit | `12.8` / `13.x` | `XianzhiKrylov 26.48` 需 CUDA 13 runtime |
| LinglongAMG | `26.46` | 压力子系统 GPU AMG 加速库 |
| XianzhiKrylov | `26.48`（公开 C API v3） | 全耦合 GPU BSR3 Krylov/ILU 加速库 |
| NVIDIA AMGX（对照） | `2.5.0` | 同平台压力子系统对照后端 |
| Ginkgo（对照） | `1.11.0` | 同平台全耦合标量 CSR ILU(0) 对照后端 |
| GPU 基准平台 | NVIDIA GeForce RTX 5090 (32 GB) | 基准对比测试卡 |

---

## 3. 目录结构

本目录包含 overlay 源码、统一补丁、示例配置、构建测试脚本、文档与 Norne N64 算例数据（不包含上游整套 OPM 源码、闭源动态库或授权文件）：

```text
OPM-with-Solved-Tech/
  README.md
  .gitignore
  docs/
    BUILDING.md                         # 各路线依赖与详细编译说明
    ROUTES.md                           # 两条路线数学路径、源码边界与接口契约
    TESTING.md                          # 静态检查、组件烟雾测试与验证口径
    RESULTS.md                          # 公开基准性能测试结果
  overlay/
    COPYING                             # 上游 OPM GPL-3.0 许可正文
    linglongamg/                        # LinglongAMG 压力预条件器与兼容层
    xianzhikrylov/                      # 公开 API 显式含井全耦合 BSR3 适配器 (NorneXianzhi.hpp)
    patches/opm-2026.04-solved-tech.patch
  examples/
    configs/                            # 示例 JSON 配置（仅使用相对路径与环境变量）
    norne/                              # Norne N64 模型说明与参考基准摘要
  norne-n64-grid-v1/                    # Norne N64 Eclipse 输入与角点网格下载
  manifests/
    source.json                         # 上游 OPM 与基准组件版本清单
    sdk-contract.json                   # 两条路线所需的公开 SDK 头文件、库符号与隔离契约
  scripts/
    common.sh                           # 公共校验与路线辅助函数
    check-env.sh                        # 按路线检查工具链、SDK 头文件、符号与授权
    prepare-source.sh                   # 复制上游源码、应用统一补丁并注入所选路线文件
    configure.sh                        # 统一配置入口 (--route=linglongamg|xianzhikrylov)
    build.sh                            # 统一构建入口并校验动态库隔离
    run-case.sh                         # 统一运行入口
    run-tests.sh                        # 宿主静态校验与环境检查入口
    verify-public-tree.sh               # 公开仓库合规扫描
  tests/
    config/validate_configs.py          # JSON 配置与清单校验
    smoke/                              # 单组件 GPU 烟雾测试
    norne-day1/                         # Norne 摘要结构校验与验证说明
```

---

## 4. 环境准备与分模式依赖

**SDK 与授权获取**：LinglongAMG 与 XianzhiKrylov 的 SDK 头文件、动态库（`.so`）及设备授权文件可联系解得科技（solved-tech.com）或通过官方渠道获取，本仓库不分发闭源库与授权文件。

脚本与配置通过标准环境变量引用外部路径：

- `OPM_SOURCE_ROOT`：上游 `opm-simulators` 2026.04 源码根目录
- `CUDA_HOME`：CUDA Toolkit 安装根目录（`26.48` 需 CUDA 13 runtime）
- `LINGLONGAMG_ROOT`：LinglongAMG SDK 根目录（含 `include/` 或 `lib/linglongamg.h`、`lib/liblinglongamg.so` 及授权）
- `XIANZHIKRYLOV_ROOT`：XianzhiKrylov `26.48` SDK 根目录（含 `include/XianzhiKrylovAPI.h`、`include/XianzhiKrylovData.h`、`lib/libXianzhiKrylov.so` 及授权）

### 各模式依赖矩阵

| 使用模式 | 必需依赖 | 不需要 / 不链接 |
| --- | --- | --- |
| **宿主静态检查** | Bash、Python 3、Git/`patch` | 无需 CUDA、无需任何加速库、无需授权 |
| **LinglongAMG 路线** | `OPM_SOURCE_ROOT`、OPM 依赖（Dune 2.11、MPI 等）、`CUDA_HOME`、`LINGLONGAMG_ROOT` + 授权 | 不需要 `XIANZHIKRYLOV_ROOT`，不链接 `libXianzhiKrylov` |
| **XianzhiKrylov 路线** | `OPM_SOURCE_ROOT`、OPM 依赖、`CUDA_HOME`、`XIANZHIKRYLOV_ROOT` (`26.48`) + 授权 | 不需要 `LINGLONGAMG_ROOT`，不链接 `liblinglongamg` |

---

## 5. 集成方式与 Overlay 说明

本仓库不复制完整的 `opm-simulators` 代码树。`scripts/prepare-source.sh` 执行以下步骤：

1. 将 `$OPM_SOURCE_ROOT` 复制到 `work/source-<route>/`（保持用户原始 OPM 源码只读）。
2. 应用唯一最小补丁 `overlay/patches/opm-2026.04-solved-tech.patch`，仅改动 4 个上游文件：
   - `CMakeLists.txt`：增加 `OPM_SOLVED_TECH_ROUTE`（`none|linglongamg|xianzhikrylov`）选项，按路线隔离编译宏、CUDA 优先头文件路径与链接库（无需开启 `USE_GPU_BRIDGE`）。
   - `opm/simulators/linalg/StandardPreconditioners_serial.hpp`：仅在 `linglongamg` 路线下注册 `"linglongamg"` 压力预条件器。
   - `opm/simulators/linalg/ISTLSolver.hpp` 与 `ISTLSolver.cpp`：提供可选的性能与独立残差校验统计，并在 `xianzhikrylov` 路线下显式路由至 `NorneXianzhi`（单 MPI、`--matrix-add-well-contributions=true`）；未启用时保持 OPM 原有默认路径与公共 ABI。
3. 仅注入所选路线对应的 `overlay/<route>/` 文件，保证产物目录纯净隔离。

---

## 6. 配置、构建与运行

运行前请先[下载并解压 Norne N64 输入包](norne-n64-grid-v1/README.md#downloads)，并将 `NORNE_N64.DATA` 与 `GRID_PROPS.INC` 放在同一目录。短周期测试可按算例说明裁剪模拟结束日期。

### 宿主静态检查（无需 GPU 与 SDK）

```bash
bash scripts/verify-public-tree.sh
bash scripts/run-tests.sh --static-only
```

若本机已准备干净的 `OPM_SOURCE_ROOT`，可附加 `--opm-source="$OPM_SOURCE_ROOT"` 同时校验两条路线的补丁应用与文件隔离。

### 路线一：LinglongAMG

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda
export LINGLONGAMG_ROOT=/path/to/LinglongAMG

bash scripts/configure.sh --route=linglongamg
bash scripts/build.sh
bash scripts/run-case.sh --route=linglongamg \
  --deck=/path/to/Norne_N64/NORNE_N64.DATA \
  --config=examples/configs/linglongamg.json
```

### 路线二：XianzhiKrylov

```bash
export OPM_SOURCE_ROOT=/path/to/opm-simulators-2026.04
export CUDA_HOME=/path/to/cuda-13.x
export XIANZHIKRYLOV_ROOT=/path/to/XianzhiKrylov-26.48

bash scripts/configure.sh --route=xianzhikrylov
bash scripts/build.sh
bash scripts/run-case.sh --route=xianzhikrylov \
  --deck=/path/to/Norne_N64/NORNE_N64.DATA \
  --config=examples/configs/xianzhikrylov.json \
  --xk-config=examples/configs/xk-bsr3-bicgstab-ilu0.json \
  --no-output --profile --true-residual
```

详细构建参数（`CMAKE_PREFIX_PATH`、CUDA 架构 `--cuda-arch=120`、授权路径覆盖等）见 [docs/BUILDING.md](docs/BUILDING.md)。

---

## 7. 基准测试结果摘要（NVIDIA GeForce RTX 5090）

基准算例为 **Norne N64（4×4×4 网格加细）**：网格维度 `184×448×88`，`2,843,196` 个活跃单元，`8,529,588` 标量未知量，单张 **NVIDIA GeForce RTX 5090 (32 GB)** 实测。网格与算例详情见 [Norne N64 算例说明](examples/norne/README.md#corner-point-grid-files)，完整分段数据见 [docs/RESULTS.md](docs/RESULTS.md)。

### A. LinglongAMG（Norne 176 天压力子系统求解：GMRES + AMG）

在 CPRW 压力子系统求解阶段（`2,843,232` 行、`20,012,068` CSR 非零元），统计 176 天模拟全过程所有压力子系统建立与迭代求解的累计耗时：

- **AMGX `2.5.0`**：累计 **`63.46 s`**（Setup `33.09 s` + Solve `30.38 s`，`675` 次调用 / `5,043` 次迭代）
- **LinglongAMG `26.46`**：累计 **`33.73 s`**（Setup `6.05 s` + Solve `27.67 s`，`669` 次调用 / `5,611` 次迭代）
- **加速比**：相对同卡 AMGX 实现 **`1.88×`** 加速。

### B. XianzhiKrylov（Norne 首日第 1 步显式含井全耦合：BiCGStab + ILU(0)）

在 OPM Flow 显式含井全耦合口径下（`--matrix-add-well-contributions=true`，`2,843,196` 块行、`8,529,588` 阶、`55,065,494` 个 `3×3` 块 / `495,589,446` 标量非零元，相对残差容差 `5e-3`），对比同卡标量 CSR ILU(0) 与 BSR3 块 ILU(0)：

- **Ginkgo `1.11.0`（标量 CSR ILU(0)）**：总耗时 **`158.85 s`**（预处理准备 `153.61 s` + 迭代求解 `5.24 s`，`13` 步）
- **XianzhiKrylov `26.48`（BSR3 块 ILU(0)）**：总耗时 **`74.77 s`**（结构分析 `0.35 s` + 数值分解 `71.30 s` + 迭代求解 `3.13 s`，`14` 步）
- **加速比**：相对同卡 Ginkgo 实现 **`2.12×`** 加速。

---

## 8. 适用范围

- `linglongamg` 面向单进程标量 FP64 CPRW 压力子系统（32 位 CSR 索引）。
- `xianzhikrylov` 面向单 MPI 进程、FP64 `BSR3` 三相黑油显式含井全耦合系统（需开启 `--matrix-add-well-contributions=true`）。

---

## 9. 上游组件来源与许可说明

- **OPM (`opm-simulators` 2026.04)**：基于 GNU General Public License v3.0 or later（GPL-3.0-or-later）发布；本仓库 `overlay/` 中基于 OPM 的修改保留 GPL-3.0 声明，许可正文见 [overlay/COPYING](overlay/COPYING)。[OPM 项目仓库](https://github.com/OPM/opm-simulators)
- **Dune-ISTL (`2.11.0`)**：GPL-2.0 with runtime exception。[Dune 项目](https://dune-project.org/)
- **NVIDIA AMGX (`2.5.0`)**：BSD-3-Clause（用于压力子系统基准对照）。[AMGX](https://github.com/NVIDIA/AMGX)
- **Ginkgo (`1.11.0`)**：BSD-3-Clause（用于全耦合线性求解基准对照）。[Ginkgo](https://github.com/ginkgo-project/ginkgo)
- **LinglongAMG / XianzhiKrylov**：解得科技专有 GPU 加速求解库及其授权协议，不包含在本仓库中，需由用户通过正式交付渠道获取并按对应许可使用。
