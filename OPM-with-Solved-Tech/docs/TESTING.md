# 测试规范与验算指南 (TESTING)

本文定义本集成仓库的完整检查体系，明确无 SDK 宿主静态检查与有 GPU/SDK 运行期测试的边界。

---

## 1. 依赖预检

使用 `scripts/check-env.sh` 检查目标路线的工具链、公开 SDK 头文件、库导出符号及授权文件：

```bash
bash scripts/check-env.sh --route=linglongamg
bash scripts/check-env.sh --route=xianzhikrylov
```

- `linglongamg` 仅校验 `CUDA_HOME` 与 `LINGLONGAMG_ROOT`；即使未设置 `XIANZHIKRYLOV_ROOT` 也应通过。
- `xianzhikrylov` 仅校验 `CUDA_HOME` 与 `XIANZHIKRYLOV_ROOT`（要求 `XianzhiKrylovAPI.h`、`XianzhiKrylovData.h` 与 9 个公开 API 符号）；即使未设置 `LINGLONGAMG_ROOT` 也应通过。

---

## 2. 纯宿主静态检查（无 GPU / 无 SDK 环境）

在未安装任何闭源 SDK、未连接 GPU 的宿主机上，运行：

```bash
bash scripts/run-tests.sh --static-only
```

该命令依次完成：
1. **脚本语法检查**：对 `scripts/*.sh` 执行 `bash -n`。
2. **配置与清单解析**：运行 `tests/config/validate_configs.py` 与 `tests/norne-day1/check-day1-summary.py`，检查所有 JSON 语法、字段完整性以及无私有路径。
3. **入口帮助检查**：校验所有 `scripts/*.sh --help` 返回 0。
4. **非法路线与无 SDK 报错检查**：确认 `configure.sh --route=invalid` 与 `configure.sh --route=cpr` 明确报错退出，确认缺失环境变量时 `check-env.sh` 明确报错并指明缺失变量。
5. **公开目录扫描**：调用 `scripts/verify-public-tree.sh`。
6. **补丁与隔离覆盖检查**（当传入 `--opm-source="$OPM_SOURCE_ROOT"` 时）：在临时目录中对两条路线分别测试 `prepare-source.sh`，验证补丁零冲突应用以及各路线目录严格隔离。

---

## 3. LinglongAMG-only 组件预检

源码位于 [tests/smoke/check-linglong-cycle.cu](../tests/smoke/check-linglong-cycle.cu)（编译方法见 [tests/smoke/README.md](../tests/smoke/README.md)）：

- 动态加载 `liblinglongamg.so` 并验证授权初始化。
- 构造 4×4 稀疏测试矩阵，执行 `linglongamg_setup` 与单次 V-cycle（`linglongamg_apply_cycle` 或 `max_iters=1` 的 `linglongamg_solve` + `linglongamg_get_solve_info`）。
- 在主机端独立复算相对残差 $\|b - A x\|_2 / \|b\|_2$ 并确认残差下降（`< 1.0`）。
- 扫描 `/proc/self/maps`，确认进程未加载 `libXianzhiKrylov`。

---

## 4. XianzhiKrylov-only 组件预检（公开 API）

源码位于 [tests/smoke/check-xianzhikrylov-bsr3.cpp](../tests/smoke/check-xianzhikrylov-bsr3.cpp)：

- 纯 C++ 编译链接 `libXianzhiKrylov.so`（不依赖内部组件头文件与内部 `cprw` 符号）。
- 对 2 块非对称对角占优 `BSR3`（6×6 标量）系统测试分阶段接口：`PHASE_ANALYSIS` → 两次 `PHASE_NUMERIC_PREPARE` / `PHASE_SOLVE` 数值复用。
- 在主机端独立复算真实相对 L2 残差（要求 `< 0.005`）。
- 扫描 `/proc/self/maps`，确认进程未加载 `liblinglongamg`。

---

## 5. Norne 首日验算口径

从[发布的 Norne N64 输入包](../examples/norne/README.md#simulation-input)复制一份 `NORNE_N64.DATA`，并保持 `GRID_PROPS.INC` 与副本在同一目录。在副本的 `SCHEDULE` 中保留至 `7 'NOV' 1997` 的报告日期，在该 `DATES` 段结束后加入 `END`，将其作为首日测试的输入。用 `--profile --true-residual` 运行并进行以下核对。

- 核对退出码（0）、收敛状态、物理天数（`1.0` 天）、时间步数、Newton 迭代次数与线性迭代次数。
- 核对 `OPM_SOLVED_TECH_TRUE_RESIDUAL` 独立真实相对残差均满足设定容差（`< 0.005`），且与 SDK 报告残差一致。
- 若完整完成首日，核对 6 项场级关键物性汇总量：`FPR`、`FOPT`、`FWPT`、`FGPT`、`FWIT`、`FGIT`。若在显式井 Schur 补全耦合（`55,065,494` 块）上受限于作业墙钟上限超时（如退出码 124），须如实记录已完成的步数与逐系统收敛残差，不得外推。

---

## 6. 计时口径规范

| 计时项目 | 口径说明 |
| --- | --- |
| 端到端（无模拟输出） | 使用 `--enable-ecl-output=false`（`--no-output`）禁用 `.UNSMRY/.SMSPEC/.UNRST` 磁盘写入，记录进程墙钟 |
| 端到端（含输出） | 进程墙钟 `WALL_SECONDS`，包含读卡、网格初始化、装配、求解及摘要输出 |
| OPM 线性阶段总计 | OPM 报告的 `Linear solve time`，已包含 `Linear setup`、传输、封装与失败步耗时 |
| `OPM_XK_ANALYSIS` | 结构分析阶段（含 `diagonal_positions` 构建、设备缓冲区分配与 `PHASE_ANALYSIS`） |
| `OPM_XK_PREPARE` | 数值准备拆解：`pack_s`（主机 `BSR3` 打包）+ `values_h2d_s`（系数上传）+ `numeric_s`（`PHASE_NUMERIC_PREPARE`） |
| `OPM_XK_LINEAR` | 单次求解拆解：`rhs_h2d_s` + `sdk_solve_s`（`PHASE_SOLVE`）+ `x_d2h_s` |

---

## 7. 失败处理与分类记录

运行日志与校验报告须区分并明确记录以下失败类型：

1. **授权失败**：LinglongAMG 返回 `LINGLONGAMG_STATUS_LICENSE_ERROR (6)` 或 XianzhiKrylov 返回 `XIANZHIKRYLOV_LICENSE_ERROR (11)`，须记录授权路径定位（不得打印或保存授权密钥正文）。
2. **资源不足**：CUDA OOM（`cudaErrorMemoryAllocation`）或主机内存超限，记录网格/块规模与可用显存。
3. **迭代 breakdown / 未收敛**：达到 `max_iterations` 未收敛或产生 NaN/Inf，独立真实残差校验抛出异常中止。
4. **调度器超时**：如 Slurm 20 分钟分区墙钟限制（例如退出码 124），须如实记录已推进的物理天数与单步性能特征。
5. **缺失依赖**：`check-env.sh` 或 `configure.sh` 在预检期明确指出缺失变量、头文件或库，不静默回退到其他路线。

---

## 8. 无 SDK 环境下的测试完成边界

- **无 SDK 可完成的检查**：`bash -n` 语法校验、`validate_configs.py`、`check-day1-summary.py`、`--help` 入口校验、非法 `--route`（含已移除的 `cpr`）报错测试、无 SDK 环境变量时的 `check-env.sh` 报错测试、`verify-public-tree.sh` 公开树扫描，以及（当提供 `OPM_SOURCE_ROOT` 时）补丁应用与按路线隔离覆盖测试。
- **需 GPU 与真实 SDK 才可运行的检查**：CMake/NVCC 完整 OPM 编译、`ldd` 链接产物校验、`check-linglong-cycle.cu`、`check-xianzhikrylov-bsr3.cpp` 烟雾测试，以及 Norne N64 算例实机运行。

---

## 9. 提交与上传前必跑命令

每次修改文件后、提交或推送公开仓库前，必须执行：

```bash
bash scripts/verify-public-tree.sh
```

确认输出 `0 violations`，无任何 `.so`、`.a`、可执行文件、授权文件、私钥、私有绝对路径或多余大文件。
