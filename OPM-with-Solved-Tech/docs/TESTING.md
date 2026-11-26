# 测试规范与验算指南 (TESTING)

本文定义本集成仓库的完整检查体系，明确无 SDK 宿主静态检查与有 GPU/SDK 运行期测试的边界。

---

## 1. 依赖预检

使用 `scripts/check-env.sh` 检查目标路线的工具链、头文件、库导出符号及授权文件：

```bash
bash scripts/check-env.sh --route=linglongamg
bash scripts/check-env.sh --route=xianzhikrylov
bash scripts/check-env.sh --route=cpr
```

- `linglongamg` 仅校验 `CUDA_HOME` 与 `LINGLONGAMG_ROOT`；即使未设置 `XIANZHIKRYLOV_ROOT` 也应通过。
- `xianzhikrylov` 仅校验 `CUDA_HOME` 与 `XIANZHIKRYLOV_ROOT`；即使未设置 `LINGLONGAMG_ROOT` 也应通过。
- `cpr` 校验 `CUDA_HOME`、`LINGLONGAMG_ROOT` 与 `XIANZHIKRYLOV_ROOT`；任一缺失即非零退出，不允许静默降级。

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
4. **非法路线与无 SDK 报错检查**：确认 `configure.sh --route=invalid` 明确报错退出，确认缺失环境变量时 `check-env.sh` 明确报错并指明缺失变量。
5. **公开目录扫描**：调用 `scripts/verify-public-tree.sh`。
6. **补丁与隔离覆盖检查**（当传入 `--opm-source="$OPM_SOURCE_ROOT"` 时）：在临时目录中对三条路线分别测试 `prepare-source.sh`，验证补丁零冲突应用以及各路线目录严格隔离。

---

## 3. LinglongAMG-only 组件预检

源码位于 [tests/smoke/check-linglong-cycle.cu](../tests/smoke/check-linglong-cycle.cu)（编译方法见 [tests/smoke/README.md](../tests/smoke/README.md)）：

- 动态加载 `liblinglongamg.so` 并验证授权初始化。
- 构造 4×4 稀疏测试矩阵，执行 `linglongamg_setup` 与单次 V-cycle（`linglongamg_apply_cycle` 或 `max_iters=1` 的 `linglongamg_solve` + `linglongamg_get_solve_info`）。
- 在主机端独立复算相对残差 $\|b - A x\|_2 / \|b\|_2$ 并确认残差下降（`< 1.0`）。
- 扫描 `/proc/self/maps`，确认进程未加载 `libXianzhiKrylov`。

---

## 4. XianzhiKrylov-only 组件预检

源码位于 [tests/smoke/check-xianzhikrylov-ilu.cu](../tests/smoke/check-xianzhikrylov-ilu.cu)：

- 测试 `xk_internal_cprw_krylov_create/solve/destroy` 对 BiCGStab（`method=3`）与 GMRES（`method=7`）设备算子回调的求解与真实残差收敛。
- 测试 `xk_internal_cprw_ilu_create`、两次 `xk_internal_cprw_ilu_prepare` 数值刷新及多次 `xk_internal_cprw_ilu_apply`，核对精确块分解解向量精度（`< 1e-12`）。
- 扫描 `/proc/self/maps`，确认进程未加载 `liblinglongamg`。

---

## 5. CPR 双后端组件与子步骤预检

在 `cpr` 路线的 `run-case.sh --true-residual`（即设置 `NORNE_GPU_CPRW_AUDIT=1`）首次线性系统求解前，适配器 `NorneCprwSolver.hpp` 自动完成：

1. **含井算子校验** (`NORNE_CPRW_OPERATOR_AUDIT`)：对比 GPU `fineSpmv + WellContributionsCuda` 与 OPM CPU `WellModelMatrixAdapter::apply`，要求相对误差 `< 1e-10`。
2. **压力 CSR 与 True-IMPES 限制校验** (`NORNE_CPRW_CORRECTION_AUDIT`)：对比 GPU 限制向量与 OPM `PressureBhpTransferPolicy::moveToCoarseLevel`，要求相对误差 `< 1e-12`；检验压力系统（含 36 个井压力方程）GMRES+LinglongAMG V-cycle 残差 `<= 1.01e-6`。
3. **完整 CPRW 修正校验**：用独立 CPU `ParallelOverlappingILU0` 复算后平滑及压力注入叠加，对比 GPU 完整 CPRW 输出，要求相对误差 `< 1e-9`。
4. **外层 BiCGStab 与真实残差**：每次求解结束独立计算原始系统 $\|b - A x\|_2 / \|b\|_2$ 并核对与 SDK 报告残差一致（差值 `<= 1e-9` 且收敛时 `<= tol + 1e-10`）。

---

## 6. Norne 首日验算口径

在获取 Norne N64 首日 Deck（`NORNE1.DATA` + `GRID_PROPS.INC`）后，使用 `--profile --true-residual` 运行首日：

- 核对退出码（0）、收敛状态、物理天数（`1.0` 天）、时间步数、Newton 迭代次数与线性迭代次数。
- 核对独立真实相对残差均满足设定容差（外层 BiCGStab `< 0.005`，压力 GMRES `< 1e-6`）。
- 核对 6 项场级关键物性汇总量：`FPR`、`FOPT`、`FWPT`、`FGPT`、`FWIT`、`FGIT` 与同算例参考（[reference-day1-summary.json](../examples/norne/reference-day1-summary.json)）在摘要输出精度下一致。

---

## 7. 计时口径规范

| 计时项目 | 口径说明 |
| --- | --- |
| 端到端（无模拟输出） | 使用 `--enable-ecl-output=false`（`--no-output`）禁用 `.UNSMRY/.SMSPEC/.UNRST` 磁盘写入，记录进程墙钟或 `Setup time - Deck input + Simulation time - Output write time` |
| 端到端（含输出） | 进程墙钟 `WALL_SECONDS`，包含读卡、网格初始化、装配、求解及摘要输出 |
| OPM 线性阶段总计 | OPM 报告的 `Linear solve time`，已包含 `Linear setup`、传输、封装与失败步耗时 |
| 线性准备 (`Linear setup`) | 权重计算、压力 CSR 构建、矩阵打包（`pack_s`）、井数据准备（`well_prepare_s`）、H2D 上传（`values_h2d_s`）、AMG setup 与 ILU 分解 |
| 纯外层求解 | 外层 Krylov 迭代调用（包含其内部调用的预条件器/压力 GMRES/AMG V-cycle/细层 ILU apply） |
| 嵌套子阶段（仅诊断） | `amg_cycle_s` 包含在 `pressure_gmres_s` 中，`pressure_gmres_s`、`ilu_apply_s`、`fine_spmv_s`、`well_s` 包含在外层求解中，**禁止把嵌套子项直接相加** |

---

## 8. 失败处理与分类记录

运行日志与校验报告须区分并明确记录以下失败类型：

1. **授权失败**：LinglongAMG 返回 `LINGLONGAMG_STATUS_LICENSE_ERROR (6)` 或 XianzhiKrylov 返回授权相关状态，须记录设备识别码格式与授权路径权限（不得打印或保存授权密钥正文）。
2. **资源不足**：CUDA OOM（`cudaErrorMemoryAllocation`）或主机内存超限，记录网格/块规模与可用显存。
3. **迭代 breakdown / 未收敛**：达到 `maxiter` 未收敛或产生 NaN/Inf，独立真实残差校验抛出异常中止。
4. **调度器超时**：如 Slurm 20 分钟分区墙钟限制（例如退出码 124），须如实记录已推进的物理天数，不得外推为完整 176 天结果。
5. **缺失依赖**：`check-env.sh` 或 `configure.sh` 在预检期明确指出缺失变量、头文件或库，不静默回退到其他路线。

---

## 9. 无 SDK 环境下的测试完成边界

- **无 SDK 可完成的检查**：`bash -n` 语法校验、`validate_configs.py`、`check-day1-summary.py`、`--help` 入口校验、非法 `--route` 报错测试、无 SDK 环境变量时的 `check-env.sh` 报错测试、`verify-public-tree.sh` 公开树扫描，以及（当提供 `OPM_SOURCE_ROOT` 时）补丁应用与按路线隔离覆盖测试。
- **必须标记“未运行 (Not Run)”的检查**：CMake/NVCC 完整 OPM 编译、`ldd` 链接产物校验、`check-linglong-cycle.cu`、`check-xianzhikrylov-ilu.cu` GPU 烟雾测试，以及 Norne N64 首日/176 天实机运行（需 CUDA GPU、闭源 SDK 与授权）。

---

## 10. 提交与上传前必跑命令

每次修改文件后、提交或推送公开仓库前，必须执行：

```bash
bash scripts/verify-public-tree.sh
```

确认输出 `0 violations`，无任何 `.so`、`.a`、可执行文件、授权文件、私钥、私有绝对路径或大体积算例输出。
