# 测试规范与验证指南 (TESTING)

本文说明本集成仓库的静态检查与 GPU 运行期验证方法。

---

## 1. 依赖预检

使用 `scripts/check-env.sh` 检查目标路线的工具链、公开 SDK 头文件、导出符号及授权文件：

```bash
bash scripts/check-env.sh --route=linglongamg
bash scripts/check-env.sh --route=xianzhikrylov
```

- `linglongamg` 仅校验 `CUDA_HOME` 与 `LINGLONGAMG_ROOT`。
- `xianzhikrylov` 仅校验 `CUDA_HOME` 与 `XIANZHIKRYLOV_ROOT`（要求 `XianzhiKrylovAPI.h`、`XianzhiKrylovData.h` 与公开 C API 符号）。

---

## 2. 宿主静态检查（无 GPU / 无 SDK 环境）

在未安装闭源 SDK 或未连接 GPU 的宿主机上，可执行：

```bash
bash scripts/run-tests.sh --static-only
```

该命令依次完成：
1. **脚本语法检查**：对 `scripts/*.sh` 执行 `bash -n`。
2. **配置与清单解析**：运行 `tests/config/validate_configs.py` 与 `tests/norne-day1/check-day1-summary.py`，校验 JSON 语法、字段完整性及无私有路径。
3. **入口帮助检查**：校验所有 `scripts/*.sh --help` 正常返回。
4. **参数校验与无 SDK 预检检查**：确认非法 `--route` 参数与缺失 SDK 环境变量时明确报错提示。
5. **公开目录扫描**：调用 `scripts/verify-public-tree.sh`。
6. **补丁与隔离覆盖检查**（当传入 `--opm-source="$OPM_SOURCE_ROOT"` 时）：在临时目录中分别验证两条路线的补丁干净应用与目录严格隔离。

---

## 3. LinglongAMG 组件烟雾测试

源码位于 [tests/smoke/check-linglong-cycle.cu](../tests/smoke/check-linglong-cycle.cu)（编译方法见 [tests/smoke/README.md](../tests/smoke/README.md)）：

- 动态加载 `liblinglongamg.so` 并验证授权初始化。
- 构造 4×4 稀疏测试矩阵，执行 `linglongamg_setup` 与单次 V-cycle。
- 在主机端独立复算相对残差并确认残差下降（`< 1.0`）。
- 扫描 `/proc/self/maps`，确认进程未加载 `libXianzhiKrylov`。

---

## 4. XianzhiKrylov 组件烟雾测试

源码位于 [tests/smoke/check-xianzhikrylov-bsr3.cpp](../tests/smoke/check-xianzhikrylov-bsr3.cpp)：

- 纯 C++ 编译链接 `libXianzhiKrylov.so`（基于 `XianzhiKrylovAPI.h` 与 `XianzhiKrylovData.h`）。
- 对 2 块非对称对角占优 `BSR3`（6×6 标量）系统测试分阶段接口：`PHASE_ANALYSIS` → 两次 `PHASE_NUMERIC_PREPARE` / `PHASE_SOLVE` 数值复用。
- 在主机端独立复算真实相对 L2 残差（要求 `< 0.005`）。
- 扫描 `/proc/self/maps`，确认进程未加载 `liblinglongamg`。

---

## 5. Norne 算例验算与计时口径

从[发布的 Norne N64 输入包](../examples/norne/README.md#simulation-input)获取 `NORNE_N64.DATA` 与 `GRID_PROPS.INC`，并保持两者位于同一目录；短周期验证可在副本的 `SCHEDULE` 中保留目标日期段并添加 `END`。使用 `--profile --true-residual` 运行并核对以下口径：

| 计时项目 | 口径说明 |
| --- | --- |
| `OPM_AMG_TIMING` | `linglongamg` 路线纯压力子系统 GPU API 耗时（`setup` + `solve`，已同步 GPU，不含宿主传输与 CPR/ILU） |
| `OPM_XK_ANALYSIS` | `xianzhikrylov` 结构分析阶段（含 `diagonal_positions` 构建、设备缓冲区分配与 `PHASE_ANALYSIS`） |
| `OPM_XK_PREPARE` | `xianzhikrylov` 数值准备拆解：`pack_s`（主机 `BSR3` 打包）+ `values_h2d_s`（系数上传）+ `numeric_s`（`PHASE_NUMERIC_PREPARE`） |
| `OPM_XK_LINEAR` | `xianzhikrylov` 单次求解拆解：`rhs_h2d_s` + `sdk_solve_s`（`PHASE_SOLVE`）+ `x_d2h_s` |
| `OPM_SOLVED_TECH_TRUE_RESIDUAL` | 主机端基于原始系数与右端向量独立复算的真实相对 L2 残差 |

---

## 6. 运行期常见错误排查

1. **授权失败**：LinglongAMG 返回 `LINGLONGAMG_STATUS_LICENSE_ERROR (6)` 或 XianzhiKrylov 返回 `XIANZHIKRYLOV_LICENSE_ERROR (11)`，请检查授权文件路径与设备绑定权限。
2. **显存不足**：遇到 `cudaErrorMemoryAllocation` 时请检查矩阵块规模与当前卡可用显存。
3. **未收敛**：检查 `max_iterations` 与配置容差，独立真残差校验会在未达设定容差时中止并报错。
4. **缺失依赖**：`check-env.sh` 或 `configure.sh` 会在预检阶段明确提示缺失的环境变量、头文件或动态库路径。

---

## 7. 公开树合规校验

在修改仓库文件后，可运行以下命令确认目录中不包含二进制产物、授权文件或私有路径：

```bash
bash scripts/verify-public-tree.sh
```
