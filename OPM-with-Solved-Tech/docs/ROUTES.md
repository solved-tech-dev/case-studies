# 两条求解路线与代码边界说明 (ROUTES)

> [!NOTE]
> 关于预条件器单独 `apply` 接口与双后端 `cpr` 组合路线的说明：此前版本中曾出现基于内部实验符号的单独 `apply` 接口，但该接口并未包含在已发布的 XianzhiKrylov `26.48` 发行包中。当前公开仓库已移除 `cpr` 路线，**OPM Flow + XianzhiKrylov 26.48 仅支持通过公开 API（`XianzhiKrylovAPI.h`）直接求解显式全耦合矩阵**；预条件器单独 `apply` 功能请等待 `26.49` 版本更新后再行接入。

## 1. 路线架构对比

| 维度 | `linglongamg` | `xianzhikrylov` |
| --- | --- | --- |
| 数学定位 | 混合 CPRW（仅替换压力 AMG） | 显式全耦合单级块 ILU(0)+BiCGStab |
| 外层迭代 | OPM/Dune CPU BiCGStab | XianzhiKrylov GPU BiCGStab (`XIANZHIKRYLOV_BSR`, `block_size=3`) |
| 细层平滑/预条件 | OPM/Dune CPU `paroverilu0` | XianzhiKrylov GPU 原生 BSR3 块 ILU(0) `level=0` |
| 粗层/压力转移 | OPM `PressureBhpTransferPolicy` (True-IMPES) | 无 |
| 压力内层 Krylov / AMG | LinglongAMG SDK 内置 GMRES(10) + GPU AMG | 无 |
| 井贡献处理 | OPM CPU `WellModelOperator` (`--matrix-add-well-contributions=false`) | OPM 显式组装井 Schur 补至主矩阵（`--matrix-add-well-contributions=true`） |
| 公开 SDK 接口 | `linglongamg.h` (`dlopen` 动态加载) | `XianzhiKrylovAPI.h` + `XianzhiKrylovData.h`（分阶段 C API v3） |
| 编译期宏 | `OPM_SOLVED_TECH_ENABLE_LINGLONGAMG` | `OPM_SOLVED_TECH_ENABLE_XIANZHIKRYLOV` |

---

## 2. 源码模块与调用链

### 2.1 `linglongamg` 路线

- **注册入口**：`opm/simulators/linalg/StandardPreconditioners_serial.hpp` 在 `#if defined(OPM_SOLVED_TECH_ENABLE_LINGLONGAMG) && HAVE_CUDA` 下注册 `"linglongamg"`。
- **适配器**：`overlay/linglongamg/opm/simulators/linalg/LinglongAMGPreconditioner.hpp` 继承 `Dune::PreconditionerWithUpdate<Vector, Vector>`，通过 `dlopen`/`dlsym` 动态加载 `liblinglongamg.so`，将 `1x1` 标量 FP64 CSR 矩阵上传至设备并调用 `linglongamg_setup` 与 `linglongamg_solve`。
- **计时记录**：`AMGCallTiming.hpp` 在设备同步前后记录纯 SDK `setup` 与 `solve` 耗时（排除数据传输与宿主 CPR/ILU）。

### 2.2 `xianzhikrylov` 路线（公开 API 显式全耦合）

- **调度入口**：`ISTLSolver.cpp` 中的 `FlexibleSolverInfo::create` 在 `OPM_SOLVED_TECH_ROUTE=xianzhikrylov`（或 `NORNE_XK_CONFIG`）时校验单 MPI 进程与 `!wellOperator_`（要求 `--matrix-add-well-contributions=true`），构造 `Dune::MatrixAdapter` 与 `NorneXianzhi<Matrix, Vector>`。
- **公开 C API 分阶段调用**（`overlay/xianzhikrylov/opm/simulators/linalg/NorneXianzhi.hpp`）：
  1. `init` / `rebuildStructure()`：提取 `BSR3` 的 `row_offsets`、`column_indices` 与必填 `diagonal_positions`（`uint64_t`），上传设备并通过 `xianzhikrylov_matrix_create`、`xianzhikrylov_vector_create` 注册设备缓冲区，执行 `XIANZHIKRYLOV_PHASE_ANALYSIS`。
  2. `update()`：如检测到结构变化自动重建描述符并重新 `ANALYSIS`，否则仅打包主机 `BSR3` 数值块上传至 `dvalues_`，调用 `XIANZHIKRYLOV_PHASE_NUMERIC_PREPARE`。
  3. `apply()`：打包并上传右端向量 `drhs_`，调用 `XIANZHIKRYLOV_PHASE_SOLVE` 原位求解设备向量 `dx_` 并下载；若运行时通过 `dlsym` 查得可选统计符号 `xk_internal_get_last_stats` 则读取迭代步与相对残差，否则在主机侧独立复算残差；编译期严格只依赖 `XianzhiKrylovAPI.h` 公开符号。

---

## 3. 路线隔离保证

1. `configure.sh` 与 `CMakeLists.txt` 对 `linglongamg` 不查找也不链接 `libXianzhiKrylov`；对 `xianzhikrylov` 不查找也不包含 `LinglongAMG`（也不需要 `USE_GPU_BRIDGE=ON`）。
2. `build.sh` 在链接结束后检查 `ldd bin/flow`，防止任何构建选项误引入无关后端。
3. 未设置 `OPM_SOLVED_TECH_ROUTE`（或使用 `CMake` 默认 `OPM_SOLVED_TECH_ROUTE=none`）时，OPM 原生 `Dune::FlexibleSolver` 路径、默认配置键与公共 ABI 保持完全不变。
