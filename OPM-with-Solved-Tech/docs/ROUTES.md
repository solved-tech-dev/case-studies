# 三条求解路线与代码边界说明 (ROUTES)

## 1. 路线架构对比

| 维度 | `linglongamg` | `xianzhikrylov` | `cpr` |
| --- | --- | --- | --- |
| 数学定位 | 混合 CPRW（仅替换压力 AMG） | 单级全耦合 ILU(0)+BiCGStab | 双后端全 GPU CPRW |
| 外层迭代 | OPM/Dune CPU BiCGStab | XianzhiKrylov GPU BiCGStab | XianzhiKrylov GPU BiCGStab |
| 细层平滑/预条件 | OPM/Dune CPU `paroverilu0` | XianzhiKrylov GPU 原序 BSR3 块 ILU(0) | XianzhiKrylov GPU 原序 BSR3 块 ILU(0) |
| 粗层/压力转移 | OPM `PressureBhpTransferPolicy` (True-IMPES) | 无 | OPM `PressureBhpTransferPolicy` (True-IMPES) + 井压力方程 |
| 压力内层 Krylov | LinglongAMG SDK 内置 GMRES(10) | 无 | XianzhiKrylov GPU GMRES(10) |
| 压力 AMG | LinglongAMG GPU AMG | 无 | LinglongAMG GPU 单次 V-cycle |
| 井贡献算子 | OPM CPU `WellModelOperator` | `NorneWellOperator` + GPU `WellContributionsCuda` | `NorneWellOperator` + GPU `WellContributionsCuda` |
| 编译期宏 | `OPM_SOLVED_TECH_ENABLE_LINGLONGAMG` | `OPM_SOLVED_TECH_ENABLE_XIANZHIKRYLOV` | `OPM_SOLVED_TECH_ENABLE_CPR` |

---

## 2. 源码模块与调用链

### 2.1 `linglongamg` 路线

- **注册入口**：`opm/simulators/linalg/StandardPreconditioners_serial.hpp` 在 `#if defined(OPM_SOLVED_TECH_ENABLE_LINGLONGAMG) && HAVE_CUDA` 下注册 `"linglongamg"`。
- **适配器**：`overlay/linglongamg/opm/simulators/linalg/LinglongAMGPreconditioner.hpp` 继承 `Dune::PreconditionerWithUpdate<Vector, Vector>`，通过 `dlopen`/`dlsym` 动态加载 `liblinglongamg.so`，将 `1x1` 标量 FP64 CSR 矩阵上传至设备并调用 `linglongamg_setup` 与 `linglongamg_solve`。
- **计时记录**：`AMGCallTiming.hpp` 在设备同步前后记录纯 SDK `setup` 与 `solve` 耗时（排除数据传输与宿主 CPR/ILU）。

### 2.2 `xianzhikrylov` 路线

- **调度入口**：`ISTLSolver.hpp` 在 `OPM_SOLVED_TECH_ROUTE=xianzhikrylov`（或 `NORNE_GPU_ILU_ONLY`）时构造 `detail::NorneWellOperator`；`ISTLSolver.cpp` 中的 `FlexibleSolverInfo::create` 构造 `NorneIluSolver<Matrix, Vector>`。
- **井算子修复**：`cuWellContributions.cu` 中的 `apply_well_contributions` 使用 `for (entry = idx_t; entry < dim * val_size; entry += blockDim.x)` 覆盖全部射孔并采用 `atomicAdd` 累加同一单元的多射孔贡献。
- **设备侧求解**：`NorneIluDevice.cu` 调用 `XianzhiKrylov_CprComponents.h` 的 `xk_internal_cprw_ilu_create/prepare/apply` 与 `xk_internal_cprw_krylov_create/solve/destroy`，外层 `fullOperator` 回调串联储层 `BSR3` SpMV 与独立井 Schur 补。

### 2.3 `cpr` 路线

- **调度入口**：`ISTLSolver.hpp` 与 `ISTLSolver.cpp` 在 `OPM_SOLVED_TECH_ROUTE=cpr`（或 `NORNE_GPU_CPRW`）时构造 `NorneCprwSolver<Matrix, Vector>`。
- **压力系统构建**：复用 OPM `PressureBhpTransferPolicy<Operator, SequentialInformation, double, false>` 在 True-IMPES 权重下构造含井压力方程的压力 CSR。
- **CPRW 两级组合**（`NorneCprwDevice.cu`）：
  1. `restrictPressure`：仅用 True-IMPES 权重限制储层残差（井压力右端初始化为零）。
  2. `inner`（XianzhiKrylov GPU GMRES(10)，`tol=1e-6, maxiter=100`）：调用 `pressureSpmv` 与 LinglongAMG 单次 V-cycle (`linglongamg_apply_cycle` 或单步 `linglongamg_solve` + `linglongamg_get_solve_info` 校验 `iterations == 1`)。
  3. `prolongPressure` + `fullOperator` + `combine`：将压力校正回填至储层压力分量并计算更新残差。
  4. `xk_internal_cprw_ilu_apply`：仅对储层 `BSR3` 执行一次后置块 ILU(0) 平滑并叠加压力校正。

---

## 3. 路线隔离保证

1. `configure.sh` 与 `CMakeLists.txt` 对 `linglongamg` 不查找也不链接 `libXianzhiKrylov`；对 `xianzhikrylov` 不查找也不包含 `LinglongAMG`；对 `cpr` 必须同时检出两套 SDK 头文件与库。
2. `build.sh` 在链接结束后检查 `ldd bin/flow`，防止任何构建选项误引入无关后端。
3. 未设置 `OPM_SOLVED_TECH_ROUTE`（或使用 `CMake` 默认 `OPM_SOLVED_TECH_ROUTE=none`）时，OPM 原生 `Dune::FlexibleSolver` 路径、默认配置键与公共 ABI 保持完全不变。
