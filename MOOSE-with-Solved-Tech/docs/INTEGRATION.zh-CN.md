# PanguLU 接入入口与复用

## 完整应用调用链

`MOOSE → libMesh → PETSc KSP/PC LU → integration/petsc/solver_pangulu.cpp → adapter/pangulu_solver.hpp → PanguLU 26.46`

`app/src/base/MooseSolvedApp.C` 注册原有模块，不新增或替换物理内核。通过 `-dll_prepend` 动态加载 `libmoose_pangulu.so` 后，`PetscDLLibraryRegister_moose_pangulu()` 为顺序 AIJ 矩阵注册名称为 `pangulu` 的通用 LU 解法器。输入中选择 `preonly lu pangulu`，求解结果通过 PETSc 向量返回 MOOSE，后续迭代和输出仍由原应用完成。

## 主要文件

| 文件 | 职责 |
|---|---|
| `integration/petsc/solver_pangulu.cpp` | 注册解法器，接收 PETSc 矩阵、单个或成块 RHS，返回解和错误状态 |
| `adapter/sparse.hpp` | 保留完整 CSR 图及数值、建立 CSC 索引映射、独立计算原 CSR 残差 |
| `adapter/pangulu_solver.hpp` | 26.46 ABI 2 API、设备传输、分析/分解/求解、复用与原生计时 |
| `app/` | 注册原有 Solid Mechanics、Heat Transfer 及所需 Ray Tracing 模块 |
| `tools/run_case.py` | 网格准备、附加完整 Jacobian 配置、完整应用运行及结果检查 |

## 矩阵与复用

插件使用实数 FP64、通用 LU，不转为 Cholesky，不只提取上三角或下三角。CSR 转换为 CSC 时保留所有存储项，包括显式零；不增加微小扰动、不删零。每次数值回调都核对维数、完整行指针和列索引，并逐项比较全部系数。

- 新图：释放原求解上下文，重新分析并数值分解。
- 同图、系数改变：保留分析结果，写入全部当前系数，调用数值重分解。
- 同图、全部系数相同：复用数值分解，仅处理后续 RHS。

PanguLU 配置为 `PANGULU_REORDERING_METHOD=L1`、`PANGULU_REUSE_LEVEL=L2`，与此次六例最新测试使用的版本和复用路径一致。状态附着在源 PETSc 矩阵上，PETSc 重建因子对象时仍能访问该缓存；源矩阵销毁后释放。不同源矩阵分别维护状态，不把不同物理系统混在一起。

## RHS 与物理过程

`MatSolve` 使用单列多向量描述符；`MatMatSolve` 将 PETSc 提供的稠密独立 RHS 块按列打包，调用一次原生多右端项 SOLVE。插件不会预先制造下一轮 Newton RHS，也不会把依赖先前解的顺序迭代改成批量独立求解。

每个 RHS 求解后，都以原 CSR 和原 RHS 在 CPU 上独立计算 `||Ax-b||₂/||b||₂`，默认要求有限且不大于 `1e-8`；零 RHS 使用绝对残差。验证时间单独输出。非零 API 状态、非有限结果和超限残差通过 PETSc 错误机制返回，不静默改用其他求解器。

## 计时

`MOOSE_PANGULU_NATIVE` 行记录同步后的非重叠原生阶段：`analysis_order_symbolic`（包含 API 的分析/预处理工作）、`numeric`、`refactor_numeric`、`solve`。不将合并的分析阶段虚拆成无法独立测量的重排与符号时间。原生阶段之和是求解器核心时间。

格式转换、内存分配、传输、独立残差、MOOSE 组装、非线性控制和文件输出不属于该核心时间。`application_wall_seconds` 是完整应用子进程的墙钟时间；它包含上述应用内开销，不能与核心时间再次相加。没有依据历史数据设置运行时间门限。

## 支持范围

Linux x86_64、PanguLU 26.46 ABI 2、兼容 CUDA 12.9 的运行环境、实数双精度 PETSc。当前接入只支持单 MPI 进程的顺序 AIJ 矩阵和单 GPU。插件使用 PETSc 因子回调及内部矩阵结构，必须用与 MOOSE 相同的 PETSc 头文件、库和 MPI 构建；不同 PETSc 安装不能混用。当前开发核验基线为 PETSc 3.25.4。分布式矩阵、复数和转置求解未实现。
