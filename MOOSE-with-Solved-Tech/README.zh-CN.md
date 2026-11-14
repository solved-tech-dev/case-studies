# MOOSE with Solved Tech

[English](README.md)

本项目将 PanguLU 26.46 通过 PETSc 直接解法器接口接入 MOOSE，提供六个有限元应用案例。MOOSE 继续负责网格、材料、边界条件、方程组装、非线性迭代及结果输出；新增插件接收应用组装的完整矩阵和右端项，调用 PanguLU 后将解送回原应用流程。

| 案例目录 | 计算内容 | 百万阶配置的矩阵阶数 |
|---|---|---:|
| `thermomech` | 五个时间步中的温度、热膨胀和位移 | 1,002,501 |
| `annular-torsion` | 空心圆柱受扭后的位移与应力 | 1,005,312 |
| `mug` | 杯形区域内的稳态标量扩散 | 1,349,545 |
| `pinched-cylinder` | 点载荷作用下薄圆柱壳的位移和转角 | 1,017,005 |
| `gap-heat` | 两块固体通过间隙传导热量 | 982,802 |
| `hex20-heat` | 二次六面体单元的导热补片问题 | 988,281 |

每个案例都提供 `small` 小网格和 `million` 百万阶配置，便于先验证安装，再运行大网格。两种配置沿用相同的材料、载荷与时间控制。二次单元导热的加密配置采用结构化单位立方体网格替代原畸变种子网格，具体变更见案例说明。

## 使用入口

- [构建与运行](docs/BUILD-AND-RUN.zh-CN.md)：依赖、构建步骤、六例运行命令。
- [案例说明](docs/CASES.zh-CN.md)：各案例在计算什么，网格如何加密，哪些配置发生了变化。
- [接入入口说明](docs/INTEGRATION.zh-CN.md)：原生调用路径、矩阵转换、复用和计时口径。
- [结果验证](docs/VALIDATION.zh-CN.md)：原矩阵残差、完整物理场及参考结果对照。

主要接入代码为 `integration/petsc/solver_pangulu.cpp`，PanguLU API 封装为 `adapter/pangulu_solver.hpp`。`app/` 仅注册原有的 Solid Mechanics、Heat Transfer 模块；`tools/run_case.py` 运行完整 MOOSE 应用。

当前支持范围：**单 MPI 进程、单 NVIDIA GPU、实数双精度顺序 AIJ 矩阵、通用 LU**。若 PETSc 提供独立的多右端项块，插件通过原生多向量 API 一次求解；Newton 迭代中依赖上一次结果的右端项仍按原顺序产生和求解。

代码包不包含实测性能表、性能图、历史计算结果、SDK 或设备许可证。运行时可以输出本次耗时，验证只检查计算结果，不要求复现某个时间数值。上游输入原样保留，来源和版本见 `upstream.lock.json`；[许可与版权说明](LICENSE-NOTICE.md)。
