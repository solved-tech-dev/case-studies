# ngspice 方程采集与 PanguLU 接入

| 文件 | 作用 |
| --- | --- |
| [`integration/ngspice/capture.patch`](../integration/ngspice/capture.patch) | 在 ngspice 的 `src/maths/ni/niiter.c` 实数求解位置加入可选的采集调用，并把采集源码加入启用 KLU 的构建。 |
| [`integration/ngspice/capture.c`](../integration/ngspice/capture.c) | `SolvedTechCaptureSolve()` 读取已组装的 KLU CSC 和右端项，将完整原系统导出为 CSR，再调用原生 `SMPsolve()`，保存完整参考解。 |
| [`adapter/replay.cpp`](../adapter/replay.cpp) | 将完整 CSR 转为 CSC，设置 L0 重排与 L1 因子复用，通过 `pangulu_execute()` 执行分析、数值 LU 分解和三次求解。 |
| [`tools/run_case.py`](../tools/run_case.py) | 核对网表、原生工作点输出、矩阵结构与系数、阶段顺序、原方程残差及参考解误差。 |

网表读取、元件方程、矩阵组装和工作点迭代都由 ngspice 处理。设置 `NGSPICE_CAPTURE` 后，采集代码在原生 KLU 求解前后记录方程和解。独立回放程序加载一次矩阵，再按顺序读入每次调用的原始右端项，逐次执行求解并验证结果。三次线性求解属于同一个直流工作点。

采集保留全部实数双精度系数及显式零元。Linux x86-64 下文件使用小端字节序。行指针 `.rp.i64` 与列索引 `.ci.i32` 都从零开始，系数、右端项、KLU 全解分别为 `.a.f64`、`.b.f64` 和 `.x.f64`。当矩阵结构或系数逐字节一致时，采集日志引用已有文件。每次调用仍单独保留右端项和解。

本包回放入口针对清单中的固定矩阵序列。程序校验完整原始 CSR 相对残差不超过 1e-8、与 KLU 完整解向量的相对 L2 差不超过 1e-6，并核对一次分析、一次分解和三次求解。
