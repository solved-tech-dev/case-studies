# PanguLU 如何接入 FEBio

[English](INTEGRATION.md)

FEBio 读取模型中的 `Control/solver/linear_solver`。取值 `febio_gpu_lu` 对应 `FEBio-with-Solved-Tech/adapter/febio_gpu.cpp` 中通过 `PluginInitialize` 注册的 `FebioGpuLU` 类。运行时用 `-import` 加载 `FEBio-with-Solved-Tech/build/pangulu.so`。

| FEBio 入口 | 插件行为 |
| --- | --- |
| CreateSparseMatrix / SetSparseMatrix | 提供或接收零基完整 CSR，由 FEBio 写入当前切线矩阵。 |
| PreProcess / Factor | 检查行指针、列索引与上一次是否完全相同。新结构建立 CSC 映射并执行分析和首次 LU 分解。结构相同但系数变化时更新全部系数并重新数值分解；系数也逐字节相同时沿用已有因子。 |
| BackSolve | 每次使用 FEBio 当次提供的右端项，调用原生求解，将解传回 FEBio。以原始 CSR 和留存的右端项逐次验证残差。 |
| Destroy | 标记现有系统尚未就绪。FEBio 重建矩阵后由 Factor 再核对结构；求解器实例析构时释放后端资源。 |

`FEBio-with-Solved-Tech/adapter/gpu_backend.inc` 实现格式转换及 PanguLU 调用。CSC 保留显式零项和原始系数位置映射，采用双精度通用 LU、L0 重排和 L2 复用。原生阶段均通过 `pangulu_execute` 执行。FEBio 的 BackSolve 每次交付一个右端项，后端用一列向量描述符接收。每次只用当前已经组装好的矩阵和右端项。
