# PanguLU 如何接入 FEALPy

[English](INTEGRATION.md) · [构建与运行](BUILD-AND-RUN.zh-CN.md)

## 从有限元组装到物理场

1. 各案例在 NumPy 后端创建 FEALPy 网格和有限元空间，组装矩阵与载荷并处理边界条件。入口为 [elastic3d/run.py](../cases/elastic3d/run.py)、[maxwell2d/run.py](../cases/maxwell2d/run.py) 和 [advdiff2d/run.py](../cases/advdiff2d/run.py)。电磁案例在组装时计入复系数 Robin 边界项。
2. [adapter/driver.py](../adapter/driver.py) 将 FEALPy 组装结果转换为完整的 SciPy CSR，合并同一坐标上的重复项并对列号排序，保留显式零元。Python 按矩阵的数据类型选择实数库或复数库，再通过 `ctypes` 传入连续数组。
3. [adapter/gpu.cpp](../adapter/gpu.cpp) 通过 C 接口导出 `solve(n, row_ptr, columns, values, rhs, solution)`。下标从零开始，行指针采用 64 位整数，列号采用 32 位整数，右端向量和解向量长度均为 `n`。线弹性和输运模型采用双精度实数，电磁模型的矩阵、右端项和解均为双精度复数。
4. [adapter/gpu_backend.inc](../adapter/gpu_backend.inc) 保持算子方向，将同一矩阵由 CSR 转为 CSC，设置 L0 重排和通用 LU，通过 `pangulu_matrix_create()` 建立稀疏矩阵，通过 `pangulu_multiple_vector_create()` 建立一列、主维为 `n` 的设备向量。随后调用 `pangulu_execute()` 依次执行分析、数值分解和求解，同步 GPU 后将解传回主机。三例均为一个系统、一个右端项，故各执行一次分析和一次分解。调用结束或出现异常时，适配层释放本次创建的 PanguLU 对象和 CUDA 缓冲区。
5. Python 在原始 CSR 方程上核对解向量的相对残差。案例脚本将解写回 FEALPy 的场对象，随后对照解析解验证。电磁案例使用棱元基函数还原电场后计算积分误差。校验通过后，可用 `--output` 保存自由度系数。

三个入口均调用独立安装的 FEALPy 包。[案例对应的 FEALPy 上游源码](https://github.com/weihuayi/fealpy/tree/e033533bd8e24f796f9c5284685b20a98b4d11cb)
