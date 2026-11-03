# CalculiX与PanguLU的接入原理

[English](INTEGRATION.md) | [简体中文](INTEGRATION.zh-CN.md) · [首页](../README.zh-CN.md)

## 求解入口

```text
CalculiX：组装A、b，处理约束与边界条件
  → spooles(...) 或 spooles_factor(...) / spooles_solve(...)
  → 适配器：转换存储格式，管理PanguLU对象，求解
  → 将x写入CalculiX解缓冲区
  → CalculiX：更新物理场，判断收敛，输出结果
```

[ccx_pangulu.cpp](../adapter/ccx_pangulu.cpp)提供CalculiX的SPOOLES求解分支调用的C ABI函数。构建时，该对象替代CalculiX的`spooles.c`包装层。有限元组装、材料计算、边界条件、Newton控制和物理场输出由CalculiX完成。

[spooles_entry.h](../adapter/spooles_entry.h)声明应用侧接口。CalculiX约束处理所需的SPOOLES辅助例程由外部SPOOLES库提供，与主线性方程求解入口分开。

## 矩阵与向量的数据流

1. 读取组装后的`ad/au/icol/irow`数组，将对称系统展开为完整CSR存储，包含显式零元。
2. 转换为CSC：64位列指针、32位行号、FP64系数，并保存CSR到CSC的索引映射以更新系数。
3. 创建PanguLU句柄、矩阵和向量，设置排序L0、复用级别L2。
4. 首次调用依次执行`ANALYSIS`、`NUMERIC_FACTORISATION`和`SOLVE`。
5. 后续调用比较完整图结构数组。图不变时刷新全部矩阵系数，执行`REFACTORISATION`和`SOLVE`；图改变时重新分析。
6. 将本次右端项传到GPU求解，再将解复制到CalculiX缓冲区。
7. 使用完整CSR矩阵和本次右端项计算相对残差。

残差条件为`||Ax−b||₂ / ||b||₂ ≤ 1e−8`；右端项为零时分母取1。适配器可使用当前分解结果额外进行至多两次修正求解，开销计入原生求解时间。

## 对象生命周期与复用

图结构不变的调用之间保留PanguLU句柄、矩阵结构和系数数组。每次分解调用均更新系数，并执行数值分解或数值重分解。

CalculiX的cleanup入口恢复应用线程设置；图改变或进程退出时释放PanguLU对象。适配器使用进程内状态，按每个进程运行一个模型的方式使用。

五个稳态模型均使用一个增量。Capacitor有一次线性系统调用，其余四例各两次。实测迭代轨迹见[性能与验证](RESULTS.zh-CN.md)。

## 计时输出

| 层级 | 内容 |
|---|---|
| `native` | GPU同步后的PanguLU分析、数值分解/重分解及求解调用 |
| `service` | 适配器工作，包含格式转换、对象设置、系数更新、传输及原生调用 |
| `evidence` | 矩阵/右端项保存及独立CSR残差计算 |
| `lifecycle` | PanguLU对象和GPU缓冲区释放 |

核心求解时间是互不重叠的`native`条目之和。分析接口合并了重排、符号分析和预处理。`service`包含`native`，两层分别报告。完整应用时间还包含CalculiX组装与物理场输出。

保存的数组描述每次调用的实际输入。跨解法器比较还需要在对应物理状态下比较图结构、系数和右端项；后续Newton右端项可能因浮点轨迹出现差异。

## 接口范围

适配器支持文档中的五个实数、对称、无移位稳态标量案例，使用单进程、单GPU和默认32位CalculiX整数ABI。其他分析类型及并发模型不在已验证接口范围内。

编译与运行命令见[构建与运行](BUILD-AND-RUN.zh-CN.md)。
