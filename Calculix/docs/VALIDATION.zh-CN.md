# 计算结果检查

[English](VALIDATION.md) | [简体中文](VALIDATION.zh-CN.md) · [首页](../README.zh-CN.md)

`tools/run_case.py`运行CalculiX后，会检查求解日志和原生物理场。检查使用当前输入模型及本次运行产生的输出。

应用可以输出运行时间作为诊断信息。时间数值不作为正确性条件，日志检查工具不评判或汇总运行时间。

## 线性方程组检查

`tools/check_log.py`检查应用是否正常完成、方程规模、求解调用次数、PanguLU库信息以及分析／重分解调用顺序。每个线性系统均使用完整原始CSR矩阵和对应右端项进行独立残差检查：

```text
||Ax - b||₂ / ||b||₂ <= 1e-8
```

右端项为零时分母取1。矩阵保存过程保留显式零元。图结构不变的后续调用会更新全部系数，并重新进行数值分解。

## 物理场检查

`tools/validate_fields.py`读取原生`model.frd`和`model.dat`输出，检查给定边界值、物理场有限性、平衡关系，以及各案例适用的解析或物理条件。这些检查与线性方程残差检查相互补充。

如需与独立提供的同网格结果进行完整场比较，可为运行工具添加`--reference-frd /path/to/reference.frd.gz`。原生打印场的精度与线性系统残差容限分别处理。

## 输出文件

在所选运行目录中查看`application-exit.json`、`log-validation.json`、`field-validation.json`和`validation.json`。排查失败时保留应用日志及原生物理场。命令见[构建与运行](BUILD-AND-RUN.zh-CN.md)，求解数据流见[接入原理](INTEGRATION.zh-CN.md)。

## 源码级检查

输入准备工具可用以下命令检查：

```bash
python3 -m unittest discover -s tests -v
```

`tests/changed_values.cpp`检查图结构不变、矩阵系数与右端项变化时的复用行为，构建与运行命令见[构建与运行](BUILD-AND-RUN.zh-CN.md)。
