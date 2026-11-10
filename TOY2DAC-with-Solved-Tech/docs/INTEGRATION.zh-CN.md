# PanguLU 接入入口

[English](INTEGRATION.md) · [简体中文](INTEGRATION.zh-CN.md)

- [apps/marmousi3d/model.c](../apps/marmousi3d/model.c) 从上游速度剖面生成三维真值模型，并按当前速度组装复数 COO 矩阵。
- [apps/marmousi3d/fwi.cpp](../apps/marmousi3d/fwi.cpp) 生成合成观测，计算正演和共轭转置伴随波场，并用 L-BFGS 与线搜索更新纵波速度。正演和伴随共用 `case_factor()`、`case_solve_many()` 和 `case_cleanup()`。
- [adapter/pangulu_fwi_adapter.cpp](../adapter/pangulu_fwi_adapter.cpp) 将矩阵从零基 COO 转为按行索引排序的 CSC，保留全部存储项。它以 L0 重排和复数双精度通用 LU 配置 PanguLU，通过 `pangulu_execute()` 调用原生分析、数值分解和求解。适配器用 `pangulu_multiple_vector_create()` 描述一列 GPU 右手边和解向量，再将解传回反演程序。

每个系统独立创建并释放求解器状态，因此每次求解都执行完整的原生分析和分解。反演程序对返回的波场计算原始 COO 残差。[案例配置](../cases/fwi-grids.json)记录了网格、频率和矩阵规模。
