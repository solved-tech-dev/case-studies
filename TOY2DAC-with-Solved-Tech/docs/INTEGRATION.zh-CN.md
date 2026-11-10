# PanguLU 接入入口

[English](INTEGRATION.md) · [简体中文](INTEGRATION.zh-CN.md)

- [apps/marmousi3d/frequency_elastic.c](../apps/marmousi3d/frequency_elastic.c) 读取上游速度剖面，生成三维模型，组装复数 COO 矩阵和点力源右手边，输出位移波场。
- [adapter/pangulu_adapter.h](../adapter/pangulu_adapter.h) 定义 C 接口。应用用 `case_factor` 传入原始矩阵，经 `case_solve` 传入一个复数右手边，再由 `case_cleanup` 释放求解器状态。
- [adapter/pangulu_adapter.cpp](../adapter/pangulu_adapter.cpp) 将零基 COO 转为按列排序的 CSC，保留原始存储项。适配器以 L0 重排、L1 复用配置复数双精度通用 LU，调用原生分析和数值分解。随后用 `pangulu_multiple_vector_create` 描述单列向量并求解。解的实部和虚部传回应用，用于原始 COO 残差计算和波场输出。

GPU 向量和 PanguLU 句柄由适配器管理。模型、原始矩阵、震源和输出由正演应用管理。[案例配置](../cases/grids.json)定义每档的网格间距、单元数、频率和预期矩阵尺寸。
