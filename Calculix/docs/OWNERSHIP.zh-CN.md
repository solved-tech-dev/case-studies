# 代码结构与归属

[English](OWNERSHIP.md) | [简体中文](OWNERSHIP.zh-CN.md) · [首页](../README.zh-CN.md)

## 组件说明

| 组件 | 实现与位置 |
|---|---|
| 有限元单元、材料、组装、约束、Newton控制、场输出 | CalculiX2.23，源码包位于`provenance/upstream/` |
| 上游案例输入 | 四个CalculiX验证案例和社区Pan案例，位于`cases/<case>/source/` |
| 网格加密与边界映射 | 项目网格准备工具`tools/generate_mesh.py` |
| PanguLU求解入口实现 | `adapter/ccx_pangulu.cpp` |
| CalculiX入口声明 | `adapter/spooles_entry.h` |
| 重排、分析、数值分解和求解算法 | 外部PanguLU库 |
| 构建、运行、残差及物理场检查 | `tools/`与`tests/` |
| Newton阶段计时 | 可选补丁`patches/nonlingeo-timing.patch` |

## 构建配置

构建工具编译CalculiX2.23与PanguLU入口实现，采用32位ITG和LP64 BLAS/LAPACK，并启用CalculiX的SPOOLES调用分支、ARPACK、MATRIXSTORAGE及NETWORKOUT。

应用侧的`spooles.h`声明由适配器头文件提供，`spooles.c`求解包装层由适配对象替代。约束处理辅助例程使用外部SPOOLES库，`dgesv`由外部LAPACK提供。可执行文件由主程序对象和CalculiX源码静态库链接生成。

可选计时补丁在Newton阶段周围增加时间测量。

## 模型准备

`source/`保存上游输入模型，`input/model.inp.gz`保存运行工具使用的加密模型。细分过程划分父单元，将集合和边界面映射到子单元，并启用完整物理场输出。

各例输入及`physics-cards.inp`列出材料参数、物理载荷、边界值和步骤控制。Ax6ht的`STEP,NLGEOM`以及节点1、46上的温度约束属于模型设置。

## 署名与许可

CalculiX及上游案例保留各自的作者与许可声明。适配器和模型准备工具是解得科技提供的项目组件。PanguLU单独安装并授权。

参见[许可说明](../LICENSE-NOTICE.zh-CN.md)及各例`parameters.json`中的来源链接。

解得科技新增部分采用AGPL-3.0-only，见[许可说明](../LICENSE-NOTICE.zh-CN.md)。
