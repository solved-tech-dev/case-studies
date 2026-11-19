# FEBio-with-Solved-Tech

[English](README.md) · [构建、运行与结果校验](docs/BUILD-AND-RUN.zh-CN.md) · [求解器接入说明](docs/INTEGRATION.zh-CN.md)

本项目将 PanguLU 26.46 作为线性求解器接入 FEBio。FEBio 负责材料响应、方程组装、非线性迭代和物理场输出；插件在每次线性求解时接收矩阵和右端项，将解向量交还 FEBio。构建过程沿用原版 FEBio 的插件接口。

| 案例 | 模型 | HEX8 单元数 | 实测方程数 |
| --- | --- | ---: | ---: |
| [bp10](cases/bp10/README.zh-CN.md) | 双相四分之一圆筒受压 | 276,480 | 904,176 |
| [tu01](cases/tu01/README.zh-CN.md) | 超弹性长方体单轴拉伸 | 256,000 | 803,518 |
| [mi06](cases/mi06/README.zh-CN.md) | neo-Hookean 薄板面内加载 | 202,500 | 1,017,000 |

从构建指南开始，可从固定版本的官方模型生成推送所用网格，核对输入文件，编译插件并运行案例。计算完成后，可在 FEBio Studio 中打开原生 XPLT 文件查看物理场。[依赖与许可说明](NOTICE.md)
