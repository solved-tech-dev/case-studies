# FEALPy 接入 PanguLU

[English](README.md) · [构建与运行](docs/BUILD-AND-RUN.zh-CN.md) · [求解器接入说明](docs/INTEGRATION.zh-CN.md)

本包用 FEALPy 和 PanguLU 26.46 计算三种物理场：三维固体位移、二维频域电场和二维稳态标量场。每个案例从 FEALPy 的原模型生成网格、组装方程，将一个右端项交给 GPU 完成一次通用 LU 分解与求解，并核对原方程残差和有限元解场。

| 案例与运行入口说明 | 离散方式 | 矩阵阶数 | CSR 存储项数 |
| --- | --- | ---: | ---: |
| [三维线弹性](cases/elastic3d/README.zh-CN.md) | 一次四面体，每节点三个位移分量 | 996,300 | 33,228,054 |
| [二维频域电磁场](cases/maxwell2d/README.zh-CN.md) | 第一类 Nédélec 棱元，复数系数 | 999,824 | 11,585,824 |
| [二维对流–扩散–反应](cases/advdiff2d/README.zh-CN.md) | 二次 Lagrange 三角形 | 1,002,001 | 11,456,053 |

先按构建文档安装指定版本的 FEALPy，编译实数和复数适配库。运行时加 `--small` 可验证小网格，使用默认参数即可复现表中规模。无需修改 FEALPy 上游源码。[来源与许可](NOTICE.md)
