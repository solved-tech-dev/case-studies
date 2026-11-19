# mi06 · 薄板在局部位移约束下的面内变形

[English](README.md)

在 neo-Hookean 薄板的指定节点施加面内位移，通过原约束控制厚度方向运动。加载只作用于原模型指定的位置，因此加密时必须保留角点约束的空间位置。输出可用于观察局部加载在板内产生的位移和应力分布。

10 × 10 × 1 网格在面内两个方向各细分 45 份，厚度仍为一层，得到 450 × 450 × 1 个单元。 生成模型含 202,500 个 HEX8 单元，实测线性系统为 1,017,000 阶。

物理场可在原生 XPLT 中查看，包括 displacement, relative volume, stress。

运行输入为 `FEBio-with-Solved-Tech/cases/mi06/mi06-refined.feb`。[上游原模型](https://github.com/febiosoftware/TestSuite/blob/986005fc1865f5de7f2c789b8153d6f8fe2c0acd/Verify/mi06.feb) · [生成与运行](../../docs/BUILD-AND-RUN.zh-CN.md)
