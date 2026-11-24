# bp10 · 受压的双相四分之一圆筒

[English](README.md)

圆筒截面利用两个对称面，并约束厚度方向位移。内壁施加随时间变化的压力，内外壁分别设置流体压力条件。Holmes–Mow 固相与具有应变依赖性的渗透率共同决定变形、孔压和流体通量。

每个父单元沿三个方向各细分 12 份，按原节点插值，保留内壁面的方向及原压力节点集合。 生成模型含 276,480 个 HEX8 单元，实测线性系统为 904,176 阶。

物理场可在原生 XPLT 中查看，包括 displacement, stress, fluid pressure, fluid flux。

运行输入为 `FEBio-with-Solved-Tech/cases/bp10/bp10-12.feb`。[上游原模型](https://github.com/febiosoftware/TestSuite/blob/986005fc1865f5de7f2c789b8153d6f8fe2c0acd/Verify/bp10.feb) · [生成与运行](../../docs/BUILD-AND-RUN.zh-CN.md)
