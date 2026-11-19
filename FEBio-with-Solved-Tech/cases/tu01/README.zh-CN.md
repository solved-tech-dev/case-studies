# tu01 · 超弹性长方体的位移控制拉伸

[English](README.md)

一端固定，另一端沿轴向施加位移，两端横向位移受限。Mooney–Rivlin 模型描述试件在拉伸过程中的非线性响应，可由反力、应力和相对体积观察受力与变形。

原 20 × 5 × 5 网格在三个方向各细分 8 份，得到 160 × 40 × 40 个六面体；端面约束随网格映射。 生成模型含 256,000 个 HEX8 单元，实测线性系统为 803,518 阶。

物理场可在原生 XPLT 中查看，包括 displacement, reaction forces, relative volume, stress。

运行输入为 `FEBio-with-Solved-Tech/cases/tu01/tu01-8.feb`。[上游原模型](https://github.com/febiosoftware/TestSuite/blob/986005fc1865f5de7f2c789b8153d6f8fe2c0acd/Verify/tu01.feb) · [生成与运行](../../docs/BUILD-AND-RUN.zh-CN.md)
