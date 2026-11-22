# 二维稳态对流–扩散–反应

[English](README.md) · [构建与运行](../../docs/BUILD-AND-RUN.zh-CN.md) · [案例入口](run.py)

在单位正方形内，扩散张量为 diag(2, 3)，对流速度为 (1, −1)，反应系数为 4。四边施加齐次 Dirichlet 条件，源项按解析场 u = sin(πx) sin(πy) 构造，对流项使矩阵非对称。 [上游模型](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/diffusion_convection_reaction/exp0001.py)

500 × 500 个方格剖分为 500,000 个三角形。二次 Lagrange 单元在顶点和边中点设置自由度，共得到 1,002,001 阶矩阵和 11,456,053 个 CSR 存储项。小网格采用 8 × 8 个方格，共 289 阶。

## 运行与结果

在源码包根目录执行：

```bash
python cases/advdiff2d/run.py --small
python cases/advdiff2d/run.py --output output/advdiff2d.npy
```

FEALPy 采用积分参数 q = 5 组装扩散项、质量项、对流项及载荷，并施加 Dirichlet 条件。保存的 float64 向量按 P2 空间插值点顺序给出场值，与解析标量场对照的全插值点相对误差须 ≤ 0.005。 每个网格均求解一个右端项，先独立验证原始 CSR 的相对残差 ‖Ax − b‖₂ / ‖b‖₂ ≤ 10⁻⁸，再检验解场并保存。
