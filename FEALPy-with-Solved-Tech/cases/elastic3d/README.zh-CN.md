# 三维线弹性

[English](README.md) · [构建与运行](../../docs/BUILD-AND-RUN.zh-CN.md) · [案例入口](run.py)

均匀各向同性的单位立方体内部体力为零，六个面施加随坐标线性变化的位移，模型同时产生拉伸和剪切。FEALPy 线弹性例 5 提供位移边界条件、三分量解析位移及材料参数 λ = μ = 4 × 10⁵。 [上游模型](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/linear_elasticity/exp0005.py)

单位立方体沿三个方向划分为 1,024 × 17 × 17 份，最终剖分为 1,775,616 个一次四面体。332,100 个节点各含三个位移自由度，对应 996,300 阶矩阵、33,228,054 个 CSR 存储项。网格沿一个方向加密。小网格采用 8 × 8 × 8 份，共 2,187 阶。

## 运行与结果

在源码包根目录执行：

```bash
python cases/elastic3d/run.py --small
python cases/elastic3d/run.py --output output/elastic3d.npy
```

模型通过 linear_system() 组装刚度矩阵和载荷，再调用 apply_bc() 施加位移边界。保存的一维 float64 数组按节点交错存储 x、y、z 位移分量，可重排为 (节点数, 3) 读取。与解析位移对照的全节点相对误差须 ≤ 0.005。 每个网格均求解一个右端项，先独立验证原始 CSR 的相对残差 ‖Ax − b‖₂ / ‖b‖₂ ≤ 10⁻⁸，再检验解场并保存。
