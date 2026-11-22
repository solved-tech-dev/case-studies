# 二维频域电磁场

[English](README.md) · [构建与运行](../../docs/BUILD-AND-RUN.zh-CN.md) · [案例入口](run.py)

在单位正方形上求解 curl curl E − k²E = f，波数 k = 1。边界采用含复系数的 Robin 阻抗条件。体源项和边界激励取自 FEALPy curlcurl 例 1，解析电场的两个分量分别为 x y (1 − x)(1 − y) 和 sin(πx) sin(πy)。 [上游模型](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/curlcurl/exp0001.py)

316 × 316 个方格剖分为 199,712 个三角形。标准第一类 Nédélec 棱元格式（FEALPy 参数 p = 1）产生 999,824 个自由度，CSR 存储项数为 11,585,824。小网格采用 6 × 6 个方格，共 384 阶。

## 运行与结果

在源码包根目录执行：

```bash
python cases/maxwell2d/run.py --small
python cases/maxwell2d/run.py --output output/maxwell2d.npy
```

矩阵和右端项均为 complex128。保存的复数向量是有限元自由度系数。可通过 model.space.function(dtype=bm.complex128) 建立同一空间的场对象，写回系数后用 value 在单元内求值。FEALPy 对重建电场与解析矢量场之差计算绝对 L² 误差，门限为 0.05。 每个网格均求解一个右端项，先独立验证原始 CSR 的相对残差 ‖Ax − b‖₂ / ‖b‖₂ ≤ 10⁻⁸，再检验解场并保存。
