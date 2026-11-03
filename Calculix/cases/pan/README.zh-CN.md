# Pan

[English](README.md) | [简体中文](README.zh-CN.md)

采用CAX6单元的铁锅轴对称模型，导热系数为50，底部通量为31500，空气侧/食物侧换热系数为10/150，环境温度为293.15。热平衡按完整旋转体的2πr权重计算。

## 模型参数

| 项目 | 数值 |
|---|---|
| 单元类型 | CAX6 |
| 网格维度 | 2维轴对称 |
| 细分参数 | 16 |
| 节点数 | 890,630 |
| 单元数 | 442,880 |
| 方程阶数n | 890,630 |
| 完整nnz | 10,206,762 |
| 每个完整应用的线性系统数 | 2 |

## 文件说明

- [加密输入](input/model.inp.gz)：gzip压缩的完整CalculiX模型，由运行工具解压。
- [参数文件](parameters.json)：网格、解法器配置与来源信息。
- [物理参数卡](physics-cards.inp)：供查看的材料、边界、载荷、步骤和输出卡片。
- [上游输入文件](source/)：CalculiX社区Launcher案例，[来源](https://github.com/calculix/examples/tree/a079253398c05f96b019e8c07ce076e370ddeb75/launcher/Therm/pan)。

## 运行

按[构建说明](../../docs/BUILD-AND-RUN.zh-CN.md)完成编译与环境配置后，在`Calculix/`运行：

```bash
mkdir -p _work/runs
python3 tools/run_case.py pan \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/pan
```

添加`--profile small`并选择另一个输出目录可运行上游网格。适配器采用L0排序和L2结构复用，图不变时更新系数并执行数值重分解。

物理含义见[案例介绍](../../docs/CASES.zh-CN.md)，数值检查见[计算结果检查](../../docs/VALIDATION.zh-CN.md)。
