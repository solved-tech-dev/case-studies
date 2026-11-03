# Dam

[English](README.md) | [简体中文](README.zh-CN.md)

采用CPS8单元的二维正交各向异性扩散模型，对应地下水/标量扩散类比，上下游指定标量值为8.5344和3.9624。

## 模型参数

| 项目 | 数值 |
|---|---|
| 单元类型 | CPS8 |
| 网格维度 | 2维 |
| 细分参数 | 32 |
| 节点数 | 985,756 |
| 单元数 | 327,680 |
| 方程阶数n | 984,474 |
| 完整nnz | 15,390,186 |
| 每个完整应用的线性系统数 | 2 |

## 文件说明

- [加密输入](input/model.inp.gz)：gzip压缩的完整CalculiX模型，由运行工具解压。
- [参数文件](parameters.json)：网格、解法器配置与来源信息。
- [物理参数卡](physics-cards.inp)：供查看的材料、边界、载荷、步骤和输出卡片。
- [上游输入文件](source/)：CalculiX2.23验证案例，[来源](https://www.dhondt.de/ccx_2.23.test.tar.bz2)。

## 运行

按[构建说明](../../docs/BUILD-AND-RUN.zh-CN.md)完成编译与环境配置后，在`Calculix/`运行：

```bash
mkdir -p _work/runs
python3 tools/run_case.py dam \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam
```

添加`--profile small`并选择另一个输出目录可运行上游网格。适配器采用L0排序和L2结构复用，图不变时更新系数并执行数值重分解。

物理含义见[案例介绍](../../docs/CASES.zh-CN.md)，数值检查见[计算结果检查](../../docs/VALIDATION.zh-CN.md)。
