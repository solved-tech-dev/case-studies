# Oneel20fi

[English](README.md) | [简体中文](README.zh-CN.md)

采用C3D20R单元的单位立方体导热模型，导热系数为50，一面温度为0，对面向环境温度−12、换热系数10的介质散热。稳态解析解为T=−2x，热流为(100,0,0)。

## 模型参数

| 项目 | 数值 |
|---|---|
| 单元类型 | C3D20R |
| 网格维度 | 3维 |
| 细分参数 | 44 |
| 节点数 | 358,425 |
| 单元数 | 85,184 |
| 方程阶数n | 352,440 |
| 完整nnz | 19,849,804 |
| 每个完整应用的线性系统数 | 2 |

## 文件说明

- [加密输入](input/model.inp.gz)：gzip压缩的完整CalculiX模型，由运行工具解压。
- [参数文件](parameters.json)：网格、解法器配置与来源信息。
- [物理参数卡](physics-cards.inp)：供查看的材料、边界、载荷、步骤和输出卡片。
- [上游输入文件](source/)：CalculiX2.23验证案例，[来源](https://www.dhondt.de/ccx_2.23.test.tar.bz2)。

## 运行

按[构建说明](../../docs/BUILD-AND-RUN.zh-CN.md)完成编译与环境配置后，在`CalculiX-with-Solved-Tech/`运行：

```bash
mkdir -p _work/runs
python3 tools/run_case.py oneel20fi \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/oneel20fi
```

添加`--profile small`并选择另一个输出目录可运行上游网格。适配器采用L0排序和L2结构复用，图不变时更新系数并执行数值重分解。

物理含义见[案例介绍](../../docs/CASES.zh-CN.md)，数值检查见[计算结果检查](../../docs/VALIDATION.zh-CN.md)。
