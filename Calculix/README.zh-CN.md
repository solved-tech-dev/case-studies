# CalculiX + PanguLU

[English](README.md) | [简体中文](README.zh-CN.md)

使用PanguLU运行五个CalculiX有限元案例：**Dam、Ax6ht、Capacitor、Oneel20fi和Pan**。

本项目在CalculiX 2.23的线性方程求解入口接入PanguLU。CalculiX负责有限元组装、材料、边界条件、Newton迭代和结果输出；适配器将组装后的方程组交给PanguLU求解，再将解返回CalculiX。

PanguLU是解得科技的解法器产品。如需购买PanguLU并获取许可证，请访问[解得科技官网](https://solved-tech.com/)。PanguLU SDK和许可证需单独安装配置。

解得科技编写的适配器、项目工具、测试和文档采用[AGPL-3.0-only](LICENSE)。CalculiX及上游案例保留各自版权和许可，见[许可说明](LICENSE-NOTICE.zh-CN.md)。

## 文档导航

- [构建与运行](docs/BUILD-AND-RUN.zh-CN.md)：环境依赖、编译和案例运行命令。
- [接入原理](docs/INTEGRATION.zh-CN.md)：求解入口、数据转换与结构复用。
- [案例介绍](docs/CASES.zh-CN.md)：物理模型与网格参数。
- [计算结果检查](docs/VALIDATION.zh-CN.md)：残差检查与物理场验证。
- [代码结构与归属](docs/OWNERSHIP.zh-CN.md)：CalculiX原有部分和接入组件。
- [许可说明](LICENSE-NOTICE.zh-CN.md)。

## 目录结构

```text
adapter/                  CalculiX与PanguLU适配器
cases/<case>/source/      上游案例输入及include文件
cases/<case>/input/       gzip压缩的加密网格输入
cases/<case>/parameters.json  网格参数与方程组规模
provenance/upstream/      CalculiX 2.23源码包及许可证
patches/                 可选的Newton阶段计时补丁
tools/                   构建、运行、验证及网格生成工具
tests/                   输入准备与系数更新检查
docs/                    中英文使用文档
```

## 快速开始

从`case-studies`仓库根目录进入`Calculix/`，按[构建说明](docs/BUILD-AND-RUN.zh-CN.md)完成编译和运行环境配置后，运行：

```bash
mkdir -p _work/runs
python3 tools/run_case.py dam \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam
```

将`dam`替换为`ax6ht`、`capacitor`、`oneel20fi`或`pan`即可运行其他案例。计算结果写入指定输出目录。

## 支持配置

案例在Linux环境中采用单进程、单NVIDIA GPU，求解实数、对称、无移位的稳态标量系统。适配器使用PanguLU 26.42、ABI2、CUDA12.9，排序设置为L0，结构复用级别为L2。图结构不变时，更新矩阵系数并执行数值重分解。

五例包含平面扩散、轴对称导热及三维标量问题。各模型的几何维度和物理含义见案例介绍。
