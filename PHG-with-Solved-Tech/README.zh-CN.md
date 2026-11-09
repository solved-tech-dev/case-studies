# PHG × PanguLU

[English](README.md) · [简体中文](README.zh-CN.md)

本项目提供三个 PHG 有限元案例在 Linux 上调用 PanguLU 26.44 的运行入口。PHG 和 PHG-Solid 继续负责网格、材料、边界、方程组装和结果输出。集成层在 PHG 原有的求解器注册表中增加 `pangulu`，由适配器求解组装后的实数双精度方程组，再把解送回原程序。

| 案例 ID | 上游入口 | 矩阵阶数 | 线性求解次数 |
| --- | --- | ---: | ---: |
| `solid-adaptive-stress` | PHG-Solid `examples/adaptive.xml` | 421,443 → 1,268,499 | 6 |
| `solid-static-p2` | PHG-Solid `examples/static.xml` | 833,091 | 1 |
| `interface-sphere-dg1` | PHG `examples/interface.c` | 786,432 | 1 |

案例参数和每次组装的矩阵规模记录在 `cases/<case>/case.json`。准备脚本从官方归档取得原始代码和网格，运行工具在输出目录生成本次配置。三个案例均沿用上游的物理模型和组装程序。

按[构建与运行](docs/BUILD-AND-RUN.zh-CN.md)配置好 Linux、CUDA 和 PanguLU 后，在项目根目录运行：

```bash
python3 tools/run_case.py solid-static-p2 --output _work/runs/static-p2 --vtk
```

运行目录包含应用日志、原生物理场和核验报告。PanguLU SDK 和目标 GPU 的许可证由用户另行取得。 [解得科技官网](https://solved-tech.com/)

[案例说明](docs/CASES.zh-CN.md) · [接口与输出](docs/INTEGRATION.zh-CN.md) · [许可与来源](LICENSE-NOTICE.md)
