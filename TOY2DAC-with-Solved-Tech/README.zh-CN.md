# Marmousi 3D 全波形反演与 PanguLU

[English](README.md) · [简体中文](README.zh-CN.md)

本项目提供三档 Marmousi 三维弹性频域全波形反演（FWI）案例。程序用 3 Hz 和 5 Hz 的单炮无噪合成观测，从平滑模型出发执行 12 次纵波速度更新。每个正演或伴随系统有一个复数双精度右手边，由单张 NVIDIA GPU 上的 PanguLU 26.46 求解。

独立的三维程序使用 TOY2DAC 2.6 随附的 Marmousi 纵波速度剖面。原始上游文件保存在 [upstream/TOY2DAC-2.6](upstream/TOY2DAC-2.6)。本项目的三维反演入口和 PanguLU 适配器分别位于 [apps/marmousi3d/fwi.cpp](apps/marmousi3d/fwi.cpp) 和 [adapter/pangulu_fwi_adapter.cpp](adapter/pangulu_fwi_adapter.cpp)。

| 案例 ID | 网格间距 | 单元数（x × y × z） | 矩阵阶数 | 原始 COO 存储项 |
| --- | ---: | ---: | ---: | ---: |
| `marmousi3d-h85` | 85 m | 99 × 14 × 35 | 162,000 | 1,981,272 |
| `marmousi3d-h75` | 75 m | 112 × 16 × 40 | 236,283 | 2,963,553 |
| `marmousi3d-h65` | 65 m | 129 × 18 × 46 | 348,270 | 4,460,910 |

按[构建与运行](docs/BUILD-AND-RUN.zh-CN.md)配置 SDK、CUDA 和 GPU 许可证，再从项目目录运行。

```bash
python3 tools/run_case.py marmousi3d-h75 --output _work/runs/marmousi3d-h75
```

[案例定义与模型来源](provenance/README.zh-CN.md) · [PanguLU 接入入口](docs/INTEGRATION.zh-CN.md) · [许可与署名](LICENSE-NOTICE.md)
