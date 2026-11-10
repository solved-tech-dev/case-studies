# Marmousi 3D 与 PanguLU

[English](README.md) · [简体中文](README.zh-CN.md)

本项目提供三档 Marmousi 三维弹性波频域正演案例。每档网格施加一个 5 Hz 垂向点力源，形成一个复数双精度右手边，在单张 NVIDIA GPU 上调用 PanguLU 26.46 求解。

独立的三维正演程序使用 TOY2DAC 2.6 随附的 Marmousi 纵波速度剖面构建三维模型。原始 TOY2DAC 文件保存在 [upstream/TOY2DAC-2.6](upstream/TOY2DAC-2.6)。三维正演程序与 PanguLU 适配代码分别位于 `apps/` 和 `adapter/`。

| 案例 ID | 网格间距 | 单元数（x × y × z） | 矩阵阶数 | 原始 COO 存储项 |
| --- | ---: | ---: | ---: | ---: |
| `marmousi3d-h85` | 85 m | 99 × 14 × 35 | 162,000 | 1,981,272 |
| `marmousi3d-h75` | 75 m | 112 × 16 × 40 | 236,283 | 2,963,553 |
| `marmousi3d-h65` | 65 m | 129 × 18 × 46 | 348,270 | 4,460,910 |

按[构建与运行](docs/BUILD-AND-RUN.zh-CN.md)配置 SDK、CUDA 运行库及 GPU 许可证，再从项目根目录运行。

```bash
python3 tools/run_case.py marmousi3d-h75 --output _work/runs/marmousi3d-h75
```

[案例定义与模型来源](provenance/README.zh-CN.md) · [接入入口](docs/INTEGRATION.zh-CN.md) · [许可与署名](LICENSE-NOTICE.md)
