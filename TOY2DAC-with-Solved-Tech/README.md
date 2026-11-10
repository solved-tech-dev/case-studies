# Marmousi 3D full-waveform inversion with PanguLU

[English](README.md) · [简体中文](README.zh-CN.md)

These cases run 3D elastic frequency-domain full-waveform inversion (FWI) on three Marmousi-derived grids. The application generates noise-free, single-shot observations at 3 and 5 Hz and makes 12 P-wave velocity updates from a smoothed starting model. Each forward or adjoint system has one double-precision complex right-hand side and is solved by PanguLU 26.46 on one NVIDIA GPU.

The standalone 3D application uses the Marmousi velocity profile supplied with TOY2DAC 2.6. The upstream files are preserved in [upstream/TOY2DAC-2.6](upstream/TOY2DAC-2.6). The inversion driver and PanguLU adapter are in [apps/marmousi3d/fwi.cpp](apps/marmousi3d/fwi.cpp) and [adapter/pangulu_fwi_adapter.cpp](adapter/pangulu_fwi_adapter.cpp).

| Case ID | Spacing | Elements (x × y × z) | Matrix order | Stored COO entries |
| --- | ---: | ---: | ---: | ---: |
| `marmousi3d-h85` | 85 m | 99 × 14 × 35 | 162,000 | 1,981,272 |
| `marmousi3d-h75` | 75 m | 112 × 16 × 40 | 236,283 | 2,963,553 |
| `marmousi3d-h65` | 65 m | 129 × 18 × 46 | 348,270 | 4,460,910 |

Follow [Build and run](docs/BUILD-AND-RUN.md) to configure the SDK, CUDA, and GPU license, then run from the project directory.

```bash
python3 tools/run_case.py marmousi3d-h75 --output _work/runs/marmousi3d-h75
```

[Cases and model provenance](provenance/README.md) · [PanguLU integration](docs/INTEGRATION.md) · [License and attribution](LICENSE-NOTICE.md)
