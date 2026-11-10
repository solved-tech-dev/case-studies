# Marmousi 3D with PanguLU

[English](README.md) · [简体中文](README.zh-CN.md)

Three frequency-domain elastic forward-modeling cases use PanguLU 26.46 to solve the complex sparse system on one NVIDIA GPU. Each case uses one 5 Hz vertical point source and one right-hand side in complex double precision.

The standalone 3D application builds its velocity field from the Marmousi P-wave profile shipped with TOY2DAC 2.6. The original TOY2DAC files are kept in [upstream/TOY2DAC-2.6](upstream/TOY2DAC-2.6). The 3D forward application and PanguLU adapter live separately in `apps/` and `adapter/`.

| Case ID | Spacing | Elements (x × y × z) | Matrix order | Stored COO entries |
| --- | ---: | ---: | ---: | ---: |
| `marmousi3d-h85` | 85 m | 99 × 14 × 35 | 162,000 | 1,981,272 |
| `marmousi3d-h75` | 75 m | 112 × 16 × 40 | 236,283 | 2,963,553 |
| `marmousi3d-h65` | 65 m | 129 × 18 × 46 | 348,270 | 4,460,910 |

Follow [Build and run](docs/BUILD-AND-RUN.md) to configure the SDK, CUDA runtime and GPU license, then run from the project root.

```bash
python3 tools/run_case.py marmousi3d-h75 --output _work/runs/marmousi3d-h75
```

[Case definitions and model provenance](provenance/README.md) · [Integration entry points](docs/INTEGRATION.md) · [License and attribution](LICENSE-NOTICE.md)
