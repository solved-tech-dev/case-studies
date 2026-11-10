# Cases and model provenance

[English](README.md) · [简体中文](README.zh-CN.md)

The original TOY2DAC 2.6 files and Marmousi template are preserved in [upstream/TOY2DAC-2.6](../upstream/TOY2DAC-2.6), alongside their original [notice and references](../upstream/TOY2DAC-2.6/0LEGAL_STATEMENT). The package uses the template's `vp_Marmousi_exact` P-wave velocity array as the input to the standalone 3D forward application. [SEISCOPE TOY2DAC 2.6 release archive](https://seiscope2.osug.fr/IMG/tgz/toy2dac_v2.6_2019_05_24-2.tgz)

The 3D domain is approximately 8.4 × 1.2 × 3.0 km. The application samples the central part of the 2D profile, extends it obliquely at 45° along the crossline, and adds a smooth crossline-and-depth shift. It bilinearly interpolates the velocity onto the 3D grid. The S-wave speed is `0.67 × Vp` and density is 2,000 kg/m³. The cited study provides context for the domain dimensions. [Liu and Chen, 2026](https://doi.org/10.1007/s11430-025-1867-y)

The solver receives a second-order finite-difference discretization of the isotropic elastic equation at 5 Hz. A vertical point force is placed near the horizontal center, two grid steps below the top. The outermost nodes have zero displacement, with damping over the four inner boundary layers. Each grid node contributes three complex displacement components.

| Spacing | Domain extent x × y × z (m) | Nodes x × y × z |
| ---: | ---: | ---: |
| 85 m | 8,415 × 1,190 × 2,975 | 100 × 15 × 36 |
| 75 m | 8,400 × 1,200 × 3,000 | 113 × 17 × 41 |
| 65 m | 8,385 × 1,170 × 2,990 | 130 × 19 × 47 |

The small differences in domain extent follow from the integer number of elements in each direction. The grid definitions and original COO storage counts are recorded in [cases/grids.json](../cases/grids.json).
