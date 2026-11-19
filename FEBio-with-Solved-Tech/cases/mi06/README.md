# mi06 · In-plane deformation of a locally constrained plate

[简体中文](README.zh-CN.md)

Prescribed in-plane displacement acts at selected nodes of a neo-Hookean plate. The original constraints control out-of-plane motion. Keeping the corner constraints at their original coordinates preserves the local loading mechanism as the grid is refined.

The original 10 × 10 × 1 grid is subdivided 45 ways along both in-plane axes. One element spans the thickness, giving 450 × 450 × 1 hexahedra. The refined model has 202,500 HEX8 elements and 1,017,000 equations in the assembled system.

The native XPLT includes displacement, relative volume, stress.

Generated input: `FEBio-with-Solved-Tech/cases/mi06/mi06-refined.feb`. [Original model](https://github.com/febiosoftware/TestSuite/blob/986005fc1865f5de7f2c789b8153d6f8fe2c0acd/Verify/mi06.feb) · [Build and run](../../docs/BUILD-AND-RUN.md)
