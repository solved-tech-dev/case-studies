# tu01 · Displacement-controlled extension of a hyperelastic block

[简体中文](README.zh-CN.md)

One end is fixed, the opposite end undergoes prescribed axial displacement, and transverse motion is restrained at both ends. The Mooney–Rivlin material gives a nonlinear response that can be inspected through reaction forces, stress and relative volume.

Eight subdivisions in each direction turn the original 20 × 5 × 5 grid into 160 × 40 × 40 hexahedra. The end-face constraints are mapped onto the refined faces. The refined model has 256,000 HEX8 elements and 803,518 equations in the assembled system.

The native XPLT includes displacement, reaction forces, relative volume, stress.

Generated input: `FEBio-with-Solved-Tech/cases/tu01/tu01-8.feb`. [Original model](https://github.com/febiosoftware/TestSuite/blob/986005fc1865f5de7f2c789b8153d6f8fe2c0acd/Verify/tu01.feb) · [Build and run](../../docs/BUILD-AND-RUN.md)
