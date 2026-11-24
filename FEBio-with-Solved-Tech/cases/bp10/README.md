# bp10 · Pressurized biphasic quarter-cylinder

[简体中文](README.zh-CN.md)

The quarter-cylinder uses two symmetry planes and suppressed out-of-plane displacement. A time-dependent pressure acts on the inner wall. Fluid pressure conditions on the inner and outer walls drive flow through the Holmes–Mow solid skeleton with strain-dependent permeability.

Each parent element is split 12 ways along all three axes. Interpolation preserves the original nodal geometry, pressure sets and orientation of loaded faces. The refined model has 276,480 HEX8 elements and 904,176 equations in the assembled system.

The native XPLT includes displacement, stress, fluid pressure, fluid flux.

Generated input: `FEBio-with-Solved-Tech/cases/bp10/bp10-12.feb`. [Original model](https://github.com/febiosoftware/TestSuite/blob/986005fc1865f5de7f2c789b8153d6f8fe2c0acd/Verify/bp10.feb) · [Build and run](../../docs/BUILD-AND-RUN.md)
