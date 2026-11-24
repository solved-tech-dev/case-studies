# NGSPICE-with-Solved-Tech

[简体中文](README.zh-CN.md)

Five generated power-distribution circuits for real double-precision DC operating-point analysis in ngspice. Each produces a sparse system with roughly one million unknowns. Use the supplied generator to obtain the exact netlists from the benchmark.

ngspice assembles the circuit equations and completes its native KLU operating-point run. The capture hook records the original equation sequence. The separate PanguLU 26.46 replay executable runs analysis, factorization and all three original solves. Since the matrix structure and values match across calls, it retains the same numeric factors for the second and third solves and checks every solution against both the original equation and KLU.

| Circuit | Grid | Matrix order | CSR entries |
| --- | ---: | ---: | ---: |
| Interleaved feeders | 1,120 × 900 | 1,008,002 | 3,792,813 |
| Diagonal straps | 1,100 × 910 | 1,001,002 | 4,320,777 |
| Slotted power plane | 1,200 × 840 | 1,008,002 | 3,973,723 |
| Dual power rails | 2 × (1,000 × 500) | 1,000,002 | 4,837,391 |
| Bus-strap network | 1,250 × 800 | 1,000,002 | 3,653,469 |

- [Build and run](docs/BUILD-AND-RUN.md)
- [Capture hook and solver entry points](docs/INTEGRATION.md)
- [Circuit descriptions and original inputs](cases/README.md)
- [Native capture validation](validation/README.md)
- [Licenses and attribution](LICENSE-NOTICE.md)
