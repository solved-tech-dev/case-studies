# Reference benchmark

The table records the original measurements from the case study. Values are cumulative native solver seconds. GPU values are medians from three paired tests with alternating execution order on one NVIDIA A100 PCIe 40 GB. CPU configurations are single measurements on two Xeon CPU Max 9462 processors with 64 physical cores available. The PARDISO run uses 64 threads. MUMPS entries are MPI ranks × OpenMP threads per rank.

| Case | PanguLU | cuDSS | PARDISO | MUMPS 1×64 | MUMPS 8×8 | MUMPS 64×1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Interleaved feeders | 1.1618 | 2.4548 | 3.8038 | 6.8899 | 6.0671 | 6.0867 |
| Diagonal straps | 1.3826 | 2.8361 | 4.9675 | 9.8812 | 7.9120 | 8.2552 |
| Slotted power plane | 1.2764 | 2.4359 | 3.8649 | 7.3736 | 5.8363 | 5.8751 |
| Dual power rails | 1.4672 | 2.7103 | 4.6163 | 8.6359 | 6.8675 | 7.0958 |
| Bus-strap network | 1.1759 | 2.3390 | 3.7392 | 6.9179 | 5.7515 | 5.8528 |

[benchmark.json](benchmark.json) contains the original phase timings, ratios for each GPU pair and complete-vector validation metrics. Each solver analyzed and factored the matrix once and called solve three times. PanguLU uses 26.46 with L0/L1, cuDSS 1.0.0 uses 40 host reordering threads, PARDISO is from oneMKL 2026.0.2, and MUMPS is 5.9.1. PanguLU analysis time includes reordering and symbolic analysis as one native phase. Assembly, external format conversion, external transfers, validation and file I/O are outside these solver totals. The release scripts reproduce the ngspice capture and PanguLU replay. The other solver figures are the recorded comparisons.

The release build completed the native ngspice run for all five cases. The captured CSR, all RHS and KLU solution vectors, and the full physical plot match the original runs byte for byte, excluding the raw file title, date and build identification. See [capture-validation.json](capture-validation.json).
