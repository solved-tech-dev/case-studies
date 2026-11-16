# ngspice capture and PanguLU replay

| File | Role |
| --- | --- |
| [`integration/ngspice/capture.patch`](../integration/ngspice/capture.patch) | Adds an optional capture call at the real solve site in ngspice's `src/maths/ni/niiter.c` and adds the capture source to the KLU-enabled build. |
| [`integration/ngspice/capture.c`](../integration/ngspice/capture.c) | `SolvedTechCaptureSolve()` reads the assembled KLU CSC and RHS, records the full original system as CSR, calls the native `SMPsolve()` and records the KLU reference vector. |
| [`adapter/replay.cpp`](../adapter/replay.cpp) | Converts the full captured CSR to CSC, sets L0 ordering and L1 factor reuse, and calls `pangulu_execute()` for analysis, numeric LU and the three solves. |
| [`tools/run_case.py`](../tools/run_case.py) | Checks the deck, native operating-point output, unchanged structure and values, solve sequence, residuals and reference-vector agreement. |

The upstream schematic/netlist flow, device equations, matrix assembly and operating-point iteration all run inside ngspice. The hook is active when `NGSPICE_CAPTURE` names a capture directory. For these published circuits, the separate replay then verifies the linear solver against the captured trajectory. The replay takes the matrix once and reads each captured RHS in order. The three solves belong to one operating point.

Capture uses real double-precision entries and preserves explicitly stored zeros. Files are little-endian on the supported x86-64 platform. Row pointers and column indices are zero based, stored as `.rp.i64` and `.ci.i32`. Coefficients, RHS and KLU solution vectors are `.a.f64`, `.b.f64` and `.x.f64`. The capture log refers back to a prior graph or coefficients when the bytes are unchanged. Each solve gets its own RHS and solution files.

The release driver expects the fixed-matrix trajectories in the case manifest. Test qualification requires every full original-CSR relative residual to be at most 1e-8 and every full-vector relative L2 deviation from KLU to be at most 1e-6. It checks one native analysis, one factorization and all three solves.
