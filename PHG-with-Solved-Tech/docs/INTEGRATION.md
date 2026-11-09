# Integration and outputs

[简体中文](INTEGRATION.zh-CN.md)

`tools/prepare.py` registers the `pangulu` solver in PHG and adds an optional hook for exporting native solution fields. The official PHG and PHG-Solid archives provide the case sources, materials, mesh reader, assembly routines, and boundary-condition handling. `tools/run_case.py` writes the selected case configuration to the run directory.

`integration/phg/solver-pangulu.c` implements the PHG solver callbacks. It collects all entries of the assembled matrix in CSR format and passes the matrix and right-hand side to `phg_pangulu_solve_csr()` in `adapter/pangulu_bridge.cpp`. The adapter converts CSR to CSC while preserving every stored entry, calls `pangulu_execute()` for analysis, numerical factorization, and solve, and writes the solution back to the PHG vector. All three cases use real double-precision general LU with one right-hand side per application call. Each adaptive step analyzes and factors the system assembled on its current mesh.

The application log records the native solver stages for diagnostics. The run report checks stage return codes, matrix dimensions, the relative residual against the original CSR matrix, and whether every value in the exported solution fields is finite. The relative residual is `||Ax − b||₂ / ||b||₂`, with a threshold of `1e-8`. The solid cases export `Displacement`, and the spherical interface case exports `u1_h` and `u2_h`. Use `--capture-linear` to save the original CSR matrix, right-hand side, and full solution for each call. These files require additional disk space on large meshes.

For development, run `python3 tools/validate_connections.py --candidate pangulu --output _work/check-gpu` to compare the full GPU solution fields for the three small cases against the project's CPU reference results. Use the `--reference-only` build to check the integration on a CPU.
