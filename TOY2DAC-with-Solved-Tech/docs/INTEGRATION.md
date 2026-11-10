# PanguLU integration

[English](INTEGRATION.md) · [简体中文](INTEGRATION.zh-CN.md)

- [apps/marmousi3d/model.c](../apps/marmousi3d/model.c) generates the 3D truth model from the upstream velocity profile and assembles a complex COO matrix from the current velocity.
- [apps/marmousi3d/fwi.cpp](../apps/marmousi3d/fwi.cpp) generates synthetic observations, computes forward and conjugate-transpose adjoint wavefields, and updates P-wave velocity with L-BFGS and a line search. Both systems use `case_factor()`, `case_solve_many()`, and `case_cleanup()`.
- [adapter/pangulu_fwi_adapter.cpp](../adapter/pangulu_fwi_adapter.cpp) converts zero-based COO to CSC sorted by row within each column and preserves every stored entry. It configures L0 ordering and double-precision complex general LU, then calls `pangulu_execute()` for native analysis, factorization, and solve. It uses `pangulu_multiple_vector_create()` to describe one GPU right-hand-side column and one solution column, and returns the solution to the application.

Each system creates and releases its own solver state, so native analysis and factorization run for every system. The application checks each returned wavefield against the original COO system. [Case definitions](../cases/fwi-grids.json) specify the grids, frequencies, and matrix sizes.
