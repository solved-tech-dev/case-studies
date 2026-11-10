# Integration entry points

[English](INTEGRATION.md) · [简体中文](INTEGRATION.zh-CN.md)

- [apps/marmousi3d/frequency_elastic.c](../apps/marmousi3d/frequency_elastic.c) reads the upstream velocity profile, builds the 3D model, assembles the complex COO matrix and point-force right-hand side, and writes the displacement field.
- [adapter/pangulu_adapter.h](../adapter/pangulu_adapter.h) defines the C interface. The application passes the original matrix to `case_factor`, sends one complex right-hand side to `case_solve`, and releases solver state with `case_cleanup`.
- [adapter/pangulu_adapter.cpp](../adapter/pangulu_adapter.cpp) converts zero-based COO to CSC with sorted row indices, retaining all stored entries. It creates a complex FP64 general LU system with L0 ordering and L1 reuse, calls the native analysis and factorization phases, then solves one column using `pangulu_multiple_vector_create`. The solved real and imaginary components return to the application for the original COO residual calculation and output.

The adapter owns the GPU vectors and the PanguLU handles. The forward application owns the model, assembled matrix, source and output. [Case parameters](../cases/grids.json) supply the spacing, element counts, frequency and expected matrix sizes for each run.
