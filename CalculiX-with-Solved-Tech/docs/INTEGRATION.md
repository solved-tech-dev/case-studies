# CalculiX–PanguLU integration

[English](INTEGRATION.md) | [简体中文](INTEGRATION.zh-CN.md) · [Home](../README.md)

## Solver entry points

```text
CalculiX: assemble A and b; apply constraints and boundary conditions
  → spooles(...) or spooles_factor(...) / spooles_solve(...)
  → adapter: convert storage, manage the PanguLU objects, solve
  → write x into the CalculiX solution buffer
  → CalculiX: update fields, evaluate convergence, write results
```

[ccx_pangulu.cpp](../adapter/ccx_pangulu.cpp) provides the C ABI functions called by the CalculiX SPOOLES solver branch. The build links this object in place of the CalculiX `spooles.c` wrapper. Finite element assembly, material evaluation, boundary conditions, Newton control, and field output remain in CalculiX.

[spooles_entry.h](../adapter/spooles_entry.h) declares the application-side interface. External SPOOLES routines are also linked for CalculiX constraint-processing helpers. These helpers are separate from the main linear-system solver entry.

## Matrix and vector flow

1. Read the assembled `ad/au/icol/irow` arrays and expand the symmetric system into complete CSR storage, including explicit zeros.
2. Convert CSR to CSC with 64-bit column pointers, 32-bit row indices, and FP64 coefficients. Keep the CSR-to-CSC index map for coefficient updates.
3. Create a PanguLU handle, matrix, and vectors. Set ordering to L0 and reuse level to L2.
4. On the first call, execute `ANALYSIS`, `NUMERIC_FACTORISATION`, and `SOLVE`.
5. On subsequent calls, compare the complete graph arrays. For an unchanged graph, refresh all matrix coefficients and execute `REFACTORISATION` and `SOLVE`. A changed graph starts a new analysis.
6. Transfer the current right-hand side to the GPU, solve, and copy the solution into the CalculiX buffer.
7. Evaluate the relative residual using the complete CSR matrix and the right-hand side supplied to that call.

The residual criterion is `||Ax−b||₂ / ||b||₂ ≤ 1e−8`; for a zero right-hand side, the denominator is 1. The adapter can perform up to two additional correction solves with the current factors. These solves are included in the native solve timings.

## Object lifetime and reuse

The PanguLU handle, matrix pattern, and coefficient arrays persist across graph-unchanged calls. Each factorization call updates the coefficient array and performs numerical factorization or refactorization.

The CalculiX cleanup entry restores the application thread setting. PanguLU objects are released when the graph changes or the process exits. The adapter uses process-local state and is intended for one model per process.

The supplied steady-state models each use one increment. Capacitor has one linear-system call; the other four examples have two. See [Checking a calculation](VALIDATION.md) for the log and field checks.

## Timing output

| Layer | Contents |
|---|---|
| `native` | Synchronized PanguLU analysis, numerical factorization/refactorization, and solve calls |
| `service` | Adapter work including storage conversion, object setup, coefficient updates, transfers, and native calls |
| `evidence` | Matrix/RHS capture and independent CSR residual evaluation |
| `lifecycle` | PanguLU object and device-buffer destruction |

Core solver time is the sum of nonoverlapping `native` entries. The analysis API combines ordering, symbolic analysis, and preprocessing. The `service` layer contains `native`, so the two totals are reported separately. Application elapsed time also includes CalculiX assembly and field output.

Captured arrays describe the actual inputs of each call. Cross-solver comparisons additionally require graph, coefficient, and RHS comparisons at corresponding physical states; later Newton right-hand sides may differ through floating-point trajectories.

## Interface scope

The adapter supports the five documented real, symmetric, unshifted steady scalar examples using one process and one GPU. It uses the default 32-bit CalculiX integer ABI. Other analysis types and concurrent models are outside the validated interface.

For compilation and execution, see [Build and run](BUILD-AND-RUN.md).
