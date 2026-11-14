# Validation

Validation checks numerical and physical results; it never enforces a target execution time.

## Per application run

`tools/run_case.py` checks the application exit code. `tools/validate.py` rejects reported nonlinear/linear divergence, incomplete physical time histories, missing required fields and nonfinite nodal or element fields. The PanguLU path also requires the expected provider, successful native phases, and finite original-CSR residuals no greater than `1e-8` for every reported RHS. Each RHS is checked independently inside the adapter before returning to PETSc. Zero RHS uses absolute residual.

For `hex20-heat`, the complete final temperature is compared with `200*x+100*y+200*z` at every node, with relative L2 tolerance `1e-6`.

## Reference comparison

Run the same case and profile with `--solver reference`, then pass that output directory to `tools/validate.py --reference`. The checker requires identical output times, coordinates and field sets, and compares every saved nodal/element/global field over the full time history. The acceptance bound is `||difference||₂ <= 1e-12 * sqrt(number of values) + 1e-6 * ||reference||₂`. The absolute RMS floor is in the model’s units and handles physically zero components represented by roundoff. Reports retain relative L2, absolute L2/RMS and maximum absolute errors, including large relative errors for near-zero components. This is an additional cross-solver check, not a substitute for the original linear-system residual.

The default `petsc` reference solver is intended for small grids. Exit code zero and field finiteness alone do not establish cross-solver agreement. A report without `reference_field_relative_l2` has not performed that comparison.

The annular torsion and shell cases also check every saved physical postprocessor and its complete time history in `physics.csv`. Reference comparison uses `||difference||₂ <= 1e-12 + 1e-6 * ||reference||₂`, with the absolute term expressed in the model's units. This accommodates zero-valued symmetry components; both absolute and relative errors are reported.

## Adapter tests

- `tests/test_sparse.cpp`: CSR/CSC correspondence, explicit zeros, value transitions, malformed-graph rejection and a known original-system residual. Runs without CUDA.
- `tests/test_native.cpp`: actual PanguLU factorization and block solve, exact-matrix factor reuse, changed-value refactorization and changed-graph reanalysis. Requires a licensed GPU.
- `tests/test_petsc.cpp`: native PETSc `MatSolve`/`MatMatSolve`, factor-object recreation with cache retention, and full coefficient updates. Requires a licensed GPU and the built plugin.
- `tools/verify_sources.py`: the nine preserved upstream inputs/seed meshes match `upstream.lock.json`, and bilingual documentation exists.

The native test counts actual analysis and numerical-factor calls. It does not infer reuse from a configuration label or timing reduction. The application check logs also include analysis/numeric counters.

## Scope of claims

Successful small-grid checks establish that these native application entries, fields and solver callbacks work together on the checked software stack. They are not new million-grid benchmark results. The large-grid profiles retain the previously prepared physical configurations; results should be validated when run on a new installation. This source package does not ship historical timings or certify a runtime, acceleration ratio, worker affinity, NUMA page placement or reproducibility of a particular performance measurement.
