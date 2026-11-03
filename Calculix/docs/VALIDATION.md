# Checking a calculation

[English](VALIDATION.md) | [简体中文](VALIDATION.zh-CN.md) · [Home](../README.md)

`tools/run_case.py` runs CalculiX and checks its solver log and native physical fields. The checks use the input model and the output produced by that run.

The application may print elapsed-time diagnostics. Their values are not correctness criteria, and the log checker does not evaluate or aggregate running times.

## Linear-system checks

`tools/check_log.py` checks application completion, equation counts, solver-call counts, PanguLU provider information, and the analysis/refactorization sequence. Every linear system must pass an independent residual check using its complete original CSR matrix and right-hand side:

```text
||Ax - b||₂ / ||b||₂ <= 1e-8
```

For a zero right-hand side, the denominator is 1. The matrix capture retains explicit zeros. For graph-unchanged calls, the adapter refreshes all coefficients and performs numerical refactorization.

## Physical-field checks

`tools/validate_fields.py` reads the native `model.frd` and `model.dat` output. It checks prescribed boundary values, finite fields, balance, and the analytic or physical conditions applicable to each example. These checks complement, rather than replace, the linear-system residual check.

For a full-field comparison against an independently supplied result on the same mesh, pass `--reference-frd /path/to/reference.frd.gz` to the runner. Printed-field precision is distinct from the linear residual tolerance.

## Output files

Inspect `application-exit.json`, `log-validation.json`, `field-validation.json`, and `validation.json` in the selected run directory. Keep the application log and native fields when investigating a failed check. See [Build and run](BUILD-AND-RUN.md) for commands and [Integration](INTEGRATION.md) for the solver data flow.

## Source-level checks

Run input-preparation checks with:

```bash
python3 -m unittest discover -s tests -v
```

`tests/changed_values.cpp` checks reuse with an unchanged graph, modified matrix coefficients, and a modified right-hand side. Build and run instructions are included in [Build and run](BUILD-AND-RUN.md).
