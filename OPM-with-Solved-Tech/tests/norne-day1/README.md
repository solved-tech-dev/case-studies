# Norne Day-1 Validation

After building a route (`linglongamg` or `xianzhikrylov`) and preparing the Norne N64 day-1 deck (`NORNE1.DATA` + `GRID_PROPS.INC`, see [Testing Guide](../../docs/TESTING.md)), run day-1 verification with independent true-residual checking enabled:

```bash
bash scripts/run-case.sh --route=xianzhikrylov \
  --deck=/path/to/NORNE1.DATA \
  --output-dir=work/norne-day1-xk \
  --no-output --profile --true-residual
```

Compare converged solver residuals, timestep/Newton statistics, and (when completed within your scheduler wall limit) field summaries against `examples/norne/reference-day1-summary.json`.

For offline summary validation without GPU hardware:

```bash
python3 tests/norne-day1/check-day1-summary.py
```
