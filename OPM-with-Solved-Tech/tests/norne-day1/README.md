# Norne N64 Validation

After building a route (`linglongamg` or `xianzhikrylov`) and preparing the Norne N64 input deck (`NORNE_N64.DATA` + `GRID_PROPS.INC`, see [Testing Guide](../../docs/TESTING.md)), run verification with independent true-residual checking enabled:

```bash
bash scripts/run-case.sh --route=xianzhikrylov \
  --deck=/path/to/NORNE_N64.DATA \
  --output-dir=work/norne-xk \
  --no-output --profile --true-residual
```

Compare converged residuals and timing breakdowns against `examples/norne/reference-day1-summary.json`.

For offline summary validation without GPU hardware:

```bash
python3 tests/norne-day1/check-day1-summary.py
```
