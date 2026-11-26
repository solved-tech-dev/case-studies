# Norne Day-1 Validation

After building a route (`linglongamg`, `xianzhikrylov`, or `cpr`) and obtaining the Norne N64 day-1 deck (`NORNE1.DATA` + `GRID_PROPS.INC`), run day-1 verification with independent true-residual checking enabled:

```bash
bash scripts/run-case.sh --route=cpr \
  --deck=/path/to/NORNE1.DATA \
  --output-dir=work/norne-day1-cpr \
  --profile --true-residual
```

Then compare the run metrics (`FPR`, `FOPT`, `FWPT`, `FGPT`, `FWIT`, `FGIT`, timestep count, Newton iterations, linear iterations, and independent residual bounds) against `examples/norne/reference-day1-summary.json`.

For offline summary validation without GPU hardware:

```bash
python3 tests/norne-day1/check-day1-summary.py
```
