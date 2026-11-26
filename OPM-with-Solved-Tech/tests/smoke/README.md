# Component Smoke Tests (Requires GPU + SDK + License)

These two standalone CUDA smoke tests verify each backend in isolation before running full OPM Norne N64 simulations:

1. `check-linglong-cycle.cu`:
   - Loads `liblinglongamg.so` via `dlopen` (validating SDK license).
   - Solves a 4x4 CSR single V-cycle step (`linglongamg_apply_cycle` or single-step `linglongamg_solve` + `linglongamg_get_solve_info` checking `iterations == 1`).
   - Verifies residual reduction (`||b - A x||_2 / ||b||_2 < 1.0`).
   - Checks `/proc/self/maps` to confirm `libXianzhiKrylov` is **not** loaded.

   Build & run:
   ```bash
   "$CUDA_HOME/bin/nvcc" -std=c++17 -O2 \
     -I"$LINGLONGAMG_INCLUDE_DIR" \
     tests/smoke/check-linglong-cycle.cu -ldl -o work/check-linglong-cycle
   ./work/check-linglong-cycle "$LINGLONGAMG_ROOT/lib/liblinglongamg.so" examples/configs/linglong-vcycle.json
   ```

2. `check-xianzhikrylov-ilu.cu`:
   - Exercises `xk_internal_cprw_krylov_create` / `solve` / `destroy` for both BiCGStab (`method=3`) and GMRES (`method=7`) with device callbacks.
   - Exercises `xk_internal_cprw_ilu_create`, two `xk_internal_cprw_ilu_prepare` numeric refreshes, and repeated `xk_internal_cprw_ilu_apply` calls on a 2-block `BSR3` matrix.
   - Checks `/proc/self/maps` to confirm `liblinglongamg` is **not** loaded.

   Build & run:
   ```bash
   "$CUDA_HOME/bin/nvcc" -std=c++17 -O2 \
     -I"$XIANZHIKRYLOV_ROOT/include" \
     tests/smoke/check-xianzhikrylov-ilu.cu \
     "$XIANZHIKRYLOV_ROOT/lib/libXianzhiKrylov.so" \
     -Wl,-rpath,"$XIANZHIKRYLOV_ROOT/lib" \
     -o work/check-xianzhikrylov-ilu
   ./work/check-xianzhikrylov-ilu
   ```
