# Component Smoke Tests (Requires GPU + SDK + License)

These two standalone tests verify each backend in isolation before running full OPM Norne N64 simulations:

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

2. `check-xianzhikrylov-bsr3.cpp`:
   - Uses ONLY public `XianzhiKrylovAPI.h` and `XianzhiKrylovData.h` (`xianzhikrylov_create`, `xianzhikrylov_config_set`, `xianzhikrylov_matrix_create`, `xianzhikrylov_vector_create`, `xianzhikrylov_execute`, `xianzhikrylov_print_info`, and destroy APIs).
   - Exercises staged `PHASE_ANALYSIS` + repeated `PHASE_NUMERIC_PREPARE` / `PHASE_SOLVE` on a 2-block `BSR3` (6x6 scalar) nonsymmetric matrix with `examples/configs/xk-bsr3-bicgstab-ilu0.json`.
   - Independently computes host true relative L2 residual (`< 0.005`, observed `~1.4e-16` / `1.6e-16`) and checks `/proc/self/maps` to confirm `liblinglongamg` is **not** loaded.

   Build & run:
   ```bash
   g++ -std=c++17 -O2 \
     -I"$XIANZHIKRYLOV_ROOT/include" -I"$CUDA_HOME/include" \
     tests/smoke/check-xianzhikrylov-bsr3.cpp \
     -L"$XIANZHIKRYLOV_ROOT/lib" -lXianzhiKrylov \
     -L"$CUDA_HOME/lib64" -lcudart \
     -Wl,--enable-new-dtags \
     -Wl,-rpath,"$XIANZHIKRYLOV_ROOT/lib" \
     -Wl,-rpath,"$CUDA_HOME/lib64" \
     -o work/check-xianzhikrylov-bsr3
   ./work/check-xianzhikrylov-bsr3 examples/configs/xk-bsr3-bicgstab-ilu0.json
   ```
