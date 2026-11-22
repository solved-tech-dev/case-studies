# How PanguLU is called from FEALPy

[简体中文](INTEGRATION.zh-CN.md) · [Build and run](BUILD-AND-RUN.md)

## From model to solution

1. The case script selects the NumPy backend, creates the upstream mesh and space, assembles the matrix and load, and applies the model's boundary terms. Entry points: [elastic3d/run.py](../cases/elastic3d/run.py), [maxwell2d/run.py](../cases/maxwell2d/run.py), [advdiff2d/run.py](../cases/advdiff2d/run.py). Maxwell incorporates the complex Robin term during assembly.
2. [adapter/driver.py](../adapter/driver.py) converts the assembled FEALPy object to a full SciPy CSR, sums duplicate coordinates and sorts column indices. Explicit zeros are retained. It selects the real or complex adapter from the matrix dtype and passes contiguous arrays through Python's `ctypes`.
3. [adapter/gpu.cpp](../adapter/gpu.cpp) exports `solve(n, row_ptr, columns, values, rhs, solution)` with C linkage. Indices start at zero. The row pointers use 64-bit integers, the column indices use 32-bit integers and both vectors have length `n`. Elasticity and transport use 64-bit real scalars. Maxwell uses 128-bit complex scalars for the matrix, load and solution.
4. [adapter/gpu_backend.inc](../adapter/gpu_backend.inc) converts CSR to CSC without transposing the mathematical operator. It configures L0 ordering and general LU, supplies the full CSC arrays to `pangulu_matrix_create()`, and constructs the device vectors with `pangulu_multiple_vector_create()` using one column and leading dimension `n`. It calls `pangulu_execute()` for analysis, numeric factorization and solve, synchronizes the GPU, and returns the result to host memory. Each case uses one system and one RHS, so it performs one analysis and one factorization. PanguLU and CUDA resources are released when the call returns, including on exceptions.
5. The driver checks the returned solution against the original CSR system. The case then assigns it to the FEALPy field and compares against its analytical solution. Maxwell evaluates the vector field with its edge-element basis before integrating the error. Successful checks allow the coefficients to be saved with `--output`.

The Python entry points use the installed upstream FEALPy package. [Upstream source used by the cases](https://github.com/weihuayi/fealpy/tree/e033533bd8e24f796f9c5284685b20a98b4d11cb)
