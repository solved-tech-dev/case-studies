# Results and validation

[English](RESULTS.md) | [简体中文](RESULTS.zh-CN.md) · [Home](../README.md)

## Benchmark hardware and software

| Configuration | Hardware / allocation |
|---|---|
| CPU solvers | 2 × Intel Xeon CPU Max 9462; 64 physical cores total, no SMT; 110 GiB host memory |
| GPU solvers | 1 × NVIDIA A100-PCIE-40GB; 40 host CPU cores and 100 GiB host memory |
| Application assembly and results | 4 threads |

The software configuration is CalculiX 2.23, PanguLU 26.42 ABI2 CUDA 12.9, cuDSS 0.8.0.10, oneMKL PARDISO 2024.2.2, MUMPS 5.9.1 with sequential MPI and oneMKL 2024.2.2, and SPOOLES 2.2.

PARDISO uses a 64-thread limit. MUMPS uses coordinated 8-thread execution with spread/cores binding, selecting the fastest numerically qualified tested configuration for each case. SPOOLES uses its default serial configuration. GPU backends have a 40-host-thread limit. PanguLU uses L0 ordering and L2 reuse. cuDSS uses MC64 maximum-diagonal-product matching for Capacitor and no MC64 for the other four cases.

## Core solver times

Values are seconds summed over all native solver phases and all linear systems in one complete application.

| Case | PanguLU | cuDSS | PARDISO | MUMPS | SPOOLES |
|---|---:|---:|---:|---:|---:|
| Dam | 6.068 | 10.496 | 7.059 | 12.049 | 49.583 |
| Ax6ht | 4.615 | 8.984 | 5.936 | 7.244 | 37.198 |
| Capacitor | 9.790 | 25.197 | 15.552 | 14.109 | 94.379 |
| Oneel20fi | 9.912 | 10.528 | 13.763 | 31.627 | 610.639 |
| Pan | 3.975 | 7.891 | 6.035 | 6.853 | 29.172 |

Full precision and configurations are available in [five-methods.json](../results/five-methods.json).

**Statistics.** PanguLU reports the minimum of three complete applications. cuDSS/PARDISO also report the minimum of three, except that Pan uses one sample and Oneel20fi PARDISO uses one sample. MUMPS and SPOOLES each use one numerically qualified complete application.

Core time includes the native analysis, numerical factorization/refactorization, and solve work exposed by each interface. Combined stages are reported together. Adapter service time and evidence capture are separate from core time; service contains native calls. Application elapsed time additionally includes assembly, output, and evidence work.

The solvers retain their respective factorization algorithms and object lifecycles. SPOOLES reanalyzes on each call; the adapters reuse analysis when the graph is unchanged. The compared inputs were checked for graph/coefficient agreement and first-RHS agreement. Later Newton right-hand sides can differ through floating-point trajectories.

## Numerical validation

The [validation record](../results/release-validation.json) contains one full-application validation run for each supplied refined model. Its per-phase timings describe those individual runs, independently of the benchmark statistics above.

All five runs completed with exit status 0 and the documented equation counts and linear-call trajectories. Every system passed the independent complete-CSR residual criterion 1e-8.

| Check | Observed result |
|---|---|
| Maximum independent relative residual | 1.285454117085e-10 |
| Complete printed temperature fields vs same-mesh references | Relative L2 difference 0 for all five cases |
| Maximum printed heat-flux relative L2 difference | 5.5985671437086086e-9 |
| Physical checks | Boundary values, balance, and case-specific analytic checks passed |

Field comparisons use the precision of CalculiX's native printed output. The small upstream meshes and the changed-coefficient test also pass their checks. The coefficient-update test retains the graph, changes three diagonal entries and the RHS, and recovers the known solution (1,2,3) after refactorization.
