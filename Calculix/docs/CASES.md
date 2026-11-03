# Examples

[English](CASES.md) | [简体中文](CASES.zh-CN.md) · [Home](../README.md)

The first four examples come from the CalculiX 2.23 verification inputs. Pan comes from the CalculiX community Launcher collection. Each example includes its upstream input, a refined input, and parameters.

The table lists measured equation counts after constraints. Complete-storage nnz includes explicit zeros.

| Case | Element / mesh | Subdivision | Equations n | Complete nnz | Parameters |
|---|---|---|---:|---:|---|
| Dam | CPS8 / planar 2D | 32 | 984,474 | 15,390,186 | [JSON](../cases/dam/parameters.json) |
| Ax6ht | CAX6 / 2D axisymmetric | 125 | 1,002,499 | 11,509,973 | [JSON](../cases/ax6ht/parameters.json) |
| Capacitor | C3D20 / 3D | 3×16×400 | 952,073 | 47,134,131 | [JSON](../cases/capacitor/parameters.json) |
| Oneel20fi | C3D20R / 3D | 44 | 352,440 | 19,849,804 | [JSON](../cases/oneel20fi/parameters.json) |
| Pan | CAX6 / 2D axisymmetric | 16 | 890,630 | 10,206,762 | [JSON](../cases/pan/parameters.json) |

## Dam: steady diffusion in layered media

This planar example represents a groundwater/scalar-diffusion problem through the CalculiX heat-transfer formulation. It contains two orthotropic materials, with scalar values 8.5344 and 3.9624 prescribed at the upstream and downstream boundaries. The CPS8 parent-element geometry and material regions are retained during refinement.

Validation checks the complete scalar and flux fields, boundary values, and global balance.

## Ax6ht: axisymmetric heat transfer in a disk

Quadratic CAX6 triangles describe the axisymmetric region. The conductivity is 50, and the model includes density, specific heat, and `STEP,NLGEOM`. Temperatures 0 and 100 are prescribed at nodes 1 and 46.

These are point-temperature conditions. The localized field near the prescribed nodes changes with refinement, so this model is best interpreted as a numerical verification problem. Validation checks complete temperature/heat-flux output, point constraints, and global balance.

## Capacitor: scalar potential in a coaxial sector

The three-dimensional quadratic-hexahedron model represents a finite-angle sector between coaxial surfaces. It uses the heat-transfer scalar-diffusion formulation as an electrostatic-potential analogy, with coefficient 8.8541878176e-12 and inner/outer values 2 and 1. The temperature output represents the potential analogue.

Validation compares the full scalar/flux fields, boundary values, balance, logarithmic radial solution, and section flux. The geometric approximation of the quadratic sector is included in the comparison tolerances.

## Oneel20fi: a cube with convection

A unit cube is refined using C3D20R elements. Its conductivity is 50. One face has prescribed temperature 0; the opposite face has convection to ambient temperature −12 with heat-transfer coefficient 10. The step has time increment 1 and end time 1.

The steady solution is T=−2x with heat flux (100,0,0), enabling checks over all output nodes. The supplied model has 352,440 equations.

## Pan: heating and cooling of an iron pan

The axisymmetric pan profile uses conductivity 50 and a bottom flux of 31500. Convection coefficients are 10 on the air side and 150 on the food side, both with ambient temperature 293.15. Refinement maps each boundary condition to the corresponding child faces.

Validation checks the complete temperature and heat-flux fields and evaluates convection balance using the full-revolution 2πr weighting. CalculiX's native expanded-sector volume is reported separately.

## Inputs and execution

Each case stores `source/`, `input/model.inp.gz`, `parameters.json`, and `physics-cards.inp`. The runner decompresses the complete model automatically. `physics-cards.inp` is a compact reference for materials, loads, boundary conditions, and step controls; use the complete model for execution.

All five models use one steady-state increment. Capacitor has one main linear-system call; the other four have two. For commands, see [Build and run](BUILD-AND-RUN.md). Numerical checks are described in [Checking a calculation](VALIDATION.md).
