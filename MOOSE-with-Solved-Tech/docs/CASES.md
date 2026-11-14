# Cases and grids

All six applications come from MOOSE tutorials or module tests. `cases/<case>/upstream/` preserves the original input and required seed mesh. `upstream.lock.json` records original paths, a pinned revision and file digests. Profiles in `case.json` and additive options in `tools/run_case.py` configure the solver and grid without rewriting upstream files.

## 1. Thermomechanical coupling — `thermomech`

A 2-by-1 rectangular body heats internally while its two ends follow prescribed temperatures. The initial temperature is 300 K, the left boundary stays at 300 K, the right boundary follows `300+5*t`, and the volume heat source is `5e4`. Vertical displacement is constrained at the bottom and horizontal displacement at the origin. MOOSE solves temperature and two displacements together, including thermal expansion, finite-strain elasticity and auxiliary stress outputs.

Original thermal, density, elastic and expansion parameters remain intact. The time step is 1 and the end time is 5. The small grid is 10×10; the enlarged grid is 666×500 quadrilaterals, with 334,167 nodes and 1,002,501 primary unknowns.

The runner explicitly connects `temperature=T`, `eigenstrain_names=thermal_expansion` and `use_finite_deform_jacobian=true`, matching the fully coupled Jacobian configuration used for these cases. These are disclosed application-level additions, not hidden solver changes. Analysis is reused for an unchanged graph; changed coefficients trigger numerical refactorization. Newton iteration counts are determined by application convergence, not hard-coded.

## 2. Annular torsion — `annular-torsion`

A hollow cylinder has inner radius 0.95, outer radius 1 and length 2. All bottom displacements are fixed. The original `ADTorque` conditions apply axial torque at the top. The application computes displacement, rotation angle and stresses using the upstream small-strain elastic model, polar-moment calculation and single load step.

The small grid retains 1 radial layer, 30 circumferential divisions and 5 axial layers. The enlarged grid uses 16×352×55 divisions: 309,760 HEX8 elements, 335,104 nodes and 1,005,312 displacement unknowns. Additional stress and strain component output supports full-field validation.

## 3. Cup-shaped diffusion — `mug`

This is the `ex01_inputfile` MOOSE example on the supplied `mug.e` mesh. A steady diffusion equation is solved through the cup wall and base, with the scalar fixed to 1 at the bottom and 0 at the top. Other boundaries retain the original natural conditions. The field `diffused` is a scalar diffusion variable; it is not relabeled as fluid velocity or thermodynamic temperature.

The small profile uses the seed directly. Three uniform refinements produce 1,267,712 HEX8 elements and 1,349,545 unknowns. The actual mesh contains the wall and base; no handle or new geometry is added.

## 4. Pinched cylindrical shell — `pinched-cylinder`

Opposing point loads pinch a thin cylindrical shell at mid-length. Symmetry reduces the calculation to one eighth of the geometry. The reference cylinder has radius 1, length 2, thickness 0.01, Young's modulus `1e6` and Poisson ratio 0.3. The reduced-model `ConstantPointSource` value is −2.5. The five primary fields are three displacements and two local rotations.

Original displacement/rotation constraints and named boundaries are retained. The upstream comments discuss end idealization; the actual input constraints must not be replaced by fully clamped ends. The small mesh is the original 10×10 surface mesh. Dividing every parent QUAD4 edge into 45 segments produces 202,500 quadrilaterals and 1,017,005 unknowns.

`tools/refine_quad.py` performs bilinear interpolation within each original parent, preserves blocks and boundary sets, and does not reproject nodes onto an analytic cylinder. This reproduces the enlarged grid's geometry policy.

## 5. Thermal contact across a gap — `gap-heat`

Two conducting solids exchange heat across a gap. Left and right external temperatures are 1 and 2. Solid conductivity is 1, gap conductance is 2.5, and both emissivities are zero. Temperature unknowns exist only on blocks 1 and 3; the inactive gap geometry therefore contributes nodes without equivalent temperature equations.

The small grid uses four subdivisions per parent edge. The enlarged grid uses 70, yielding 980,000 active quadrilaterals and 982,802 temperature unknowns. The full geometric mesh additionally includes the inactive gap block.

The coarse grid needs more iterations to reach the original convergence criterion. Only the `small` profile sets `nl_max_its=100`; convergence tolerances are unchanged. The `million` profile retains the original iteration limit.

Gmsh imports shared edges on both adjacent elements. `tools/fix_gap_boundary.py` restores the original active-material side of each contact surface while preserving geometry, connectivity and blocks. Omitting this step changes the contact definition. RHS vectors remain in the original application order. Exactly equal full graphs and coefficients reuse the existing numerical factors.

## 6. Quadratic heat-conduction patch — `hex20-heat`

A unit cube has exterior temperature `T=200*x+100*y+200*z` and conductivity `4.85e-4`. HEX20 geometry and second-order Lagrange temperature interpolation are used. The same linear expression gives the exact interior temperature and an independent full-field check.

The upstream input uses a distorted coarse seed. This package preserves its materials, temperature function, quadrature and solve controls, while generating a structured unit cube: 4³ HEX20 elements for the small profile or 62³ for the enlarged profile. The latter has 238,328 elements and 988,281 unknowns. All exterior faces receive boundary ID 10. This is a disclosed mesh replacement, not subdivision of the distorted seed.

## Common execution path

All cases use a full assembled Jacobian, `KSPPREONLY` and general LU. Inputs originally using PJFNK are switched to the assembled-Jacobian Newton path through additive options, equally for PanguLU and the reference solver. Equations, materials, loads and time intervals are preserved.

Matrix dimensions describe these grids and equation-numbering settings; they are not timing requirements. Execution reports non-overlapping native phases and complete physical output. See [validation](VALIDATION.md).
