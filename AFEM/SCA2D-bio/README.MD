# Summary of User-Defined Material Fortran Files

This folder contains four user-defined Fortran source files:

- `usca2d.for`: Abaqus `UMAT` entry point. It passes Abaqus strain, stress, state variable, property, and characteristic-length data into `sca2d(...)` and includes the other source files at compile time.
- `sca2d.for`: Main constitutive model. This is where the 2D smeared crack / crack-band material response is implemented.
- `myaux.for`: Linear-algebra utilities, mainly matrix inversion (`inverse`, `qinverse`) and the pivoting routines they depend on.
- `myStandardSupport.for`: Thin Abaqus/Standard compatibility helpers, including `myExit()` and a wrapper around `SPRIND`.

## What `sca2d.for` is doing

`sca2d.for` implements a 2D orthotropic damage model with crack-band regularization for an Abaqus material point. In engineering strain ordering, it works with:

- `eps11`
- `eps22`
- `gamma12`

The routine starts from an undamaged orthotropic elastic law, checks whether fiber or matrix cracking has initiated, and then, once cracking exists, solves for crack strains and reduces the effective stiffness accordingly.

## Inputs and interpretation

### Material properties

`sca2d` reads the material constants from `props` in this order:

1. `E11`: elastic modulus in the fiber direction
2. `E22`: elastic modulus transverse to the fiber
3. `G12`: in-plane shear modulus
4. `nu12`: Poisson ratio
5. `sigcrf0`: fiber-direction tensile cracking stress
6. `sigcrm0`: matrix-direction tensile cracking stress
7. `taucrm0`: matrix shear cracking stress
8. `GICf`: fiber mode-I fracture energy
9. `GICm`: matrix mode-I fracture energy
10. `GIICm`: matrix mode-II fracture energy
11. `GGICf`: residual shear modulus used after fiber breakage
12. `alpha`: scaling factor used in post-failure stiffness reduction

`damp` is hard-coded to `0.01d0` and used to build a diagonal viscous regularization matrix for the crack-strain solve.

### State variables

The code documents the state variables as:

1. characteristic length in direction 1
2. characteristic length in direction 2
3. maximum crack strain 11
4. maximum crack strain 21 (engineering), used here as fiber shear crack history
5. maximum crack strain 22
6. maximum crack strain 12 (engineering), used here as matrix shear crack history
7. normalized dissipated fiber mode-I energy
8. normalized dissipated matrix mode-I energy
9. normalized dissipated matrix mode-II energy
10. deletion flag
11. last crack strain 11
12. last crack strain 21
13. last crack strain 22
14. last crack strain 12

In the current version, `statev(1:2)` are not actually used for the characteristic length because the routine sets:

- `l1 = L`
- `l2 = L`

where `L` comes from the caller and is `CELENT` in the Abaqus `UMAT`.

## High-level algorithm inside `sca2d`

### 1. Build the undamaged orthotropic elasticity

The routine forms the compliance matrix

\[
S_{co} =
\begin{bmatrix}
1/E11 & -\nu12/E11 & 0 \\
-\nu12/E11 & 1/E22 & 0 \\
0 & 0 & 1/G12
\end{bmatrix}
\]

and then inverts it with `inverse(...)` to get the continuum stiffness `Dco`.

This `Dco` is the intact material stiffness used before cracking and also the background continuum stiffness used after cracking.

### 2. Define crack transformation operators

Two fixed transformation matrices are defined:

- `N1`: used for matrix-crack kinematics
- `N2`: used for fiber-crack kinematics

These map the four crack-strain internal variables into the 2D engineering strain/stress space. The branch taken later decides which transformation is active.

### 3. Check crack-band regularization limits

Before evolving damage, the routine checks whether the element characteristic length is too large for the specified fracture energies and strengths:

- `2*GICf*E11/sigcrf0^2 < l1`
- `2*GICm*E22/sigcrm0^2 < l2`
- `2*GIICm*G12/taucrm0^2 < l2`

If any of these hold, the code prints diagnostics and calls `myExit()`. This is a mesh-objectivity safeguard: the crack-band softening law requires an element length small enough to represent the chosen fracture energies without becoming nonphysical.

### 4. Recover crack history

The code reads the historical maxima:

- `epscrfmax = statev(3)`
- `gamcrfmax = statev(4)`
- `epscrmmax = statev(5)`
- `gamcrmmax = statev(6)`

These history variables make the damage evolution irreversible by always using the maximum crack strain seen so far.

### 5. If no crack has formed yet, use pure elasticity

If all crack-history variables are essentially zero, the code computes:

- `sig = Dco * eps`
- `ddsdde = Dco` when `isImplicit = 1`

It then checks whether any initiation criterion has been exceeded:

- `sig(1) > sigcrf0` for fiber tension
- `sig(2) > sigcrm0` for matrix tension
- `abs(sig(3)) > taucrm0` for matrix shear

If any is exceeded, it does not immediately compute a finite crack opening. Instead, it seeds all crack variables with a tiny value `1.0d-14`. That tiny seed is used to switch the model from the intact elastic branch into the cracked branch on the next pass.

### 6. Matrix-crack branch

If matrix cracking has started (`epscrmmax` or `gamcrmmax` is nonzero) and fiber cracking has not yet reached its threshold (`epscrfmax <= 1e-8`), the code takes the matrix-crack branch.

In this branch it:

1. Reads the old crack strain vector from `statev(11:14)`.
2. Builds a diagonal damping matrix `Dda = damp * I`.
3. Calls `calcEpscr(...)` with transformation matrix `N1` to solve for the current crack strain `epscr`.
4. Calls `calcDcr(...)` to compute the crack constitutive matrix `Dcr`.
5. Forms the coupled continuum-crack system:
   - `term11 = N1^T * Dco * N1`
   - `term22 = N1^T * Dco`
   - `invDcrDco = inverse(Dcr + term11 + (1/dt) * Dda)`
6. Builds the degraded effective tangent:
   - `Dcocr = Dco - term22^T * invDcrDco * term22`
7. Computes the stress:
   - `sig = Dcocr * eps - term22^T * invDcrDco * (1/dt) * Dda * epscrold`
8. Returns `Dcocr` as the Jacobian if Abaqus requested an implicit tangent.

So, mechanically, the code is treating the crack strain as an internal variable coupled to the continuum strain and solving for a reduced tangent stiffness after matrix cracking.

### 7. Fiber-crack branch

If `epscrfmax > 1e-8`, the code switches to the fiber-crack branch. The structure is almost identical to the matrix-crack branch, but it uses transformation matrix `N2` instead of `N1`.

That means:

- crack strain is solved with different kinematic coupling
- the same `calcDcr(...)` softening law is reused
- the same coupled reduction formula is used to obtain the effective stiffness and stress

This branch is intended to represent fiber tensile failure.

## How crack strain is solved: `calcEpscr`

`calcEpscr(...)` is one of the most important parts of the file.

It solves for the crack strain vector `epscr(4)` using Newton's method. The comments note that Broyden's method was attempted but did not work because the Jacobian varied too much.

### Key features of `calcEpscr`

- It performs the internal calculations in quadruple precision (`real*16`) because the crack initiation seed is extremely small.
- It starts from `x = 0` as the initial guess for the crack strain vector.
- At each Newton iteration, it:
  1. Recomputes the crack stiffness `qDcr` using `qcalcDcr(...)`.
  2. Forms the residual
     \[
     F = (Dcr + N^T Dco N + (1/dt)Dda)x - N^T Dco \, \varepsilon - (1/dt)Dda\,\varepsilon_{cr}^{old}
     \]
  3. Approximates the Jacobian numerically by finite differences, perturbing one crack-strain component at a time.
  4. Inverts that Jacobian with `qinverse(...)`.
  5. Updates the solution using the Newton step.
- It stops when the sum of absolute residual components is below `1e-4`, or after `50` iterations.

If convergence fails, the code prints a detailed diagnostic dump of the current matrices, strain, history variables, material parameters, and final residual, then exits Abaqus via `myExit()`.

So this routine is effectively the local constitutive solver for the internal crack variables.

## How the crack constitutive matrix is built: `calcDcr`

`calcDcr(...)` constructs a diagonal `4 x 4` matrix `Dcr` representing the crack-law tangent in four internal crack modes. It uses the current crack strain and the maximum historical crack strain to enforce irreversibility:

- `epscr_temp(i) = max(abs(epscr(i)), abs(epscrmax(i)))`

That means unloading does not restore stiffness; the softening state depends on the largest crack strain attained.

### Individual entries

- `Dcr(1,1)`: fiber mode-I softening
- `Dcr(2,2)`: a special "fiber mode-II" shear-related stiffness, not a true mode-II law according to the author’s comment
- `Dcr(3,3)`: matrix mode-I softening
- `Dcr(4,4)`: matrix mode-II softening

For the normal opening/shear softening laws, the code uses a linear descending relation in traction-separation slope form, capped by a very small residual stress (`1e-10` times the original strength). This keeps the stiffness finite instead of dropping exactly to zero.

### Special DCB-specific modification

Near the end of `calcDcr(...)`, there is a hard-coded block marked:

`This is for DCB`

Once matrix crack opening history `epscrmax(3)` exceeds `1e-8`, the code overrides part of the crack stiffness:

- `Dcr(4,4) = alpha * Dcr(3,3) / (2*(1+0.25))`
- `Dcr(1,1) = alpha^2 * Dcr(3,3)`

This is a problem-specific coupling for the Double Cantilever Beam setup, making shear and fiber-direction crack stiffness depend on the matrix mode-I degradation. That block is important because it means the routine is not a purely general orthotropic damage law; it contains DCB-specific assumptions.

`qcalcDcr(...)` is simply the quadruple-precision version used by the Newton solver.

## Stress update and tangent update

Once the crack strain and crack stiffness are known, the cracked material response is assembled as:

- continuum stiffness `Dco`
- crack stiffness `Dcr`
- transformation matrix `N1` or `N2`
- damping contribution `(1/dt)Dda`

The resulting effective tangent `Dcocr` is the continuum stiffness reduced by the crack compliance contribution. Stress is then computed using both the current total strain and the previous crack strain.

For Abaqus/Standard, `usca2d.for` passes `isImplicit = 1`, so `sca2d` always returns `ddsdde` as either:

- the intact stiffness `Dco`, or
- the degraded cracked tangent `Dcocr`

## State-variable updates at the end

At the end of `sca2d`, the routine stores:

- updated maximum crack histories in `statev(3:6)`
- current crack strain in `statev(11:14)`

The update uses `max(...)` or `max(abs(...), ...)`, which again enforces damage irreversibility.

There is also a deletion check:

- if `statev(7)`, `statev(8)`, or `statev(9)` exceed `0.95`, then `statev(10)=0`

However, in this version the code that accumulates the dissipated energies into `statev(7:9)` is commented out. So the deletion flag logic remains present, but the energy variables are not actively updated in the shown implementation.

## Important modeling assumptions and quirks

- The file is written for 2D engineering strain/stress ordering only.
- The routine assumes tensile cracking and does not contain a developed compression or crack-closing treatment. The very first comments explicitly list "compression" and "crack closing" as unresolved items.
- Characteristic length is currently taken directly from the Abaqus element length `L`, even though the state-variable list suggests directional lengths could be stored.
- Damage initiation is triggered by direct comparison of stresses against strengths.
- Crack evolution is irreversible and regularized by fracture energy and characteristic length.
- A small damping matrix is included in the crack solve, which likely helps numerical stability.
- The DCB-specific override in `calcDcr` makes the model specialized to this application.
- The code uses both double precision and quadruple precision, with the latter reserved for the local Newton solve.

## How the files connect together

The execution path in Abaqus/Standard is:

1. Abaqus calls `UMAT` in `usca2d.for`.
2. `UMAT` passes total strain `STRAN + DSTRAN`, stress, state variables, properties, time step, and `CELENT` into `sca2d(...)`.
3. `sca2d(...)` computes the updated stress and tangent.
4. `myaux.for` supplies the matrix inversion routines used by `sca2d`, `calcEpscr`, and `calcDcr`.
5. `myStandardSupport.for` supplies `myExit()` used when the routine encounters unrecoverable conditions.

## Bottom line

`sca2d.for` is a 2D orthotropic crack-band constitutive law for Abaqus that:

- starts from linear orthotropic elasticity,
- triggers damage when fiber tension, matrix tension, or matrix shear strengths are exceeded,
- solves for internal crack strains through a local Newton iteration,
- degrades the effective stiffness using fracture-energy-based softening,
- uses characteristic-length regularization to maintain mesh-objective energy dissipation,
- and includes a DCB-specific post-cracking stiffness coupling that tailors the model to this fracture problem.

In short, this file is the heart of the user material: it computes how the material transitions from intact orthotropic elasticity to a softened cracked response with history dependence and fracture-energy control.
