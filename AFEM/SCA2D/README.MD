# SCA2D - 2D Incremental Smeared Crack Approach (Isotropic)

An Abaqus/Standard UMAT implementation of the **Incremental Smeared Crack Approach (SCA)** for modeling crack initiation and propagation in 2D isotropic materials under plane stress conditions. This formulation avoids the need for a Newton-Raphson iterative scheme by using an incremental (AFEM) approach.

**Last Updated:** February 25, 2026
**Contact:** meka1@purdue.edu

---

## Overview

The SCA2D code models fracture in isotropic materials by:

1. Detecting crack initiation using a failure criterion (Max Principal Stress or Max Stress)
2. Tracking the crack orientation via principal stress directions at the point of initiation
3. Computing incremental crack strains in a local coordinate system aligned with the crack
4. Degrading material stiffness through a linear traction-separation law
5. Supporting unloading/reloading behavior with history-dependent crack stiffness

The implementation enforces the **Bazant size-effect limit** to ensure mesh-objectivity by checking that the element characteristic length is small enough relative to the fracture process zone.

---

## File Descriptions

### 1. `main_SCA2D_inc.for` - Abaqus UMAT Interface

**Purpose:** Entry point that bridges Abaqus/Standard with the SCA formulation.

- Implements the standard `UMAT` subroutine signature required by Abaqus for user-defined material models
- Designed for **CPS4 plane stress elements**
- Passes all Abaqus-provided variables (stress, strain, time, material properties, element geometry) to the core formulation
- Includes the two dependency files: `formulation_AFEM2D.for` and `matrix_inverse.for`
- Sets the implicit integration flag (`isIMPLICIT=1`) for consistent tangent stiffness computation
- Uses `CELENT` (characteristic element length provided by Abaqus) directly as the crack band width

**Call Hierarchy:**
```
Abaqus/Standard
  └── UMAT()                          [main_SCA2D_inc.for]
        └── sca2diso()                 [formulation_AFEM2D.for]
              ├── STIFFNESS_MATRIX_2D_PS()
              ├── mysprind() → SPRIND()  [Abaqus utility]
              ├── calcDcr()
              ├── calcEcr()
              │     ├── calcDcr()
              │     └── matrixInverse()  [matrix_inverse.for]
              └── calcDcocr()
                    └── matrixInverse()  [matrix_inverse.for]
```

---

### 2. `formulation_AFEM2D.for` - Core SCA Formulation

**Purpose:** Contains the constitutive model and all mechanics subroutines for the incremental smeared crack approach.

#### Material Properties (Input via Abaqus `PROPS` array)

| Index | Property | Symbol | Description |
|-------|----------|--------|-------------|
| 1 | `E` | *E* | Young's Modulus |
| 2 | `nu` | *v* | Poisson's Ratio |
| 3 | `GIc` | *G_Ic* | Mode-I Fracture Toughness |
| 4 | `XT` | *X_T* | Tensile Strength |
| 5 | `INDEX` | - | Failure Criterion (0 = Max Principal Stress, 1 = Max Stress) |

#### State Variables (STATEV)

| Index | Description |
|-------|-------------|
| 1 | Crack flag (0 = intact, 1 = cracked) |
| 2 | Element characteristic length *L* |
| 3 | Failure mode (+1 = tension, -1 = compression) |
| 4 | Maximum failure index at initiation |
| 5 | Crack normal x-component (cos *theta*) |
| 6 | Crack normal y-component (sin *theta*) |
| 7-8 | Maximum local crack strains encountered (normal, shear) |
| 9-10 | Previous increment local crack strains (normal, shear) |
| 11-12 | Current increment in local crack strains (normal, shear) |
| 13 | Local crack stress at initiation |
| 14 | Current local crack stress |
| 15 | Current Dcr(1,1) value |

#### Main Subroutine: `sca2diso()`

Executes the following steps each increment:

1. **Elastic stiffness computation** - Builds the 3x3 plane stress isotropic stiffness matrix *D_co*
2. **Strain update** - Accumulates the total strain: *eps = eps + d_eps*
3. **Bazant limit check** - Verifies element size satisfies *L_eff < 2 G_Ic E / sigma_cr0^2*; aborts if violated
4. **Pre-crack phase** - If no crack exists:
   - Computes elastic stress: *sigma = sigma + D_co * d_eps*
   - Evaluates principal stresses and directions using `SPRIND` (Abaqus utility)
   - Checks tensile and compressive failure indices against the strength *X_T*
   - If failure index >= 1.0, sets the crack flag and stores crack orientation
5. **Post-crack phase** - If a crack exists:
   - Retrieves crack history (max crack strain, old crack strain, crack orientation)
   - Computes the crack stiffness matrix *D_cr* via `calcDcr()`
   - Solves for the incremental local crack strain *d_ecr* via `calcEcr()`
   - Computes the degraded tangent stiffness *D_cocr* via `calcDcocr()`
   - Updates stress: *sigma = sigma + D_cocr * d_eps*

#### Helper Subroutines

**`STIFFNESS_MATRIX_2D_PS(E, nu, Dco)`**
Builds the 3x3 isotropic plane stress elasticity matrix:
```
         E        [ 1   v   0           ]
Dco = --------  * [ v   1   0           ]
      (1 - v^2)   [ 0   0   (1-v)/2    ]
```

**`getN(c, s, N)`**
Constructs the 3x2 transformation matrix *N* that maps local crack tractions (normal, shear) to global stress components, given the crack normal direction (cos *theta*, sin *theta*):
```
      [ c^2        -c*s       ]
N  =  [ s^2         c*s       ]
      [ 2*s*c    c^2 - s^2    ]
```

**`mysprind(S, PS, AN, LSTR, NDI, NSHR)`**
Wrapper around Abaqus's built-in `SPRIND` utility to compute principal stresses and their directions from the stress vector.

**`myExit()`**
Wrapper around Abaqus's `XIT` subroutine for graceful termination.

**`calcEcr(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N, dEcr, Ecr, GIc, sigcr0, Leff, nu, scr_old, scr)`**
Solves for the increment of local crack strain. The crack strain increment is found from:

```
d_ecr = [D_cr + N^T * D_co * N]^(-1) * N^T * D_co * d_eps
```

After solving, it updates the crack strain history:
- *Ecr = Ecr_OLD + dEcr*
- *Ecr_MAX = max(Ecr, Ecr_MAX)*

**`calcDcr(Ecr_OLD, Ecr_MAX, GIc, sigcr0, Leff, nu, Dcr, scr_old, scr, dEcr)`**
Computes the 2x2 local crack stiffness matrix using a **linear traction-separation law**:

- **Loading** (crack strain advancing beyond maximum): Secant stiffness from origin to current point on the softening curve
- **Unloading/Reloading** (crack strain below maximum): Tangent stiffness of the softening law
- **Complete failure** (crack strain exceeds critical value *2 G_Ic / (L_eff * sigma_cr0)*): Near-zero residual stiffness (1e-8)
- Shear crack stiffness is set proportional to the normal stiffness, scaled down by a factor of *1e-5 / (2(1+v))*

**`calcDcocr(Dco, Dcr, N, Dcocr)`**
Computes the 3x3 degraded (cracked) tangent stiffness using the standard SCA formula:

```
D_cocr = D_co - D_co * N * [D_cr + N^T * D_co * N]^(-1) * N^T * D_co
```

This represents the homogenized stiffness of the cracked material element.

---

### 3. `matrix_inverse.for` - Linear Algebra Utilities

**Purpose:** Provides matrix inversion and linear system solvers based on LU decomposition (Crout's method with partial pivoting). Adapted from *Numerical Recipes in Fortran*.

#### Subroutines

| Subroutine | Precision | Description |
|------------|-----------|-------------|
| `matrixInverse(a, y, n, np)` | REAL\*8 | Inverts matrix *a* into *y* using LU decomposition. Solves *A * X_i = e_i* for each column of the identity matrix. |
| `qmatrixInverse(a, y, n, np)` | REAL\*16 | Quad-precision version of `matrixInverse`. |
| `solveLinSysLU(A, b, n, flag)` | REAL\*8 | Solves *Ax = b* in-place (solution returned in *b*). Returns `flag=1` on success, `flag=0` if matrix is singular. |
| `qSolveLinSysLU(A, b, n, flag)` | REAL\*16 | Quad-precision version of `solveLinSysLU`. |
| `ludcmp(a, n, np, indx, d, flag)` | REAL\*8 | LU decomposition using Crout's method with partial pivoting. Decomposes *A = P * L * U*. |
| `qludcmp(a, n, np, indx, d, flag)` | REAL\*16 | Quad-precision version of `ludcmp`. |
| `lubksb(a, n, np, indx, b)` | REAL\*8 | Forward/back substitution to solve *L * U * x = P * b* given an LU-decomposed matrix. Optimized to skip leading zeros in *b*. |
| `qlubksb(a, n, np, indx, b)` | REAL\*16 | Quad-precision version of `lubksb`. |

**Key Implementation Details:**
- Maximum supported matrix size: **NMAX = 10** (hardcoded parameter in `ludcmp`/`qludcmp`)
- Near-zero pivot tolerance: **TINY = 1.0e-20** (substituted for exact zeros to handle near-singular matrices)
- The LU decomposition modifies the input matrix *a* in-place
- Partial pivoting is used for numerical stability

---

## Smart Changes

Potential improvements to make these codes more efficient, maintainable, and robust:

### 1. Replace General LU Inversion with Analytical 2x2 Inverse

The `calcEcr` and `calcDcocr` subroutines only ever invert **2x2 matrices**, yet they call the general-purpose `matrixInverse` routine which performs full LU decomposition with pivoting. A direct analytical formula would be significantly faster and simpler:

```fortran
! Direct 2x2 inverse: inv(A) = (1/det) * [d, -b; -c, a]
det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
invA(1,1) =  A(2,2) / det
invA(1,2) = -A(1,2) / det
invA(2,1) = -A(2,1) / det
invA(2,2) =  A(1,1) / det
```

This eliminates the overhead of LU decomposition, pivoting, and back-substitution for a trivially small system.

### 2. Cache the Elastic Stiffness Matrix

`STIFFNESS_MATRIX_2D_PS` is called every increment for every integration point, but the elastic constants (*E*, *v*) never change. The 3x3 elastic stiffness matrix `Dco` could be computed once and stored as a state variable or computed only on the first call, avoiding redundant work across potentially millions of increments.

### 3. Pre-compute Repeated Matrix Products

In `calcEcr`, the product `N^T * D_co * N` is computed, and the same product is computed again inside `calcDcocr`. Since `N`, `D_co` are unchanged between these calls, this product could be computed once in `sca2diso` and passed to both subroutines, eliminating one redundant matrix multiplication chain.

### 4. Remove Duplicated Quad-Precision Subroutines

The file `matrix_inverse.for` contains four quad-precision (`REAL*16`) duplicates (`qmatrixInverse`, `qSolveLinSysLU`, `qludcmp`, `qlubksb`) that are **never called** anywhere in the SCA2D code. All actual calls use the double-precision versions. These duplicates double the code maintenance burden for no benefit. They should either be removed or, if needed in the future, consolidated using Fortran generic interfaces or preprocessor macros.

### 5. Use Named Constants Instead of Magic Numbers

Several hardcoded values lack explanation:
- `sigcr0_r = sigcr0 * 5D-5` (line 386 of `formulation_AFEM2D.for`) - the residual stress fraction
- `Dcr(2,2) = Dcr(1,1) * term * 1.0D-5` (line 412) - the shear retention factor
- `1.0d-8` used repeatedly for initial crack strains and post-failure stiffness

These should be defined as named `PARAMETER` constants at the top of the file with comments explaining their physical meaning, making the code self-documenting and easier to tune.

### 6. Add Input Validation in the UMAT Wrapper

`main_SCA2D_inc.for` passes data directly to `sca2diso` with no validation. Adding checks for:
- `NTENS == 3` (confirming plane stress)
- `NPROPS >= 5` (sufficient material properties provided)
- `NSTATV >= 15` (sufficient state variables allocated)
- Material properties are physically valid (E > 0, 0 < nu < 0.5, GIc > 0, XT > 0)

would catch common user setup errors early with clear messages, rather than producing cryptic failures deep in the formulation.

### 7. Eliminate Redundant `calcDcr` Call

`calcDcr` is called once in Step 7 of `sca2diso` (line 221) and then called again inside `calcEcr` (line 346). The first call's results (`Dcr`, `scr`) are not used before `calcEcr` overwrites them. The call at line 221 can be removed entirely, and the `scr_old`/`scr` update at line 223 should be moved to after the `calcEcr` call.

### 8. Use Fortran 90+ Features for Clarity

The code mixes Fortran 77 fixed-format conventions with some Fortran 90 features. Fully adopting modern Fortran would improve readability:
- Use `MODULE` to encapsulate the SCA formulation and its helper routines, providing explicit interfaces and avoiding implicit typing issues
- Replace `DIMENSION` statements with explicit array declarations in the argument list
- Use `INTENT(IN/OUT/INOUT)` on all subroutine arguments to make data flow explicit and enable compiler optimization
- Use derived types (`TYPE`) for state variables instead of raw array indexing (e.g., `state%crack_flag` instead of `STATEV(1)`)

### 9. Improve the Bazant Limit Check Behavior

Currently, violating the Bazant size-effect limit causes an immediate hard abort via `XIT`. A more user-friendly approach would be to:
- Use Abaqus's `PNEWDT` mechanism to request a smaller time increment before aborting
- Write the diagnostic to the Abaqus `.msg` file (unit 7) rather than stdout, so it appears in the standard Abaqus log
- Include the element number (`NOEL`) in the error message to help the user identify which elements need refinement

### 10. Raise `NMAX` or Make It Dynamic in LU Decomposition

The `ludcmp` and `qludcmp` routines have `NMAX=10` hardcoded. While the current SCA2D code only inverts 2x2 matrices, this library is a general utility. If reused for larger systems (e.g., extending to 3D or multi-crack models), it would silently produce wrong results for matrices larger than 10x10. Using allocatable arrays or at minimum raising `NMAX` with a bounds check would prevent this silent failure.
