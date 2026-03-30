# claude_AFEM2D_PEv1.for — SCA 2D Plane Strain UMAT (Progressive Damage)

## Overview

`claude_AFEM2D_PEv1.for` implements the **Smeared Crack Approach (SCA)** for isotropic
linear-elastic/softening materials under **2D plane strain** conditions within the ABAQUS
UMAT framework. It is designed for progressive damage analysis, including DCB (Double
Cantilever Beam) delamination modeling.

The model tracks **one fixed crack** per integration point and covers:

- **Pre-peak**: Linear elastic response until failure criterion is met
- **Initiation overshoot prevention**: PNEWDT-based step control at crack onset
- **Post-peak softening**: Linear cohesive (Bazant crack-band) with total-form traction
- **Unloading / Reloading**: Secant-based crack stiffness
- **Complete failure**: Residual near-zero stiffness after full crack opening
- **Viscous regularization**: Optional `eta/DT` term for convergence during softening

Called from `main_SCA2D_inc.for`:

```fortran
include 'claude_AFEM2D_PEv1.for'
include 'matrix_inverse.for'
```

---

## Element Type

**Plane Strain** (`CPE4` elements in ABAQUS):

| Variable | Value | Meaning                                      |
|----------|-------|----------------------------------------------|
| `NTENS`  | 4     | sig_xx, sig_yy, sig_zz, tau_xy               |
| `NDI`    | 3     | Three normal stress/strain components        |
| `NSHR`   | 1     | One engineering shear component              |

---

## Material Properties (`PROPS`)

| Index | Symbol | Description                        | Units    | Required |
|-------|--------|------------------------------------|----------|----------|
| 1     | E      | Young's modulus                    | Pa       | Yes      |
| 2     | nu     | Poisson's ratio                    | --       | Yes      |
| 3     | GIc    | Mode I fracture toughness          | J/m^2    | Yes      |
| 4     | XT     | Tensile (& compressive) strength   | Pa       | Yes      |
| 5     | INDEX  | Failure criteria selector          | --       | Yes      |
| 6     | eta    | Viscous regularization parameter   | s        | Optional |

**eta guidance**: Start with `eta = 1e-5`. Increase to `1e-4` or `1e-3` if convergence
is poor. Decrease if results deviate from reference. Set to 0 (or omit PROPS(6))
to disable viscous regularization.

ABAQUS input example:
```
*USER MATERIAL, CONSTANTS=6
 E, nu, GIc, XT, INDEX, eta
```

---

## State Variables (`STATEV`)

| Index | Symbol             | Description                                        |
|-------|--------------------|----------------------------------------------------|
| 1     | isCRACKED          | Crack flag: 0 = intact, 1 = cracked               |
| 2     | MODE               | Failure mode: +1 = tensile, -1 = compressive       |
| 3     | Leff               | Effective element size (crack band width)          |
| 4     | FI_MAX             | Maximum failure index attained                     |
| 5     | cos(theta)         | Cosine of crack-normal angle                       |
| 6     | sin(theta)         | Sine of crack-normal angle                         |
| 7     | Ecr_MAX(1)         | Maximum normal crack strain ever reached           |
| 8     | Ecr_MAX(2)         | Maximum shear crack strain ever reached            |
| 9     | Ecr_OLD(1)         | Normal crack strain at end of current increment    |
| 10    | Ecr_OLD(2)         | Shear crack strain at end of current increment     |
| 11    | dEcr(1)            | Normal crack strain increment (this step)          |
| 12    | dEcr(2)            | Shear crack strain increment (this step)           |
| 13    | sigcr0 * FI_MAX    | Peak local crack normal stress stored              |
| 14    | scr(1)             | Normal crack traction (consistent with STATEV(9))  |
| 15    | scr(2)             | Shear crack traction (consistent with STATEV(10))  |
| 16    | sigcr0             | Peak local normal crack stress at initiation       |
| 17    | taucr0             | Peak local shear crack stress at initiation        |

**Key consistency guarantee**: STATEV(14) and STATEV(9) are always evaluated at the
same point. The FIX-GAP second `calcDcr_PE` call uses Ecr_new (= STATEV(9:10)) to
recompute scr (= STATEV(14:15)), so plotting STATEV(14) vs STATEV(9) produces points
that lie exactly on the softening curve (or secant line during unloading).

---

## Mathematical Formulation

### 1. Elastic Stiffness Matrix (Plane Strain)

For isotropic material under plane strain:

```
         E
D_co = ----------- * [ 1-nu  nu    nu    0            ]
       (1+nu)(1-2nu)  [ nu    1-nu  nu    0            ]
                       [ nu    nu    1-nu  0            ]
                       [ 0     0     0     (1-2nu)/2    ]
```

Incremental stress update: `d{sig} = D_co * d{eps}`

---

### 2. Failure Criterion -- Maximum Principal Stress

Principal stresses (sig_1, sig_2, sig_3) obtained via ABAQUS `SPRIND`.

```
FI_T = max(sig_i) / XT        (tensile)
FI_C = |min(sig_i)| / XT      (compressive)
FI_MAX = max(FI_T, FI_C)
```

Crack initiation: `FI_MAX >= 1.0`

The crack normal is perpendicular to the controlling principal stress direction.

---

### 3. Initiation Overshoot Prevention

**Problem**: The elastic predictor `SIG = SIG_prev + D_co * DSTRAN` can overshoot
the failure surface significantly in a single increment, producing `FI_MAX >> 1.0`.
This creates a large initial crack strain jump and an inconsistent starting point
for softening.

**Solution**: If `FI_MAX > 1.05` (overshoot > 5%):

1. Compute FI at the **start** of the increment from `SIG_prev` (saved before update):
   ```
   FI_prev = max_principal(SIG_prev) / XT
   ```

2. Linear interpolation for the fraction `alpha` at which FI = 1.0:
   ```
   alpha = (1.0 - FI_prev) / (FI_MAX - FI_prev)
   alpha = clamp(alpha, 0.10, 0.95)
   ```

3. Set `PNEWDT = alpha` and RETURN. ABAQUS abandons the increment and retries
   with `DT_new = alpha * DT`, placing the stress within 5% of XT at initiation.

**Result**: The crack always initiates with stress approximately equal to XT,
and the first post-peak increment starts from a consistent state at the top
of the softening curve.

---

### 4. Crack Kinematics -- N Matrix

Crack unit normal: `n = [cos(theta), sin(theta)]^T`

The N matrix (4x2) maps local crack strains to global Voigt strains:

```
     [ c^2      -cs       ]
N =  [ s^2       cs       ]
     [ 0         0        ]
     [ 2sc       c^2-s^2  ]
```

---

### 5. Crack Softening Law -- Linear Cohesive (Bazant Crack Band)

The normal crack traction follows a linear softening law:

```
scr(Ecr_n) = sigcr0 + slope_soft * Ecr_n

where: slope_soft = -sigcr0^2 / (2 * GIc / Leff)
```

**Ultimate crack strain** (zero traction):

```
Ecr_n_ult = 2 * GIc / (Leff * sigcr0)
```

**Fracture energy equivalence**:

```
0.5 * sigcr0 * Ecr_n_ult * Leff = GIc
```

**Bazant size limit** (prevents snap-back):

```
Leff <= 2 * GIc * E / sigcr0^2
```

---

### 6. Crack Constitutive Matrix D_cr (2x2)

The normal crack stiffness has three regimes:

#### (a) Softening branch: `Ecr_OLD >= Ecr_MAX` and `Ecr < Ecr_n_ult`

Crack is on the damage envelope, propagating forward.

```
Dcr(1,1) = slope_soft + eta/DT          (viscous-regularized slope)

scr(1) = sigcr0 + slope_soft * Ecr_OLD  (TOTAL form, analytical)
       + (eta/DT) * dEcr                (viscous rate term)
scr(1) = max(scr(1), sigcr0_r)          (clamp to residual)
```

**Why total form**: The traction `scr = sigcr0 + slope_soft * Ecr` lies exactly on
the analytical softening line regardless of accumulated numerical errors from previous
increments. This is the key to eliminating the one-increment lag between STATEV(9) and
STATEV(14). The incremental form `scr = scr_old + slope*dEcr` drifts over many
increments and causes stresses to not reach zero at full failure.

#### (b) Secant branch: `Ecr_OLD < Ecr_MAX` and `Ecr < Ecr_n_ult`

Crack is unloading or reloading within the damage surface.

```
scr_max = sigcr0 + slope_soft * Ecr_MAX    (traction at envelope)
scr_max = max(scr_max, 0)

Dcr(1,1) = scr_max / Ecr_MAX               (secant slope)
scr(1)   = Dcr(1,1) * Ecr_OLD              (total form)
```

Special case: if `Ecr_MAX ~ 0` (just initiated), use softening slope instead
to avoid division by near-zero.

#### (c) Residual branch: `Ecr >= Ecr_n_ult`

Element is fully failed.

```
Dcr(1,1) = RES_K = 1e-8    (small residual stiffness)
scr(1)   = sigcr0_r         (near-zero residual traction)
```

#### Shear component

```
Dcr(2,2) = Dcr(1,1) * 1/(2(1+nu)) * 1e-3    (pre-failure)
Dcr(2,2) = RES_K / (2(1+nu))                 (post-failure)
```

---

### 7. Viscous Regularization

Viscous regularization adds a rate-dependent term to the crack stiffness:

```
Dcr_eff = slope_soft + eta / DT
```

This makes the effective tangent less negative during softening, improving
Newton-Raphson convergence. The viscous contribution to traction:

```
scr_visc = (eta / DT) * dEcr
```

This mirrors the viscous regularization in ABAQUS cohesive elements.
When `eta = 0` or `DT -> infinity`, the term vanishes and the formulation
reduces to standard rate-independent SCA.

The viscous term only applies to the **softening branch**. The secant and
residual branches are rate-independent.

---

### 8. Crack Strain Increment

Stress continuity at the crack face:

```
N^T * D_co * (d_eps - N * d_Ecr) = D_cr * d_Ecr
```

Solving for `d_Ecr`:

```
Y = D_cr + N^T * D_co * N       (2x2 matrix)

d_Ecr = Y^{-1} * N^T * D_co * DSTRAN
```

Y is inverted via LU decomposition. If Y is near-singular, `PNEWDT = 0.5`
is requested and the increment is abandoned.

---

### 9. Failure Overshoot Detection (Post-Peak)

If crack strain would jump past `Ecr_n_ult` by more than 5%:

```
IF (Ecr_OLD + dEcr > 1.05 * Ecr_n_ult):
    alpha = (Ecr_n_ult - Ecr_OLD) / dEcr
    PNEWDT = min(PNEWDT, alpha)
```

This ensures the transition from softening to residual is gradual. Small
overshoots (<5%) are absorbed by the residual branch which clamps `scr`
to `sigcr0_r`.

**Secant-to-softening overshoot** (Ecr crosses Ecr_MAX during reloading)
does NOT trigger a step cut. The total-form traction and the FIX-GAP second
`calcDcr_PE` call handle this transition naturally.

---

### 10. Consistent Tangent Operator

The cracked-continuum tangent for ABAQUS Newton-Raphson:

```
D_cocr = D_co - D_co * N * Y^{-1} * N^T * D_co
```

Returned in `DDSDDE`. Reduces to `D_co` when no crack; approaches
rank-deficient as `D_cr -> 0`.

---

### 11. Interpenetration Prevention (OS-4)

For tensile-mode cracks, the normal crack strain must not go negative
(crack faces cannot interpenetrate):

```
IF (MODE > 0 AND Ecr(1) < 0):
    dEcr(1)    = 0 - STATEV(9)     (true previous Ecr)
    Ecr(1)     = 0
    Ecr_OLD(1) = 0
```

Uses `STATEV(9)` (not the local `Ecr_OLD`, which has already been updated
at this point) to compute the correct `dEcr`.

---

### 12. FIX-GAP: Ecr/scr Consistency

After computing `dEcr` and updating `Ecr_OLD := Ecr_new`, the code calls
`calcDcr_PE` a **second time** with the updated `Ecr_OLD`:

```
Step 7:   calcDcr_PE(Ecr_OLD_old, ...) -> scr at Ecr_OLD_old   [first call]
Step 8:   calcEcr_PE(...)              -> dEcr, Ecr_new
          STATEV(9:10) = Ecr_new
FIX-GAP:  calcDcr_PE(Ecr_new, ...)    -> scr at Ecr_new        [second call]
          STATEV(14:15) = scr
```

This ensures that:
- **Softening branch**: `scr(1) = sigcr0 + slope_soft * Ecr_new` (total form
  evaluated at the new Ecr)
- **Secant branch**: `scr(1) = D_sec * Ecr_new`

So STATEV(14) and STATEV(9) are always consistent -- they both refer to the
**same point** on the softening curve or secant line.

### 12a. Zero Seeds at Initiation

At crack initiation, all crack strains are initialized to **zero** (not 1e-8):

```
STATEV(7:12)  = 0.0      (Ecr_MAX, Ecr_OLD, dEcr all zero)
STATEV(14:15) = sigcr0    (traction at the top of softening curve)
```

This places the initial state exactly at `(Ecr=0, scr=sigcr0)` -- the origin
of the softening curve. The old approach of seeding with `TINY_ = 1e-8` created
a small but visible offset between STATEV(9) and the point where the softening
curve begins.

The secant branch handles `Ecr_MAX = 0` safely: if `Ecr_MAX < ECR_MIN (1e-14)`,
it falls through to the softening slope instead of dividing by near-zero.

---

## Algorithm / Pseudo-code

```
SUBROUTINE sca2diso_PE(NTENS, ..., PROPS, EPS, DSTRAN, SIG,
                       STATEV, DDSDDE, ...)

  --- INITIALISATION ---
  READ  E, nu, GIc, XT, INDEX       from PROPS(1:5)
  READ  eta                          from PROPS(6) [optional, default 0]
  SET   Leff = CELENT                (element characteristic length)
  READ  sigcr0, taucr0               from STATEV(16:17)

  STEP 1: Compute D_co (4x4 plane-strain elastic stiffness)
  STEP 2: Update total strain: EPS += DSTRAN

  STEP 3: Bazant size check
          IF Leff > 2*GIc*E/sigcr0^2 -> ABORT

  --- PRE-PEAK BLOCK [STATEV(1) < 0.5] ---

  Save SIG_prev = SIG                  (for overshoot detection)
  SIG += D_co * DSTRAN                 (elastic predictor)
  DDSDDE = D_co

  Compute principal stresses -> FI_MAX
  Determine failure mode and crack orientation

  IF (FI_MAX >= 1.0):
      IF (FI_MAX > 1.05):
          Compute FI_prev from SIG_prev
          alpha = (1.0 - FI_prev) / (FI_MAX - FI_prev)
          PNEWDT = alpha
          RETURN                          (ABAQUS retries with smaller step)

      CRACK ACTIVATION:
          STATEV(1) = 1
          Compute N matrix
          sigcr0 = N^T * SIG             (local stress at initiation)
          STATEV(14:15) = [sigcr0, taucr0]
          STATEV(7:12)  = 0.0             (zero crack strain seeds)

  --- POST-PEAK BLOCK [STATEV(1) >= 1] ---

  STEP 6: Retrieve history: Ecr_MAX, Ecr_OLD, dEcr, scr_old
  STEP 7: calcDcr_PE(Ecr_OLD, ...) -> D_cr, scr     [first call]
  STEP 8: calcEcr_PE(...) -> dEcr, Ecr
          Interpenetration clamp (Ecr(1) >= 0)
          Update STATEV(7:12)

  FIX-GAP: calcDcr_PE(Ecr_new, ...) -> D_cr, scr    [second call]
           STATEV(14:15) = scr                        [consistent with STATEV(9)]

  STEP 9: calcDcocr_PE(D_co, D_cr, N) -> D_cocr
          SIG += D_cocr * DSTRAN
          DDSDDE = D_cocr

  RETURN
END
```

---

## PNEWDT Strategy (Step Size Control)

The formulation uses a **minimal-intervention** approach to step size control.
ABAQUS's automatic time stepping handles step growth; the UMAT only cuts
when strictly necessary.

| Condition | Action | Rationale |
|-----------|--------|-----------|
| FI_MAX > 1.05 at initiation | `PNEWDT = alpha` | Prevent stress overshoot at crack onset |
| Y matrix singular in calcEcr_PE | `PNEWDT = 0.5` | dEcr would be garbage; must retry |
| Ecr jumps past 1.05*Ecr_n_ult | `PNEWDT = alpha` | Prevent massive overshoot at full failure |
| Secant-to-softening transition | No cut | Total-form traction + FIX-GAP handles it |
| calcDcocr_PE near-singular | No cut (user-commented) | Elastic fallback as tangent; NR iterates |

**Key design principle**: Every abandoned increment wastes all iteration work.
Accepting a slightly imperfect increment and letting Newton-Raphson correct
in the next iteration is far cheaper than restarting.

---

## Subroutine Reference

| Subroutine              | Role                                                     |
|-------------------------|----------------------------------------------------------|
| `sca2diso_PE`           | Main SCA driver; coordinates all steps                   |
| `STIFFNESS_MATRIX_2D_PE`| Assembles 4x4 plane-strain isotropic D_co               |
| `getN_PE`               | Builds 4x2 crack transformation matrix N                 |
| `mysprind_PE`           | Wrapper around ABAQUS SPRIND                             |
| `myExit_PE`             | Wrapper around ABAQUS XIT                                |
| `calcEcr_PE`            | Computes crack strain increment via stress continuity    |
| `calcDcr_PE`            | Computes crack constitutive matrix D_cr + traction scr   |
| `calcDcocr_PE`          | Assembles cracked-continuum tangent D_cocr               |
| `matrixInverse`         | LU-based 2x2 matrix inverse (in matrix_inverse.for)     |

---

## Change Log

### v1 (03/03/2026) -- Overshoot Fixes
- OS-1: Total-form scr in secant branch
- OS-2: Reloading overshoot detection (later removed, see DCB fixes)
- OS-3: Removed (caused infinite step-cut loop)
- OS-4: Interpenetration clamp using STATEV(9) for correct dEcr
- FIX-GAP: Second calcDcr_PE call for Ecr/scr consistency

### DCB Fixes (03/10/2026)
- Added PROPS(6) = eta for viscous regularization
- Removed PNEWDT=1.0 overrides that defeated step-cutting
- Fixed OS-2new failure overshoot: checked Ecr_n_ult (not Ecr_MAX)
- Removed PNEWDT cut from calcDcocr_PE singularity (elastic fallback)
- Removed OS-2 step cut entirely (total-form + FIX-GAP handles it)

### Progressive Damage Fixes (03/11/2026)
- **Initiation overshoot prevention**: PNEWDT cut when FI_MAX > 1.05,
  using linear interpolation from SIG_prev to find alpha at FI = 1.0
- **Zero seeds**: Crack strains initialized to 0.0 (not 1e-8) at initiation,
  eliminating the offset between STATEV(9) and STATEV(14) at crack onset
- **Total-form softening traction**: `scr = sigcr0 + slope_soft * Ecr_OLD`
  replaces incremental `scr = scr_old + slope*dEcr`. Guarantees scr is always
  on the analytical softening curve, reaching exactly zero at Ecr_n_ult
- **Safe secant at Ecr_MAX~0**: Falls through to softening slope when
  Ecr_MAX < 1e-14, avoiding division by near-zero at first post-peak increment

---

## References

- Bazant, Z.P. & Oh, B.H. (1983). Crack band theory for fracture of concrete.
  *Materiaux et Construction*, **16**(3), 155-177.
- de Borst, R., Crisfield, M.A., Remmers, J.J.C., & Verhoosel, C.V. (2012).
  *Nonlinear Finite Element Analysis of Solids and Structures* (2nd ed.). Wiley.
- Rots, J.G. (1988). *Computational Modeling of Concrete Fracture*. PhD Thesis, TU Delft.
- ABAQUS Theory Guide -- Section on User-Defined Material Behaviors (UMAT).
