# formulation_AFEM2D_PE.for — SCA 2D Plane Strain UMAT

## Overview

`formulation_AFEM2D_PE.for` implements the **Smeared Crack Approach (SCA)** for isotropic
linear-elastic/softening materials under **2D plane strain** conditions within the ABAQUS UMAT
framework. It is `#include`d by `main_SCA2D_inc.for` and provides the subroutine `sca2diso`
plus all required helper subroutines.

The model tracks **one fixed crack** per integration point and covers:

- **Pre-peak**: Linear elastic response until a failure criterion is met
- **Post-peak**: Crack propagation with linear cohesive (Bazant crack-band) softening
- **Unloading / Reloading**: Secant-based crack stiffness on unloading
- **Complete failure**: Residual near-zero stiffness after full crack opening

Called from `main_SCA2D_inc.for` via:

```fortran
include 'formulation_AFEM2D_PE.for'
```

---

## Element Type

**Plane Strain** (`CPE` elements in ABAQUS):

| Variable | Value | Meaning                                      |
|----------|-------|----------------------------------------------|
| `NTENS`  | 4     | σ\_xx, σ\_yy, σ\_zz, τ\_xy                  |
| `NDI`    | 3     | Three normal stress/strain components        |
| `NSHR`   | 1     | One engineering shear component              |

---

## Material Properties (`PROPS`)

| Index | Symbol   | Description                        | Units    |
|-------|----------|------------------------------------|----------|
| 1     | E        | Young's modulus                    | MPa      |
| 2     | ν        | Poisson's ratio                    | —        |
| 3     | G\_Ic    | Mode I fracture toughness          | N/mm     |
| 4     | X\_T     | Tensile (& compressive) strength   | MPa      |
| 5     | INDEX    | Failure criteria selector (0 or 1) | —        |

---

## State Variables (`STATEV`)

| Index | Symbol                   | Description                                      |
|-------|--------------------------|--------------------------------------------------|
| 1     | isCRACKED                | Crack flag: 0 = intact, 1 = cracked              |
| 2     | MODE                     | Failure mode: +1 = tensile, −1 = compressive     |
| 3     | L                        | Element size (crack band width)                  |
| 4     | FI\_MAX                  | Maximum failure index attained                   |
| 5     | cos θ                    | Cosine of crack-normal angle                     |
| 6     | sin θ                    | Sine of crack-normal angle                       |
| 7     | ε\_cr,n^MAX              | Maximum normal local crack strain                |
| 8     | ε\_cr,t^MAX              | Maximum shear local crack strain                 |
| 9     | ε\_cr,n^OLD              | Normal local crack strain (previous increment)   |
| 10    | ε\_cr,t^OLD              | Shear local crack strain (previous increment)    |
| 11    | Δε\_cr,n                 | Normal crack strain increment (current step)     |
| 12    | Δε\_cr,t                 | Shear crack strain increment (current step)      |
| 13    | σ\_cr0 × FI\_MAX         | Peak local crack normal stress stored            |
| 14    | σ\_cr^OLD                | Local normal crack stress (previous increment)   |
| 15    | τ\_cr^OLD                | Local shear crack stress (previous increment)    |
| 16    | σ\_cr0                   | Peak local normal crack stress at initiation     |
| 17    | τ\_cr0                   | Peak local shear crack stress at initiation      |

---

## Mathematical Formulation

### 1. Elastic Stiffness Matrix (Plane Strain)

For an isotropic material under plane strain the 4×4 stiffness matrix is:

$$
\mathbf{D}_{co} = \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
  1-\nu & \nu   & \nu   & 0 \\
  \nu   & 1-\nu & \nu   & 0 \\
  \nu   & \nu   & 1-\nu & 0 \\
  0     & 0     & 0     & \dfrac{1-2\nu}{2}
\end{bmatrix}
$$

where the stress and strain vectors follow Voigt notation:

$$\{\sigma\} = [\sigma_{xx},\;\sigma_{yy},\;\sigma_{zz},\;\tau_{xy}]^T, \qquad
\{\varepsilon\} = [\varepsilon_{xx},\;\varepsilon_{yy},\;\varepsilon_{zz},\;\gamma_{xy}]^T$$

The incremental linear-elastic stress update is:

$$\Delta\{\sigma\} = \mathbf{D}_{co}\,\Delta\{\varepsilon\}$$

---

### 2. Failure Criterion — Maximum Principal Stress

Principal stresses $\sigma_1 \geq \sigma_2 \geq \sigma_3$ (and their eigenvectors) are
obtained from the ABAQUS built-in routine `SPRIND`.

**Tensile failure index:**

$$\text{FI}_T = \frac{\sigma_{\max}}{X_T}, \qquad \sigma_{\max} = \max(\sigma_1,\,\sigma_2,\,\sigma_3)$$

**Compressive failure index:**

$$\text{FI}_C = \frac{|\sigma_{\min}|}{X_T}, \qquad \sigma_{\min} = \min(\sigma_1,\,\sigma_2,\,\sigma_3)$$

**Overall failure index and dominant mode:**

$$\text{FI}_{\max} = \max(\text{FI}_T,\;\text{FI}_C)$$

- If $\text{FI}_T \geq \text{FI}_C$: tensile mode (MODE = +1), crack normal ∥ direction of σ\_max
- Otherwise: compressive mode (MODE = −1), crack normal ∥ direction of σ\_min

**Crack initiation condition:**

$$\text{FI}_{\max} \geq 0.99$$

The crack plane is **perpendicular** to the controlling principal stress direction.
The crack angle θ is stored as (cos θ, sin θ) = first two components of the controlling
eigenvector.

---

### 3. Crack Kinematics — N Matrix

Let the crack unit normal be $\mathbf{n} = [c,\;s]^T$ with $c = \cos\theta$, $s = \sin\theta$.

The **N matrix** (4×2) maps local crack strain increments
$$\Delta\boldsymbol{\varepsilon}_{cr} = [\Delta\varepsilon_{cr,n},\;\Delta\gamma_{cr,t}]^T$$
to the corresponding contribution in global Voigt strain:

$$\mathbf{N} =
\begin{bmatrix}
  c^2  & -cs       \\
  s^2  & \phantom{-}cs   \\
  0    & 0         \\
  2sc  & c^2 - s^2
\end{bmatrix} \in \mathbb{R}^{4\times 2}$$

Its transpose $\mathbf{N}^T$ recovers local crack strains from global strains:

$$\Delta\boldsymbol{\varepsilon}_{cr} = \mathbf{N}^T\,\Delta\boldsymbol{\varepsilon}_{global}$$

**Verification:** for a horizontal crack ($\theta = 0$, $c=1$, $s=0$):

$$\mathbf{N} = \begin{bmatrix}1&0\\0&0\\0&0\\0&1\end{bmatrix},
\quad \mathbf{N}^T\{\varepsilon\} = [\varepsilon_{xx},\;\gamma_{xy}]^T \;\checkmark$$

---

### 4. Crack Softening Law — Linear Cohesive (Bazant Crack Band)

The normal crack traction evolves with crack opening strain following a **linear softening law**:

$$\sigma_{cr}(\varepsilon_{cr,n}) = \sigma_{cr0} - \frac{\sigma_{cr0}^{\,2}}{2\,G_{Ic}/L_{eff}}\,\varepsilon_{cr,n}$$

where:

| Symbol         | Meaning                                          |
|----------------|--------------------------------------------------|
| $\sigma_{cr0}$ | Peak crack traction (= normal stress at failure) |
| $G_{Ic}$       | Mode I fracture toughness [J/m²]                 |
| $L_{eff}$      | Effective element size (crack band width)        |

**Ultimate crack strain** (zero traction):

$$\varepsilon_{cr,n}^{u} = \frac{2\,G_{Ic}}{L_{eff}\;\sigma_{cr0}}$$

**Fracture energy equivalence** (area under the traction–separation curve × element size = G\_Ic):

$$\frac{1}{2}\,\sigma_{cr0}\,\varepsilon_{cr,n}^{u}\,L_{eff} = G_{Ic} \;\checkmark$$

---

### 5. Bazant Size Limit (Mesh Objectivity)

To prevent snap-back (negative energy dissipation) and ensure mesh-size-independent fracture
energy, the element size must satisfy:

$$\boxed{L_{eff} \leq \frac{2\,G_{Ic}\,E}{\sigma_{cr0}^{\,2}}}$$

If this condition is violated, the analysis aborts with a diagnostic message recommending a
smaller element size.

---

### 6. Crack Constitutive Matrix D\_cr

The 2×2 crack stiffness matrix $\mathbf{D}_{cr}$ (in local crack coordinates) governs the
traction–crack-strain relationship: $\Delta\boldsymbol{\sigma}_{cr} = \mathbf{D}_{cr}\,\Delta\boldsymbol{\varepsilon}_{cr}$.

The normal and shear components are **decoupled** ($D_{cr,12} = D_{cr,21} = 0$).

**Normal component** — three regimes:

| Condition | Regime | $D_{cr,nn}$ |
|-----------|--------|-------------|
| $\varepsilon_{cr,n}^{temp} < \varepsilon_{cr,n}^{u}$ AND $\varepsilon_{cr,n}^{OLD} < \varepsilon_{cr,n}^{MAX}$ | Loading (secant) | $\dfrac{\sigma_{cr}^{MAX}}{\|\varepsilon_{cr,n}^{MAX}\|}$ |
| $\varepsilon_{cr,n}^{temp} < \varepsilon_{cr,n}^{u}$ AND $\varepsilon_{cr,n}^{OLD} \geq \varepsilon_{cr,n}^{MAX}$ | Unloading/Reloading | $-\dfrac{(\sigma_{cr0,r} - \sigma_{cr0})^2}{2\,G_{Ic}/L_{eff}} \approx -\dfrac{\sigma_{cr0}^2}{2\,G_{Ic}/L_{eff}}$ |
| $\varepsilon_{cr,n}^{temp} \geq \varepsilon_{cr,n}^{u}$ | Full failure | $10^{-5}$ |

where:

$$\sigma_{cr}^{MAX} = \sigma_{cr0} - \frac{\sigma_{cr0}^2}{2\,G_{Ic}/L_{eff}}\;\varepsilon_{cr,n}^{MAX}
\quad \text{(traction at current crack strain maximum)}$$

$$\sigma_{cr0,r} = 5\times10^{-5}\;\sigma_{cr0}
\quad \text{(residual traction for numerical stability)}$$

$$\varepsilon_{cr,n}^{temp} = \max\!\left(|\varepsilon_{cr,n}^{OLD}|,\;|\varepsilon_{cr,n}^{MAX}|\right)$$

**Shear component:**

$$D_{cr,tt} = D_{cr,nn} \times \frac{1}{2(1+\nu)} \times 10^{-8} \approx 0$$

(Shear crack stiffness is set to nearly zero — no frictional sliding resistance.)

---

### 7. Local Crack Strain Increment

**Stress continuity** at the crack face (traction continuity):

$$\mathbf{N}^T\,\mathbf{D}_{co}\bigl(\Delta\boldsymbol{\varepsilon} - \mathbf{N}\,\Delta\boldsymbol{\varepsilon}_{cr}\bigr)
= \mathbf{D}_{cr}\,\Delta\boldsymbol{\varepsilon}_{cr}$$

Rearranging to solve for $\Delta\boldsymbol{\varepsilon}_{cr}$:

$$\left(\mathbf{D}_{cr} + \mathbf{N}^T\mathbf{D}_{co}\mathbf{N}\right)\Delta\boldsymbol{\varepsilon}_{cr}
= \mathbf{N}^T\mathbf{D}_{co}\,\Delta\boldsymbol{\varepsilon}$$

$$\boxed{\Delta\boldsymbol{\varepsilon}_{cr} =
\underbrace{\left(\mathbf{D}_{cr} + \mathbf{N}^T\mathbf{D}_{co}\mathbf{N}\right)^{-1}}_{\mathbf{Y}^{-1}}
\mathbf{N}^T\mathbf{D}_{co}\,\Delta\boldsymbol{\varepsilon}}$$

where $\mathbf{Y} = \mathbf{D}_{cr} + \mathbf{N}^T\mathbf{D}_{co}\mathbf{N} \in \mathbb{R}^{2\times2}$
is inverted via LU decomposition (`matrixInverse` from `matrix_inverse.for`).

---

### 8. Consistent Tangent Operator (Algorithmic Jacobian)

Substituting the expression for $\Delta\boldsymbol{\varepsilon}_{cr}$ back into the stress
increment gives the **cracked-continuum consistent tangent**:

$$\Delta\boldsymbol{\sigma} = \mathbf{D}_{cocr}\,\Delta\boldsymbol{\varepsilon}$$

$$\boxed{\mathbf{D}_{cocr} = \mathbf{D}_{co}
- \mathbf{D}_{co}\,\mathbf{N}
\underbrace{\left(\mathbf{D}_{cr} + \mathbf{N}^T\mathbf{D}_{co}\mathbf{N}\right)^{-1}}_{\mathbf{Y}^{-1}}
\mathbf{N}^T\mathbf{D}_{co}}$$

This is returned in `DDSDDE` for the ABAQUS Newton–Raphson global solver.
It reduces to $\mathbf{D}_{co}$ when there is no crack and approaches a rank-deficient matrix
as the crack stiffness $\mathbf{D}_{cr} \to \mathbf{0}$.

---

## Algorithm / Pseudo-code

```
SUBROUTINE sca2diso(NTENS, NDI, NSHR, ..., PROPS, EPS, DSTRAN, SIG,
                    STATEV, DDSDDE, ...)

  ─── INITIALISATION ──────────────────────────────────────────────────
  READ  E, ν, GIc, XT, INDEX  from PROPS(1:5)
  READ  MODE, L               from STATEV(2), STATEV(3)
  SET   Leff ← L              (1st-order element: Leff = L)
  READ  sigcr0 ← STATEV(16)   (peak crack stress at initiation)
  READ  taucr0 ← STATEV(17)

  STEP 1: Compute D_co (4×4 plane-strain elastic stiffness)
          via STIFFNESS_MATRIX_2D(E, ν, D_co)

  STEP 2: Update total strain tensor
          EPS(1:NTENS) ← EPS(1:NTENS) + DSTRAN(1:NTENS)

  STEP 3: Bazant size check [only if sigcr0 > 0]
          BL ← 2·GIc·E / sigcr0²
          IF (Leff > BL):
            PRINT warning with recommended element size
            CALL XIT  [abort analysis]

  ─── PRE-PEAK BLOCK  [STATEV(1) < 0.5  →  no crack yet] ─────────────

  STEP 4: Linear elastic stress update
          SIG   ← SIG + D_co · DSTRAN
          DDSDDE ← D_co

  STEP 5: Compute principal stresses and directions
          CALL SPRIND(SIG, {σ₁,σ₂,σ₃}, T, ...)

          FI_T ← max(σ_i) / XT           [tensile  failure index]
          FI_C ← |min(σ_i)| / XT         [compressive failure index]
          FI_MAX ← max(FI_T, FI_C)

          IF (FI_T ≥ FI_C):
            MODE ← +1;  J_FAIL ← index of max tensile principal stress
          ELSE:
            MODE ← −1;  J_FAIL ← index of max compressive principal stress

          STATEV(2)   ← MODE
          STATEV(4)   ← FI_MAX
          STATEV(5:6) ← T(1:2, J_FAIL)   [cos θ, sin θ]

  STEP 5a: IF (FI_MAX ≥ 0.99):  CRACK INITIATION
            STATEV(1)   ← 1.0             [set crack flag]
            N           ← getN(cos θ, sin θ)
            sigcr0      ← (N^T · SIG)(1)  [local normal stress at failure]
            taucr0      ← (N^T · SIG)(2)  [local shear  stress at failure]
            STATEV(16)  ← sigcr0
            STATEV(17)  ← taucr0
            STATEV(13)  ← FI_MAX · XT
            STATEV(14:15) ← [sigcr0, taucr0]  [scr_old]
            STATEV(7:12)  ← 1×10⁻⁸           [seed crack strains]

  ─── POST-PEAK BLOCK  [STATEV(1) ≥ 1  →  crack active] ──────────────

  STEP 6: Retrieve history variables
          Ecr_MAX ← STATEV(7:8)
          Ecr_OLD ← STATEV(9:10)
          dEcr    ← STATEV(11:12)
          scr_old ← STATEV(14:15)
          N       ← getN(STATEV(5), STATEV(6))

  STEP 7: Compute D_cr (2×2) and predicted crack stress
          CALL calcDcr(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff, ν,
                       D_cr, scr_old, scr, dEcr)
          STATEV(14:15) ← scr

  STEP 8: Compute crack strain increment Δε_cr
          CALL calcEcr(Ecr_OLD, Ecr_MAX, DSTRAN, D_co, N, dEcr, Ecr,
                       GIc, sigcr0, taucr0, Leff, ν, scr_old, scr)

          [Inside calcEcr:]
            Y    ← D_cr + N^T · D_co · N
            Δε_cr ← Y⁻¹ · (N^T · D_co · DSTRAN)
            Ecr  ← Ecr_OLD + Δε_cr
            Ecr_MAX ← component-wise max(Ecr, Ecr_MAX)
            Ecr_OLD ← Ecr

          STATEV(7:8)   ← Ecr_MAX
          STATEV(9:10)  ← Ecr_OLD  (= Ecr)
          STATEV(11:12) ← Δε_cr

  STEP 9: Compute cracked tangent D_cocr
          CALL calcDcocr(D_co, D_cr, N, D_cocr)

          [Inside calcDcocr:]
            Y      ← D_cr + N^T · D_co · N
            D_cocr ← D_co − D_co · N · Y⁻¹ · N^T · D_co

          SIG    ← SIG + D_cocr · DSTRAN
          DDSDDE ← D_cocr

  RETURN
END SUBROUTINE
```

---

## Subroutine Reference

| Subroutine              | File                       | Role                                                                |
|-------------------------|----------------------------|---------------------------------------------------------------------|
| `sca2diso`              | `formulation_AFEM2D_PE.for`| Main SCA driver; coordinates all steps                             |
| `STIFFNESS_MATRIX_2D`   | `formulation_AFEM2D_PE.for`| Assembles 4×4 plane-strain isotropic **D**\_co                     |
| `getN`                  | `formulation_AFEM2D_PE.for`| Builds 4×2 crack transformation matrix **N**                        |
| `mysprind`              | `formulation_AFEM2D_PE.for`| Thin wrapper around ABAQUS `SPRIND` (principal stresses/directions) |
| `myExit`                | `formulation_AFEM2D_PE.for`| Thin wrapper around ABAQUS `XIT` (terminate analysis)              |
| `calcEcr`               | `formulation_AFEM2D_PE.for`| Computes crack strain increment Δε\_cr via stress continuity       |
| `calcDcr`               | `formulation_AFEM2D_PE.for`| Computes crack constitutive matrix **D**\_cr (loading/unloading)   |
| `calcDcocr`             | `formulation_AFEM2D_PE.for`| Assembles cracked-continuum tangent **D**\_cocr                    |
| `matrixInverse`         | `matrix_inverse.for`       | LU-based 2×2 matrix inverse (double precision)                     |

---

## Identified Bugs in the Original Code

| # | Location                     | Bug Description                                                                  |
|---|------------------------------|----------------------------------------------------------------------------------|
| 1 | `sca2diso`, line 23          | `DFGRD1` declared `INTEGER` — should be `REAL*8` (deformation gradient)         |
| 2 | `sca2diso`, lines 66–84      | `sigcr0` (= STATEV(16)) used before being assigned (line 189) → divide-by-zero on first call |
| 3 | `calcDcr`, line 385          | Quad-precision literal `2.0Q0` mixed into `REAL*8` context — should be `2.0D0` |
| 4 | `sca2diso`, lines 193–194    | `write(*,*)` debug statements execute every increment, flooding the `.msg` file |
| 5 | `calcDcr` / `calcEcr`        | No guard against `sigcr0 = 0` (first call before crack, causes NaN/Inf)         |
| 6 | `calcEcr`                    | No crack interpenetration check: normal crack strain can go negative              |
| 7 | `calcDcr`, line 397          | Shear stiffness `Dcr(2,2) = Dcr(1,1)*term*1e-8` is purely ad-hoc               |
| 8 | `calcEcr` / `calcDcocr`      | No near-singular check on **Y** before inversion (silent failure if rank-deficient) |

---

## References

- Bazant, Z.P. & Oh, B.H. (1983). Crack band theory for fracture of concrete.
  *Matériaux et Construction*, **16**(3), 155–177.
- de Borst, R., Crisfield, M.A., Remmers, J.J.C., & Verhoosel, C.V. (2012).
  *Nonlinear Finite Element Analysis of Solids and Structures* (2nd ed.). Wiley.
- Rots, J.G. (1988). *Computational Modeling of Concrete Fracture*. PhD Thesis, TU Delft.
- ABAQUS Theory Guide — Section on User-Defined Material Behaviors (UMAT).
