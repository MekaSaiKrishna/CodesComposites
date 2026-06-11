# AFEM3D — Abaqus UMAT: Smeared Crack Approach for Progressive Damage

## Overview

UMAT subroutines implementing **incremental progressive damage** using the **Smeared Crack Approach (SCA)** for composite materials. Currently implemented for **isotropic** materials; orthotropic (composite) extension is in progress via `meka_formulation_MAXS.for`.

Compatible element type: **C3D8** (3D 8-node brick).

---

## Active Version

> **`meka_formulation_AFEM7.for`** — current active formulation (last debugged: 02/17/2026)
>
> Entry point: **`main_UD_inc.for`** (incremental UMAT wrapper, calls AFEM7)

---

## File Summary

### Entry Points (Main UMAT Wrappers)

| File | Description | Last Debugged |
|---|---|---|
| `main_UD_inc.for` | **Active entry point.** Incremental UMAT wrapper; calls `meka_formulation_AFEM7.for` | 01/29/2026 |
| `main_UD_tot.for` | Total-strain UMAT wrapper; calls `meka_formulation_AFEM5.for` | 01/29/2026 |
| `main_UD_elas.for` | Elasticity-only UMAT wrapper; calls `meka_formulation_elas.for` | 01/29/2026 |

---

### Formulation Files (Isotropic SCA — Evolutionary Sequence)

| File | Description | Status | Last Debugged |
|---|---|---|---|
| `meka_formulation_AFEM7.for` | **CURRENT ACTIVE.** Incremental SCA with Newton-Raphson local crack strain update, quadruple precision, unloading/reloading, Bazant regularization | Active | 02/17/2026 |
| `meka_formulation_AFEM5.for` | Incremental SCA; fixed crack orientation, multiple failure criteria (Max Principal Stress / Max Stress), full Jacobian | Archived | 01/29/2026 |
| `meka_formulation_AFEM5v1.for` | Cleaner-formatted version of AFEM5 | Archived | — |
| `meka_formulation_AFEM5 - Copy.for` | Mid-development backup of AFEM5 | Archived | — |
| `meka_formulation_AFEM6.for` | Incremental SCA with overshoot cutback to peak; global softening slope | Archived | 11/12/2025 |
| `meka_formulation_AFEM4.for` | Total SCA; adds nondimensional energy ratio parameter (STATEV 9) and element length storage (STATEV 10) | Archived | — |
| `meka_formulation_AFEM3.for` | Total SCA; adds secant stiffness tracking (STATEV 8) | Archived | 03/09/2025 |
| `meka_formulation_AFEM2.for` | Total SCA; improved Jacobian and element size regularization via characteristic length | Archived | 02/26/2025 |
| `meka_formulation_AFEM.for` | **Original.** Total SCA isotropic (tension only), 3-zone traction-separation law | Archived | 01/31/2025 |

---

### Orthotropic Formulation

| File | Description | Status | Last Debugged |
|---|---|---|---|
| `meka_formulation_MAXS.for` | Orthotropic total SCA with 9-mode Max Stress failure criteria (fiber/matrix tension, compression, 3 shear planes). Quad precision Newton-Raphson. | In development | 01/31/2025 |
| `meka_formulation_elas.for` | Pure isotropic/orthotropic elasticity (no damage) | Support | 11/17/2025 |

---

### Numerical Utilities

| File | Description |
|---|---|
| `matrix_inverse.for` | LU decomposition, matrix inversion, and linear system solver in both double (`REAL*8`) and quadruple (`REAL*16`) precision. Used by all formulation files. Based on Numerical Recipes. |
| `standard_support.for` | Abaqus wrapper utilities: `mysprind()` (principal stress eigenvalues), `myExit()` |

---

## Test Cases

See [CHECKLIST.MD](CHECKLIST.MD) for the full test matrix.

Metrics tracked per test:
- Simulation run time
- Peak stress (must match defined material strength)
- Failure strain
- Area under stress-strain curve (fracture energy)

---

## Key Technical Features

| Feature | Detail |
|---|---|
| Failure criteria | Max Principal Stress (isotropic), Max Stress (orthotropic, 9 modes) |
| Crack orientation | Fixed at initiation; tracked via transformation matrix |
| Regularization | Bazant element size check; characteristic length `L` stored in STATEV |
| Numerical stability | Newton-Raphson iteration for consistent crack strain update (AFEM7) |
| Precision | Quadruple precision (`REAL*16`) in critical NR solver subroutines |
| Unloading/reloading | Fully tracked via zone flags and stored max crack strain (AFEM7) |



# AFEM3D — ABAQUS UMAT: Incremental Smeared Crack Approach

> **Active formulation:** `meka_formulation_AFEM8.for`  
> **Entry point:** `main_UD_inc8.for`  
> **Last debugged:** 06/11/2026  
> **Contact:** meka1@purdue.edu

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [File Index](#2-file-index)
3. [AFEM8 — Theory and Code Walkthrough](#3-afem8--theory-and-code-walkthrough)
4. [State Variables and Properties](#4-state-variables-and-properties)
5. [Test Cases](#5-test-cases)
6. [Known Issues and Plan of Action](#6-known-issues-and-plan-of-action)
7. [Debug Log](#7-debug-log)

---

## 1. Project Overview

This project implements **incremental progressive damage** for isotropic materials in ABAQUS using the **Smeared Crack Approach (SCA)**. Rather than inserting a discrete crack surface, the crack is smeared over the element volume using a characteristic length `Leff`, and its effect on the stress-strain response is captured through a modified continuum tangent `Dcocr`.

**Compatible element:** C3D8 (3D 8-node brick, linear).  
**Failure criterion (current):** Max Principal Stress.  
**Material scope (current):** Isotropic only.

---

## 2. File Index

### Entry Points

| File | Description | Status |
|---|---|---|
| `main_UD_inc8.for` | **ACTIVE.** UMAT wrapper for AFEM8. Passes `PNEWDT` for overshoot cutbacks. | Active |
| `main_UD_inc.for` | UMAT wrapper for AFEM7 (no `PNEWDT`). | Archived |
| `main_UD_tot.for` | UMAT wrapper for AFEM5 (total-strain form). | Archived |
| `main_UD_elas.for` | Elasticity-only UMAT (no damage). | Archived |

### Isotropic Formulations (evolutionary sequence)

| File | Description | Status | Last Debugged |
|---|---|---|---|
| `meka_formulation_AFEM8.for` | **ACTIVE.** AFEM7 + PNEWDT overshoot, J_FAIL tracking, double `calcDcr` call, interpenetration clamp, near-singular guards | Active | 06/11/2026 |
| `meka_formulation_AFEM7.for` | Incremental SCA, unloading/reloading, Bazant regularization. No overshoot detection. | Archived | 02/17/2026 |
| `meka_formulation_AFEM6.for` | Incremental SCA with global (not local) overshoot cutback | Archived | 11/12/2025 |
| `meka_formulation_AFEM5.for` | Incremental SCA, fixed crack, multiple criteria, full Jacobian | Archived | 01/29/2026 |
| `meka_formulation_AFEM5v1.for` | Clean reformat of AFEM5 | Archived | — |
| `meka_formulation_AFEM4.for` | Total SCA + energy ratio tracking | Archived | — |
| `meka_formulation_AFEM3.for` | Total SCA + secant stiffness tracking | Archived | 03/09/2025 |
| `meka_formulation_AFEM2.for` | Total SCA + Bazant regularization | Archived | 02/26/2025 |
| `meka_formulation_AFEM.for` | **Original.** Total SCA, 3-zone TSL | Archived | 01/31/2025 |

### Orthotropic / Utility

| File | Description | Status |
|---|---|---|
| `meka_formulation_MAXS.for` | Orthotropic total SCA, 9-mode Max Stress, quad-precision NR | In development |
| `meka_formulation_elas.for` | Pure elasticity (iso/ortho), no damage | Support |
| `matrix_inverse.for` | LU decomposition + inversion, double and quad precision | Support |
| `standard_support.for` | ABAQUS wrappers: `mysprind()`, `myExit()` | Support |

---

## 3. AFEM8 — Theory and Code Walkthrough

### 3.1 What is the Smeared Crack Approach?

A discrete crack (discontinuity) is represented as additional strain smeared over an element of width `Leff`. The total strain is split:

```
   ε_total = ε_elastic + ε_crack
                                          ┌──────────────────┐
  Applied                                 │                  │
  displacement                        ────┤   CRACK ZONE     ├────
      │                                   │  (smeared over   │
      ▼                                   │    Leff = CELENT)│
  ┌───┬───┐        ┌───┬╔═══╗───┐         └──────────────────┘
  │   │   │  ──►  │   │║///║   │
  │   │   │        │   │║///║   │   ε_crack smeared
  └───┴───┘        └───┴╚═══╝───┘   over the element width
  intact            cracked
```

The **characteristic length** `Leff = CELENT` (from ABAQUS) regularises the fracture energy so the result is mesh-independent (Bazant scaling).

---

### 3.2 Traction-Separation Law (TSL)

The crack normal traction `σ_cr` is a linear function of crack opening strain `Ecr_n`:

```
  σ_cr
   │
XT ┼────────────╮  ◄── peak (crack initiates here)
   │            │\
   │            │  \   slope_soft = -XT² / (2·GIc/Leff)
   │            │    \
   │            │      \
 0 ┼────────────┼────────╯──────── Ecr_n
   0                 Ecr_n_ult
                = 2·GIc / (Leff·XT)

  Fracture energy GIc = area under triangle = ½·XT·Ecr_n_ult·Leff
```

**Loading** (Ecr_n increasing):  `σ_cr = XT + slope_soft · Ecr_n`  
**Unloading** (Ecr_n < Ecr_MAX): secant from origin to (Ecr_MAX, σ_cr_at_max)  
**Complete failure** (Ecr_n > Ecr_n_ult): residual stiffness `RES_K = 1e-8` to prevent singularity

---

### 3.3 Crack Transformation Matrix N

ABAQUS `SPRIND` returns three principal stress eigenvectors `T(1,:)`, `T(2,:)`, `T(3,:)`. The index `J_FAIL` (1, 2, or 3) identifies which eigenvector is the crack **normal** direction — the one whose principal stress first reached `XT`.

The **N matrix** (6×3) maps the local crack strain vector `[Ecr_n, Ecr_s1, Ecr_s2]` to global Voigt strain:

```
  LOCAL CRACK FRAME              GLOBAL VOIGT
                                 ┌────┐
   Ecr_n  (opening,              │ ε₁₁│
           along T(J_FAIL,:))    │ ε₂₂│
   Ecr_s1 (sliding,         [N]  │ ε₃₃│
           along T(j≠J_FAIL))──► │ γ₁₂│
   Ecr_s2 (sliding,              │ γ₁₃│
           along T(k≠J_FAIL))    │ γ₂₃│
                                 └────┘
  Column J of N = sym(n ⊗ d_J) in Voigt,
  where n = T(J_FAIL,:), d_J = T(J,:)
```

Because ABAQUS orders `PRINS(1) ≥ PRINS(2) ≥ PRINS(3)` for stress, `J_FAIL = 1` in almost all tensile scenarios. It is tracked explicitly for robustness.

---

### 3.4 The Cracked Continuum Tangent Dcocr

Once the crack strain increment `dEcr` is known, the cracked tangent (the Woodbury identity applied to the composite material system) is:

```
  Dcocr = Dco  −  Dco·N·[Dcr + Nᵀ·Dco·N]⁻¹·Nᵀ·Dco
  ─────   ────     ────────────────────────────────
  6×6     elastic  correction (rank-3 update)

  Meaning: the overall stiffness = elastic stiffness
           MINUS the stiffness "lost" to crack compliance
```

The 3×3 matrix `Dcr` is diagonal. Only `Dcr(J_FAIL, J_FAIL)` follows the TSL. The two shear components get a near-zero residual stiffness scaled by `term·1e-5` where `term = 1/(2(1+ν))`.

---

### 3.5 Increment Flow (Post-Peak Block)

Each cracked increment steps through:

```
  ┌── READ ──────────────────────────────────────────────┐
  │  Ecr, Ecr_MAX, J_FAIL, T rows  ← STATEV             │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP A ─────────────────────────────────────────────┐
  │  calcDcr(Ecr_old, Ecr_MAX, J_FAIL) → Dcr_old         │
  │  (used only to build system matrix Y; scr discarded)  │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP B ─────────────────────────────────────────────┐
  │  Y = Dcr_old + Nᵀ·Dco·N                              │
  │  dEcr = Y⁻¹ · Nᵀ·Dco·DSTRAN         (solve)         │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP C ─────────────────────────────────────────────┐
  │  POST-PEAK OVERSHOOT CHECK:                           │
  │  if Ecr(J_FAIL)+dEcr(J_FAIL) > 1.05·Ecr_n_ult       │
  │     → compute fraction α, set PNEWDT = α, continue   │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP D ─────────────────────────────────────────────┐
  │  Ecr      = Ecr + dEcr                               │
  │  Ecr_MAX  = max(Ecr, Ecr_MAX)   component-wise       │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP E ─────────────────────────────────────────────┐
  │  INTERPENETRATION CLAMP:                              │
  │  if Ecr(J_FAIL) < 0 → force Ecr(J_FAIL) = 0         │
  │  (tensile crack cannot close past its reference)     │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP F ─────────────────────────────────────────────┐
  │  calcDcr(Ecr_new, Ecr_MAX, J_FAIL) → Dcr_new, scr_n  │
  │  (consistent update: eliminates 1-increment lag)     │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP G ─────────────────────────────────────────────┐
  │  calcDcocr(Dco, Dcr_new, N) → Dcocr                  │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP H ─────────────────────────────────────────────┐
  │  WRITE → STATEV (Ecr, Ecr_MAX, dEcr, Dcr_n, scr_n)  │
  └──────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP I ─────────────────────────────────────────────┐
  │  SIG    = SIG + Dcocr · DSTRAN                       │
  │  DDSDDE = Dcocr           (implicit solver)          │
  └──────────────────────────────────────────────────────┘
```

---

### 3.6 Pre-Peak Overshoot Handling

When the elastic trial stress overshoots the failure surface by more than 5%, the code interpolates the fraction `α` of the increment that lands exactly on the surface and requests ABAQUS cut the step:

```
  FI
   │
1.05────────────────── ← overshoot detected
   │                 ╱
1.0─────────────────╱── ← failure surface (XT)
   │           ╱   ╱
   │          ╱   ╱
   │         ╱
0  ┼─────────────────── increment
   t_n   α·Δt  t_{n+1}

  α = (1.0 - FI_prev) / (FI_new - FI_prev)
  PNEWDT = α  →  ABAQUS retries with Δt_new = α·Δt
```

---

## 4. State Variables and Properties

### State Variables (NSTATV = 26)

| Slot | Variable | Description |
|---|---|---|
| 1 | `isCRACKED` | 0 = intact, 1 = cracked |
| 2 | `MODE` | Failure mode (1 = tensile) |
| 5 | `L` | Characteristic length `CELENT` at crack initiation |
| 6–8 | `Ecr_MAX(1:3)` | Max crack strain reached (component J_FAIL = crack opening) |
| 9–11 | `T(1,:)` | SPRIND eigenvector row 1 (stored as-is, global frame) |
| 12–14 | `T(2,:)` | SPRIND eigenvector row 2 |
| 15–17 | `T(3,:)` | SPRIND eigenvector row 3 |
| 18–20 | `Ecr(1:3)` | Current crack strain |
| 21–23 | `dEcr(1:3)` | Crack strain increment this step |
| 24 | `Dcr_n` | `Dcr(J_FAIL,J_FAIL)` — normal crack stiffness |
| 25 | `scr_n` | Normal crack traction (post-processing only) |
| 26 | `J_FAIL` | Index of max principal stress eigenvector (1, 2, or 3) |

### Material Properties (NPROPS = 5)

| Slot | Symbol | Description |
|---|---|---|
| 1 | `E` | Young's modulus |
| 2 | `ν` | Poisson's ratio |
| 3 | `GIc` | Mode-I fracture energy |
| 4 | `XT` | Tensile strength (peak stress) |
| 5 | `INDEX` | Failure criterion: 0 = Max Principal Stress |

### `.inp` file keyword

```
*Depvar
     26,
*User Material, constants=5
E,  nu,  GIc,  XT,  INDEX
```

---

## 5. Test Cases

All test cubes are unit cubes (1×1×1 m) with a single C3D8 element. Material: E=10000, ν=0, GIc=1, XT=100, INDEX=0 (unless noted). Applied displacement = 0.025 m.

**Expected results for uniaxial tension (analytical):**
- Failure strain: `εf = XT/E = 100/10000 = 0.01`
- Ecr_n_ult: `2·GIc/(Leff·XT) = 2·1/(1·100) = 0.02`
- Total failure displacement: `0.01 + 0.02 = 0.03 m`

### Isotropic 1-crack tests (no unloading)

| # | Input File | Load Direction | Expected Failure Stress | Status |
|---|---|---|---|---|
| T-1 | `plus11_iso.inp` | +X tension (S11) | 100 MPa | ✅ Running |
| T-2 | `plus22_iso.inp` | +Y tension (S22) | 100 MPa | ✅ BCs fixed 06/11/2026 |
| T-3 | `tension33_iso.inp` | +Z tension (S33) | 100 MPa | Not tested |
| S-1 | `simpshear12_iso.inp` | Simple shear 12 (γ12) | τ = XT·G/E | ⚠ Shear not degrading |
| S-2 | `pureshear12_iso.inp` | Pure shear 12 (S12) | 100 MPa (principal) | ⚠ Shear not degrading |

### Metrics tracked per test

| Metric | How measured |
|---|---|
| Run time | Wall clock, noted in this table |
| Peak stress | Must match `XT` from material definition |
| Failure strain | SDV history output at failure increment |
| Fracture energy | Area under S-ε curve ≈ `GIc / Leff` |

---

## 6. Known Issues and Plan of Action

### Issue 1 — Spurious twisting / convergence failure at tail end of softening

**Symptom:** Near complete failure (Ecr_n → Ecr_n_ult), the simulation becomes unstable. The element develops spurious off-axis stresses and the job aborts.

**Root cause:** When `Dcr(J_FAIL,J_FAIL) → RES_K = 1e-8` (residual), the shear crack stiffness is set to `Dcr(j,j) = RES_K × term × 1e-5 ≈ 0`. The N matrix has non-zero entries coupling all 6 Voigt components. At near-zero `Dcr`, the Woodbury formula gives `Dcocr ≈ Dco - Dco·N·(Nᵀ·Dco·N)⁻¹·Nᵀ·Dco`, which removes the stress component along the crack normal but still has full off-diagonal coupling. The residual terms drive spurious shear stresses that the implicit solver cannot resolve, causing hourglass-type twisting.

**Proposed fix:** After reaching complete failure (`Ecr_n > Ecr_n_ult`), zero out the rows and columns of `Dcocr` associated with all directions coupled by N, rather than relying on the near-zero `Dcr` to do it numerically.

---

### Issue 2 — Shear stress not degrading under shear loading

**Symptom:** Under simple shear (γ12 applied) or pure shear, the crack initiates at 45° correctly, but `S12` does not soften — it either remains constant or continues growing.

**Root cause:** Under shear loading, the max principal stress is at 45°, so the crack opens in mode-I at 45°. Analysing the Dcocr formula for this geometry:

```
  Shear loading (γ12):                  Crack at 45°:
                                          n = [1/√2, 1/√2, 0]
     ┌─────────┐
  ───┤  γ12→   ├───     Principal stresses: σ₁ = +τ (at 45°)
     │         │                            σ₂ = -τ (at 135°)
     └─────────┘

  N(:,1) = [0.5, 0.5, 0, 1, 0, 0]   (maps crack opening to ε11, ε22, γ12)

  Dcocr(4,4) = G - G² / (Dcr₁₁ + Nᵀ·Dco·N)₁₁
             = G - G² / (slope_soft + 13461)

  During LOADING  (slope_soft = -5000): Dcocr(4,4) = 3846 - 1749 = 2097
  At RESIDUAL     (Dcr → 0):           Dcocr(4,4) = 3846 - 1099 = 2748
```

The shear stiffness **does not reach zero** because the crack opens at 45° (mode-I), not parallel to the shear plane (mode-II). The current single-crack model can only reduce S12 by approximately 45% — it cannot fully degrade shear to zero.

**Proposed fix — Mode-II shear tracking:**  
Activate the crack-sliding components `Ecr_s1, Ecr_s2` with their own TSL (mode-II fracture energy `GIIc`). For the 45° crack under shear: the shear sliding component `dEcr(2)` (along T(2,:) = [-1/√2, 1/√2, 0]) also becomes active. Adding a mode-II TSL for `Dcr(j,j)` (j ≠ J_FAIL) would allow the shear to degrade fully.

**Alternative fix — Direct shear failure criterion:**  
For pure shear cases, detect shear failure directly (`τ12 ≥ XS`) and orient the crack normal perpendicular to the shear plane (not at 45°). This requires an additional failure criterion index (INDEX = 1 for Max Shear).

---

### Issue 3 — plus22_iso.inp had wrong boundary conditions (FIXED 06/11/2026)

The original file applied x-displacement (DOF 1) to the x=1 face while locking the x=0 face in x — effectively a second x-tension test, not y-tension. Fixed by:
- Loading: y=1 face (nodes 1,3,5,7) displaced in DOF 2 by 0.025
- Fixed end: y=0 face (nodes 2,4,6,8) fully fixed in DOF 2
- Anti-rigid-body: nodes 6,8 in DOF 1; node 8 in DOF 3

---

## 7. Debug Log

| Date | File(s) Changed | Change | Effect |
|---|---|---|---|
| 06/10/2026 | `meka_formulation_AFEM8.for` (new), `main_UD_inc8.for` (new) | Created AFEM8: ported overshoot logic from `logic-2d.for`, added `PNEWDT`, double `calcDcr` call, interpenetration clamp | AFEM7 archived |
| 06/10/2026 | `meka_formulation_AFEM8.for` | Added `J_FAIL` tracking: parametric crack-normal index from SPRIND output, stored in `STATEV(26)`. `calcDcr` updated to accept `J_FAIL` and set `Dcr(J_FAIL,J_FAIL)` as normal stiffness. Off-diagonal shear: `Dcr(j,j) = Dcr(J_FAIL,J_FAIL)·term·1e-5` | Fixes hardcoded index-1 assumption |
| 06/10/2026 | `meka_formulation_AFEM8.for` | `STATEV(24)` (Dcr normal) uncommented and updated to `Dcr(J_FAIL,J_FAIL)` | Correct SDV output |
| 06/10/2026 | `plus11_iso.inp`, `plus22_iso.inp`, `simpshear12_iso.inp`, `pureshear12_iso.inp` | `*Depvar` increased from 25 to 26 | Matches new STATEV layout |
| 06/11/2026 | `plus22_iso.inp` | Corrected boundary conditions for y-tension. Fixed DOFs, loading face (nodes 1,3,5,7), load direction DOF 2. Removed spurious x-face constraint. | Eliminated twisting/abort |
| 02/17/2026 | `meka_formulation_AFEM7.for` | Fixed `Leff = 1.0*L` (was commented out → Leff=0 → infinite softening slope → stress plateau). Fixed `*Depvar 12→25`. Removed dead STATEV reads. | Stress now softens in 11-tension |

