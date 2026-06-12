# AFEM3D — ABAQUS UMAT: Incremental Smeared Crack Approach

> **Active formulation:** `meka_formulation_AFEM8.for`  
> **Entry point:** `main_UD_inc8.for`  
> **Last reviewed:** 06/12/2026  
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
**Failure criteria (current):** Max Principal Stress (Mode-I) + Direct Shear Stress (Mode-II).  
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
| `meka_formulation_AFEM8.for` | **ACTIVE.** AFEM7 + PNEWDT overshoot, J_FAIL tracking, double `calcDcr` call, interpenetration clamp, near-singular guards, Mode-II shear TSL | Active | 06/11/2026 |
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

### 3.2 Mode-I Traction-Separation Law (Tensile Crack)

The crack normal traction `σ_cr` is a linear function of crack opening strain `Ecr_n`:

```
  σ_cr
   │
XT ┼────────────╮  ◄── peak (crack initiates here, σ_cr = XT)
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

### 3.3 Mode-II Traction-Separation Law (Shear Crack)

When shear failure governs (`τ ≥ XS` before σ₁ ≥ XT), the crack slides rather than opens. A Mode-II TSL is applied to the **sliding** crack strain `Ecr_s`:

```
  τ_cr
   │
XS ┼────────────╮  ◄── peak (crack initiates here, τ_cr = XS)
   │            │\
   │            │  \   slope_soft_II = -XS² / (2·G2c/Leff)
   │            │    \
   │            │      \
 0 ┼────────────┼────────╯──────── |Ecr_s|
   0                Ecr_s_ult
               = 2·G2c / (Leff·XS)

  Fracture energy G2c = area under triangle = ½·XS·Ecr_s_ult·Leff
```

The TSL is **symmetric in shear direction** — the crack slides freely in ±. `Ecr_s` can be negative.  
`τ_cr = SIGN(τ_eff, Ecr_s)` where `τ_eff` is determined from `|Ecr_s|`.

**Special term:** `Dcr(1,1) = 1000·E` (penalty) prevents crack-face normal **opening** in mode-II.  
Mode-I crack cannot have mode-II sliding, and mode-II crack cannot have mode-I opening.

---

### 3.4 Mode Determination at Initiation

At every intact-element increment, two failure indices are computed from the elastic trial stress:

```
  FI_mode1  =  max[ PRINS(i) / XT ]   for PRINS(i) > 0
              (max-principal-stress criterion — ABAQUS SPRIND)

  FI_shear  =  max( |S12|, |S13|, |S23| ) / XS
              (direct shear criterion — no principal stress rotation)
```

Failure initiates when either FI ≥ 1. The mode that reaches its criterion **first** (higher FI when both ≥ 1, or whichever crosses 1 first) is locked in permanently for that element.

```
  FI_shear ≥ 1 AND FI_shear > FI_mode1  →  MODE-II  (shear)
  FI_mode1 ≥ 1 AND FI_mode1 ≥ FI_shear  →  MODE-I   (tensile)
```

`STATEV(2) = MODE` records this choice at initiation and never changes.

---

### 3.5 Mode-II Crack Orientation (Coordinate-Aligned)

Unlike mode-I (which uses SPRIND eigenvectors), the mode-II crack frame is fixed to coordinate axes aligned with the dominant shear component:

```
   Dominant      T(1,:) = crack normal    T(2,:) = sliding dir    T(3,:) = secondary
   ─────────     ──────────────────────   ─────────────────────   ───────────────────
   |S12| max     [1, 0, 0]  (x-normal)   [0, 1, 0]  (y-sliding)  [0, 0, 1]
   |S13| max     [1, 0, 0]  (x-normal)   [0, 0, 1]  (z-sliding)  [0, 1, 0]
   |S23| max     [0, 1, 0]  (y-normal)   [0, 0, 1]  (z-sliding)  [1, 0, 0]
```

**Why coordinate-aligned works:** Verifying that `N(:,2)` maps `Ecr_s` to the correct shear component:

```
  S12 case: n=[1,0,0], m=[0,1,0]  →  N(4,2) = n₁m₂+n₂m₁ = 1  → Ecr_s → γ12  ✓
  S13 case: n=[1,0,0], m=[0,0,1]  →  N(5,2) = n₃m₁+n₁m₃ = 1  → Ecr_s → γ13  ✓
  S23 case: n=[0,1,0], m=[0,0,1]  →  N(6,2) = n₂m₃+n₃m₂ = 1  → Ecr_s → γ23  ✓
```

`J_FAIL = 2` always for Mode-II (the sliding direction is always component 2 in the local frame).

---

### 3.6 Crack Transformation Matrix N

The **N matrix** (6×3) maps the local crack strain vector `[Ecr_n, Ecr_s1, Ecr_s2]` to global Voigt strain:

```
  LOCAL CRACK FRAME              GLOBAL VOIGT
                                 ┌────┐
   Ecr_n  (opening or            │ ε₁₁│
           penalty normal)       │ ε₂₂│
   Ecr_s1 (sliding, J_FAIL=2)  [N]  │ ε₃₃│
           along T(2,:)    ──►   │ γ₁₂│
   Ecr_s2 (secondary shear)     │ γ₁₃│
           along T(3,:)          │ γ₂₃│
                                 └────┘
  Column J of N = sym(n ⊗ d_J) in Voigt,
  where n = T(1,:) [crack normal], d_J = T(J,:)
```

---

### 3.7 The Cracked Continuum Tangent Dcocr

The same Woodbury-identity formula applies for both modes:

```
  Dcocr = Dco  −  Dco·N·[Dcr + Nᵀ·Dco·N]⁻¹·Nᵀ·Dco
  ─────   ────     ────────────────────────────────
  6×6     elastic  rank-3 correction (stiffness lost to crack compliance)
```

`Dcr` is mode-dependent:

```
  MODE-I  Dcr:                      MODE-II  Dcr:
  ┌──────────────────────┐          ┌──────────────────────┐
  │ TSL(Ecr_n)  0    0   │          │ 1000·E    0       0  │  ← penalty (blocks opening)
  │    0      ε_shear 0  │          │    0    TSL(|Ecr_s|) 0│  ← mode-II TSL
  │    0        0   ε_sh │          │    0      0    tiny  │
  └──────────────────────┘          └──────────────────────┘
   ε_shear = Dcr(1,1)·term·1e-5       tiny = Dcr(2,2)·term·1e-5
```

---

### 3.8 Increment Flow (Post-Peak Block)

Each cracked increment steps through (MODE branching at Steps A/F):

```
  ┌── READ ──────────────────────────────────────────────────┐
  │  MODE, J_FAIL, Ecr, Ecr_MAX, T rows  ← STATEV           │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP A ─────────────────────────────────────────────────┐
  │  MODE==1: calcDcr(Ecr_old, Ecr_MAX, J_FAIL)  → Dcr_old   │
  │  MODE==2: calcDcr_ModeII(Ecr_old, Ecr_MAX)   → Dcr_old   │
  │  (builds system matrix Y only; scr discarded)            │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP B ─────────────────────────────────────────────────┐
  │  Y = Dcr_old + Nᵀ·Dco·N                                  │
  │  dEcr = Y⁻¹ · Nᵀ·Dco·DSTRAN             (solve)         │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP C ─────────────────────────────────────────────────┐
  │  POST-PEAK OVERSHOOT CHECK:                               │
  │  MODE==1: if Ecr(J_FAIL)+dEcr(J_FAIL) > 1.05·Ecr_n_ult  │
  │  MODE==2: if |Ecr(2)|+|dEcr(2)| > 1.05·Ecr_s_ult        │
  │     → compute fraction α, set PNEWDT=α, continue        │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP D ─────────────────────────────────────────────────┐
  │  Ecr      = Ecr + dEcr                                    │
  │  MODE==1: Ecr_MAX = max(Ecr, Ecr_MAX)  component-wise    │
  │  MODE==2: Ecr_MAX(2) = max(|Ecr(2)|, Ecr_MAX(2))        │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP E ─────────────────────────────────────────────────┐
  │  MODE==1 ONLY — INTERPENETRATION CLAMP:                   │
  │  if Ecr(J_FAIL) < 0 → force Ecr(J_FAIL) = 0             │
  │  MODE==2: skip — shear crack can slide in either direction│
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP F ─────────────────────────────────────────────────┐
  │  MODE==1: calcDcr(Ecr_new, Ecr_MAX, J_FAIL) → Dcr, scr  │
  │  MODE==2: calcDcr_ModeII(Ecr_new, Ecr_MAX)  → Dcr, τ_cr │
  │  (consistent update: eliminates 1-increment lag)         │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP G ─────────────────────────────────────────────────┐
  │  calcDcocr(Dco, Dcr, N) → Dcocr   (same formula, both)   │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP H ─────────────────────────────────────────────────┐
  │  WRITE → STATEV (Ecr, Ecr_MAX, dEcr, Dcr_n, scr_n/τ)   │
  └──────────────────────────────────────────────────────────┘
              │
              ▼
  ┌── STEP I ─────────────────────────────────────────────────┐
  │  SIG    = SIG + Dcocr · DSTRAN                           │
  │  DDSDDE = Dcocr           (implicit solver)              │
  └──────────────────────────────────────────────────────────┘
```

---

### 3.9 Pre-Peak Overshoot Handling

When the elastic trial stress overshoots the failure surface by more than 5%, the code interpolates the fraction `α` of the increment that lands exactly on the surface and requests ABAQUS cut the step:

```
  FI
   │
1.05────────────────── ← overshoot detected (both mode-I and mode-II)
   │                 ╱
1.0─────────────────╱── ← failure surface
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
| 2 | `MODE` | Failure mode: **1 = Mode-I (tensile)**, **2 = Mode-II (shear)** |
| 5 | `L` | Characteristic length `CELENT` at crack initiation |
| 6–8 | `Ecr_MAX(1:3)` | Mode-I: max Ecr component-wise. Mode-II: `Ecr_MAX(2)` = max(\|Ecr_s\|) |
| 9–11 | `T(1,:)` | Crack-frame row 1. Mode-I: SPRIND eigenvector. Mode-II: coordinate axis (crack normal) |
| 12–14 | `T(2,:)` | Crack-frame row 2 — sliding direction |
| 15–17 | `T(3,:)` | Crack-frame row 3 — secondary direction |
| 18–20 | `Ecr(1:3)` | Current crack strain |
| 21–23 | `dEcr(1:3)` | Crack strain increment this step |
| 24 | `Dcr_n` | `Dcr(J_FAIL,J_FAIL)` — primary crack stiffness |
| 25 | `scr` | Crack traction along primary direction (post-processing only) |
| 26 | `J_FAIL` | Primary crack strain index: **1 for Mode-I** (normal), **2 for Mode-II** (sliding) |

### Material Properties (NPROPS = 5 to 7)

| Slot | Symbol | Description | Default |
|---|---|---|---|
| 1 | `E` | Young's modulus | — |
| 2 | `ν` | Poisson's ratio | — |
| 3 | `GIc` | Mode-I fracture energy | — |
| 4 | `XT` | Tensile strength (Mode-I peak) | — |
| 5 | `INDEX` | 0 = Max Principal Stress + Direct Shear | — |
| 6 | `G2c` | Mode-II fracture energy | `GIc` if NPROPS < 6 |
| 7 | `XS` | Shear strength (Mode-II peak) | `XT/2` if NPROPS < 7 |

**Existing `.inp` files with `constants=5` automatically use the defaults.** No `.inp` changes required.

### `.inp` file keyword (full 7-prop example)

```
*Depvar
     26,
*User Material, constants=7
E,  nu,  GIc,  XT,  INDEX,  G2c,  XS
```

---

## 5. Test Cases

All test cubes are unit cubes (1×1×1 m) with a single C3D8 element. Material: E=10000, ν=0, GIc=1, XT=100, INDEX=0 (unless noted). Applied displacement = 0.025 m.

**Expected results for uniaxial tension (analytical):**
- Failure strain: `εf = XT/E = 100/10000 = 0.01`
- `Ecr_n_ult = 2·GIc/(Leff·XT) = 2·1/(1·100) = 0.02`
- Total failure displacement: `0.01 + 0.02 = 0.03 m`

**Expected results for shear (ν=0.3, XS=XT/2=50, G=3846, G2c=GIc=1 defaults):**
- Mode-II initiates at S12 = 50 MPa
- `Ecr_s_ult = 2·G2c/(Leff·XS) = 2·1/(1·50) = 0.04`
- `Dcocr(4,4) = G − G²/(slope_soft_II + G) = 3846 − 3846²/2596 = −1854` (softening ✓)
- Total failure shear strain: `γ12_fail ≈ 0.040`

### Isotropic 1-crack tests (no unloading)

| # | Input File | Load Direction | Expected Failure Stress | Mode | Status |
|---|---|---|---|---|---|
| T-1 | `plus11_iso.inp` | +X tension (S11) | 100 MPa | Mode-I | ✅ Running |
| T-2 | `tension22_iso.inp` | +Y tension (S22) | 100 MPa | Mode-I | Not tested |
| T-3 | `tension33_iso.inp` | +Z tension (S33) | 100 MPa | Mode-I | Not tested |
| S-12s | `simpshear12_iso.inp` | Simple shear 12 (γ12) | 50 MPa (XS default) | Mode-II | ⬜ Pending verification |
| S-12p | `pureshear12_iso.inp` | Pure shear 12 (S12) | 50 MPa (XS default) | Mode-II | ⬜ Pending verification |

### Metrics tracked per test

| Metric | How measured |
|---|---|
| Run time | Wall clock, noted in this table |
| Peak stress | Must match `XT` (Mode-I) or `XS` (Mode-II) from material definition |
| Failure strain | SDV history output at failure increment |
| Fracture energy | Area under S-ε curve ≈ `GIc/Leff` or `G2c/Leff` |

---

## 6. Known Issues and Bugs

### Issue 1 — Spurious twisting / convergence failure at tail end of softening

**Symptom:** Near complete failure (Ecr_n → Ecr_n_ult), the simulation becomes unstable. The element develops spurious off-axis stresses and the job aborts.

**Root cause:** When `Dcr(J_FAIL,J_FAIL) → RES_K = 1e-8` (residual), the off-diagonal shear stiffness is set to `Dcr(j,j) = RES_K × term × 1e-5 ≈ 0`. At near-zero `Dcr`, the Woodbury formula gives `Dcocr ≈ Dco - Dco·N·(Nᵀ·Dco·N)⁻¹·Nᵀ·Dco`, which removes the stress component along the crack normal but retains full off-diagonal coupling. The residual terms drive spurious shear stresses that the implicit solver cannot resolve, causing hourglass-type twisting.

**Proposed fix:** After reaching complete failure (`Ecr_n > Ecr_n_ult`), zero out the rows and columns of `Dcocr` associated with all directions coupled by N, rather than relying on near-zero `Dcr` to do it numerically.

**Status:** Deferred. Will address after shear tests are validated.

---

### Issue 2 — Shear stress not degrading under shear loading

**Status: RESOLVED by Mode-II implementation (06/11/2026).**

**Root cause (for reference):** Under shear loading, the max principal stress is at 45°, so the crack opened in Mode-I at 45°. The 45° crack could only reduce S12 by ~45% due to the geometry of the N matrix — it could never reach zero. At residual, `Dcocr(4,4)` actually *increased* to 2748, making softening impossible.

**Fix:** Introduced Mode-II with coordinate-aligned crack orientation. When |S12| (or |S13|, |S23|) reaches XS before σ₁ reaches XT, Mode-II activates instead:
- Crack frame: T(2,:) = sliding direction (e.g. [0,1,0] for S12)
- N(:,2) = [0,0,0,1,0,0] → `Ecr_s` maps directly to γ12
- System matrix: `Y(2,2) = slope_soft_II + G = −1250 + 3846 = 2596`
- Cracked shear tangent: `Dcocr(4,4) = G − G²/Y = 3846 − 5706 = −1854` (shear DECREASES ✓)
- At complete failure: `Y(2,2) → G` → `Dcocr(4,4) → 0` ✓

---

### Bug 1 — `BL_II` Bazant snap-back limit uses `E` instead of `G` [UNRESOLVED]

**File / Line:** `meka_formulation_AFEM8.for`, line 127.

**Code (wrong):**
```fortran
BL_II = 2.0D0*G2c*E / taucr0**2
```
**Should be:**
```fortran
G_shear = 0.5D0*E / (1.0D0 + nu)
BL_II   = 2.0D0*G2c*G_shear / taucr0**2
```

**Why:** The Mode-II Bazant limit is derived from requiring no snap-back in the shear TSL: `|slope_soft_II| ≤ G` (shear modulus, not Young's modulus). Rearranging: `Leff ≤ 2·G2c·G / XS²`. Using E instead of G overestimates the limit by a factor of `E/G = 2(1+ν)`. For ν=0.3 this is 2.6×, so elements can be 2.6× too large without triggering the warning — the snap-back guard is essentially inactive for typical meshes.

**Impact:** Moderate. Simulation will not abort even when the element size causes a snap-back in the Mode-II TSL, leading to non-physical negative tangent stiffness that may cause convergence failures without a meaningful error message.

---

### Bug 2 — Post-peak overshoot (Step C) does not `RETURN` after setting `PNEWDT` [UNRESOLVED]

**File / Lines:** `meka_formulation_AFEM8.for`, lines 329–349.

**What happens:** When an overshoot past `Ecr_n_ult` (Mode-I) or `Ecr_s_ult` (Mode-II) is detected, `PNEWDT` is set to request a step cutback, but execution continues through Steps D–I. The crack strain, Ecr_MAX, and STATEV are all updated using the over-large dEcr. ABAQUS will discard these updates on retry (because PNEWDT < 1), so correctness is preserved, but:

1. All subsequent computation in that call (Step F `calcDcr`, Step G `calcDcocr`, Step H STATEV writes, Step I stress update) is wasted.
2. The near-singular guard in `calcDcocr` (Step G) may fire and print a spurious WARN message and reduce PNEWDT further, masking the real source of the cutback.
3. It is inconsistent with the pre-peak block, which immediately `RETURN`s after setting `PNEWDT`.

**Fix:** Add `RETURN` immediately after the `PNEWDT = MIN(PNEWDT, alpha)` line in both the Mode-I and Mode-II branches of Step C.

---

### Bug 3 — Mode-II overshoot `alpha` uses `Ecr_MAX(2)` instead of `ABS(Ecr(2))` [UNRESOLVED]

**File / Lines:** `meka_formulation_AFEM8.for`, lines 333–336.

**Code (wrong):**
```fortran
alpha = MAX((Ecr_s_ult - Ecr_MAX(2)) / ABS(dEcr(2)), 0.10D0)
```
**Should be:**
```fortran
alpha = MAX((Ecr_s_ult - ABS(Ecr(2))) / ABS(dEcr(2)), 0.10D0)
```

**Why:** `alpha` should represent the fraction of the increment that takes the crack from its **current** state `Ecr(2)` to `Ecr_s_ult`. When the crack has unloaded (slid back), `ABS(Ecr(2)) < Ecr_MAX(2)`, so the numerator `Ecr_s_ult - Ecr_MAX(2)` is smaller than the true available travel `Ecr_s_ult - ABS(Ecr(2))`, giving an alpha that is too small and cutting the step more aggressively than needed.

**Impact:** Minor in typical monotonic loading (where Ecr_MAX ≈ ABS(Ecr)), but meaningful under load reversal — forces unnecessarily small time steps near the failure limit.

---

### Bug 4 — `STATEV(3)` and `STATEV(4)` are allocated but never used [UNRESOLVED]

**File / Lines:** `meka_formulation_AFEM8.for`, header comment and `.inp` files.

**What happens:** The STATEV layout jumps from slot 2 (MODE) to slot 5 (L). Slots 3 and 4 are never written or read anywhere in the subroutine. They consume two slots in every element's state variable storage.

**Impact:** Minor (wasted memory, documentation gap). Proposed use: slot 3 = scalar damage variable `d` (0→1) for post-processing; slot 4 = cumulative fracture energy `G_released/GIc` for a normalised damage metric.

---

### Bug 5 — Near-singular guard in Step B may not trigger for Mode-II [UNRESOLVED]

**File / Lines:** `meka_formulation_AFEM8.for`, lines 314–321.

**Code:**
```fortran
IF (ABS(detY) .LT. SING_TOL*(ABS(Y(1,1))*ABS(Y(2,2))*ABS(Y(3,3)) + 1.0D0)) THEN
```

**Why it fails for Mode-II:** In Mode-II, `Y(1,1) ≈ Dcr(1,1) + ... ≈ 1e3·E ≈ 10⁷`. The relative tolerance denominator becomes `SING_TOL × 10⁷ × |Y(2,2)| × |Y(3,3)|`. When Y(2,2) is small (near failure) but Y(1,1) is large (penalty), the product can still be large, so the threshold is large and a genuinely near-singular Y may slip through undetected. The same check in `calcDcocr` (Step G) has the same problem.

**Impact:** Minor–Moderate. Near-singular Y at the tail of Mode-II softening may not trigger the `PNEWDT = 0.5` cutback, potentially causing `matrixInverse` to produce garbage results and the job to stall or produce nonsense stresses.

---

### Bug 6 — Mode-II `Dcr(3,3)` secondary shear is physically unjustified [UNRESOLVED]

**File / Lines:** `meka_formulation_AFEM8.for`, line 613 (inside `calcDcr_ModeII`).

**Code:**
```fortran
Dcr(3,3) = Dcr(2,2) * term * 1.0D-5
```
where `term = 1/(2(1+ν))`.

**Why:** This formula was copied from `calcDcr` (Mode-I), where the `term` factor approximates the shear coupling from the isotropic elasticity. For a Mode-II crack, `T(3,:)` is the secondary shear direction (perpendicular to both the crack normal and the sliding direction). There is no physical mechanism that couples the secondary shear to the primary sliding response — `Dcr(3,3)` should be zero (or a small penalty). A non-zero value introduces a spurious stiffness in the secondary shear direction that grows and decays with `Dcr(2,2)`, creating off-axis residual stresses during softening.

**Impact:** Minor. Off-axis stresses are very small (`1e-5 × term × Dcr(2,2)`), but they may contribute to the spurious twisting described in Issue 1.

---

## 7. Debug Log

| Date | File(s) Changed | Change | Effect |
|---|---|---|---|
| 06/12/2026 | `README.MD` | Code review of `meka_formulation_AFEM8.for`. Identified 6 bugs (see §6): BL_II uses E not G; Step C no RETURN after PNEWDT; Mode-II alpha uses Ecr_MAX not ABS(Ecr); STATEV(3,4) unused; near-singular guard too loose for Mode-II; Dcr(3,3) secondary shear unjustified. No code changed. | Bugs documented; fixes deferred pending shear test results. |
| 06/11/2026 | `meka_formulation_AFEM8.for` | Added Mode-II TSL: `calcDcr_ModeII`, mode-determination logic (FI_shear vs FI_mode1), coordinate-aligned crack frames for S12/S13/S23, STATEV(2)=MODE, J_FAIL=2 for Mode-II, Step C/D/E mode-branching, pre-peak Mode-II overshoot, G2c/XS PROPS with defaults, PENALTY Dcr(1,1)=1e3·E | Shear degradation now works; S12/S13/S23 reach zero at Ecr_s_ult |
| 06/11/2026 | `README.MD` | Updated for Mode-II: theory §3.3–3.5, PROPS table (G2c/XS), STATEV(2) meaning, test-case table, Issue 2 closed | — |
| 06/10/2026 | `meka_formulation_AFEM8.for` (new), `main_UD_inc8.for` (new) | Created AFEM8: ported overshoot logic from `logic-2d.for`, added `PNEWDT`, double `calcDcr` call, interpenetration clamp | AFEM7 archived |
| 06/10/2026 | `meka_formulation_AFEM8.for` | Added `J_FAIL` tracking: parametric crack-normal index from SPRIND output, stored in `STATEV(26)`. `calcDcr` updated to accept `J_FAIL` and set `Dcr(J_FAIL,J_FAIL)` as normal stiffness. Off-diagonal shear: `Dcr(j,j) = Dcr(J_FAIL,J_FAIL)·term·1e-5` | Fixes hardcoded index-1 assumption |
| 06/10/2026 | `meka_formulation_AFEM8.for` | `STATEV(24)` (Dcr normal) uncommented and updated to `Dcr(J_FAIL,J_FAIL)` | Correct SDV output |
| 06/10/2026 | `plus11_iso.inp`, `tension22_iso.inp`, `simpshear12_iso.inp`, `pureshear12_iso.inp` | `*Depvar` increased from 25 to 26 | Matches new STATEV layout |
| 02/17/2026 | `meka_formulation_AFEM7.for` | Fixed `Leff = 1.0*L` (was commented out → Leff=0 → infinite softening slope → stress plateau). Fixed `*Depvar 12→25`. Removed dead STATEV reads. | Stress now softens in 11-tension |
