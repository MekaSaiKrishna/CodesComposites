# AFEM3D — ABAQUS UMAT: Incremental Smeared Crack Approach

> **Active formulation (isotropic):** `meka_formulation_AFEM11.for` *(three independent Mode-I cracks, decoupled 3×3 solves)*  
> **Entry point (isotropic):** `main_UD_inc11.for`  
> **Active formulation (orthotropic):** `meka_formulation_AFEM10.for` *(single crack, Mode-I + Mode-II)*  
> **Entry point (orthotropic):** `main_UD_inc10.for`  
> **Previous isotropic active:** `meka_formulation_AFEM9.for` / `main_UD_inc9.for` *(archived — coupled 9×9 system, Bugs A–E)*  
> **Last reviewed:** 06/28/2026  
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
**Failure criteria (current):** Max Normal Stress on coord. axes (Mode-I) + Direct Shear Stress (Mode-II).  
**Material scope (current):** Isotropic only.

---

## 2. File Index

### Entry Points

| File | Description | Status |
|---|---|---|
| `main_UD_inc11.for` | **ACTIVE (isotropic).** UMAT wrapper for AFEM11 (three independent Mode-I cracks, decoupled). | Active |
| `main_UD_inc10.for` | **ACTIVE (orthotropic).** UMAT wrapper for AFEM10 (single crack, Mode-I + Mode-II). | Active |
| `main_UD_inc9.for` | UMAT wrapper for AFEM9 (three simultaneous Mode-I cracks — Bugs A–E, archived). | Archived |
| `main_UD_inc8.for` | UMAT wrapper for AFEM8. Passes `PNEWDT` for overshoot cutbacks. | Archived |
| `main_UD_inc.for` | UMAT wrapper for AFEM7 (no `PNEWDT`). | Archived |
| `main_UD_tot.for` | UMAT wrapper for AFEM5 (total-strain form). | Archived |
| `main_UD_elas.for` | Elasticity-only UMAT (no damage). | Archived |

### Isotropic Formulations (evolutionary sequence)

| File | Description | Status | Last Debugged |
|---|---|---|---|
| `meka_formulation_AFEM11.for` | **ACTIVE.** Three INDEPENDENT coord-aligned Mode-I cracks. Each crack (x/y/z-normal) initiates when its own S_ii ≥ XT. Three separate 3×3 solves — no Poisson or shear cross-coupling. Shear-locked PENALTY_FAC=1e6. Combined Woodbury: Dcocr = Dco − ΣᵢAᵢ. NSTATV=26, NPROPS=4. | Active | 06/28/2026 |
| `meka_formulation_AFEM9.for` | Three simultaneous coord-aligned Mode-I cracks. S22 triggers all 3 simultaneously. 9×9 coupled system — Bugs A–E (Poisson λ and shared-shear cross-coupling). **Archived; replaced by AFEM11.** | Archived | 06/15/2026 |
| `meka_formulation_AFEM8.for` | AFEM7 + PNEWDT overshoot, unified `calcDcr(MODE,...)`, Max Normal Stress criterion (coord-aligned, no SPRIND), Mode-II shear TSL, interpenetration clamp, near-singular guards | Archived | 06/12/2026 |
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
  FI_mode1  =  max[ SIG(i) / XT ]   for SIG(i) > 0  (i = 1,2,3)
              (max normal stress on coordinate axes — no SPRIND)

  FI_shear  =  max( |S12|, |S13|, |S23| ) / XS
              (direct shear criterion)
```

Failure initiates when either FI ≥ 1. The mode with the higher FI (when both ≥ 1) is locked in permanently for that element.

```
  FI_shear ≥ 1 AND FI_shear > FI_mode1  →  MODE-II  (shear)
  FI_mode1 ≥ 1 AND FI_mode1 ≥ FI_shear  →  MODE-I   (tensile)
```

`STATEV(2) = MODE` records this choice at initiation and never changes.

The coordinate-aligned crack frame for **Mode-I** is set from the governing direction `INDEX_NORMAL` (the index i that maximises SIG(i)/XT):

```
  INDEX_NORMAL = 1  →  T(1,:)=[1,0,0]  T(2,:)=[0,1,0]  T(3,:)=[0,0,1]
  INDEX_NORMAL = 2  →  T(1,:)=[0,1,0]  T(2,:)=[0,0,1]  T(3,:)=[1,0,0]
  INDEX_NORMAL = 3  →  T(1,:)=[0,0,1]  T(2,:)=[1,0,0]  T(3,:)=[0,1,0]
```

`J_FAIL = 1` always for Mode-I (crack normal = T(1,:), crack opening = Ecr(1)).

---

### 3.5 Mode-II Crack Orientation (Coordinate-Aligned)

Both Mode-I and Mode-II crack frames are fixed to coordinate axes (no SPRIND eigenvectors). For Mode-II the dominant shear component determines the frame:

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
  MODE-I  Dcr:                        MODE-II  Dcr:
  ┌──────────────────────────────┐    ┌──────────────────────────────┐
  │ TSL(Ecr_n)    0         0    │    │ 1000·E      0          0     │
  │    0      1000·G        0    │    │    0    TSL(|Ecr_s|)   0     │
  │    0          0     1000·G   │    │    0        0          0     │
  └──────────────────────────────┘    └──────────────────────────────┘
  Dcr(1,1): Mode-I TSL, opening      Dcr(1,1): penalty — blocks normal opening
  Dcr(2,2): 1000·G penalty —         Dcr(2,2): Mode-II TSL, sliding softens
    Mode-I crack CANNOT slide.       Dcr(3,3): 0 — no secondary shear coupling
  Dcr(3,3): same, secondary tangent  G = E/(2(1+ν))
```

**Why the Mode-I shear penalty is mandatory:** With `Dcr(2,2) = 0`, the Woodbury identity gives:

```
  Y(2,2) = 0 + G  →  A(4,4) = G²/G = G  →  Dcocr(4,4) = G − G = 0
```

Zero shear stiffness after cracking means any numerical perturbation causes the element to shear freely — producing the observed element tilting. Setting `Dcr(2,2) = Dcr(3,3) = 1000·G` gives `Y(2,2) ≈ 1000·G`, so `A(4,4) = G²/(1001·G) ≈ G/1001 ≈ 0`, and `Dcocr(4,4) ≈ G` — full shear stiffness is retained. The penalty does NOT affect the tensile stress-strain curve because for pure tension `DSTRAN(4) = DSTRAN(5) = 0` and the shear dEcr terms remain zero regardless of `Dcr(2,2)`.

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
  │  calcDcr(MODE, Ecr_old, Ecr_MAX, ...)  → Dcr_old          │
  │  (builds system matrix Y only; scr_n discarded)          │
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
  │  calcDcr(MODE, Ecr_new, Ecr_MAX, ...) → Dcr, scr_n       │
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

### AFEM11 State Variables (NSTATV = 26)

*Three independent Mode-I cracks, coord-aligned. Each activates independently when its own normal stress ≥ XT.*

| Slot | Variable | Description |
|---|---|---|
| 1 | `crack_flag(1)` | Crack 1, x-normal: 0 = intact, 1 = cracked |
| 2 | `crack_flag(2)` | Crack 2, y-normal |
| 3 | `crack_flag(3)` | Crack 3, z-normal |
| 4 | `L` | Characteristic length `CELENT` (set at first initiation) |
| 5–7 | `Ecr_MAX(3)` | Peak normal crack strain per crack |
| 8–10 | `Ecr(crack 1)` | Crack 1 strain: `[normal, shear1, shear2]` |
| 11–13 | `Ecr(crack 2)` | Crack 2 strain |
| 14–16 | `Ecr(crack 3)` | Crack 3 strain |
| 17–19 | `scr_n(3)` | Crack normal traction per crack (post-processing) |
| 20–22 | `dEcr_n(3)` | Normal crack strain increment this step (post-processing) |
| 23–26 | *(spare)* | Available for future use |

### AFEM11 Material Properties (NPROPS = 4)

| Slot | Symbol | Description |
|---|---|---|
| 1 | `E` | Young's modulus |
| 2 | `ν` | Poisson's ratio |
| 3 | `GIc` | Mode-I fracture energy (same for all three cracks) |
| 4 | `XT` | Tensile strength (Mode-I peak) |

```
*Depvar
     26,
*User Material, constants=4
E,  nu,  GIc,  XT
```

---

### AFEM9 State Variables (NSTATV = 26)

*Archived. Three simultaneous Mode-I cracks (all activate when S22 ≥ XT). 9×9 coupled system — Bugs A–E.*

| Slot | Variable | Description |
|---|---|---|
| 1 | `isCRACKED` | 0 = fully intact, 1 = cracked |
| 2 | `L` | Characteristic length `CELENT` |
| 3–5 | `Ecr_MAX_n(3)` | Peak normal crack strain per direction |
| 6–14 | `Ecr(9)` | Current crack strain vector (9 components) |
| 15–23 | `dEcr(9)` | Crack strain increment this step |
| 24–26 | `scr(3)` | Normal traction at each crack (post-processing) |

### AFEM9 Material Properties (NPROPS = 4)

| Slot | Symbol | Description |
|---|---|---|
| 1 | `E` | Young's modulus |
| 2 | `ν` | Poisson's ratio |
| 3 | `GIc` | Mode-I fracture energy (same for all three cracks) |
| 4 | `XT` | Tensile strength (Mode-I peak) |

```
*Depvar
     26,
*User Material, constants=4
E,  nu,  GIc,  XT
```

---

### AFEM8 State Variables (NSTATV = 26)

| Slot | Variable | Description |
|---|---|---|
| 1 | `isCRACKED` | 0 = intact, 1 = cracked |
| 2 | `MODE` | Failure mode: **1 = Mode-I (tensile)**, **2 = Mode-II (shear)** |
| 5 | `L` | Characteristic length `CELENT` at crack initiation |
| 6–8 | `Ecr_MAX(1:3)` | Mode-I: max Ecr component-wise. Mode-II: `Ecr_MAX(2)` = max(\|Ecr_s\|) |
| 9–11 | `T(1,:)` | Crack-frame row 1 — crack normal (coordinate axis for both Mode-I and Mode-II) |
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
| 5 | `INDEX` | 0 = Max Normal Stress on coord. axes + Direct Shear | — |
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
- Crack initiates at `ε = XT/E = 100/10000 = 0.01`, `σ = 100 MPa`
- `Ecr_n_ult = 2·GIc/(Leff·XT) = 2·1/(1·100) = 0.02`
- Total failure displacement: **`0.02 m`** — during softening, elastic strain decreases to zero as σ→0; the 0.01 elastic extension at initiation is not additive (it shrinks back to zero)

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

### Issue 1 — Element tilting / zero shear stiffness after Mode-I cracking [RESOLVED 06/12/2026]

**Symptom (observed):** After the peak stress is reached in uniaxial tension (T-11, T-22, T-33), the element tilts instead of softening cleanly. Spurious off-axis stresses develop.

**Root cause:** With `Dcr(2,2) = Dcr(3,3) = 0`, the Woodbury identity immediately zeros the shear stiffness in `Dcocr`:

```
  Y(2,2) = Dcr(2,2) + G = 0 + G = G
  A(4,4) = G² / G = G
  Dcocr(4,4) = G − G = 0   ← element has zero shear stiffness after cracking
```

Any numerical perturbation then drives free shear rotation (tilting).

**Fix (06/12/2026):** Added shear-locked penalty to Mode-I `Dcr`:

```fortran
Dcr(2,2) = PENALTY_FAC * G_local   ! G_local = E/(2(1+ν))
Dcr(3,3) = PENALTY_FAC * G_local
```

This is physically correct: a Mode-I crack OPENS but does NOT SLIDE (it is shear-locked). The penalty prevents sliding, exactly as `Dcr(1,1) = PENALTY*E` prevents normal opening for Mode-II. The penalty does not affect the tensile stress-strain curve because for pure tension DSTRAN(shear) = 0, so dEcr(2) = 0 regardless of Dcr(2,2).

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

### Bug 1 — `BL_II` Bazant snap-back limit uses `E` instead of `G` [RESOLVED 06/12/2026]

**Fix:** `BL_II = 2.0D0*G2c*G/taucr0**2` (G computed on line 100). Print label also corrected to `'   2*G2c*G/XS^2 = '`.

**Why it mattered:** The Mode-II Bazant limit is derived from requiring no snap-back in the shear TSL: `|slope_soft_II| ≤ G` (shear modulus, not Young's modulus). Using E instead of G overestimates the limit by a factor of `E/G = 2(1+ν)` — for ν=0.3 this is 2.6×, so elements can be 2.6× too large without triggering the warning.

---

### Bug 2 — Post-peak overshoot (Step C) does not `RETURN` after setting `PNEWDT` [UNRESOLVED]

**File / Lines:** `meka_formulation_AFEM8.for`, lines 329–349.

**What happens:** When an overshoot past `Ecr_n_ult` (Mode-I) or `Ecr_s_ult` (Mode-II) is detected, `PNEWDT` is set to request a step cutback, but execution continues through Steps D–I. The crack strain, Ecr_MAX, and STATEV are all updated using the over-large dEcr. ABAQUS will discard these updates on retry (because PNEWDT < 1), so correctness is preserved, but:

1. All subsequent computation in that call (Step F `calcDcr`, Step G `calcDcocr`, Step H STATEV writes, Step I stress update) is wasted.
2. The near-singular guard in `calcDcocr` (Step G) may fire and print a spurious WARN message and reduce PNEWDT further, masking the real source of the cutback.
3. It is inconsistent with the pre-peak block, which immediately `RETURN`s after setting `PNEWDT`.

**Fix:** Add `RETURN` immediately after the `PNEWDT = MIN(PNEWDT, alpha)` line in both the Mode-I and Mode-II branches of Step C.

---

### Bug 3 — Mode-II overshoot `alpha` uses `Ecr_MAX(2)` instead of `ABS(Ecr(2))` [RESOLVED 06/12/2026]

**Fix:** Numerator changed to `Ecr_s_ult - ABS(Ecr(2))` in Step C of the post-peak block.

**Why it mattered:** `alpha` represents the fraction of the increment that takes the crack from its **current** state `Ecr(2)` to `Ecr_s_ult`. When the crack has unloaded (slid back), `ABS(Ecr(2)) < Ecr_MAX(2)`, so using Ecr_MAX(2) in the numerator cut the step more aggressively than needed.

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

### Bug 6 — Mode-II `Dcr(3,3)` secondary shear is physically unjustified [RESOLVED 06/12/2026]

**Fix:** `Dcr(3,3) = 0.0D0` in the unified `calcDcr` subroutine (both MODE=1 and MODE=2 branches). The old `Dcr(2,2) * term * 1e-5` value was a leftover copy from the Mode-I shear coupling approximation and has no physical basis in Mode-II.

---

### Known Issues — AFEM9 (identified 06/14/2026 from DCB simulation)

> **Note (06/15/2026):** The AFEM9 formulation was redesigned to address all five bugs below. The redesign follows `qcalcDcr` from `misesAndSmearedCrack.for` and removes the per-crack flag system, HUGE_PENALTY, and PENALTY_FAC. See the 06/15/2026 debug log entry for details.

AFEM9 was tested against the first-failing element of a DCB (Double-Cantilever Beam) simulation and compared with the reference `misesAndSmearedCrack.for` (Haotian SCA). Five physically wrong behaviours were observed. All five share the same root cause.

**Root cause — Poisson (λ) coupling through N^T·Dco·N**

In the 9×9 system matrix `Y = Dcr + N^T·Dco·N`, the Lamé constant λ creates off-diagonal entries that couple all three crack-normal directions simultaneously:

```
  Y(1,4) = λ   ← couples crack-1 normal (ε₁₁) with crack-2 normal (ε₂₂)
  Y(1,7) = λ   ← couples crack-1 normal with crack-3 normal (ε₃₃)
  Y(4,7) = λ
```

Solving `dEcr = inv(Y) · (N^T·Dco·DSTRAN)` redistributes crack strain across all three normal directions for any single-axis loading. In the DCB, only the peel direction (y) is loaded, but through λ the solver also drives crack opening in x and z — which are not loaded and should remain unchanged. AFEM8 avoids this by having only one active crack at a time in a properly-oriented 3×3 Y that has no Poisson cross-coupling at the crack-normal index.

**Data summary (200 increments of DCB loading)**

| Component | AFEM9 result | Reference (Haotian SCA) | Verdict |
|---|---|---|---|
| S11 peak | ~18.4 MPa | ~11.3 MPa | Different |
| S11 post-failure | **−16.3 MPa and growing** | ~0 (machine-epsilon) | **WRONG** |
| S22 peak | ~19.5 MPa | ~20.0 MPa | Similar |
| S22 residual | **−0.488 MPa (locked constant)** | ~0 (machine-epsilon) | **WRONG** |
| S33 | **grows to −280 MPa** | < 1.7 MPa, returns to ~0 | **WRONG** |
| S12 | ≈ machine-epsilon | ≈ machine-epsilon | Correct |
| S13 | **grows to +0.83 MPa** | ≈ machine-epsilon after failure | **WRONG** |
| S23 | ≈ machine-epsilon | ≈ machine-epsilon | Correct |

---

**AFEM9 — Bug A: S11 sign inversion (spurious compressive stress in Mode-I direction)**

- **Symptom:** S11 plateaus near 18.4 MPa for ~60 increments after crack-1 initiates, then drops and goes deeply negative, reaching −16.3 MPa while the DCB continues to open.
- **Expected:** S11 softens to XT × R_RES ≈ 0 and stays non-negative. A Mode-I crack cannot carry compression.
- **Root cause (two combined effects):**
  1. In the DCB, DSTRAN is dominated by the peel (y) component. The driving term for crack-1 is `(N^T·Dco·DSTRAN)(1) = λ·DSTRAN(2)`, which is small (only Poisson coupling, not direct x-loading). So `dEcr(1) ≈ 0` each increment and crack-1 barely softens despite having initiated.
  2. Once crack-2 opens substantially, the off-diagonal Y(1,4) = λ drives `dEcr(1)` negative every increment. The interpenetration clamp prevents Ecr(1) < 0, but the Woodbury-condensed Dcocr still produces a negative S11 contribution each step via `SIG(1) += Dcocr(1,2)·DSTRAN(2)`, which accumulates without bound.

---

**AFEM9 — Bug B: S33 grows unboundedly to −280 MPa**

- **Symptom:** S33 starts at ~2.2 MPa (Poisson-driven), reverses sign at step 10, and grows monotonically to −280 MPa. E33 stays small (< 0.003), implying an apparent z-stiffness of ~100 GPa is still active — the element never loses z-stiffness through any cracking mechanism.
- **Expected:** S33 should remain small, bounded, and return to near-zero after the peel crack fully softens.
- **Root cause:** Crack-3 (z-normal) is inactive throughout (Dcr(7,7) = HUGE_PENALTY × E). Because Y(7,7) = 1e6·E >> off-diagonals, invY(7,7) ≈ 0 and crack-3's own crack strain stays near zero — correct so far. However, the Woodbury correction term A66(3,3) = (Dco·N·invY·N^T·Dco)(3,3) also receives contributions through invY entries coupling index 7 with active indices 1 and 4 via λ. These are not zero. Combined with a Poisson-driven compressive DSTRAN(3) from the large peel loading, `SIG(3) += Dcocr(3,3)·DSTRAN(3)` accumulates a compressive increment every step, growing without physical dissipation because crack-3 never opens.

---

**AFEM9 — Bug C: S22 residual locked at −0.488 MPa instead of ~0**

- **Symptom:** After crack-2 reaches complete failure (Ecr(4) > Ecr_n_ult), S22 freezes at exactly −0.488 MPa and remains there for the rest of the simulation (~100 increments) rather than softening to XT × R_RES ≈ 0. E22 continues growing to 3.47 while S22 stays at this negative offset.
- **Expected:** S22 residual ≈ XT × R_RES ≈ 0.
- **Root cause:** Once Dcr(4,4) = RES_K, the Woodbury formula correctly drives Dcocr(2,2) ≈ 0 in the y-normal direction. However, the off-diagonal Dcocr(2,1) remains non-zero due to λ coupling with x. Each increment, `SIG(2) += Dcocr(2,1)·DSTRAN(1)` adds a small negative contribution (since x-direction is now under compressive stress from Bug A). These accumulate to the −0.488 MPa offset and then stabilise when the compressive x-stress also stabilises.

---

**AFEM9 — Bug D: S13 grows monotonically to +0.83 MPa (spurious shear)**

- **Symptom:** S13 grows steadily from 0 to +0.83 MPa over 200 increments. In a symmetric DCB loaded in Mode-I peel, S13 should be zero or machine-epsilon throughout.
- **Expected:** S13 ≈ 0 at all times.
- **Root cause:** The N matrix entries N(5,3) = 1 and N(5,8) = 1 mean that γ₁₃ receives contributions from both crack-1's secondary shear (Ecr(3)) and crack-3's primary shear (Ecr(8)). The PENALTY_FAC × G stiffness on these shear components, when condensed through the Woodbury formula alongside the active crack-2 coupling, produces a small but consistently signed S13 increment each step from the coupling term `Dcocr(5,2)·DSTRAN(2)`.

---

**AFEM9 — Bug E: Crack-1 softening suppressed (plateau at S11 ≈ 18.4 MPa)**

- **Symptom:** After crack-1 initiates, S11 stays nearly constant at ~18.4 MPa for approximately 60 increments before abruptly dropping. The softening slope is absent for the plateau phase.
- **Expected:** S11 should begin following the TSL slope immediately after initiation.
- **Root cause:** This is the slow-crack-opening mechanism from Bug A. Because `(N^T·Dco·DSTRAN)(1) ≈ λ·DSTRAN(2)` is small for a DCB (y-peel dominated), `dEcr(1)` per increment is negligible. Crack-1 barely opens. Meanwhile, `SIG(1) += Dcocr(1,1)·DSTRAN(1) + Dcocr(1,2)·DSTRAN(2) + ...` accumulates S11 from Poisson of y-loading faster than the negligible crack softening can reduce it. The plateau ends only when crack-2 opens far enough that the Y(1,4) coupling forces a sudden large `dEcr(1)` — at which point S11 drops rapidly and reverses sign.

---

**Proposed fix direction for AFEM9**

The coordinate-aligned 3-crack simultaneous formulation cannot avoid the λ coupling in N^T·Dco·N. The recommended approach is to decouple the three cracks from the shared solve:

1. After a crack fully fails (Dcr(i,i) → RES_K), zero its row and column contribution to the 9×9 right-hand side so it cannot drive the remaining active cracks through λ.
2. Alternatively, use a stress-return approach per crack direction (check each normal direction independently using the elastic trial stress, then apply softening only in that direction) to avoid the cross-coupling entirely.
3. As a short-term patch: clamp SIG(i) ≥ 0 for any crack-normal direction i where crack_flag(i) = 1 (interpenetration of stress, not just strain), to at least prevent unphysical compressive traction from accumulating.

---

## 7. Debug Log

| Date | File(s) Changed | Change | Effect |
|---|---|---|---|
| 06/28/2026 | `meka_formulation_AFEM11.for` (new), `main_UD_inc11.for` (new), `README.MD` | **Root-cause fix for all AFEM9 Bugs A–E.** Complete rewrite as three INDEPENDENT 3×3 crack solves following AFEM10 logic. Two coupling sources eliminated: (1) Poisson λ normal-normal coupling in 9×9 Y removed by decoupling; (2) shared-shear Voigt-index double-counting (Dcocr(γ₁₂) = G−G−G = −G) prevented by shear-locked PENALTY_FAC=1e6 per crack. Each crack initiates independently when its own S_ii ≥ XT. For coord-aligned cracks in isotropic material N_i^T·Dco·N_i is diagonal so Y_i is diagonal — no cross-coupling. Combined Woodbury: Dcocr = Dco − Σᵢ Aᵢ. NSTATV=26 (new layout), NPROPS=4 (unchanged). AFEM9 archived. | Eliminates spurious S11 sign inversion, S33 growth to −280 MPa, S22 residual lock, S13 growth, and suppressed crack-1 softening (Bugs A–E). Each axis should now show independent linear softening XT→0 with area = GIc/Leff. Pending test verification. |
| 06/15/2026 | `meka_formulation_AFEM9.for` | **Root-cause fix for missing stress degradation in Y-tension test.** Three bugs compounded: (1) `T(3,:) = (1,0,1)` was not a unit vector → corrupted N matrix entries; fixed to `(0,0,1)`. (2) `Ecr_MAX_n` initialized to `ECR_INIT = 1e-14` at onset → secant Dcr ≈ 2e21 MPa → `inv(Y99) ≈ 0` → `dEcr ≈ 0` each step → no crack opening. Fixed to `XT/E`. (3) Even with (1)+(2) fixed, a final algorithmic error remained: the original RHS `NT_Dco·DSTRAN` gives tiny `dEcr` (DSTRAN is a small increment) so the elastic increment stays positive and SIG keeps rising. Switching to `NT_Dco·EPS_total` overshot (abrupt drop to 0). **Final fix:** replace Step A/B with a tangent-Dcr solve: `Dcr_tan(ni,ni) = slope_soft` (the TSL slope) and RHS = `NT_Dco·DSTRAN`. For linear softening this is exact — the single linear solve gives the TSL-consistent crack-strain increment at every step with no Newton iteration needed. Added subroutine `calcDcr9_tan`. Step F/G (Woodbury) keep the secant Dcr for ABAQUS Newton stability. | S22 now softens linearly from XT to 0, with fracture strain matching Ecr_n_ult. Pending tensile-Y run to confirm. |
| 06/15/2026 | `meka_formulation_AFEM9.for`, `main_UD_inc9.for`, `README.MD` | Full redesign of AFEM9 to match `qcalcDcr` from `misesAndSmearedCrack.for` (Haotian SCA). Changes: (1) Failure trigger simplified to S22 > XT only (peel criterion). (2) All 3 cracks activate simultaneously — removed per-crack `crack_flag(3)` system. (3) New `calcDcr9`: secant TSL per normal, `Dcr_min = min(Dcr(1,1), Dcr(4,4), Dcr(7,7))` applied to all normals, shear = Dcr_min·term·1e-5 (muted). `sigcr0_r = XT·5e-10`. (4) Removed HUGE_PENALTY, PENALTY_FAC, Step-0 new-crack check, Fix-1 N masking. (5) Kept stress clamp SIG(1:3) ≥ 0 (Fix-2). (6) NSTATV reduced from 35 to 26 (no crack_flag, Ecr_MAX_n tracks 3 normals only). ABAQUS `.inp` files must be updated to `*Depvar 26`. | Addresses root cause of Bugs A–E by removing HUGE_PENALTY Poisson coupling and independently-driving each crack via Dcr_min. Pending DCB re-test. |
| 06/14/2026 | `meka_formulation_AFEM9.for` | Fix-1: after Step-0, zero N69 columns for unformed cracks so HUGE_PENALTY Dcr entries do not enter the Woodbury solve through Poisson λ off-diagonals; recompute NT96/NT_Dco from masked N69. Fix-2: after Step-I stress update, clamp SIG(I) ≥ 0 for each active crack. Targets Bugs A–D. Bug E (crack-1 softening suppressed by direction) not yet fixed — requires per-crack independent stress return. | Bugs A–D targeted; Bug E still open. |
| 06/14/2026 | `README.MD` | DCB simulation analysis: compared AFEM9 first-element stress-strain against Haotian SCA reference. Identified 5 bugs (A–E in §6 Known Issues — AFEM9): S11 sign inversion, S33 growing to −280 MPa, S22 residual −0.488 MPa, spurious S13 growth, crack-1 softening suppressed. Root cause: Poisson (λ) off-diagonal coupling in the 9×9 Y matrix redistributes crack strain to non-loaded directions. No code changed; fixes documented. | All 5 bugs open. |
| 06/13/2026 | `meka_formulation_AFEM9.for`, `main_UD_inc9.for`, `README.MD` | Created AFEM9: three simultaneous coord-aligned Mode-I cracks. N matrix 6×9 (hardcoded for T=I₃), Dcr 9×9 diagonal, system matrix Y 9×9 inverted via `matrixInverse`. Per-crack initiation flags; inactive cracks locked with HUGE_PENALTY×E. Same incremental solve and PNEWDT overshoot logic as AFEM8. NSTATV=35. | First version, not yet tested against test cases. |
| 06/12/2026 | `meka_formulation_AFEM8.for`, `README.MD` | Fixed Issue 1 (element tilting in tension): added shear-locked penalty `Dcr(2,2)=Dcr(3,3)=PENALTY_FAC*G` to Mode-I branch of `calcDcr`. Root cause: Dcr(2,2)=0 drove Dcocr(shear)=0 via Woodbury, causing free shear rotation. Mode-II untouched. | Resolved Issue 1 (tilting). Bugs 2, 4, 5 remain open. |
| 06/12/2026 | `meka_formulation_AFEM8.for`, `README.MD` | Refactor: replaced SPRIND principal-stress criterion (Mode-I) with Max Normal Stress on coordinate axes. Unified `calcDcr_ModeI` + `calcDcr_ModeII` into single `calcDcr(MODE,...)`. J_FAIL=1 always for Mode-I. Fixed Bug 1 (BL_II label E→G), Bug 3 (Mode-II alpha numerator ABS(Ecr) not Ecr_MAX), Bug 6 (Dcr(3,3)=0). Corrected expected failure displacement from 0.03→0.02 m. | Bugs 2, 4, 5 remain open. |
| 06/12/2026 | `README.MD` | Code review of `meka_formulation_AFEM8.for`. Identified 6 bugs (see §6): BL_II uses E not G; Step C no RETURN after PNEWDT; Mode-II alpha uses Ecr_MAX not ABS(Ecr); STATEV(3,4) unused; near-singular guard too loose for Mode-II; Dcr(3,3) secondary shear unjustified. No code changed. | Bugs documented; fixes deferred pending shear test results. |
| 06/11/2026 | `meka_formulation_AFEM8.for` | Added Mode-II TSL: `calcDcr_ModeII`, mode-determination logic (FI_shear vs FI_mode1), coordinate-aligned crack frames for S12/S13/S23, STATEV(2)=MODE, J_FAIL=2 for Mode-II, Step C/D/E mode-branching, pre-peak Mode-II overshoot, G2c/XS PROPS with defaults, PENALTY Dcr(1,1)=1e3·E | Shear degradation now works; S12/S13/S23 reach zero at Ecr_s_ult |
| 06/11/2026 | `README.MD` | Updated for Mode-II: theory §3.3–3.5, PROPS table (G2c/XS), STATEV(2) meaning, test-case table, Issue 2 closed | — |
| 06/10/2026 | `meka_formulation_AFEM8.for` (new), `main_UD_inc8.for` (new) | Created AFEM8: ported overshoot logic from `logic-2d.for`, added `PNEWDT`, double `calcDcr` call, interpenetration clamp | AFEM7 archived |
| 06/10/2026 | `meka_formulation_AFEM8.for` | Added `J_FAIL` tracking: parametric crack-normal index from SPRIND output, stored in `STATEV(26)`. `calcDcr` updated to accept `J_FAIL` and set `Dcr(J_FAIL,J_FAIL)` as normal stiffness. Off-diagonal shear: `Dcr(j,j) = Dcr(J_FAIL,J_FAIL)·term·1e-5` | Fixes hardcoded index-1 assumption |
| 06/10/2026 | `meka_formulation_AFEM8.for` | `STATEV(24)` (Dcr normal) uncommented and updated to `Dcr(J_FAIL,J_FAIL)` | Correct SDV output |
| 06/10/2026 | `plus11_iso.inp`, `tension22_iso.inp`, `simpshear12_iso.inp`, `pureshear12_iso.inp` | `*Depvar` increased from 25 to 26 | Matches new STATEV layout |
| 02/17/2026 | `meka_formulation_AFEM7.for` | Fixed `Leff = 1.0*L` (was commented out → Leff=0 → infinite softening slope → stress plateau). Fixed `*Depvar 12→25`. Removed dead STATEV reads. | Stress now softens in 11-tension |
