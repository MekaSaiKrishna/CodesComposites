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
