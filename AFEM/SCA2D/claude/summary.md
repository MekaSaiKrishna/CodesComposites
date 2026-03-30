# Folder Summary — 2D SCA claude_updates

**Path:** `c:\ABAQUS+UMAT\Week123 [020526]\2D SCA\claude_updates\`
**Last updated:** 2026-03-30

This folder contains the ABAQUS UMAT implementation of the Smeared Crack Approach (SCA) for 2D plane strain, developed and debugged iteratively from the original `formulation_AFEM2D_PE.for`. All files implement the same constitutive model unless noted.

---

## File Index

### Entry Point (compile this)

| File | Description |
|------|-------------|
| `main_SCA2D_inc.for` | ABAQUS UMAT wrapper (`SUBROUTINE UMAT`). Routes ABAQUS arguments to `sca2diso_PE`. Include this file in every ABAQUS job. Also includes `matrix_inverse.for`. |

### Core Dependency

| File | Description |
|------|-------------|
| `matrix_inverse.for` | LU-based 2×2 and N×N matrix inverse (`matrixInverse`). Required by all formulation files. Not edited. |

### Formulation Files — Plane Strain (PE)

Listed in chronological order of creation. Each file defines `sca2diso_PE` and its helper subroutines.

| File | Date | Status | Description |
|------|------|--------|-------------|
| `formulation_AFEM2D_PE.for` | original | **Do not use** | Original unmodified PE formulation. Contains known bugs (DFGRD1 declared INTEGER, sigcr0 used before assignment, quad-precision literals, unconditional debug WRITE, wrong interpenetration fix). Reference only. |
| `claude_AFEM2D_PE.for` | 2026-02-27 | Archive | First enhanced version. Fixes 10 bugs from the original (see file header). Does not have overshoot fixes. Fix #9 (second Dcr call) commented out. Incremental scr form throughout. Use only to reproduce pre-overshoot baseline results. |
| `claude_AFEM2D_PEv1 - Copy.for` | 2026-03-03 | Archive (working snapshot) | Adds overshoot fixes OS-1 through OS-4 and FIX-GAP over `claude_AFEM2D_PE.for`. **This is the last known-working snapshot before DCB-specific modifications were made to `claude_AFEM2D_PEv1.for`.** Contains a minor bug in OS-2new (checks against `Ecr_MAX` instead of `Ecr_n_ult`). No viscous regularization, no initiation overshoot prevention. |
| `claude_AFEM2D_PEv1.for` | 2026-03-11 | Active development | Adds viscous regularization (PROPS(6)=η), initiation overshoot prevention, zero seeds at initiation, total-form softening traction, and corrects OS-2new threshold. Contains 4 residual issues addressed in v2 (see `claude_AFEM2D_PEv2.for` header). `Leff` toggleable via comment. |
| `claude_AFEM2D_PEv2.for` | 2026-03-30 | **Recommended for most use** | Targeted bug-fix release over v1. Fixes: (A) removes unused vars `J`, `Sco`; (B) adds `STATEV(4)=0` to pre-peak init; (C) removes `PNEWDT=1.0` override that disabled step control for large overshoots; (D) restores `PNEWDT` argument in `calcDcocr_PE` so singularities signal ABAQUS. All physics identical to v1. |
| `claude_AFEM2D_PEv3.for` | 2026-03-30 | **Recommended for new work** | Structural refactor of v2. `calcEcr_PE` eliminated (inlined). `calcDcr_PE` called exactly twice per cracked increment (pre-solve + post-update). One-increment STATEV lag eliminated by construction — no FIX-GAP patch needed. Same physics as v2. Cleaner code, fewer subroutine calls. |
| `formulation_AFEM2D_PEv2.for` | 2026-03-~ | **Experimental — do not use** | Parallel experiment with a different STATEV layout (sigcr0/taucr0 moved to STATEV(19:20), extra `Ecr_NEW` variable). FIX-GAP disabled. Contains a bug: `scr_old` initialised from STATEV(16:17) before those are written, yielding zero traction at initiation. No viscous regularization. |

### Documentation

| File | Description |
|------|-------------|
| `README_PE.MD` | Detailed mathematical derivation and pseudocode for the PE formulation (written after `claude_AFEM2D_PE.for`). |
| `README_PD.MD` | Documentation for the progressive damage features added in v1 (initiation overshoot prevention, zero seeds, total-form softening). |
| `fixes_overshooting.md` | Design notes for the OS-1 through OS-4 overshoot fixes introduced in `claude_AFEM2D_PEv1 - Copy.for`. |
| `summary_v3.md` | Technical summary of v3: structural changes, lag-elimination proof, PROPS/STATEV tables, constitutive logic, PNEWDT rules, and known limitations. |
| `summary.md` | This file. |

---

## When to Use Which File

### For a standard monotonic loading simulation (DCB, ENF, CT specimen)

**Use `claude_AFEM2D_PEv3.for` + `main_SCA2D_inc.for` + `matrix_inverse.for`.**

Update `main_SCA2D_inc.for` to `INCLUDE` v3 instead of whichever formulation file it currently points to. v3 has the cleanest implementation, no known bugs, and the lag-free STATEV layout.

### For a simulation requiring viscous regularization (convergence issues in softening)

**Use `claude_AFEM2D_PEv3.for`.** Set `PROPS(6) = η` (typical: `1e-4 × Δt_typical`). Start with `η = 0` and increase if Newton-Raphson diverges repeatedly in the post-peak regime.

### To reproduce results from before 2026-03-10 (pre-DCB work)

**Use `claude_AFEM2D_PEv1 - Copy.for`.** This is the last snapshot before DCB-specific changes (hardcoded Leff, viscous term) entered the active development file.

### To reproduce the exact results from `claude_AFEM2D_PEv1.for` runs

**Use `claude_AFEM2D_PEv1.for`.** Note: if the run used `Leff = 0.02` (hardcoded, now commented out), restore that line manually before running.

### To reproduce the original (buggy) behaviour for comparison

**Use `formulation_AFEM2D_PE.for`.** Do not use this for production.

### Never use in production

- `formulation_AFEM2D_PE.for` — known bugs
- `formulation_AFEM2D_PEv2.for` — experimental STATEV reorg with a zero-traction initiation bug and FIX-GAP disabled

---

## How to Swap Formulation Files

`main_SCA2D_inc.for` uses an INCLUDE or CALL statement to pull in the formulation. To switch versions, locate the include/call line in `main_SCA2D_inc.for` and update the filename. The subroutine signature of `sca2diso_PE` is identical across all working versions (v1-Copy, v1, v2, v3).

---

## STATEV Layout (v2 and v3)

```
STATEV(1)  isCRACKED   0 = intact, 1 = active
STATEV(2)  MODE        +1 = tensile, -1 = compressive
STATEV(3)  Leff        stored element size [m]
STATEV(4)  FI_MAX      peak failure index (pre-peak)
STATEV(5)  cos θ       crack normal cosine
STATEV(6)  sin θ       crack normal sine
STATEV(7)  Ecr_MAX(1)  max normal crack strain
STATEV(8)  Ecr_MAX(2)  max shear crack strain
STATEV(9)  Ecr(1)      current normal crack strain
STATEV(10) Ecr(2)      current shear crack strain
STATEV(11) dEcr(1)     normal crack strain increment
STATEV(12) dEcr(2)     shear crack strain increment
STATEV(13) σ_init      max principal stress at initiation
STATEV(14) scr(1)      current normal crack traction
STATEV(15) scr(2)      current shear crack traction
STATEV(16) σ_cr0       normal traction at initiation
STATEV(17) τ_cr0       shear traction at initiation
```

Minimum NSTATV = 17 for all working versions.

Note: `formulation_AFEM2D_PEv2.for` uses a different layout (sigcr0 at STATEV(19:20)) and is incompatible with the above.
