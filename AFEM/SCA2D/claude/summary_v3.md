# claude_AFEM2D_PEv3.for — Technical Summary

**Date:** 2026-03-30
**Based on:** claude_AFEM2D_PEv2.for
**Subroutine:** `sca2diso_PE`
**Model:** Isotropic incremental SCA, 2D plane strain (CPE4)

---

## 1. Purpose

v3 is a structural refactor of v2. Every fix accumulated since the original file is retained, but the code is reorganised so that the one-increment STATEV lag **cannot occur by design**, rather than being patched after the fact. One subroutine (`calcEcr_PE`) is eliminated, reducing the call graph and removing a redundant internal `calcDcr_PE` call.

---

## 2. The One-Increment Lag — Root Cause and Elimination

### What the lag was

In every pre-v3 version, the post-peak block had the following pattern:

```
1. Call calcDcr_PE(Ecr_OLD)  → Dcr_old, scr_pred   [STATEV(14:15) = scr_pred]
2. Call calcEcr_PE(...)      → internally calls calcDcr_PE(Ecr_OLD) again
                               solves dEcr, updates Ecr_OLD → Ecr_new
                               STATEV(9:10) = Ecr_new
3. FIX-GAP: Call calcDcr_PE(Ecr_new) → Dcr_new, scr_new
                               STATEV(14:15) = scr_new   (overwrite step 1)
4. calcDcocr_PE(Dcr_new)
```

The FIX-GAP (step 3) was added because without it:
- `STATEV(9:10)` held `Ecr_new` (from step 2)
- `STATEV(14:15)` held `scr_pred` computed from `Ecr_OLD` (from step 1)
- Plotting `STATEV(14)` vs `STATEV(9)` showed a one-increment horizontal offset

Even with FIX-GAP, `calcDcr_PE` was called **three times** per cracked increment:
- Once in Step 1 (for the first STATEV write and to get scr_pred)
- Once inside `calcEcr_PE` (to form the Y matrix)
- Once as FIX-GAP (final correct state)

### How v3 eliminates it

v3 replaces the three-call pattern with an explicit two-call pattern:

```
Step A: Call calcDcr_PE(Ecr)       → Dcr_old   [scr output discarded]
Step B: Solve Y·dEcr = N^T Dco DSTRAN  [Y = Dcr_old + N^T Dco N]
Step C: Overshoot check, set PNEWDT if needed
Step D: Ecr_new = Ecr + dEcr;  Ecr_MAX updated
Step E: Interpenetration clamp (Ecr_new(1) clamped to 0 if negative)
Step F: Call calcDcr_PE(Ecr_new)   → Dcr_new, scr_new
Step G: calcDcocr_PE(Dcr_new)
Step H: Write STATEV(7:12) and STATEV(14:15) — ALL from same Ecr_new
```

Because `STATEV(9:10) = Ecr_new` and `STATEV(14:15) = scr_new` are both written in Step H from the same `Ecr_new`, the lag is structurally impossible. The FIX-GAP is not a patch — it is the natural architecture.

`calcEcr_PE` is eliminated entirely. The Y-system solve (3 lines) and history update are inlined into the main subroutine, where they are clearer and produce no extra call overhead.

---

## 3. Changes from v2

| # | Location | v2 | v3 |
|---|----------|----|----|
| 1 | Subroutine count | 7 (including `calcEcr_PE`) | 6 (`calcEcr_PE` eliminated) |
| 2 | `calcDcr_PE` calls per cracked step | 3 | 2 |
| 3 | FIX-GAP as a patch | Yes (comment explains it) | Not needed (structural) |
| 4 | STATEV write location | Separate for Ecr and scr | Single pass, Step H |
| 5 | `scr_prev` naming | Implicit (passed as `scr_old`) | Explicit `scr_prev = STATEV(14:15)` |
| 6 | `Ecr_save` | Named `Ecr_OLD` (confusion with `Ecr`) | Explicit `Ecr_save` before update |
| 7 | Inline comments | Throughout all subroutines | Removed from code; only header, PROPS/STATEV block, step labels in main |

All physics from v2 is preserved identically:
- Viscous regularization (`eta`, `PROPS(6)`)
- Initiation overshoot prevention (FI > 1.05 → PNEWDT + RETURN)
- Zero seeds at initiation
- Total-form softening and secant traction
- Residual branch clamp (`scr(1) = sigcr0_r`)
- Failure overshoot check (OS-2new, corrected threshold = 1.05·Ecr_n_ult)
- Interpenetration clamp OS-4 (`dEcr(1) = -Ecr_save(1)`)
- PNEWDT in `calcDcocr_PE` with MIN rule
- Bazant size limit check with `myExit_PE`

---

## 4. PROPS Reference

| Index | Symbol | Description | Units |
|-------|--------|-------------|-------|
| PROPS(1) | E | Young's modulus | Pa |
| PROPS(2) | ν | Poisson's ratio | — |
| PROPS(3) | G_Ic | Mode-I fracture energy | J/m² |
| PROPS(4) | X_T | Tensile initiation strength | Pa |
| PROPS(5) | INDEX | Failure criterion index (integer) | — |
| PROPS(6) | η | Viscous regularization coefficient (optional, default 0) | Pa·s |

**Bazant size limit:** For mesh-objective softening, element size must satisfy
`Leff ≤ 2·G_Ic·E / X_T²`. The code checks this and calls XIT if violated.

**Viscous regularization:** Adds `η/Δt` to the softening slope, making the effective tangent less negative. Mirrors the ABAQUS cohesive element viscosity parameter. Set `η = 0` to disable. Typical range: `1e-5` to `1e-3` × characteristic time step.

---

## 5. STATEV Reference

| Index | Symbol | Description | Initial | Units |
|-------|--------|-------------|---------|-------|
| STATEV(1) | isCRACKED | Crack flag: 0 = intact, 1 = active | 0 | — |
| STATEV(2) | MODE | +1 = tensile, −1 = compressive | 0 | — |
| STATEV(3) | Leff | Effective element size (output only) | 0 | m |
| STATEV(4) | FI_MAX | Peak failure index (pre-peak tracking) | 0 | — |
| STATEV(5) | cos θ | Crack-normal cosine (global x) | 0 | — |
| STATEV(6) | sin θ | Crack-normal sine (global x) | 0 | — |
| STATEV(7) | Ecr_MAX(1) | Max normal crack strain | 0 | m/m |
| STATEV(8) | Ecr_MAX(2) | Max shear crack strain | 0 | m/m |
| STATEV(9) | Ecr(1) | Current normal crack strain | 0 | m/m |
| STATEV(10) | Ecr(2) | Current shear crack strain | 0 | m/m |
| STATEV(11) | dEcr(1) | Normal crack strain increment | 0 | m/m |
| STATEV(12) | dEcr(2) | Shear crack strain increment | 0 | m/m |
| STATEV(13) | σ_init | Max principal stress at initiation | 0 | Pa |
| STATEV(14) | scr(1) | Current normal crack traction | 0 | Pa |
| STATEV(15) | scr(2) | Current shear crack traction | 0 | Pa |
| STATEV(16) | σ_cr0 | Normal peak traction at initiation | 0 | Pa |
| STATEV(17) | τ_cr0 | Shear peak traction at initiation | 0 | Pa |

Minimum NSTATV = 17.

---

## 6. Constitutive Logic

### Softening law (linear cohesive, Bazant crack band)

```
σ_cr(ε_cr) = σ_cr0 + slope_soft · ε_cr
slope_soft = −σ_cr0² / (2·G_Ic / Leff)
ε_cr_ult   = 2·G_Ic / (Leff · σ_cr0 · (1 − R_RES))
```

### Branch selection (normal component)

```
Ecr_eff = max(|Ecr(1)|, |Ecr_MAX(1)|)

if Ecr_eff < ε_cr_ult:
    if Ecr(1) < Ecr_MAX(1):   → SECANT   D = scr_max / Ecr_MAX(1)
                                           scr = D · Ecr(1)            [total]
    else:                      → SOFTENING D = slope_soft + η/Δt
                                           scr = σ_cr0 + slope_soft·Ecr + (η/Δt)·dEcr  [total + viscous]
else:                          → RESIDUAL  D = RES_K
                                           scr = σ_cr0 · R_RES
```

### Cracked tangent

```
D_cocr = D_co − D_co·N·(D_cr + N^T·D_co·N)^{−1}·N^T·D_co
```

---

## 7. Step Control (PNEWDT) Rules

All PNEWDT modifications use `PNEWDT = MIN(PNEWDT, value)` — never a direct assignment — to avoid overriding reductions made by other checks.

| Trigger | PNEWDT |
|---------|--------|
| Initiation overshoot (FI > 1.05) | `MIN(PNEWDT, α)` where α = linear interp to FI=1.0; floor 0.10 |
| Failure overshoot (Ecr would exceed 1.05·ε_cr_ult) | `MIN(PNEWDT, α)` where α = fraction to reach ε_cr_ult; floor 0.10 |
| Y near-singular in solve | `MIN(PNEWDT, 0.5)` + RETURN |
| D_cocr near-singular | `MIN(PNEWDT, 0.5)` + fallback to D_co |

---

## 8. Known Limitations

1. **Shear traction in softening/residual branches uses incremental form.** A physically consistent total form for shear would require an explicit shear softening law independent of the normal degradation variable. The current `D_cr(2,2) = D_cr(1,1) · term · 1e-3` coupling means D_cr(2,2) is negative in the softening regime, making a total-form `scr(2) = D(2,2)·Ecr(2)` produce a wrong sign. The incremental form is acceptable since the shear traction is secondary and small.

2. **Fixed crack angle.** The crack orientation is frozen at initiation (STATEV(5:6)). No crack rotation or secondary cracking is implemented.

3. **Mode I dominated.** The residual traction `sigcr0_r = sigcr0 · 1e-8` is a numerical floor, not a physical frictional residual. Mixed-mode fracture energies are not differentiated.

4. **Linear softening only.** The Bazant crack-band model uses a linear traction-separation law. Exponential or bilinear laws would require changes to `calcDcr_PE`.
