# Overshooting in the Two-Branch D_cr Model — Analysis and Fixes

## 1. Background: The Two-Branch Crack Stiffness Model

The `calcDcr_PE` subroutine assigns the normal crack stiffness `D_cr(1,1)` using a simple
if-else based on whether the current crack strain equals the historical maximum:

```
Condition                         Branch            D_cr(1,1)
─────────────────────────────────────────────────────────────────────────
Ecr_OLD(1) >= Ecr_MAX(1)         "UNLOADING"       slope_soft   (< 0)
Ecr_OLD(1)  < Ecr_MAX(1)         "LOADING"         secant       (> 0)
```

Despite the branch names, the actual physical meaning is:

| Condition | Physical state | D_cr sign | Effect on D_cocr |
|-----------|---------------|-----------|------------------|
| `Ecr_OLD = Ecr_MAX` | **On the softening envelope** — crack is actively propagating | Negative (slope_soft) | Softening tangent |
| `Ecr_OLD < Ecr_MAX` | **Inside the damage surface** — crack is unloading or reloading | Positive (secant) | Stiffening secant |

The key values are:

$$D_{cr,nn}^{\text{softening}} = \underbrace{-\frac{\sigma\_{cr0}^2}{2G\_{Ic}/L\_{eff}}}_{\text{slope-soft  < 0}}$$

$$
D_{cr,nn}^{\text{secant}} = \frac{\sigma_{cr}(\varepsilon_{cr,n}^{MAX})}{\varepsilon_{cr,n}^{MAX}}= \frac{\sigma_{cr0} + \text{(slope-soft)}\cdot\varepsilon_{cr,n}^{MAX}}{\varepsilon_{cr,n}^{MAX}} \quad (>0)
$$

---

## 2. Tracing a Load Cycle: Where the Overshoot Lives

### Phase 1 — Monotonic loading (all increments on the softening envelope)

Each increment: `Ecr_OLD = Ecr_MAX` (they track each other) → **softening branch** always active.

```
σ_cr
 ↑                              softening curve
 σ_cr0  ●──────────────────────────────────────────── ε_cr_n
        ↘  slope_soft
         ↘
          ↘   ← each increment moves right along this line
           ↘
```

At the end of any loading increment: `Ecr_OLD = Ecr_MAX` is stored in STATEV.

---

### Phase 2 — First unloading increment (the "branch switch" increment)

DSTRAN reverses → `calcDcr_PE` still sees `Ecr_OLD = Ecr_MAX` → **softening branch** → `D_cr = slope_soft` (< 0).

The crack strain increment `dEcr` is computed as (see `calcEcr_PE`):

$$\Delta\varepsilon_{cr} = \underbrace{(D_{cr} + \mathbf{N}^T\mathbf{D}_{co}\mathbf{N})^{-1}}_{\mathbf{Y}^{-1}} \mathbf{N}^T\mathbf{D}_{co}\,\Delta\varepsilon < 0$$

Because `D_cr < 0` here, `Y` is reduced. If `|slope_soft|` is comparable to `N^T D_co N`,
`Y` can become very small → `dEcr` is large. The crack strain closes more than expected.

**End of this increment:** `Ecr = Ecr_MAX - δ` (δ > 0), so `Ecr_OLD < Ecr_MAX`.

---

### Phase 3 — Second unloading increment (THE PROBLEM INCREMENT)

`Ecr_OLD < Ecr_MAX` → **secant branch** is active → `D_cr = secant > 0`.

But `scr_old` stored in STATEV was computed using `slope_soft` in the PREVIOUS increment:

$$\sigma_{cr}^{old} = \sigma_{cr}^{env} + |\text{slope-soft}| \cdot \delta = \sigma_{cr0} + \text{slope-soft}\cdot(\varepsilon_{cr,n}^{MAX} - \delta)$$

This is the **softening curve** value at `Ecr_OLD`. However, the correct traction for the
secant branch at the same point is:

$$\sigma_{cr}^{secant} = \frac{\sigma_{cr}^{MAX}}{\varepsilon_{cr,n}^{MAX}} \cdot (\varepsilon_{cr,n}^{MAX} - \delta)$$

These are **not equal**. The softening curve traction is always ABOVE the secant line for
`ε_cr < ε_cr_MAX`:

$$\sigma_{cr}^{env}(\varepsilon) > \sigma_{cr}^{secant}(\varepsilon) \quad \forall\; 0 < \varepsilon < \varepsilon_{cr}^{MAX}$$

**Consequence:** The secant branch inherits an inflated `scr_old`, so:

$$\sigma_{cr}^{new} = \underbrace{\sigma_{cr}^{old}}_{\text{on softening curve}} + D_{cr}^{secant}\cdot\Delta\varepsilon_{cr}$$

The stress path traced on the traction–separation diagram:

```
σ_cr
  ↑    softening curve
  ●────────────────────────────── ε_cr_n
  │  ↘                          ↗  reloading
  │    ↘   (phase 1 loading)   ↗  secant
  │      ↘                   ↗
  │        ●────────────────●  ← phase 2: closure on softening curve (correct)
  │          ↘
  │            ●  ← start of phase 3: scr_old jumps UP because secant is used
  │            |    but scr_old is still the softening-curve value
  |   incorrect ↕  "jump" in traction here
  ●────────────●────────────────  ε_cr_n
               ↑
          Ecr_MAX - δ        (this is Ecr_OLD at start of phase 3)
```

The incremental stress update computes a traction that is **inconsistent** with the
traction–separation law, producing the erratic stress path the user observes.

---

### Case B — Reloading overshoot

When the crack is in the secant zone (unloaded, `Ecr_OLD < Ecr_MAX`) and the strain
increment is large enough that the crack reopens **past** `Ecr_MAX`:

$$\varepsilon_{cr,n}^{OLD} + \Delta\varepsilon_{cr,n} > \varepsilon_{cr,n}^{MAX}$$

The secant Dcr would produce a traction **above** the softening envelope, which is inadmissible.
The code then updates `Ecr_MAX` to this inflated value (see `Ecr_MAX(1) = MAX(Ecr(1), Ecr_MAX(1))`
in `calcEcr_PE`), moving the reference point on the softening envelope to the wrong place.

```
σ_cr
  ↑       softening curve
  ●─────────────────────────────────────────────── ε_cr_n
  │      ↘
  │        ●  ← correct Ecr_MAX, scr_max (on envelope)
  │         ↘
  │          |╲
  │          | ↑── spurious new Ecr_MAX after overshoot
  │          |  ↖ secant overshoots, goes ABOVE the envelope
  │          |   ↗
  ●──────────●────● ─────────────────────────────── ε_cr_n
             ↑   ↑
        Ecr_OLD Ecr_MAX
              (overshoot region)
```

---

## 3. The Overshoot Method — General Numerical Concept

In plasticity/damage mechanics, **overshooting** occurs when a trial update violates the
admissibility surface (yield or damage envelope). The standard remedy is a two-step
**predictor–corrector** (return mapping):

1. **Elastic predictor**: compute trial state assuming no new damage/plasticity.
2. **Admissibility check**: test if the trial state lies on or inside the damage surface.
3. **Return mapping** (if violated): project the trial state back onto the damage surface.

For the **linear softening** damage surface here:

$$f(\sigma_{cr},\,\varepsilon_{cr,n}^{MAX}) = \sigma_{cr} - \left[\sigma_{cr0} + \text{slope-soft}\cdot\varepsilon_{cr,n}^{MAX}\right] \leq 0$$

An overshoot occurs when the secant-computed `dEcr` carries `ε_cr_n` above `ε_cr_n^{MAX}`,
which would give `f > 0`.

---

## 4. Proposed Numerical Fixes

The following fixes address the two overshooting cases without altering the physical
constitutive law (softening curve on the envelope; secant within the envelope).

---

### Fix OS-1 — Total-form `scr` update in the secant branch

**Root cause addressed:** `scr_old` is on the softening curve when the secant branch starts,
causing an inconsistent incremental update.

**Fix:** In `calcDcr_PE`, when the secant branch is active, compute `scr` as the **total
form** (direct projection onto the secant line) rather than as an incremental update from
a potentially inconsistent `scr_old`:

$$\sigma_{cr}^{(n+1)} = D_{cr}^{secant} \cdot \varepsilon_{cr,n}^{OLD}$$

instead of

$$\sigma_{cr}^{(n+1)} = \sigma_{cr}^{old} + D_{cr}^{secant} \cdot \Delta\varepsilon_{cr,n}^{old}$$

This resets the traction to the secant line at each increment in the secant zone, eliminating
the inherited error from the softening-curve value. The incremental update is replaced by a
direct evaluation of the secant constitutive law.

**Physical justification:** The secant branch represents the relationship
`σ_cr = (σ_cr^{MAX}/ε_cr^{MAX}) · ε_cr`, which is a TOTAL (not incremental) law. Using the
total form is more faithful to the underlying damage model.

**Effect on the stress path:**
```
σ_cr
  ↑       softening curve
  ●───────────────────────────── ε_cr_n
  │      ↘
  │        ●─ (Ecr_MAX, scr_max)
  │       ╱ ↘
  │      ╱   ↘   ← still follows softening during loading
  │     ╱      ↘
  │    ╱ secant  ↘
  │   ╱ (clean,   ●──── ε_cr_n
  ●──●    no jump)
  (0,0) ← all unloading/reloading passes through origin
```

---

### Fix OS-2 — Reloading overshoot detection + PNEWDT (Case B)

**Root cause addressed:** The secant branch computes `Ecr_new > Ecr_MAX`, which would
spuriously advance the reference point on the softening envelope.

**Algorithm:**

Before updating `Ecr_MAX`, record `Ecr_MAX_in` at the entry of `calcEcr_PE`. After
computing `dEcr`:

$$\text{if } \varepsilon_{cr,n}^{OLD} < \varepsilon_{cr,n}^{MAX,in} \;\text{ AND }\; \varepsilon_{cr,n}^{OLD} + \Delta\varepsilon_{cr,n} > \varepsilon_{cr,n}^{MAX,in}$$

then the reloading front has hit the softening envelope mid-step. The PNEWDT hint is:

$$\text{PNEWDT} = \min\!\left(\text{PNEWDT},\;\underbrace{\frac{\varepsilon_{cr,n}^{MAX,in} - \varepsilon_{cr,n}^{OLD}}{\Delta\varepsilon_{cr,n}}}_{\alpha} \times 0.85\right)$$

where α is the fraction of the step at which `Ecr = Ecr_MAX`. ABAQUS will restart the
increment with `DSTRAN' = α · DSTRAN`, so the crack strain just reaches (but does not
exceed) the previous maximum. In the subsequent increment (with the remaining strain),
the softening branch will be correctly active.

---

### Fix OS-3 — Loading-reversal transition flag + PNEWDT (Case A)

**Root cause addressed:** When the crack is on the softening envelope and the global
strain reverses, the first closure increment correctly uses `slope_soft`. However, at
the END of that increment the branch will switch to secant in the NEXT increment,
causing a traction jump. Requesting a smaller step at the reversal reduces the size
of the discrepancy.

**Algorithm:** In the main `sca2diso_PE`, save `Ecr_MAX_save = STATEV(7:8)` and
`Ecr_OLD_save = STATEV(9:10)` before calling `calcEcr_PE`. After `calcEcr_PE` returns:

$$\text{if } \varepsilon_{cr,n}^{OLD,save} \geq \varepsilon_{cr,n}^{MAX,save} \;\text{ AND }\; \Delta\varepsilon_{cr,n} < 0$$

then we just left the softening envelope. Request:

$$\text{PNEWDT} = \min(\text{PNEWDT},\; 0.5)$$

This does NOT prevent the traction jump, but reduces its magnitude (smaller δ → smaller
difference between softening-curve and secant values at `Ecr_OLD`).
Together with Fix OS-1 (total-form secant), Fix OS-3 ensures smooth numerical transitions.

---

### Fix OS-4 — Interpenetration clamping (corrected formula)

**Bug found:** The original clamping in `claude_AFEM2D_PE.for`:

```fortran
dEcr(1) = -Ecr_OLD(1)   ! Sai: Need to verify this
```

is **wrong** because at this point `Ecr_OLD` has already been updated inside `calcEcr_PE`
to equal `Ecr` (the new, negative crack strain). So this sets `dEcr(1) = -Ecr(1)`, not
the physically correct value.

The correct formula uses `STATEV(9)` which still holds the **true previous** crack strain
(it has not been overwritten yet at this point in the code):

$$\Delta\varepsilon_{cr,n}^{clamp} = 0 - \varepsilon_{cr,n}^{OLD,STATEV} = -\text{STATEV(9)}$$

so that `ε_cr_n^{OLD,true} + Δε_cr_n^{clamp} = 0`, which correctly drives the crack
strain to zero from its true previous value.

---

## 5. Summary of Changes in `claude_AFEM2D_PEv1.for`

| Fix    | Location              | Change                                                         |
|--------|-----------------------|----------------------------------------------------------------|
| OS-1   | `calcDcr_PE`          | Secant branch: replace incremental `scr = scr_old + Dcr*dEcr` with total form `scr = Dcr * Ecr_OLD` |
| OS-2   | `calcEcr_PE`          | After computing `Ecr`, check if reloading overshoot occurred; set PNEWDT proportional to α = (Ecr_MAX − Ecr_OLD)/dEcr |
| OS-3   | `sca2diso_PE`         | After `calcEcr_PE`, detect loading-reversal (first closure step); set PNEWDT = 0.5 |
| OS-4   | `sca2diso_PE`         | Fix interpenetration clamp: use `STATEV(9)` (true old crack strain) instead of updated `Ecr_OLD` |

---

## 6. What Is NOT Changed

- The softening law: `σ_cr = σ_cr0 + slope_soft · ε_cr_n` (unchanged)
- The Bazant limit check (unchanged)
- The N matrix and crack kinematics (unchanged)
- The consistent tangent `D_cocr` formula (unchanged)
- The failure criterion (unchanged)
- The Dcr branch selection logic (slope_soft on envelope, secant off envelope — unchanged)

The converged, quasi-statically correct solution is identical. Only the **numerical path
to convergence** is improved.

---

## 7. Limitations and Notes

- Fix OS-2 + OS-3 use PNEWDT, which only works in ABAQUS/Standard with automatic
  time stepping enabled. They have no effect in explicit or if the user fixes the step size.
- Fix OS-1 (total form) means `scr` in the secant zone is always `D_secant × Ecr_OLD`.
  This is the exact traction predicted by damage mechanics for secant unloading.
  If the user wants a genuinely incremental (and potentially non-origin-crossing) unloading
  model, OS-1 should be disabled and only OS-2 + OS-3 used.
- A complete return-mapping algorithm (explicit analytical) is possible for the linear
  softening law and would be the most rigorous fix, but it requires sub-stepping within
  the UMAT and significantly more code.
