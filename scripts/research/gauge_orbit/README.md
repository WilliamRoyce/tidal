# Gauge-orbit chains for the localized class (GH #477 stage F1 — PARKED)

**Status: ACTIVE — stage F1 in progress.** Gate F1-a has **passed at
machine precision** (below); F1-b onward is being built. No production
code exists yet, by design: the arc's checkpoint discipline forbids
shipping any solver change before the full F1 battery closes.

## What was measured (gate F1-a — PASS, N = 16/24/32)

Hand-derived linearized-diffeomorphism orbit chains, built in the modal
solver's own basis, annihilate the **vacuum** (`Bpeak=0`) E.cal pencil at
machine precision, for every gauge-profile mode:

| generator | worst relative chain defect (N=16 / 24 / 32) | verdict |
| --- | --- | --- |
| `xi_t` | 1.2e-17 / 5.1e-17 / 4.8e-17 | EXACT |
| `xi_x` | 3.2e-17 / 7.3e-17 / 6.4e-17 | EXACT |
| `xi_y` | 3.2e-17 / 7.3e-17 / 6.4e-17 | EXACT |
| `xi_z` | 4.0e-17 / 6.7e-17 / 4.5e-17 | EXACT |
| `chi_u1` | 2.4e-01 / 2.7e-01 / 3.3e-01 | **breaks at O(1)** (intended negative control) |

Two things this establishes, independent of the parked arc:

1. **The weak directions of the E.cal pencil are exactly the linearized
   diffeomorphism orbit**, with the chain shapes written down in
   `orbit.py` — including the velocity-slot structure (`v_h_3`, `v_h_6`,
   `v_h_8`) that a field-by-field analysis cannot see. The numerically
   measured weak triplets (h_0,h_3,v_h_3), (h_1,h_6,v_h_6),
   (h_2,h_8,v_h_8) are reproduced analytically.
2. **U(1) is broken in the action, not by the background.** The `chi`
   chain fails at O(1) *even at Bpeak=0*, which is the signature of the
   Lorenz gauge-fixing term this spec carries — not of the localized
   background. Any future gauge analysis of this spec must account for
   the fixing term before attributing breaking to physics.

## The identity behind it (why the defect is what it is)

Expanding the exact Noether identity `E_i R^i_α ≡ 0` to first order about
a background gives

```
R^i_α L_ij + E_i(Φ̄)·(∂_j R^i_α) = 0        (exact)
```

with `E[Φ̄]` the **background EOM residual** and `L` the linearized
operator. Because `L` is the Hessian of a quadratic action (symmetric at
the Euler–Lagrange level), the same coefficients give both the gauge
orbit (right/column null structure) and the Noether identity (left/row
structure), each annihilated up to a defect **exactly linear in the
background-EOM violation**.

Consequence for this spec: exact gauge structure ⟺ self-consistent
background. E.cal violates that twice *by construction* — flat metric
with EM stress `T^EM ~ B₀²` (test-field approximation) and a
dual-Gaussian `B₀(z)` requiring an external current `J ~ ∂_zB₀`. Both are
the literature-standard externally-imposed treatment and the physics is
unaffected; the cost is that the gauge symmetry is *approximately* broken
(breaking ≈ B₀(z)², ~11 decades across the grid, no gap) rather than an
exact null space. That is the root reason tolerance-based identification
of the gauge sector is ill-posed at float64 on this class (GH #473
checkpoint), and it is a property of the formulation, not of TIDAL.

## Files

- `orbit.py` — pencil capture (reuses the `_pencil_deflate` monkeypatch
  pattern from `../klf_staircase/klf_gate_c.py`), the orbit chain
  construction for `xi_t`/`xi_x`/`xi_y`/`xi_z`/`chi_u1`, row
  equilibration, and the chain-identity defect / relative measure.
  Conventions verified against `tidal/solver/modal.py`: full-fftn basis,
  slot-major indexing (`slot*N + mode`), `k = 2π·fftfreq(N, dx)`, and the
  E.cal spec coordinate `"x"` being physical `z`.
- `f1a_exactness.py` — the gate above. `uv run python
  scripts/research/gauge_orbit/f1a_exactness.py [N ...]`.

## Not yet measured (the rest of stage F1)

F1-b onward: breaking-vs-closed-form and the background-EOM-violation
cross-check; the orbit-span/weak-direction principal angles; the consumer
construction (symbolic-orbit quotient with residual-graded pinning) and
its end-to-end physics gate. The full design, gates and honest-scope
statement are on GH #477.

If F1-d cannot meet its gates after honest effort, the arc stops there
and the measured evidence goes to the user — the same checkpoint
discipline that closed the #473 staircase arc rather than shipping a
partial construction.
