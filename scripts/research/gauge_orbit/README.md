# Gauge-orbit chains for the localized class (GH #477 stage F1 — PARKED)

**Status: ACTIVE — stage F1 in progress.** Gates F1-a, F1-b, F1-c and
F1-f have **passed**; F1-d (the consumer end-to-end, which carries the
arc's stop-checkpoint) is next. No production
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
- `f1a_exactness.py`, `f1b_breaking.py`, `f1c_span.py` — the gates
  above; each runs standalone, e.g. `uv run python
  scripts/research/gauge_orbit/f1c_span.py 24`.
- `orbit.joint_operator` / `orbit.grade` — the joint chain map and the
  relative grading sigma(u) = ||D(u)||/||C(u)|| that the production
  consumer would use (null(C) dropped first to avoid phantom profiles).

## F1-b (PASS) — the defect IS the background-EOM violation

The identity predicts different powers of the background amplitude for
the two violations, which is a sharp falsifiable statement no other
origin of the defect would reproduce. Measured (N=24, log-log fit over
Bpeak = 0.0025 … 0.02):

| generator | a-row slope (predict 1) | h-row slope (predict 2) |
| --- | --- | --- |
| `xi_t` | 1.000 | 1.998 |
| `xi_x` | identically 0 | 2.000 |
| `xi_y` | 0.988 | 1.988 |
| `xi_z` | 1.000 | 1.999 |
| `chi_u1` | 0 (background-independent) | identically 0 |

`xi_x` reaches no photon row at all — the a-rows read h_2, h_5 and
h_0/h_4/h_7/h_9 while the `xi_x` orbit is h_1, h_6, v_h_6 — so its
Maxwell-violation defect is *exactly* zero rather than O(Bpeak). The
defect is also spatially localized where the background lives (far-field
share 4.0e-3 for `xi_t`, 1.7e-3 for `xi_x`).

## F1-f (PASS) — the generator representative matters, by 123x

`delta a = L_xi Abar` and `delta a = xi^nu Fbar` differ by a U(1) with
chi = xi·Abar, which this spec's Lorenz gauge-fixing term breaks. Measured
on `xi_y`:

| representative | a-row defect | far-field share |
| --- | --- | --- |
| full Lie `L_xi Abar` | 8.28e-2 | 1.87e-1 |
| covariant `xi^nu Fbar` | 6.71e-4 | 1.33e-4 |

`xi_z` is the control: `Abar_z = 0`, so both representatives coincide and
the numbers match to the digit. This is why the consumer grades the
**joint** span of all candidates rather than each generator separately —
a per-generator reading would over-report `xi_y`'s breaking by 123x and
fail to pin it.

## F1-c (PASS) — the symbolic orbit IS the numerically weak sector

The experiment this issue was opened on. Two independent findings:

1. **The grading rediscovers the physics unprompted.** Told nothing about
   B0(z), the least-determined gauge profiles come out with a
   near-Gaussian energy share of **0.000** (purely far-field) while the
   most-determined sit at 0.998 — i.e. "the symmetry returns where the
   background vanishes" falls out of the operator alone.
2. **The pencil's weak directions lie in the orbit span.** The 4–8
   weakest right singular vectors of `A - lambda B` sit within 6.4e-3 rad
   of the symbolic orbit span (median 1.3e-3 over the 32 weakest), with
   an unexplained component of 0.002–0.006. The weakest direction of all
   lives on h_0=0.63, h_3=0.51, v_h_3=0.50, h_9=0.21, v_h_9=0.21 — a
   `xi_t`+`xi_z` orbit vector.

This also settles the open `xi_z` question: `xi_z` is *not* better
determined than the others (sigma median 2.4e-4 vs 1.3e-3 for `xi_t`), so
hypothesis H1 is refuted; it was missing from the 2026-08-27 weak-triplet
list because of its depth-3 (L_2) chain shape, which an L_1-shaped probe
does not resolve (H2). It is plainly present in the weakest direction
above.

## Not yet measured (the rest of stage F1)

F1-d: the consumer
construction (symbolic-orbit quotient with residual-graded pinning) and
its end-to-end physics gate; F1-e (Noether-row map); F1-g (cost). The full design, gates and honest-scope
statement are on GH #477.

If F1-d cannot meet its gates after honest effort, the arc stops there
and the measured evidence goes to the user — the same checkpoint
discipline that closed the #473 staircase arc rather than shipping a
partial construction.
