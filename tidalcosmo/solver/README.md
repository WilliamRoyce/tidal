# `tidalcosmo/solver/` — the time-dependent per-k engine (WS3)

> **Directory boundary from H4 (2026-08-31); internals settled by H3.** The H4 study
> drew this boundary *alongside* the detailed investigations rather than after them, so the
> boundary itself may still be revised — that does not require re-litigating H4. **The
> contents below are not provisional**: they record decisions settled in
> `docs/cosmology/solver_design.md`. Do not reopen them from this file.

**Responsibility.** Integrating `M(η, k)` over conformal time, per wavenumber, inside an inference
loop. Consume `M(η, k)` from `perturbations/`, return a transfer function.

**The internals are designed — see `docs/cosmology/solver_design.md`** (H3, 2026-08-31,
`c25168df`). This is no longer a reserved seam. Key decisions to build to:

- **Architecture: option (iii)** — our own solver chained to *unmodified* CAMB, both engines
  over one shared core. Decided structurally, not by cost: no Boltzmann code has per-frequency
  photon propagation, so O3 is only possible this way.
- **The binding cost is coefficient *assembly*, not the exponential** — `expm` is 0.4% of the
  legacy 82.6 ms. η-grid segmented assembly is the shared prerequisite (**#518**), and it
  gates both front-ends.
- **A stepper registry decided by bake-off**, not a fixed priority ladder: exponential
  midpoint → GL4/CF4 Magnus with Duhamel sources → adiabatic/matrix-WKB (**#519**, gated on
  the bake-off showing Magnus over budget) → piecewise-analytic transfer matrices.
- **Classifier ships soft** — logged heuristics plus a `--solver` override; hard gates deferred.
- Two prototype results worth not rediscovering: the Magnus review's printed **CF4 coefficients
  are a misprint** (correct: `(3±2√3)/12`), and **`cond(V)` guards miss a straddled Jordan
  point** (11% silent error) — use the `‖Γ‖h ≤ 0.02` trigger instead.

The legacy
`solver/modal.py` (5,525 lines) is **not ported and is not assumed to be an oracle** — it solves
`expm(M·t)` for *constant* `M` on a flat periodic grid, which is not our problem.

**Two front-ends over a shared core**, per H2 §0.1 — because O2 and O3 are different numerical
problems:

- an **oscillation-resolving mode-equation solver** for gravitational waves
  (`k ~ 10⁻⁴–1 Mpc⁻¹`, `~1–10³` oscillations over a Hubble time);
- an **eikonal amplitude engine with patch averaging** for CMB photons
  (`k ≈ 2×10²² Mpc⁻¹`, `~10²⁶` oscillations — no integrator steps through that, so the carrier is
  removed analytically and only the slowly varying amplitude is integrated).

**Workstream.** WS3 (#492). **Filled at.** M4, per `docs/cosmology/solver_design.md` — 650 lines, landed 2026-08-31. Follow-ups: **#518** (assembly, prerequisite), **#519** (matrix-WKB, bake-off-gated), **#520** (this rewrite), **#505** (O3 front-end, needs #504 + #518).

**One inherited test contract, not code.** GH #367 and #379: **every dispatch path must consume
the kinetic matrix `M` identically**, or cross-path regressions appear silently. Two such
regressions occurred in legacy. Whatever H3 designs, that property is testable and should be
tested.

**Budget.** 10× slower than CAMB is acceptable; 100× is fatal.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
