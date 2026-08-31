# `tidalcosmo/solver/` — the time-dependent per-k engine (WS3)

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

**Responsibility.** Integrating `M(η, k)` over conformal time, per wavenumber, inside an inference
loop. Consume `M(η, k)` from `perturbations/`, return a transfer function.

**The internals are H3's to design; this directory is a reserved seam.** The legacy
`solver/modal.py` (5,525 lines) is **not ported and is not assumed to be an oracle** — it solves
`expm(M·t)` for *constant* `M` on a flat periodic grid, which is not our problem.

**Two front-ends over a shared core**, per H2 §0.1 — because O2 and O3 are different numerical
problems:

- an **oscillation-resolving mode-equation solver** for gravitational waves
  (`k ~ 10⁻⁴–1 Mpc⁻¹`, `~1–10³` oscillations over a Hubble time);
- an **eikonal amplitude engine with patch averaging** for CMB photons
  (`k ≈ 2×10²² Mpc⁻¹`, `~10²⁶` oscillations — no integrator steps through that, so the carrier is
  removed analytically and only the slowly varying amplitude is integrated).

**Workstream.** WS3 (#492). **Filled at.** M4, per `docs/cosmology/solver_design.md` (H3's artifact — not yet written).

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
