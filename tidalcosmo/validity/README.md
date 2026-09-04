# `tidalcosmo/validity/` — honest per-run flags

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** Enforcing the spectator approximation numerically, per run, instead of
asserting it in a sentence the way the literature does. This is a supervisor-flagged, first-class
requirement, and a cheap methodological contribution in its own right.

- `flags.py` — the flag vocabulary and the shared **flagged-rejection** mechanism (the same one
  `spectrum/` and the `Theory`'s NaN guard use — one concept, not three).
- `spectator.py` — `ρ_new/ρ_γ` against the `ΔN_eff ≲ 0.1` bound; amplitudes `|h|, |f| ≪ 1`; and
  the **growth-impact monitor**: per mode, the ratio of new-sector to standard-sector source terms
  in the Einstein constraints — "would these perturbations have affected the growth we froze?"
- `background_eom.py` — the **background-EOM residual**. CAMB's background solves *Einstein's*
  equations while our quadratic action is PGT; consistency requires the PGT background field
  equations be satisfied on that background to sufficient accuracy. The concept is proven by the
  #477 work; **its code is not carried over.**

**Workstream.** WS2/WS4. **Filled at.** M1 (skeleton), M5 (the physics).

**Design rules.** A verdict is a flag *with a reason*, never a bare boolean. A flag is a
diagnostic, **never a gate** — GH #454: the legacy stability probe blocked runs, and the two code
paths silently drifted apart for four months as a result.

**Note.** `ρ_new/ρ_γ` is computed fresh here. It is *not* a port of `measurement/_energy.py`,
whose flat-space Hamiltonian reconstruction has no FRW counterpart.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
