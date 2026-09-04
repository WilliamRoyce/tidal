# `tidalcosmo/likelihoods/` — per-rung likelihood components

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** Cobaya `Likelihood` components for observables the ecosystem does not already
cover. Most likelihoods come free from Cobaya (Planck, ACT, SPT, DESI/BAO, supernovae, lensing) —
that is the whole payoff of chaining into it, and **nothing that already exists gets rewritten
here**.

**Workstream.** WS5 (#494). **Filled at.** M5 and later, per rung.

**Open questions.**
- No drop-in birefringence likelihood exists. The escalation recorded in the program document is:
  a Gaussian prior on the published `β`, then forking `LilleJohs/cosmic-birefringence-planck-act`
  (MIT), then SPT-3G BB lite for the anisotropic case.
- H1 §3.2's caveat applies to anything quoting evidences: stock Cobaya's PolyChord does not give
  correct `log Z` under non-uniform priors without the Ormondroyd patch.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
