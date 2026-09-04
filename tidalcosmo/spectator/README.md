# `tidalcosmo/spectator/` — the Cobaya `Theory` component

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** Holds `SpectatorTheory`, the Cobaya `Theory` class that *is* the extension:
`get_requirements() -> {"CAMBdata": None}`, `must_provide()`, `calculate()`, `get_*()`. Its
defaults live beside it as `SpectatorTheory.yaml` — Cobaya requires `<ClassName>.yaml` to sit in
the same directory as the module defining the class, which is why this is a directory rather than
a single file.

Per-channel configurations are **YAML presets over the one class**, not separate classes — the
pattern `mflike/{TT,TE,TTTEEE}.yaml` uses over a single `mflike.py`.

**Workstream.** WS5 (#494). **Filled at.** M1 (pass-through, for the O0 gate).

**Open questions.**
- Which products the class declares it provides, per rung — replaced transfer functions, or
  conversions/rotations applied to CAMB output. Depends on which of the three mechanisms
  (`repo_reshape.md` §2.4) each rung needs.
- Failure behavior is a *flagged rejection*, never a crash (§3); the mechanism is shared with
  `validity/` and `spectrum/`, not reinvented here.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
