# `tidalcosmo/presets/` — runnable ladder configurations

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** Complete, runnable Cobaya input YAMLs for each rung of the observable ladder —
O0 pass-through, O1 tabulated background, O2 tensor, and onward. Following SOLikeT's
`presets/{defaults,templates}` pattern: shipping configurations that work is a feature, not
documentation.

**Workstream.** WS5. **Filled at.** M1 onward, one per rung as it lands.

**Note.** These files name the package in dotted-path form. Until the M7 rename they will say
`tidalcosmo.…`, so **they must not be circulated outside the project** — an external user's
`theory: {tidalcosmo.SpectatorTheory: …}` breaks on the day we rename (`repo_reshape.md` §1.2).

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
