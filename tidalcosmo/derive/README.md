# `tidalcosmo/derive/` — the symbolic derivation front end

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** TOML → `.wls` generation, the `wolframscript` driver, and the ported Wolfram
modules. All symbolic processing stays in Wolfram — **never post-process equations in Python.**

**This is a port *and* a substantial extension, not a lift-and-shift.** The `.wl` modules (5,982
lines) are the asset and are largely CLI-independent. `tidal/cli/_derive.py` (6,718 lines, the
largest file in the codebase) is where TOML→`.wls` generation is *fused to `argparse`*; porting
means separating the generator from the CLI. Beyond that, WS2 must **add** what legacy never
needed:

- conformal time `η` as a first-class coordinate;
- `a(η)`, `H(η)` as **unspecified background functions** surviving decomposition and export;
- gauge projection into a CAMB-named gauge, performed in Wolfram;
- emission in CAMB/PSALTer conventions, with the gauge recorded as spec metadata;
- **per-channel source functions** (`repo_reshape.md` §2.5) — a derivation deliverable nobody
  currently owns;
- a **second export form** for H2's eikonal reduction (O3).

The known `ExportJSON.wl:1638` `t`-filtering bug (time-dependent terms dropped from the
Hamiltonian export) is fixed in the ported copy — but it is one bug inside a much larger
extension, not the extent of the work.

**Workstream.** WS2 (#491). **Filled at.** M3, as its own milestone.

**Cost basis.** The v0.33.9 measured table in `docs/tex/derivation_performance.tex`. The
per-theory timing headers in the TOMLs are declared untrustworthy there — ceilings, not estimates.

**Constraint.** One `wolframscript` at a time; the engine license is single-seat.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
