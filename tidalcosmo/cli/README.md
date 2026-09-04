# `tidalcosmo/cli/` — the thin command-line adapter

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** Parse arguments, build a typed config from `config/`, call the library. **No
physics lives here.**

That sentence is the entire point of the directory. In the legacy package the CLI *is* the config
object and holds the forward model — `tidal/cli/` is 25,698 of ~67k lines, `_simulate.py` alone is
3,126 — which is why four modules outside it (`tidal/inference/_likelihood.py`,
`tidal/measurement/{_posthoc_audit,_run_stages,_sweep_results}.py`) import **private** names back
out of a CLI module, every one of them behind a `reportPrivateUsage` suppression.

`_run_stages.py` is the instructive one: it exists *as a seam* to hide that coupling from the other
callers, and its own docstring records why relocating the wrappers would not fix anything — the
real dependency is the ~3,000-line simulation driver and the eleven private measurement
dispatchers, not the wrappers around them. A seam over a bad boundary is still a bad boundary.
Here the CLI and the Cobaya component are instead **two thin callers of one library entry point**.

**Workstream.** WS1 (#490). **Filled at.** M1 onward, growing only as needed.

**Scope.** Far smaller than legacy's eleven subcommands. `sweep`, `sample`, `analyze` and `plot`
are **dropped** — Cobaya supplies sampling and priors (D9), GetDist and anesthetic supply plots.
Expect roughly: derive, inspect, validate, and a run entry point.

**Convention that ports.** User-facing errors carry actionable hints — `error_with_hint(msg,
hints)` rather than a bare error, as used across ~60 legacy error sites.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
