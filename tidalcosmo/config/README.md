# `tidalcosmo/config/` — typed configuration objects

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.
>
> **Still genuinely undesigned** — no later handoff has settled this directory's
> internals. Treat the contents as a sketch.

**Responsibility.** Frozen dataclass settings bundles over a central config module read **at call
time, not import time** — PSALTer's idiom (`Parallelism`, `NSSettings`,
`LikelihoodHyperparameters` over `psalter/config.py`), adopted deliberately for consistency with
the supervisor's own code.

**This directory is the fix for the single worst structural defect in the legacy package.** There,
the CLI *is* the config object: the forward model takes an `argparse.Namespace`, which is why four
modules outside `tidal/cli/` import private names back out of it. Here, **the forward model is
library code taking a typed config, and the CLI and the Cobaya component are two thin callers of
the same entry point.**

**Workstream.** WS1 (#490). **Filled at.** M1.

**Useful input, not a shim.** `_backfill_simulate_args` in
`tidal/measurement/_posthoc_audit.py:278-318` enumerates the ~20 attributes the legacy forward
model direct-accesses. Read it as a *requirements list* for what a run configuration must carry —
never as something to build on.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
