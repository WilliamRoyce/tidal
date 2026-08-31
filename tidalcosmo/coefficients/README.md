# `tidalcosmo/coefficients/` — symbolic coefficient evaluation

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

**Responsibility.** Turning a symbolic coefficient from the spec into a numeric callable of
`(η, k)`, with caching appropriate to what actually varies.

**Re-implemented informed-by, not ported.** The legacy `tidal/solver/coefficients.py` has the
right *idea* — a multi-level cache keyed on what a coefficient depends on — but its axes are
wrong for us: it caches over a spatial grid and constant-in-time coefficients, whereas here the
dependence is on conformal time and wavenumber, with the background supplied as a table.

**Workstream.** WS2/WS3. **Filled at.** M3–M4.

**Port with it.** The expression sandboxing of GH #406 — coefficient strings are allowlisted **by
name**, never `eval`'d openly. That security property must travel with its docstring.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
