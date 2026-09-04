# `tidalcosmo/coefficients/` — symbolic coefficient evaluation

> **Directory boundary from H4 (2026-08-31); internals settled by H3.** The H4 study
> drew this boundary *alongside* the detailed investigations rather than after them, so the
> boundary itself may still be revised — that does not require re-litigating H4. **The
> contents below are not provisional**: they record decisions settled in
> `docs/cosmology/solver_design.md` §1 (assembly is the binding cost, #518). Do not reopen them from this file.

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

> **Amendment (coherence pass, 2026-09-04):** H3 measured that **coefficient assembly — not the matrix exponential — is the binding
> cost** of the solver: `expm` is 0.4% of the legacy 82.6 ms per call. η-grid *segmented*
> assembly is therefore the shared prerequisite for **both** solver front-ends (**#518**),
> not a supporting detail. Read this directory as the performance-critical path, not a
> low-risk port. `docs/cosmology/solver_design.md` §12 sets the implementation order.
