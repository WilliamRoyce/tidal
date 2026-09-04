# `tidalcosmo/spec/` — the equation-spec interchange contract

> **Directory boundary from H4 (2026-08-31); internals settled by H4 + H6/H8.** The H4 study
> drew this boundary *alongside* the detailed investigations rather than after them, so the
> boundary itself may still be revised — that does not require re-litigating H4. **The
> contents below are not provisional**: they record decisions settled in
> `docs/cosmology/spectrum_design.md` §6.1 (the Stage-1 contract). Do not reopen them from this file.

**Responsibility.** Loading, validating and querying the symbolic pipeline's output — the contract
between the Wolfram stage and everything numeric.

**The schema is new, not inherited.** Per the standing goal below, the spec is emitted in
CAMB/PSALTer conventions rather than legacy TIDAL notation. It also carries the **declared gauge
as first-class metadata**, so the choice travels with the equations and the CAMB seam can assert
on it.

**Workstream.** WS2 (#491). **Filled at.** M3 (though a minimal reader may be needed earlier to
consume committed specs).

**Ported here.** `spec_query.py` and `sign_algebra.py` port near-verbatim from
`tidal/symbolic/` — self-contained, high-value, #401 lineage. Their soundness-over-coverage
property must travel: a sign verdict says `unknown` rather than guessing, and every verdict names
the tactic that decided it. `json_loader.py` is *re-implemented informed-by*: keep the concepts
(operator vocabulary, kinetic handling, `implicit_dynamical_sector`, `dependency_closure`), not
the schema.

**Settled — was an open question, twice over.**

> **Amendment (coherence pass, 2026-09-04):** this read *"Whether we emit WXF alongside JSON … is unresolved."* It is resolved in both
> directions. H6 §3.2 superseded the alignment item: the WXF as it exists is *insufficient* —
> unlabeled `J`-blocks mixing parities, placeholder slots — so we emit the richer Stage-1
> contract of `spectrum_design.md` §6.1 instead. And #523 goes further: **PSALTer writes no
> WXF at all** and populates only two association keys; the committed `.wxf` files were
> produced by separate curation tooling. WXF is something we *produce*, not an interchange we
> receive. The premise "PSALTer ingests WXF" also conflated the released **JAX** package's
> `extract` (which reads WXF) with the **Wolfram** `ParticleSpectrum` (which ingests a
> Lagrangian in-session).

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
