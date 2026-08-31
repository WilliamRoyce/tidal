# `tidalcosmo/diagnostics/` — post-processing over sampler output

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

**Responsibility.** Analysis that runs *after* a chain exists, over Cobaya/anesthetic output.

**Ported from `tidal/inference/_importance.py`.** The KL-divergence and Bayesian-model-
dimensionality diagnostics — `D_KL` (total information gain), `d_G` (effective constrained
parameters), and marginal `D_KL` per parameter. These answer "which parameters does the data
actually constrain?", which no sampler gives you for free. The Handley (2019) and Handley et al.
(2015) references travel with the code.

**Workstream.** WS5. **Filled at.** M1 or later — nothing gates on it.

**Deliberately absent: plotting.** GetDist covers posterior and triangle plots (Cobaya's standard
companion) and anesthetic covers nested-sampling output. A third implementation of a solved
problem is not worth maintaining; anything genuinely bespoke belongs in a one-off figure script.

**Possibly here, possibly from PSALTer.** The cubed-sphere joint prior (`RadialAngularPrior`)
could become a Cobaya external prior. But its lineage *is* `psalter/_tile/geometry.py` — the
legacy docstring says its face and tile conventions match the supervisor's implementation — so
adopting PSALTer's directly may beat porting ours.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
