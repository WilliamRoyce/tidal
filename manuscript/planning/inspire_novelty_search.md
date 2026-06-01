# INSPIRE-HEP novelty search — PGT × Gertsenshtein

**Date executed:** 2026-05-23 (during introduction-lit-review planning session)
**Search agent:** Explore subagent with WebFetch access to `inspirehep.net/api/literature`
**Purpose:** verify the novelty of computing graviton↔photon conversion in a Poincaré Gauge / torsion / metric-affine background before any "first treatment" wording is committed in §1 prose.

## Verdict

**No DIRECT prior work found.** The calculation of graviton↔photon conversion in a
PGT / torsion / metric-affine gravity background with dynamical torsion modifying
Einstein–Maxwell is genuinely novel and defensible as a first-principles contribution.

The literature contains extensive work on each component separately
(Gertsenshtein in GR; Gertsenshtein in scalar-tensor / f(R); torsion–photon
coupling without GWs; PGT propagating modes without EM) but no paper combining
them was found.

## Methodology

25 queries executed against INSPIRE-HEP (`/api/literature?q=…&format=json&fields=titles,authors,arxiv_eprints,abstracts,publication_info&size=25`) and arXiv. Each query and its variant were checked; ~50 papers screened in detail; ~500 abstracts scanned.

Queries (sample, not exhaustive):

| # | Query | Hits | Notes |
|---|---|---|---|
| 1 | `"Poincare gauge" AND "Gertsenshtein"` | 0 | Exact concept combination absent |
| 2 | `"Poincaré gauge" AND "Gertsenshtein"` (accent) | 0 | Encoding variant — same result |
| 3 | `torsion AND "Gertsenshtein"` | 0 | Core novelty untouched |
| 4 | `torsion AND "graviton-photon"` (and variants) | 0 | Core novelty untouched |
| 5 | `"metric-affine" AND ("graviton photon" OR Gertsenshtein)` | 0 | Sister theory blank |
| 6 | `"Einstein-Cartan" AND ("Gertsenshtein" OR "graviton photon")` | 0 | EC torsion blank |
| 7 | `teleparallel AND ("Gertsenshtein" OR "graviton photon")` | 0 | Teleparallel blank |
| 8 | `"Poincare gauge gravity" AND (graviton OR "gravitational wave")` | several | All ADJACENT (no EM coupling) |
| 9 | `"high-frequency gravitational wave" AND torsion` | 0 | No HF-GW × torsion phenomenology |
| 10 | `"Poincaré gauge" AND "high-frequency gravitational"` | 0 | — |

The zero-hit pattern is itself informative: even loosely-related keyword combinations return nothing, which raises confidence that the gap is real and not a missed-query artefact.

## ADJACENT priors — to cite as "we extend / differ from" context in §1 ¶9

These are the closest prior works. They overlap one component (torsion, or
graviton-photon, or modified-gravity Gertsenshtein) but not the combination.

| Paper | arXiv | What they did | Why it's adjacent (not direct) |
|---|---|---|---|
| Kruglov 2008 | [0804.4011](https://arxiv.org/abs/0804.4011) | Light–torsion mixing in B-field from higher-derivative QG | Torsion is the mixed-with field, not a background; no GW |
| Gaete & Helayel-Neto 2008 | [0810.2989](https://arxiv.org/abs/0810.2989) | Torsion–EM coupling in path-dependent formalism, confining potential | No gravitons; treats photons classically |
| Gaete & Helayel-Neto 2009 | [0912.3767](https://arxiv.org/abs/0912.3767) | Screening / confinement from interacting EM + torsion fields | Same family as above |
| Ejlli 2020 | [2004.02714](https://arxiv.org/abs/2004.02714) *(in local cache)* | Exact GW–photon mixing in constant B (standard GR) | No torsion |
| Cembranos, Diaz & Ortiz 2018 | [1806.11020](https://arxiv.org/abs/1806.11020) | Graviton-photon oscillation in alternative theories of gravity | Scalar-tensor / f(R), not torsion |
| Cembranos et al. 2023 | [2302.08186](https://arxiv.org/abs/2302.08186) *(in local cache)* | Graviton-photon oscillation in cosmic background for general theory | Same family; leaves PGT gap explicitly open |
| Addazi & Capozziello 2024 | [2401.15965](https://arxiv.org/abs/2401.15965) | Resonant graviton–photon conversion with stochastic B in expanding universe | Standard GR; cosmological backgrounds only |
| Matsuo & Ito 2025 | [2505.08457](https://arxiv.org/abs/2505.08457) | Graviton–photon conversion in blazar jets as HF-GW probe | Standard GR; astrophysical application |
| Janssen & Jiménez-Cano 2018 | [1807.10168](https://arxiv.org/abs/1807.10168) | Metric-affine gravity projective symmetry → EM coupling | Different mechanism; sister formalism |
| PGT propagating-modes | [2109.09546](https://arxiv.org/abs/2109.09546), [1910.07506](https://arxiv.org/abs/1910.07506), [2003.00664](https://arxiv.org/abs/2003.00664) | Massive spin-1 / spin-0 torsion modes in quadratic PGT | No EM coupling, no Gertsenshtein |

## Recommended novelty-claim wording

For §1 ¶9 (adapt during drafting; this is a starting point, not committed prose):

> "To our knowledge, this is the first treatment of the Gertsenshtein effect in Poincaré Gauge Theory. While torsion–photon interactions have been studied in quantum-gravity contexts \cite{Kruglov2008, GaeteHelayelNeto2008, GaeteHelayelNeto2009}, and graviton–photon mixing has been generalised to alternative theories of gravity \cite{Cembranos2018, Cembranos2023, AddaziCapozziello2024}, the combination — GW↔photon conversion in torsion-sourced gravity — remains unexplored. The present work closes that gap by deriving the mixing kernel in a PGT background where torsion is a propagating field coupled to Maxwell electrodynamics, and surveying the parameter space of quadratic PGT for sectors where the GR Gertsenshtein result is enhanced, suppressed, or qualitatively modified."

**Important framing constraints:**

- ✅ Defensible: *"first treatment of the Gertsenshtein effect in Poincaré Gauge Theory"*
- ❌ Overclaim: *"first treatment of torsion in graviton-photon physics"* (Kruglov 2008 etc. exist)
- ❌ Overclaim: *"no prior work on torsion and electromagnetism"* (Gaete–Helayel-Neto, Hehl–Obukhov etc. exist)
- ✅ Defensible: *"the combination of GW↔photon conversion with dynamical PGT torsion remains unexplored"*

## Caveats and update protocol

1. **Supervisor-provided references must be re-tested.** If Barker or a collaborator cites a paper not in this search, run it through the title/abstract test against the table above before assuming novelty still holds.
2. **The search is INSPIRE-centric.** Conference proceedings or workshop talks not indexed in INSPIRE could in principle exist; if a relevant talk surfaces (e.g. through supervisor channels), update this artefact.
3. **If a DIRECT prior surfaces post-planning**, the §1 ¶9 framing changes from "first treatment" to "extending / refining \cite{NewPriorRef}". The rest of the lit-review plan is unaffected.

## See also

- The full lit-review plan: [planning-the-manuscript-writing-structured-goose.md](../../../home/vscode/.claude/plans/planning-the-manuscript-writing-structured-goose.md) (note: lives under `~/.claude/plans/`, not in the repo)
- Project-side narrative context: [literature_critical_analysis.md](../../../home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory/literature_critical_analysis.md)
- The three constructive paths discussed at §1 ¶8: [perturbative_reduction_v6_complete.md](../../../home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory/perturbative_reduction_v6_complete.md)
