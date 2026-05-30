# Citation Audit Tracker

> Live status for the two-part manuscript citation audit (Task 1: BibTeX entry
> verification; Task 2: citation-usage verification). Source of truth for progress
> is the inline `% AUDIT` / `% CITE` markers in `manuscript/references.bib` — run
> `uv run python -m scripts.bibaudit.audit_status` any time for the live tally.

**Last updated:** 2026-05-30 (session: Task 1 COMPLETE)
**Current phase:** Task 1 DONE (all 230 active entries audited; latexmk exit 0) → **Task 2 next**.
**Verdict tally (230 active):** MATCH 38 · MINOR 143 · ENRICHED 17 · FIXED 9 · REVIEWED 11 · NOID 4 · NOTFOUND 3 · PENDING 4 · UNCITED 1
**Compile:** latexmk exit 0, 0 undefined citations/refs; BibTeX 0 errors (12 pre-existing missing-journal warnings).

### Task 1 closeout
- **Books**: 9 enriched w/ Crossref-verified DOIs (higham, henneaux, leveque, davis, hobson, blagojevic2002, hairer ODE-I/ODE-II/geometric); 4 REVIEWED (boyd/roache/misner/MacKay — no DOI / Crossref false-match); 5 pre-arXiv historical REVIEWED (complete, predate DOIs). User did golub/trefethen×2/brenan/mathematica/randall. **Manual residual for user: NONE.**
- **xAct**: retyped `martingarcia2007xact`→@misc; added canonical **xPerm** `martingarcia2008xperm`; cited both at 3 suite sites.
- **SuperLU**: added SuperLU_MT `demmel1999superlumt` (the variant TIDAL uses), swapped into 4 cites; deleted wrong-variant SuperLU_DIST.
- **NumPy** `harris2020numpy` cited at numerical.tex. **Duplicates removed**: mekki2023skitsundae, blagojevic2013gauge (cite redirected).
- **QZ/`moler1973algorithm`** PARKed (no QZ/generalized-eigenvalue text in manuscript) — flagged for user.
**Plan:** `/home/vscode/.claude/plans/working-on-the-manuscript-glimmering-knuth.md`

## ▶ Quick handoff (read first if resuming cold)

1. **Where we are:** Phase 0 done — read-only tooling built in `scripts/bibaudit/`
   and validated (clean field extraction on all 231; classifier sanity passed).
   Task 1 network classification (`compare.py`) has been kicked off.
2. **What's running / next:** check `scripts/bibaudit/cache/report.json` (or
   `compare_run.log`). Then: review the SUBSTANTIVE list at the **STOP gate**
   *with the user* before applying any field changes.
3. **How to resume / check progress:** everything is on disk and idempotent —
   - `uv run python -m scripts.bibaudit.audit_status` → X/231 audited, Y/192 annotated
   - `uv run python -m scripts.bibaudit.compare --summary-only` → re-print verdicts
   - markers already in `references.bib` = done; `compare.py` re-run skips cached fetches.

## Method (do not deviate)

- Scripts are **read-only**; every change to `references.bib` goes through the
  **Edit tool** (exact match). Field values copied verbatim from the cache.
- **Enrichment** (sparse stub → fuller INSPIRE record) is pre-authorized.
  **SUBSTANTIVE** mismatches (author/title/year) → **STOP**, resolve with user.
- Journals standardized to **INSPIRE-abbreviated**; light normalization only
  (journal + page-dash + arXiv-id + DOI-case), preserve field order/indentation.
- Task 2 notes accrete **per-key** in the bib comments (`%@cite` block).

## Task 1 — BibTeX entry audit (231 entries)

Batched by the bib's `% === section ===` headers. Tick when markers applied +
`audit_status` confirms + diff reviewed.

| # | Section (batch) | entries | classified | markers applied |
|---|-----------------|--------:|:----------:|:---------------:|
| 1 | Gertsenshtein Effect | 12 | ☐ | ☐ |
| 2 | Torsion & Poincaré Gauge Theory | 13 | ☐ | ☐ |
| 3 | Chern-Simons Theory | 2 | ☐ | ☐ |
| 4 | Axion-Photon Mixing | 1 | ☐ | ☐ |
| 5 | Numerical Methods | 23 | ☐ | ☐ |
| 6 | Numerical Methods — Textbooks | 5 | ☐ | ☐ |
| 7 | Software & Codebases | 9 | ☐ | ☐ |
| 8 | Gravitational Wave Observatories | 5 | ☐ | ☐ |
| 9 | Pre-arXiv Foundations | 12 | ☐ | ☐ |
| 10 | Barker Software Papers | 3 | ☐ | ☐ |
| 11 | Barker–Hobson–Lasenby Joint Corpus | 14 | ☐ | ☐ |
| 12 | Barker Recent Solo/Collaborator | 9 | ☐ | ☐ |
| 13 | Key Collaborator / Field Papers | 5 | ☐ | ☐ |
| 14 | Torsion-EM Coupling Foundations | 13 | ☐ | ☐ |
| 15 | Constraint promotion / higher-derivative | 64 | ☐ | ☐ |
| 16 | Lit-review additions (intro §1) | 23 | ☐ | ☐ |
| 17 | Round-3 v3 bib closures | 7 | ☐ | ☐ |
| 18 | v8 instrument papers (LIGO/Virgo/KAGRA) | 11 | ☐ | ☐ |

### STOP / substantive-review list (resolve with user before fixing)

Classification run 2026-05-30 (Crossref-primary, hardened normalizer): 9 SUBSTANTIVE
after suppressing Crossref data-quality artifacts (curly quotes, MathML, author-list
order/truncation). Triaged:

**REAL fixes (need decision):**
1. `berlin2024axionphoton` — arXiv **2405.08865** is *Ginés, Noordhuis, Weniger, Witte*
   "Numerical analysis of resonant axion-photon mixing" (title matches). Local author
   = "Berlin, Asher and others" is WRONG (likely confused with [[Berlin:2021txa]]).
   → fix author/title to Ginés et al., OR the key intent is a different Berlin paper.
2. `shoshany2021ogre` — DOI `10.21105/joss.03125` resolves to "City2BA" (Konolige).
   **Correct DOI = `10.21105/joss.03416`** (OGRe, Shoshany 2021). → fix DOI.
3. `capozziello2022comparing` — DOI `…10052-022-10813-z` resolves to a Lake paper.
   **Correct DOI = `10.1140/epjc/s10052-022-10823-x`** (digit transposition). → fix DOI.
4. `hwangnoh2024nonlinear` — local title "Nonlinear gauge-invariant cosmological
   perturbation in Einstein-Hilbert gravity" but DOI `10.1016/j.dark.2024.101534`
   = "Graviton–photon conversions in Euler–Heisenberg nonlinear electrodynamics".
   Title/DOI mismatch → confirm intended paper, fix title or DOI.

**Title updates (review):**
5. `ruchlin2018nrpy` — DOI is the SENR/NRPy+ PRD paper; canonical title
   "SENR/NRPy+: Numerical relativity in singular curvilinear coordinate systems".
   Local "NRPy+: A Code Generation System for Numerical Relativity" → update title.
6. `barker2025irfound2` — arXiv 2507.05349 title revised on INSPIRE to
   "…II: Catalogue of all torsion-like theories including new ghost-tachyon-free
   cases". Local has older subtitle → update to published title.

**Benign (likely keep as-is):**
7. `blagojevicHehl2013` — year 2013 (ICP book) vs Crossref 2011 (WS edition). Keep 2013.
8. `lidemmel2003superlu` — Crossref returned truncated title "SuperLU_DIST"; local full
   title is correct. False positive (truncation). No change.
9. `heisenberg1936consequences` — local English title vs canonical original German
   "Folgerungen aus der Diracschen Theorie des Positrons". Translated title; keep
   (or switch to original German). No mis-attribution.

### PENDING (12) — Google Books transient during fail-fast; gentle re-run to resolve

All books w/ ISBN: golub2013matrix, trefethen2005spectra, trefethen1997numerical,
brenan1996numerical, blagojevic2002gravitation, higham2008functions,
henneaux1992quantization, leveque2007finite, davis2006sparse, hobson2006general,
misner1973gravitation, MacKay:2003. → re-run `compare.py` (gentle) to fetch via
Google Books, or verify manually.

### NOTFOUND / NOID (20) — manual verification

Pre-arXiv / no-DOI papers (gertsenshtein1962wave, sciama1962physical,
trautman1972einstein, ostrogradsky1850memoires, zakharov1970linearized,
yoshida1990symplectic, fornberg1988generation, almohy2011action), books without a
resolvable ISBN (hairer ODE I/II, hairer2006geometric, boyd2001chebyshev,
roache1998verification, blagojevic2013gauge), and software @misc (mathematica,
mekki/randall skitsundae, baratta2023dolfinx, martingarcia2007xact,
peeters2018cadabra2). → manual cross-check; most are fine as-is.

### Unmapped journals (33) — extend journal_map.json

All numerics / CS / stats journals not yet in the map (Computer Physics
Communications, ACM TOMS, SIAM Review / J. Sci. Comput. / J. Numer. Anal. / J.
Matrix Anal., Acta Numerica, J. Comput. Phys., Mathematics of Computation, JOSS,
Journal of the ACM, Physica D, Bayesian Analysis, Can. J. Math., Scholarpedia, …)
plus a few already-abbreviated ones. → extend journal_map with standard
abbreviations for consistency (these are non-HEP, so INSPIRE doesn't abbreviate
them; pick the conventional short form).

## Task 2 — citation-usage audit (192 distinct keys)

Layered verification (abstract/intro/conclusion → escalate if inconclusive).
Offline TeX available for 57 keys; 25 keys already mapped in
`literature_content_notes.md`.

| group | keys | verified |
|-------|-----:|:--------:|
| (populated from extract_cites.py / literature_content_notes groups) | — | ☐ |

### Citation MISMATCH list (paper does not support the claim)

_To be populated._

### Uncited entries (39) — candidate missing citations

39 entries defined but never `\cite`d (0 undefined refs). Codebase grep
(`tidal/`, `docs/`, `examples/`, `scripts/`) for arXiv/DOI surfaced these
**code/doc-referenced** candidates that likely belong in the manuscript:

| key | referenced in | suggested manuscript home |
|-----|---------------|---------------------------|
| `berlin2024axionphoton` | `examples/torsion_dark_photon*/theory.toml`, `docs/tex/gertsenshtein_plasma.tex`, NEXT_PHASES, supervisor mtg | dark-photon / axion-photon mixing (§1.1 or §4.2) |
| `kushwaha2024nonminimal` | `docs/tex/gertsenshtein.tex` | non-minimal coupling (§2.3 / §5.2 path 1) |
| `dunne1999chernsimons` | `docs/tex/chern_simons.tex` | Chern-Simons (if CS appears in manuscript) |
| `harris2020numpy` | `tidal/cli/_cite.py` (project depends on it) | software/methods citation (§3 / App) |
| `Ejlli:2019bqj` | `docs/references.md` only | Gertsenshtein follow-up — review |

`_cite.py` recommends citing NumPy/SciPy/SUNDIALS/xAct for any TIDAL use →
confirm these are cited in the methods/software section.

**Remaining ~34 uncited** (Barker corpus `barker2020h0/horndeski/projective`,
`hobson2021conformal`, `lin2020weylghost`; numerics/textbooks `berenger1994pml`,
`courant1928*`, `moler1973algorithm`, `trefethen2005spectra`, `davis2006sparse`,
`roache1998verification`, `maggiore2007gravitational`, `misner1973gravitation`;
software `baratta2023dolfinx`, `oskooi2010meep`, `peeters2007cadabra/2018cadabra2`,
`ruchlin2018nrpy`, `shoshany2021ogre`; instrument `abbott2023gwtc3`;
GW-frontier `Goryachev/Lasky/Paoletti/He/Campbell`): no arXiv/DOI hit in code —
**decision pending with user**: cite where relevant vs comment out `verdict=UNCITED`.

## Decision log (append-only, dated)

- **2026-05-30** — Plan approved. Phase 0 tooling built (`scripts/bibaudit/`).
  Decisions: enrichment pre-authorized, SUBSTANTIVE→STOP; journals→INSPIRE-abbrev;
  light normalization; Task 2 per-key accreting `%@cite` blocks; academic email
  `wr286@cam.ac.uk` for Crossref polite pool.
