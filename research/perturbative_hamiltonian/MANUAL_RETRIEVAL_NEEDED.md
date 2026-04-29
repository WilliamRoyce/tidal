# Papers Cited in the Investigation — Retrieval Status

**Last updated**: 2026-04-27 — after user manual retrieval + agent download pass.

**Status**: Most cited papers are now accessible. Outstanding items either truly
need university-library access (pre-arXiv journal-only) or were never on arXiv.

These need **manual retrieval** before any of the proposed publications goes out,
because:
- The agents have been transcribing claims from these papers based on second-hand
  summaries (memory files, review articles, Wikipedia, ar5iv abstracts).
- Some agent claims about these papers may not survive a direct read.
- The "verbatim quote" hallucination risk is highest for papers the agents could
  only access at second-hand.

## Tier 1 — Critical-evidence papers (status check)

| # | Reference | Status | Location |
|---|-----------|--------|----------|
| 1 | **Blagojević & Nikolić 1983** (Phys. Rev. D **28**:2455) | ✅ User-retrieved | `literature/PhysRevD.28.2455/` |
| 2 | ~~**Blagojević & Nikolić 1984** (Nuovo Cim. B **84**:25)~~ → **CITATION CORRECTED**: actual paper is M. Blagojević and I. A. Nikolić, "Hamiltonian Structure of the Theory of Gravity with $R + T^2$ Type of Lagrangian", Nuovo Cimento B **73**(2):258–273 (1983), DOI: 10.1007/BF02721794 | ✅ User-retrieved | `literature/BF02721794/BF02721794.pdf`. Paper is on $R + T^2$ PGT (linear $R$, quadratic torsion) — does NOT contain the $R^2/\tilde{R}^2$ class that TIDAL needs. |
| 3 | **Yo, Nester & Ni 1999** (gr-qc/9902032) | ✅ Downloaded | `literature/gr-qc_9902032/main.tex` |
| 4 | **Yo & Nester 2002** (gr-qc/0112030) | ✅ Downloaded | `literature/gr-qc_0112030/main.tex` |
| 5 | **Blagojević & Cvetković 2018** (arXiv:1804.05556) | ✅ Downloaded; Appendix D verified | `literature/1804.05556/main.tex`. Verbatim quote confirmed at PDF lines 2685-2688 by Meta-Review K. Appendix D contains BOTH limit-no-go AND a constructive F-matrix method (Meta-Review N). |
| 6 | **Hehl, McCrea, Mielke & Ne'eman 1995** (Phys. Rep. **258**:1) | ✅ User-retrieved + downloaded arXiv version | `literature/PhysRep258.1/` (user) and `literature/gr-qc_9402012/main.tex` (arXiv) |

## Tier 2 — Pre-arXiv classics (genuinely require university/library access)

These papers are cited as historical/contextual background but no agent claim depends
critically on the verbatim text. **Tier 2 retrieval is convenience-only** — the
investigation's verdicts don't hinge on these.

If you DO want to retrieve them: most are available via APS Digital Library (PRD,
PRL with backfile to 1970s), Springer (Nuovo Cimento, Nucl. Phys. B), or Princeton
University Press (Henneaux-Teitelboim book). Italian Phys. Soc. (Nuovo Cim. B) may
require direct journal access.

| # | Reference | Source | Access path |
|---|-----------|--------|-------------|
| 7 | **Hayashi & Shirafuji 1979** | Prog. Theor. Phys. **64**:866 (+ sequels) | Oxford Academic (PTP backfile) |
| 8 | **Sezgin & Van Nieuwenhuizen 1980** | Phys. Rev. D **21**:3269 | APS Digital Library |
| 9 | **Stelle 1977** | Phys. Rev. D **16**:953 | APS Digital Library |
| 10 | **Faddeev & Jackiw 1988** | Phys. Rev. Lett. **60**:1692 | APS Digital Library |
| 11 | **Boulware & Deser 1972** | Phys. Rev. D **6**:3368 | APS Digital Library |
| 12 | **van Dam & Veltman 1970** | Nucl. Phys. B **22**:397 | Elsevier ScienceDirect (Nucl Phys B backfile) |
| 13 | **Zakharov 1970** | JETP Lett. **12**:312 | Russian physics journal; American translation via JETP Letters archives |
| 14 | **Lee & Wald 1990** | J. Math. Phys. **31**:725 | AIP Publishing |
| 15 | **Crnković & Witten 1987** | "Three Hundred Years of Gravitation" (Hawking-Israel eds.) | Cambridge University Press book chapter |
| 16 | **Henneaux & Teitelboim 1992** | "Quantization of Gauge Systems" | Princeton University Press — physical book |
| 17 | **Jaen, Llosa & Molina 1986** (JLM) | Phys. Rev. D **34**:2302 | APS Digital Library |
| 18 | **Cheyette 1988** | Nucl. Phys. B **297**:183 | Elsevier ScienceDirect |
| 19 | **Eliezer & Woodard 1989** | Nucl. Phys. B **325**:389 | Elsevier ScienceDirect |
| 20 | **Krupka book 2015** | "The Inverse Problem of the Calculus of Variations" (Atlantis) | Atlantis Press; Springer eventually |

## Tier 3 — Mentioned but not load-bearing

These were mentioned in agent writeups but no agent claim depends critically on them.
Lower-priority retrieval.

| # | Reference | Source |
|---|-----------|--------|
| 21 | **Helmholtz 1887** | J. Reine Angew. Math. **100**:137 | Original Helmholtz conditions for the inverse problem. 1887 — physical archive only. |
| 22 | **Tonti 1969** | Italian math journal | Vainberg-Tonti formula. |
| 23 | **Vainberg 1964** | Russian math journal | Vainberg-Tonti formula. |
| 24 | **Marmo, Saletan, Simoni & Vitale 1985** | Annals Phys. **163**:204 | Helmholtz for constrained systems. Pre-arXiv. |
| 25 | **Bender & Mannheim 2008** | Phys. Rev. Lett. **100**:110402 (arXiv:0706.0207) | Pais-Uhlenbeck quantization. Available on arXiv. |
| 26 | **Hossenfelder, Mistry & Padilla** | (assorted) | Auxiliary-field higher-derivative gravity. Several papers; none load-bearing. |
| 27 | **Bopp 1940 / Podolsky 1942** | Various | Bopp-Podolsky electrodynamics origin. Pre-arXiv. |
| 28 | **Stueckelberg 1938** | Helv. Phys. Acta **11**:225 / **30**:209 | Original Stückelberg paper. Pre-arXiv. |
| 29 | **Ostrogradsky 1850** | Mem. Acad. St. Petersbourg | Original Ostrogradsky theorem. 1850 — historical archive. |

## Modern papers behind paywalls (post-arXiv but journal-only access)

| # | Reference | Source |
|---|-----------|--------|
| 30 | **Phys. Lett. B 821:136608** (Lyakhovich-Sharapov 2021) | Elsevier paywall; arXiv version is arXiv:2106.09355 (Abakumova-Lyakhovich) |
| 31 | **JHEP 11(2024)146** (Karananas et al.) | Springer; check arXiv preprint |
| 32 | **Phys. Rev. D 97 (2018) 016016** (Brambilla-Soto-Vairo) | APS paywall |

## arXiv-cached papers (now downloaded to `literature/`)

All of these are now in the local literature cache as TeX source. They had been
hard for the agents to access during the original investigation but are
now verifiable directly.

| arXiv ID | Title | Local path |
|----------|-------|-----------|
| 1804.05556 | Blagojević-Cvetković 2018 | `literature/1804.05556/main.tex` |
| 1406.6646 | Krupka-Voicu 2015 (Vainberg-Tonti) | `literature/1406.6646/main.tex` |
| 2009.05459 | Voicu 2020 (4D-GB) | `literature/2009.05459/varcompl_*.tex` |
| 2403.15564 | Ducobu-Voicu 2024 (variational bootstrapping I) | `literature/2403.15564/main.tex` |
| 2406.09540 | Ducobu-Voicu 2025 (variational bootstrapping II) | `literature/2406.09540/part2_R1.tex` |
| 2102.10579 | Lyakhovich 2021 (general Stückelberg) | `literature/2102.10579/main.tex` |
| 2106.09355 | Abakumova-Lyakhovich 2021 (reducible Stückelberg) | `literature/2106.09355/main.tex` |
| 1508.02401 | Hinterbichler-Saravani 2015 | `literature/1508.02401/curvaturesquaregravityarxiv.tex` |
| 2411.16928 | Chatzistavrakidis-Ranjbar-Zekoč 2024 | `literature/2411.16928/Stueckelberg_v2.tex` |
| 2206.00658 | Barker HiGGS 2022 | `literature/2206.00658/apstemplate.tex` |
| 2210.15980 | Chiou-Geiller-Wang 2022 | `literature/2210.15980/main.tex` |
| 2602.12114 | Chan-López et al. 2026 (modern FJ Schur-bordering) | `literature/2602.12114/draft_sfj.tex` |
| 2601.22007 | Aashish-Saif 2026 (Stückelberg antisymm tensor) | `literature/2601.22007/main.tex` |
| 2302.03545 | Hou-Cai-Li 2023 f(T) Stückelberg | `literature/2302.03545/main.tex` |
| 2410.03422 | Cai-Saridakis 2024 lowering strong coupling | `literature/2410.03422/main.tex` |
| 2311.17459 | Glavan-Zlosnik-Lin 2024 metric-affine R² (Meta-N forward-citation lead) | `literature/2311.17459/Hamiltonian_metric-affine_R2.tex` |
| 2503.18972 | Beltrán Jiménez et al. f(R,Q) | `literature/2503.18972/f_R,Q_.tex` |
| gr-qc/9902032 | Yo-Nester-Ni 1999 ("constraint bifurcation") | `literature/gr-qc_9902032/main.tex` |
| gr-qc/0112030 | Yo-Nester 2002 (higher-spin extension) | `literature/gr-qc_0112030/main.tex` |
| hep-th/0302033 | Loran 2003 (irregular constraints classification) | `literature/hep-th_0302033/main.tex` |
| hep-th/0301256 | Miskovic-Zanelli 2003 (irregular Hamiltonian systems) | `literature/hep-th_0301256/main.tex` |
| 1903.02263 | Blagojević-Cvetković 2019 (entropy in PG, BC follow-up) | `literature/1903.02263/main.tex` |

## Tier 3 (low-priority) — Lyakhovich forward-citation papers (2026-04-27 audit)

The Phase 2.3 forward-citation audit (`notes/lead_lyakhovich_forward_citations.md`)
identified the following papers that cite arXiv:2106.09355 or arXiv:2102.10579.
**None is load-bearing** — the audit verdict is (b) NO NEW LEAD, all forward
citations re-apply the original recipe to duality cases. Listed here for
completeness so future agents do not re-discover them.

| arXiv ID | Title | Status |
|----------|-------|--------|
| 2206.07891 | Abakumova-Frolovsky-Herbig-Lyakhovich 2022 (linearised Nordström dual spin-2) | Direct-read 2026-04-27, cached at `/tmp/forward_cit/2206.07891/main.tex`. Self-citation; no new technique relevant to TIDAL. |
| 2207.10634 | Guleryuz 2023 (string-landscape attractors) | Background citation only; off-topic for TIDAL. |
| 2303.02616 | Abakumova-Lyakhovich 2023 (Dualisation of free fields) | Direct-read 2026-04-27, cached at `/tmp/forward_cit/2303.02616/main.tex`. Generalises 2106.09355 with multi-layer reducible Stückelberg, but still 2nd-order EoMs at fixed parameters. Strengthens (does not weaken) Round 1 Agent D's no-go. |
| 2402.12437 | Casadio-Chataignier 2024 (relaxation of first-class constraints) | Direct-read 2026-04-27, cached at `/tmp/forward_cit/2402.12437/main.tex`. Cites Lyakhovich only to *contrast* with Fock-Stueckelberg time mechanism. Off-topic for TIDAL. |
| 2405.13706 | Delplanque-Skvortsov 2024 (chiral vs symmetric massive HS) | Direct-read 2026-04-27, cached at `/tmp/forward_cit/2405.13706/main.tex`. One-shot bibliographic mention; off-topic. |

These are NOT copied to `literature/` because no TIDAL claim depends on them.
Local cache at `/tmp/forward_cit/` is ephemeral; if a future investigation needs
to revisit, re-download via `curl -sL https://arxiv.org/e-print/<id> -o <id>.tar.gz`.

## Local literature already available (TIDAL repo)

These are in `/workspaces/torsion-gertsenshtein/literature/` and accessible directly:

- arXiv:0905.3732 (Nikiforova torsion stability)
- arXiv:1710.01562 (Glavan EFT order reduction)
- arXiv:1812.02675 (Blagojević 450 ghost-free)
- arXiv:2009.11739 (Aoki-Mukohyama bigravity-PGT — local TeX)
- arXiv:2406.12826 (Barker conformal PGT)
- arXiv:2506.17017 (Bahamonde cubic PGT)
- arXiv:gr-qc_9211002 (Parker-Simon)
- arXiv:hep-th_0103093 (Shapiro torsion phenomenology)
- ... and ~25 others (see `literature/README.md`)

If a paper's arXiv ID appears in the cited list AND in `literature/`, retrieval is
trivial (just `cat /workspaces/torsion-gertsenshtein/literature/<id>/*.tex`).

## Suggested retrieval workflow

For Tier 1 papers (#1-#6):

1. **Item #5 (BC 2018 Appendix D) is the highest-priority single retrieval.** A dedicated
   meta-review agent (Meta-Review N) is currently investigating, but if it cannot access
   the appendix via WebFetch, manual PDF retrieval may be needed. The whole "25-year
   open problem" framing depends on what this appendix actually contains.

2. **Items #1, #2 (Blagojević-Nikolić 1983, 1984) need physical-archive access** unless
   the journals have digital backfile. Phys. Rev. D 1983 is on APS Digital Library;
   Nuovo Cim. B 1984 is harder (Italian Phys. Soc. — Springer eventually digitised).

3. **Items #3, #4 (Yo-Nester-Ni)** are gr-qc preprints; arXiv has them. WebFetch should
   work — if not, try the published Int. J. Mod. Phys. D versions.

4. **Item #6 (Hehl et al. 1995)** is gr-qc/9402012, arXiv-available.

For Tier 2 papers, retrieval is convenience-only — the investigation's verdict
doesn't hinge on these.

For Tier 3 papers, retrieval is historical-curiosity-only.

## Risk mitigation

Until Tier 1 papers are manually retrieved and verified:

- **Do NOT publish** any artefact that quotes BC 2018 Appendix D verbatim (the Round 1
  Finding 2 quote is unverified — Review 2 flagged as "possibly hallucinated").
- **Do NOT cite** Blagojević-Nikolić 1983/1984 as "literally the canonical Hamiltonian
  of $R + R^2 + T^2$ PGT" without verifying. The agents inferred this from second-
  hand sources.
- **Do verify** that the "25-year history" claim is reasonable: at minimum, confirm
  Blagojević & Nikolić's 1984 Nuovo Cim. paper actually addresses $R^2 + T^2$ PGT
  (not just generic PGT).

## How to add to this list

If future agent investigations cite a paper that cannot be accessed via
WebFetch / arXiv / ar5iv / TIDAL local literature, add it here with:
- Full reference (authors, year, journal, volume, page)
- Why it matters (which agent claim depends on it)
- Tier (1 = critical, 2 = important, 3 = peripheral)
