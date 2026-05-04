# MSci Report — Phase 0 & Phase 1 Checklist

Live tracking document. Tick items as completed. See `report_plan.md` for full context.

---

## Pre-execution (done first, before downloading any papers)

- [x] Write `manuscript/planning/report_plan.md` — repo-local archive of the full plan
- [x] Write `manuscript/planning/report_checklist.md` — this file
- [x] Write `manuscript/planning/literature_content_notes.md` — per-paper notes from reading all 23 locally-available papers
- [x] Update `docs/MEMORY.md` with pointer to planning files
- [x] Commit: `chore(report): add report plan archive and Phase 0/1 checklist` (582685e)

---

## Phase 0 — Fetch references and freeze style intel

### Tier 1 — Barker software-package papers (highest priority)

- [x] Download `literature/2406.09500/` — PSALTer original (125 .tex files)
- [x] Download `literature/2506.02111/` — PSALTer v2 parity-violating (26 .tex files)
- [x] Download `literature/2512.25007/` — Hamilcar (103 .tex files)

### Tier 2 — BHL joint corpus on PGT/torsion (highest priority for main body)

- [x] Verify `literature/2406.12826/` — *Every PGT is conformal* (4 .tex files, already local)
- [x] Download `literature/2309.14783/` — *Manifestly covariant variational principle* (1 .tex, single-file gz)
- [x] Download `literature/2303.11094/` — *Gravitational confinement and galactic rotation curves* (5 .tex files)
- [x] Download `literature/2101.02645/` — *Nonlinear Hamiltonian analysis, quadratic torsion I* (4 .tex files)
- [x] Download `literature/2006.03581/` — *PGT cosmology → Horndeski* (2 .tex files)
- [x] Download `literature/2003.02690/` — *H₀ tension with emergent dark radiation* (14 .tex files)

### Tier 3 — Hobson + Lasenby joint, close-topic

- [x] Download `literature/2008.09053/` — *Fresh perspective on gauging the conformal group* (1 .tex, single-file gz)
- [x] Download `literature/2005.02228/` — *Ghost and tachyon free Weyl gauge theories* (1 .tex file)

### Tier 4 — Barker solo voice

- [x] Download `literature/2205.13534/` — solo-author multipliers paper (2 .tex files)
- [x] Download `literature/2311.11790/` — solo-author Yang-Mills gravity (2 .tex files)

### Tier 5 — Additional Barker papers (recent metric-affine and PGT work)

- [x] Download `literature/2510.17094/` — Gertsenshtein effect on curved spacetime (1 .tex, DIRECT relevance)
- [x] Download `literature/2507.09228/` — TorC torsion condensation cosmology (Handley+Barker, 2 .tex)
- [x] Verify `literature/2510.08201/` — R² Hamiltonian analysis (already local, 1 .tex)
- [x] Download `literature/2506.21662/` — Infrared foundations I: rank-3 field theories (5 .tex)
- [x] Download `literature/2507.05349/` — Infrared foundations II: torsion-like + ghost analysis (5 .tex)
- [x] Download `literature/2505.23894/` — Can metric-affine gravity be saved? (2 .tex)
- [x] Download `literature/2402.07641/` — Particle spectra of Palatini/metric-affine (2 .tex)
- [x] Download `literature/2402.14917/` — Metric-affine gravity from extended projective symmetry (4 .tex)

### Tier 6 — Key collaborator / field papers

- [x] Download `literature/2406.11956/` — Weyl-invariant Einstein-Cartan gravity (Karananas et al., 1 .tex)
- [x] Download `literature/2506.17017/` — Cosmology of Cubic PGT (Bahamonde et al., 1 .tex)

### Already-in-`literature/` papers that must be studied for CONTENT (not yet analysed)

These are local — no download needed, but must be read for §1, §2, §4, §5:
- [x] Read and note `literature/1812.02675/` — Lin, Hobson, Lasenby: ghost-free PGT classification (§1.2, §5.2)
- [x] Read and note `literature/1804.05556/` — Blagojević & Cvetković PGT Hamiltonian structure (§5.2)
- [x] Read and note `literature/2301.02072/` — Simple derivation of Gertsenshtein effect (§1.1, §4)
- [x] Read and note `literature/2310.04150/` — EM field definitions in graviton-photon conversion (§2.3)
- [x] Read and note `literature/2004.02714/` — Exact graviton-photon mixing in B-field (Ejlli, §4)
- [x] Read and note `literature/2105.04565/` — Dark Photon Limits Handbook (Caputo et al., §4.2)
- [x] Read and note `literature/hep-th_0103093/` — Shapiro: Physical aspects of torsion (§2 background)
- [x] Read and note `literature/gr-qc_0001010/` — Hehl & Obukhov: How does EM couple to gravity? (§2)
- [x] Read and note `literature/gr-qc_0307063/` — Itin & Hehl: Maxwell + quadratic torsion (§2, §5)

### Pre-arXiv must-cite papers (bibliography entries only)

These need entries in `references.bib` — no TeX to download:
- [x] `gertsenshtein1962wave` — already present in references.bib
- [x] `boccaletti1970conversion` — already present in references.bib
- [x] `raffelt1988mixing` — already present in references.bib
- [x] `kibble1961lorentz` — added (J. Math. Phys. 2, 212 (1961))
- [x] `hehl1976spin` — already present in references.bib
- [x] `sezgin1980ghostfree` — added (PRD 21, 3269 (1980))
- [x] `cillis1996graviton` — added (PRD 54, 4757 (1996))

### Style intel artefacts

- [x] Create `manuscript/style_intel/` directory
- [x] Write `manuscript/style_intel/style_analysis.md` — 117 lines, covers all 10 patterns (a–j), analytical prose
- [x] Write `manuscript/style_intel/genre_conventions.md` — 108 lines, covers (a–g) + cross-paper evolution notes
- [x] Create `manuscript/STYLE_GUIDE.md` — 96 lines, actionable checklist with 9 sections + self-review list

### Commit

- [x] Commit: `chore(report): freeze style analysis and supervisor reference corpus for MSci draft`

---

## Phase 1 — Outline approval, notation set, bibliography

### Supervisor sign-off

- [ ] Send outline (from `report_plan.md` § "Proposed report outline") to Barker for sign-off
- [ ] Receive sign-off (or incorporate feedback and re-send)

### Bibliography hygiene

- [x] Add `barker2024psalter` (2406.09500) to `manuscript/references.bib`
- [x] Add `barker2025psalter2` (2506.02111) to `manuscript/references.bib`
- [x] Add `barker2025hamilcar` (2512.25007) to `manuscript/references.bib`
- [x] Add software comparators: NRPy+, OGRe, Cadabra2 (×2 entries); Dedalus already present; xAct/xPert already present
- [x] Add dark-photon / hidden-photon review (`caputo2021darkphoton`, PRD 104, 095029 (2021))
- [x] Add ghost-instability PGT literature: `sezgin1980ghostfree` (PRD 21, 3269); `lin2019pgt` (Lin-Hobson-Lasenby 2019); `blagojevic2018hamiltonian` — covers the key no-ghost constraints
- [x] Add GW-detection context anchors: `amaro2017lisa` (LISA); `nanograv2023fifteen` (NANOGrav 15-yr); `abbott2023gwtc3` (GWTC-3 / LIGO O4 context)
- [x] Verify entry count ≥ 55 — **81 entries** total

### Notation survey and macros

- [x] Read Tier 1–2 papers focusing on: metric signature, vierbein/spin-connection notation, torsion irreducible decomposition labels, Gertsenshtein coupling symbols
- [x] Record standard choices + source traceback in `manuscript/macros.tex` comments (see §"Agreed notation set" block added to macros.tex)
- [x] Agree macro set for all report symbols: metric (-,+,+,+); indices μ,ν,ρ (spacetime) / a,b,c (frame); e^a_μ vierbein; ω^{ab}_μ spin connection; T_μ torsion trace; S^μ axial; q_{μνρ} tensor; κ²=16πG; P=sin²(κB₀D/2)
- [x] Check for conflicts with `docs/tex/preamble.tex` — no conflicts; both use `\varkappa` for κ, same torsion/contorsion names

### Commit

- [ ] Commit: `chore(report): Phase 1 — bibliography 81 entries, notation agreed in macros.tex, literature notes`

---

## Phase 0/1 complete when

- All Tier 1–4 papers present in `literature/`
- `manuscript/STYLE_GUIDE.md`, `manuscript/style_intel/style_analysis.md`, `manuscript/style_intel/genre_conventions.md` all exist
- `references.bib` ≥ 55 entries
- `manuscript/macros.tex` has agreed notation set with literature-traceable comments
- Outline sign-off received from Barker

**Phase 2 (drafting) is a separate future session — do not begin writing report prose until Phase 0/1 are complete.**
