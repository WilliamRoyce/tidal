# MSci Report — Phase 0 & Phase 1 Checklist

Live tracking document. Tick items as completed. See `report_plan.md` for full context.

---

## Pre-execution (done first, before downloading any papers)

- [x] Write `manuscript/planning/report_plan.md` — repo-local archive of the full plan
- [x] Write `manuscript/planning/report_checklist.md` — this file
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
- [ ] Read and note `literature/1812.02675/` — Nair-Park-Yoon ghost-free PGT classification (§1.2, §5.2)
- [ ] Read and note `literature/1804.05556/` — Blagojević & Cvetković PGT Hamiltonian structure (§5.2)
- [ ] Read and note `literature/2301.02072/` — Simple derivation of Gertsenshtein effect (§1.1, §4)
- [ ] Read and note `literature/2310.04150/` — EM field definitions in graviton-photon conversion (§2.3)
- [ ] Read and note `literature/2004.02714/` — Exact graviton-photon mixing in B-field (Ejlli, §4)
- [ ] Read and note `literature/2105.04565/` — Dark Photon Limits Handbook (Caputo et al., §4.2)
- [ ] Read and note `literature/hep-th_0103093/` — Shapiro: Physical aspects of torsion (§2 background)
- [ ] Read and note `literature/gr-qc_0001010/` — Hehl & Obukhov: How does EM couple to gravity? (§2)
- [ ] Read and note `literature/gr-qc_0307063/` — Itin & Hehl: Maxwell + quadratic torsion (§2, §5)

### Pre-arXiv must-cite papers (bibliography entries only)

These need entries in `references.bib` — no TeX to download:
- [ ] Add `gertsenshtein1962` — M.E. Gertsenshtein, Sov. Phys. JETP 14, 84 (1962)
- [ ] Add `boccaletti1970` — Boccaletti et al., Nuovo Cim. B (1970) [the formula we validate to 0.04%]
- [ ] Add `raffelt1988` — Raffelt & Stodolsky, PRD 37, 1237 (1988)
- [ ] Add `kibble1961` — T.W.B. Kibble, J. Math. Phys. 2, 212 (1961)
- [ ] Add `hehl1976` — Hehl et al., Rev. Mod. Phys. 48, 393 (1976)
- [ ] Add `sezgin1980` — Sezgin & Van Nieuwenhuizen, PRD 21, 3269 (1980)
- [ ] Add `cillis1996` — Cillis & Harari, PRD 54, 4757 (1996)

### Style intel artefacts

- [x] Create `manuscript/style_intel/` directory
- [x] Write `manuscript/style_intel/style_analysis.md` — 117 lines, covers all 10 patterns (a–j), analytical prose
- [x] Write `manuscript/style_intel/genre_conventions.md` — 108 lines, covers (a–g) + cross-paper evolution notes
- [x] Create `manuscript/STYLE_GUIDE.md` — 96 lines, actionable checklist with 9 sections + self-review list

### Commit

- [ ] Commit: `chore(report): freeze style analysis and supervisor reference corpus for MSci draft`

---

## Phase 1 — Outline approval, notation set, bibliography

### Supervisor sign-off

- [ ] Send outline (from `report_plan.md` § "Proposed report outline") to Barker for sign-off
- [ ] Receive sign-off (or incorporate feedback and re-send)

### Bibliography hygiene

- [ ] Add `barker2024psalter` (2406.09500) to `manuscript/references.bib`
- [ ] Add `barker2025psalter2` (2506.02111) to `manuscript/references.bib`
- [ ] Add `barker2025hamilcar` (2512.25007) to `manuscript/references.bib`
- [ ] Add software comparators: NRPy+, Dedalus, OGRe, FieldsX, Cadabra2 (check if xAct/xPert already present)
- [ ] Add dark-photon / hidden-photon review (Caputo–Millar–O'Hare–Vitagliano 2021 or equivalent)
- [ ] Add ghost-instability PGT literature (Sezgin–Van Nieuwenhuizen; Beltrán-Jiménez–Heisenberg–Koivisto; Nair–Park–Yoon)
- [ ] Add GW-detection context anchors (LISA, NANOGrav, LIGO O4 — one each)
- [ ] Verify entry count ≥ 55

### Notation survey and macros

- [ ] Read Tier 1–2 papers focusing on: metric signature, vierbein/spin-connection notation, torsion irreducible decomposition labels, Gertsenshtein coupling symbols
- [ ] Record standard choices + source traceback in `manuscript/macros.tex` comments
- [ ] Agree macro set for all report symbols (metric, indices, torsion components, contortion, irreducible decomp, Gertsenshtein constants)
- [ ] Check for conflicts with `docs/tex/preamble.tex`; update docs to match report notation where they disagree

### Commit

- [ ] Commit: `chore(report): outline approved, notation agreed, bibliography extended`

---

## Phase 0/1 complete when

- All Tier 1–4 papers present in `literature/`
- `manuscript/STYLE_GUIDE.md`, `manuscript/style_intel/style_analysis.md`, `manuscript/style_intel/genre_conventions.md` all exist
- `references.bib` ≥ 55 entries
- `manuscript/macros.tex` has agreed notation set with literature-traceable comments
- Outline sign-off received from Barker

**Phase 2 (drafting) is a separate future session — do not begin writing report prose until Phase 0/1 are complete.**
