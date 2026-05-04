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
