# MSci Report — Phase 0 & Phase 1 Checklist

Live tracking document. Tick items as completed. See `report_plan.md` for full context.

---

## Pre-execution (done first, before downloading any papers)

- [x] Write `manuscript/planning/report_plan.md` — repo-local archive of the full plan
- [x] Write `manuscript/planning/report_checklist.md` — this file
- [ ] Update `docs/MEMORY.md` with pointer to planning files
- [ ] Commit: `chore(report): add report plan archive and Phase 0/1 checklist`

---

## Phase 0 — Fetch references and freeze style intel

### Tier 1 — Barker software-package papers (highest priority)

- [ ] Download `literature/2406.09500/` — PSALTer original (Barker, Marzo, Rigouzzo 2024)
- [ ] Download `literature/2506.02111/` — PSALTer v2 parity-violating (Barker, Karananas, Tu 2025)
- [ ] Download `literature/2512.25007/` — Hamilcar (Barker sole author 2025)

### Tier 2 — BHL joint corpus on PGT/torsion (highest priority for main body)

- [x] Verify `literature/2406.12826/` — *Every PGT is conformal* (already local; check TeX completeness)
- [ ] Download `literature/2309.14783/` — *Manifestly covariant variational principle*
- [ ] Download `literature/2303.11094/` — *Gravitational confinement and galactic rotation curves*
- [ ] Download `literature/2101.02645/` — *Nonlinear Hamiltonian analysis, quadratic torsion I*
- [ ] Download `literature/2006.03581/` — *PGT cosmology → Horndeski*
- [ ] Download `literature/2003.02690/` — *H₀ tension with emergent dark radiation*

### Tier 3 — Hobson + Lasenby joint, close-topic

- [ ] Download `literature/2008.09053/` — *Fresh perspective on gauging the conformal group*
- [ ] Download `literature/2005.02228/` — *Ghost and tachyon free Weyl gauge theories*

### Tier 4 — Barker solo voice

- [ ] Download `literature/2205.13534/` — solo-author multipliers paper (2022)
- [ ] Download `literature/2311.11790/` — solo-author Yang-Mills gravity (3 pp)

### Style intel artefacts

- [ ] Create `manuscript/style_intel/` directory
- [ ] Write `manuscript/style_intel/style_analysis.md` — analytical write-up of Barker + BHL writing patterns (abstract, intro scaffold, section structure, equation handling, voice, transitions). Analytical, not extractive. Original prose only.
- [ ] Write `manuscript/style_intel/genre_conventions.md` — analytical synthesis of physics-software-paper genre conventions (from 16-paper survey). Used to plan appendices A–E.
- [ ] Create `manuscript/STYLE_GUIDE.md` — concrete checklists (section-opening rule, equation rule, label rule, voice rule, figures rule, forbidden-words list, self-review checklist)

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
