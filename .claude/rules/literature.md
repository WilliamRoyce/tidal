---
description: Literature references and citation guidance for physics papers
paths:
  - "literature/**"
  - "docs/references.md"
  - "docs/tex/**"
---

# Literature & References

When citing or referencing physics papers:
1. Check `literature/` first — 20 arXiv TeX sources stored locally (but never limit research to local papers)
2. Read the TeX source directly (faster than WebFetch)
3. Key papers available locally:
   - Gertsenshtein: 2301.02072 (Palessandro & Rothman, "A simple derivation of the Gertsenshtein effect", Phys.Dark Univ.40 101187 — note: WRONG formula P=sin²(√G B₀ D), off by √(4π); see docs/tex/gertsenshtein_formula.tex §4), 2303.07562 (Hwang-Noh, "Definition of electric and magnetic fields in curved spacetime", Annals Phys. 454 169332), 2310.04150 (Hwang-Noh, "On graviton-photon conversions in magnetic environments", Phys.Dark Univ.43 101426 — Gertsenshtein-specific application; cites Palessandro-Rothman; documents that naive EM identification is correct in TT gauge only), 2004.02714 (Ejlli, exact solution), 2405.01407 (Palessandro solo 2024, "Graviton-photon oscillations as a probe of quantum gravity", CQG 41 215011 — builds on P&R 2023, perpetuates same √(4π) error in weak-field limit), 2406.17853 (Lella, Calore, Carenza, Mirizzi, "Constraining gravitational-wave backgrounds...", Phys.Rev.D 110 083042 — agrees with TIDAL formula via Δ_{gγ}=B_T/(√2 M_P)=κ/2), 2507.16609 (Domcke-Garcia-Cely-Lee, GW scattering on magnetic fields, JCAP), 2302.08186 (modified-gravity Gertsenshtein), 2405.11786 (Hwang-Noh, EH nonlinear)
   - Torsion: 2406.12826 (Barker, PGT conformal), 2410.01355 (Obukhov, spin-torsion), hep-th/0103093 (Shapiro, spacetime torsion), 0905.3732 (Nikiforova, torsion stability), 2507.02362 (Bahamonde, torsion-EM coupling)
   - Non-minimal torsion-EM: gr-qc/0305049 (Rubilar, torsion birefringence), gr-qc/0307063 (Itin, T²F² classification), gr-qc/0001010 (Hehl-Obukhov, EM-gravity coupling)
   - Ghost-free PGT: 1812.02675 (Blagojevic, 450 ghost-free cases), 2009.11739 (Aoki, ghost-free via bigravity), 2506.17017 (Bahamonde, cubic PGT ghost elimination)
   - Related: 2405.08865 (Berlin, axion-photon numerical methods)
4. For new papers: download TeX from arXiv, extract to `literature/<id>/`
5. Update `docs/references.md` when adding new references
