# `tidalcosmo/perturbations/` — assembly of the coupled block

> **Directory boundary from H4 (2026-08-31); internals settled by H3.** The H4 study
> drew this boundary *alongside* the detailed investigations rather than after them, so the
> boundary itself may still be revised — that does not require re-litigating H4. **The
> contents below are not provisional**: they record decisions settled in
> `docs/cosmology/solver_design.md` §4. Do not reopen them from this file.

**Responsibility.** Building `M(η, k)` — the small linear system the solver integrates — from the
spec plus the background.

**What "the coupled block" means, precisely** (`repo_reshape.md` §2.2): **our fields *plus* every
standard mode they couple to.** Not just the new sector. If torsion modifies how photon
polarization evolves, then CAMB's photon polarization is *wrong for that channel*, and that mode
must be evolved here rather than read from CAMB. Getting this wrong is a silent correctness bug,
which is why it has its own section in the design document.

**Workstream.** WS2/WS3. **Filled at.** M4.

**Open question.** The state layout and `M⁻¹` handling in legacy `solver/{state,fields,_kinetic,
kinetic_matrix}.py` hold reusable ideas, but their **structure is explicitly open to redesign** —
they accreted over time, and H3's architecture decides the shape here.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.

> **Amendment (coherence pass, 2026-09-04):** the open question above — *"H3's architecture decides the shape here"* — is **answered**
> (`docs/cosmology/solver_design.md`, 2026-08-31). H3 chose a clean design with **nothing
> ported** from legacy `modal.py`: steppers consume **node-sampled arrays**, never evaluator
> objects, and the state layout follows from the segmented-assembly contract (**#518**) rather
> than from `solver/{state,fields}.py`.
