# `tidalcosmo/perturbations/` — assembly of the coupled block

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

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
