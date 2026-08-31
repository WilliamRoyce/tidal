# `tidalcosmo/observables/` — transfer functions, projection, and rotations

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

**Responsibility.** Turning solved perturbations into something a likelihood consumes. Three
files for the three mechanisms of `repo_reshape.md` §2.4:

- `transfer.py` — **mechanism 2**: we compute our own `Δ_ℓ(k)` for a channel and form `C_ℓ` from
  it, leaving the projection unchanged.
- `los.py` — **mechanism 3**: our own line-of-sight integral, `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))`,
  for when the source function itself changes shape. `adammoss/nanoCMB` (MIT) is the reference
  implementation — ~1400 readable lines, sub-percent against CAMB over `2 ≤ ℓ ≤ 2500`. Adapt with
  attribution; prefer read-and-reimplement.
- `rotation.py` — **mechanism 1**: effects expressible as a term in `S` built from quantities CAMB
  already evolves, emitted as sympy and compiled by `set_custom_scalar_sources`. Note
  `source_ell_scales = 2` is documented by CAMB as "a new polarization-like source" — the spin-2
  normalization, directly relevant since birefringence and V-modes are polarization observables.

**Workstream.** WS4 (#493). **Filled at.** M5.

**Owed, and currently unowned.** The **per-channel source functions** — given a solved
`(h, torsion)` history, what quantity is integrated along the line of sight to produce a B-mode?
Deriving `S` is real work per channel, it is not in legacy (flat space has no line of sight) and
not free from CAMB (its sources are ΛCDM's), and by project rule it should be derived
**symbolically in Wolfram**. See `repo_reshape.md` §2.5.

**Caution.** Post-processing a rotation is exact only for constant, isotropic,
frequency-independent effects; otherwise it belongs inside the LOS integral (arXiv:2209.07804).

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
