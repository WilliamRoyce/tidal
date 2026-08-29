# The spectator (test-field) route: scope investigation (H7)

**Executed:** 2026-08-29, during program planning. **Status:** complete.
**Question:** we evolve a new sector's linear perturbations on a fixed, CAMB-supplied
ΛCDM background (no backreaction), coupled to gravitons and photons, and compare imprints
to data. What is the scope of this route — who has done what, what can it deliver, and is
the niche occupied?

## 1. The methodology and its names

Same approximation, five names by subfield: **test field** (QFT in curved spacetime),
**spectator field** (inflation/early universe), **external/background field**
(axion-photon and Gertsenshtein mixing), **probe limit** (holography), **negligible
backreaction** (magnetogenesis). Canonical instances that are exactly our setup:

- **Tensor modes on FRW are the paradigm case.** `h'' + 2ℋh' + k²h = 0` is solved in
  CAMB/CLASS on a background the tensors do not influence. Everything we want to do to the
  tensor sector modifies a computation *already structured as a spectator computation*.
- **Modified GW propagation:** Amendola, Ballesteros & Pettorino arXiv:1405.7004 —
  B-mode effects of modified tensor propagation with "background and growth of matter
  perturbations standard" (the validity statement, already in print). Also
  arXiv:1408.2224 (GW friction), arXiv:1805.08731 (`d_L^gw/d_L^em`).
- **Cembranos, González Ortiz & Martín-Moruno arXiv:2302.08186** (local:
  `literature/2302.08186/conversionDef4arxiv.tex`) — the closest single template.
  Explicit spectator caveat at line 512: background evolution "assumed ΛCDM …
  discrepancies not observable at the background level." Beyond-GR graviton equation
  `h'' + 2ℋ(1+ν)h' + (c_T²k² + μ²)h = 0` with density-matrix decoherence through
  recombination/reionization. Useful negative result: graviton mass `μ ≤ 10⁻³³ eV`
  changes nothing; friction `ν` matters strongly (amplification `ν < 0`, vanishing
  `ν = −1`).
- Spectator scalars in inflation (arXiv:1302.2637), PMFs below the backreaction bound
  (arXiv:2303.04583), axion-SU(2) — the best-documented *failure* of the limit
  (arXiv:2311.07557): backreaction regions shown as excluded parameter space, the format
  we should adopt.

## 2. Validity criteria — and an open methodological hole

Two criteria used in the literature:

- **Energy budget / `ΔN_eff`:** Domcke & Garcia-Cely arXiv:2006.01161 (local), eq. (30):
  `ρ_g(T) ≤ (7/8)(4/11)^{4/3} ΔN_eff ρ_γ(T)`, `ΔN_eff ≲ 0.1` — quantitative, and
  simultaneously an observable. **This is the criterion to enforce.**
- **Amplitude:** `|h| ≪ 1`, `|f| ≪ |F̄|` (Palessandro & Rothman arXiv:2301.02072;
  Planck-suppression argument in arXiv:2405.01407).

**Hole:** in the mixing/Gertsenshtein literature the criterion is asserted in one
sentence and never numerically enforced; no paper quantifies the error of the test-field
approximation there (error quantification exists only in inflationary gauge-field lattice
work, e.g. CosmoLattice arXiv:2607.24978). Enforcing per-run validity flags is therefore
both a safeguard and a cheap methodological contribution. Program adds a fourth monitor:
the **growth-impact ratio** (new-sector vs standard-sector source terms in the Einstein
constraints, per mode).

## 3. Consistency of coupling-without-backreaction

A clean double expansion (Cembranos arXiv:2302.08186, line 114, verbatim in spirit):
order 0 in perturbations → background equations (discarded; CAMB's ΛCDM solution used —
the only entry point of the spectator assumption); order 1 → tadpoles vanish on the
background solution; order 2 → the full quadratic action including all mixing terms, kept
whole. Coupling and backreaction sit at different orders. Berlin et al. arXiv:2405.08865
(local) is the axion-photon confirmation: expansion in `g_aγγ`, produced photons never
fed back.

**Mixing with the (non-negligible) photon sector:** the photon background is *unperturbed
at zeroth order* (CAMB's equilibrium blackbody); conversion computed as a perturbation on
it; the small parameter is the conversion probability `P` (`~10⁻²⁰` cosmologically,
arXiv:2312.17636, local); the negligible distortion *is* the observable (FIRAS/ARCADE
2/EDGES). Self-consistency closes on the `ΔN_eff` bound — the same inequality is both the
validity check and the competing constraint.

**Strategic corollary:** CMB constraints from broadband conversion channels are typically
*weaker* than `N_eff` unless the signal is spectrally narrow or structurally distinctive.
V-modes (arXiv:2502.12517 — recent, thin literature) and birefringence are exactly the
distinctive channels `N_eff` is blind to — the physics case for the polarization emphasis.

**Silent-failure risk:** the CAMB background solves Einstein's equations; a PGT quadratic
action requires the PGT background equations to hold on (FRW, `T̄=0`/tracking) —
tracking/frozen solutions in arXiv:2003.02690 (local, group's own). Check = background-EOM
residual on the CAMB background, per theory.

## 4. Scope boundary

Reachable in the strict limit: **(i) propagation** effects on standard quanta (tensor
transfer, birefringence, dispersion) and **(ii) conversion** effects moving a small
energy fraction (distortions, V-modes, radio excess). **Not reachable:** gravitational
sourcing of scalar CMB anisotropies (new-sector `δρ, σ` in the Einstein constraints) —
that is perturbation-level backreaction, the axionCAMB-style full-component route
(axionCAMB/AxiECAMB arXiv:2412.15192 are *full-backreaction*, not test-field — our route
is cheaper precisely by declining what they do). Also out: anything needing `ρ_new` in
the Friedmann equation (`ΔN_eff` as a fit parameter).

## 5. Is the niche occupied? — No (near-misses enumerated)

| Framework | What it does | Why it isn't this |
|---|---|---|
| xPand / xPert / xAct (arXiv:1302.6174, 0807.0824) | symbolic cosmological perturbation equations, any theory | no numerics, no observables, no likelihood |
| CppTransport / PyTransport (arXiv:1609.00380/1) | Lagrangian → auto-generated solver → correlators | inflation only; closest in spirit — cite as precedent |
| hi_class (arXiv:1909.01828) / EFTCAMB (arXiv:1405.3590) | covariant Lagrangian/EFT → Boltzmann → likelihood | Horndeski gravity sector, full backreaction, not general, not spectator |
| SymBoltz.jl (arXiv:2509.24740) | symbolic-component differentiable Einstein–Boltzmann | takes equations, not actions; not spectator-specialized; possible future backend |
| gammaALPs (github.com/me-manu/gammaALPs) | photon-ALP transfer matrices on prescribed environments | the architectural analogue (prescribed background + small mixing system), but astrophysical LOS, no Boltzmann/likelihood |
| CosmoDS (arXiv:2603.14740) | custom sector as Cobaya Theory, background only | background only; wiring precedent |
| TorC (arXiv:2507.09228, local) | modified background + ΛCDM perturbations | the exact complement: they modified the background and kept ΛCDM perturbations; we fix the background and add new perturbations. Neither has done both. |

**No torsion-as-spectator-on-FRW paper exists** (multiple search phrasings; corroborated
by arXiv:2506.17017 and TorC both listing cosmological perturbation theory as future
work). Defensible claim shape: the spectator approximation is standard but applied one
sector at a time, by hand, with validity asserted rather than enforced; the combination
"general Lagrangian → spectator perturbations → observables → Cobaya likelihood" plus
enforced validity is new, and torsion has never been treated this way.

## 6. Economics

The spectator route replaces the Boltzmann hierarchies (~10²–10³ multipoles/`k` with
approximation switching) by a small ODE system per `k` — roughly two orders of magnitude
smaller state, no stiffness machinery. Cobaya's fast/slow blocking is designed for
exactly this split (ΛCDM slow → CAMB; couplings fast → our ODE); **dragging** recommended
for large hierarchies with fast/slow degeneracies (e.g. `r` vs `ν`, arXiv:1405.7004).
Production precedent: dark-siren pipelines fix ΛCDM and vary only propagation parameters
(arXiv:2509.04348). No published wall-clock benchmark of this split exists — measure it
ourselves.

## Caveats

- No error-quantification literature for the test-field approximation in the mixing
  context (asserted, never bounded).
- "No torsion-as-spectator paper" is absence of search hits plus corroborating
  future-work statements, not proof.
- The ~110 local papers without extractable titles were not individually checked;
  `literature/README.md` documents only 15 of 125 directories (stale index).
