# TorC pipeline audit (H1)

**Executed:** 2026-08-30. **Status:** complete.
**Subject:** Legner, Handley & Barker, *"Alleviating the Hubble tension with Torsion
Condensation (TorC)"*, [arXiv:2507.09228](https://arxiv.org/abs/2507.09228) — local at
`literature/2507.09228/paper_Qtorsion.tex` — together with its two code forks and its
Zenodo data deposit.

## What this audit is for

The TIDAL → Cosmology program (`docs/COSMOLOGY_PROGRAM.md`, umbrella #488) evolves a new
sector's **perturbations** as spectators on a ΛCDM background that unmodified CAMB
supplies. TorC did the opposite: it modified the **background** (tabulated `ρ_Λ^eff(a)`,
`P_Λ^eff(a)` fed into a patched CAMB) and kept ΛCDM perturbations.

So TorC is **not a foundation and not a template for the physics**. It is the one existing
end-to-end run by this group through CAMB + Cobaya + PolyChord against real Planck and
SH0ES likelihoods, which makes it the natural **plumbing check** for our stack — the O1
rung of the observable ladder. The tabulated-background hook that O1 needs is a *strictly
optional* feature, off by default; nothing on the package's main path modifies the
background expansion, and the package carries no TorC dependence.

Each part of this audit exists for a specific downstream reason:

| Investigation | Why the program needs it |
| --- | --- |
| The paper's `ρ_Λ^eff`, `P_Λ^eff` and background system | Defines what a "tabulated background" input actually *is*, and shows why the interface must be `(ρ, P)` rather than `w` |
| `slegner/CAMB` diff | Identifies the code we must own to feed CAMB a `(ρ, P)` table at all, and how invasive it is — the input to the inherit-vs-reapply decision |
| `slegner/cobaya` diff | Gives the Cobaya-side contract for handing a table to CAMB. The provider→consumer wiring is the same shape our own `Theory` class will use, so it is a template even where the physics is not |
| Zenodo record | Supplies acceptance targets (evidences), the exact likelihood and prior configuration, and the reference implementation that generates the tables |
| Upstream CAMB check | Decides whether any patch is needed at all |

Two decisions come out of it, recorded in [§6](#6-decisions) and [§7](#7-decision-r2--how-the-camb-patch-gets-made): **O1a** (fixed-table
pass-through) for the scope of O1, and **re-apply against latest upstream** for the CAMB
patch.

---

## 1. The paper

### 1.1 Model and parameters

The TorC Lagrangian (`eq:TorCLagrangian`, line 171) is a purely quadratic Poincaré gauge
theory, `R² + T²`, with no Einstein–Hilbert term:

- `μ`, `ν` — dimensionless curvature couplings;
- `λ` — mass-dimension-two torsion coupling, giving a torsional dark energy;
- `Λ` — the usual cosmological constant.

Tree-level ghost and tachyon freedom requires (`eq:TorCconditions`, line 181)

```text
λ ≥ 0,    μ < 0,    (ν + 2μ)(ν − μ) > 0
```

and **the analysis sets `λ = 0`** (line 187), eliminating torsional dark energy and
reducing the dark sector to `Λ`. Only `λ` and `μ` matter for the background condensation.

For cosmology the paper uses a scalar-tensor equivalent (`TorCLagrangianMA`, line 198) in
two scalars `(ϖ, φ)`, which is exact for a homogeneous isotropic torsion ansatz
(`Tsamparlis`, line 206 — the 24 torsion components reduce to two time-dependent scalars).
`φ` is algebraically eliminable (line 268), leaving `ϖ` as the single dynamical torsion
scalar. `ϖ → 1` is an attractor, and `ϖ ≡ 1` reproduces ΛCDM exactly.

Two parameters beyond ΛCDM (`tab:TorCparams`, line 239):

| Parameter | Meaning |
| --- | --- |
| `Ω_Λ` | **Bare** dark energy density. Unlike ΛCDM it is *not* fixed by closure — the torsion terms break `ΣΩ = 1` (line 301), so it is an independent sampled parameter |
| `ϖ_r` | Value of the torsion scalar in the early Universe; the initial condition of the background ODE system |

### 1.2 Background system and time variable

The modified Friedmann equations are `eq:TorCF1` and `eq:TorCF2` (lines 283–286) and the
torsion equation of motion is `eq:varpiEOM` (line 292). The system is integrated in
**conformal time**, `dt = a dτ` (line 304), from power-series initial conditions
(`eq:aPowerSeries` line 684, `eq:varpiPowerSeries` line 690) started at

```text
τ_ini = 1e-6 · 16 √Ω_r ϖ_r / (Ω_m (3 ϖ_r² + 1) H_0)          (line 703)
```

Because `Ω_Λ` shifts the late-time expansion, `H` does not reach the sampled `H_0` at
`a = 1` without a rescaling `a → α a` (`eq:rescale`, line 309), with `α` determined
numerically.

> **Note for WS2.** TorC's `τ` is the same variable as the program's working time `η`, and
> the same as CAMB's internal `tau` — `cmbant/CAMB` `fortran/equations.f90:4` reads
> *"return d tau/ d a, where tau is the conformal time"*. The program's choice of conformal
> time therefore agrees with both CAMB and this predecessor, and no conversion layer is
> needed at that seam.

### 1.3 The effective dark energy — `ρ_Λ^eff` and `P_Λ^eff`

This is the input O1 consumes, so it is transcribed exactly. The modified first Friedmann
equation is written (`eq:F1eff`, line 319) as

```text
H² = H_0² (Ω_r a⁻⁴ + Ω_m a⁻³) + ρ_Λ^eff(t) / (3 M_P²)
```

which *defines* `ρ_Λ^eff` by comparison with `eq:TorCF1`.

**Effective density** — `eq:effrho`, lines 324–327, verbatim LaTeX:

```latex
\rho_{\Lambda}^{\mathrm{eff}} \equiv \frac{3\Mp}{\varpi^2} \Biggl[
    H_0^2 \Omega_\Lambda
  - H_0^2\left(\varpi^2 -1\right)\left(\Omega_\mathrm{r} a^{-4}+\Omega_\mathrm{m} a^{-3}\right)
  - 2 H \varpi \dot{\varpi}
  - \frac{\left(1+3 \varpi^2\right) \dot{\varpi}^2}{3\left(\varpi^2-1\right)}
\Biggr]
```

that is,

```text
ρ_Λ^eff = (3 M_P² / ϖ²) · [ H_0² Ω_Λ
                          − H_0² (ϖ² − 1)(Ω_r a⁻⁴ + Ω_m a⁻³)
                          − 2 H ϖ ϖ̇
                          − (1 + 3ϖ²) ϖ̇² / (3(ϖ² − 1)) ]
```

**Effective pressure** — `eq:effP`, lines 335–341, verbatim LaTeX:

```latex
P_{\Lambda}^\mathrm{eff} = \frac{\Mp}{12 a^4 H \varpi^3(\varpi^2-1)}\Biggl[
    3 H_0^2 H \varpi \left(\varpi^2-1\right)
    \times \Bigl(\Omega_\mathrm{r} \bigl(\varpi^2 -1\bigr) - 3 a^4 \Omega_\Lambda \bigl(1 + 3\varpi^2 \bigr) \Bigr)
  - 9 a^4 H^3 \varpi^3 \bigl(\varpi^2 -1\bigr)^2
  + 6 \bigl(\varpi^2 -1\bigr) \dot{\varpi} \biggl(4 H_0^2 \bigl(\Omega_\mathrm{r} + \Omega_\mathrm{m} a + \Omega_\Lambda a^4 \bigr)
      - a^4 H^2 \varpi^2 \bigl(5 + 3 \varpi^2 \bigr)\biggr)
  - 3 a^4 H \varpi \dot{\varpi}^2 \bigl(-13 + 3\varpi^2 \bigl(6 + \varpi^2\bigr)\bigr)
  - 8 a^4 \bigl(1 + 3\varpi^2\bigr) \dot{\varpi}^3
\Biggr]
```

that is,

```text
P_Λ^eff = M_P² / (12 a⁴ H ϖ³ (ϖ² − 1)) · [
      3 H_0² H ϖ (ϖ² − 1) · ( Ω_r (ϖ² − 1) − 3 a⁴ Ω_Λ (1 + 3ϖ²) )
    − 9 a⁴ H³ ϖ³ (ϖ² − 1)²
    + 6 (ϖ² − 1) ϖ̇ · ( 4 H_0² (Ω_r + Ω_m a + Ω_Λ a⁴) − a⁴ H² ϖ² (5 + 3ϖ²) )
    − 3 a⁴ H ϖ ϖ̇² · ( −13 + 3ϖ² (6 + ϖ²) )
    − 8 a⁴ (1 + 3ϖ²) ϖ̇³ ]
```

`P` is obtained from `ρ` through the continuity equation `ρ̇ = −3H(ρ + P)` (line 330), which
the paper justifies by the Bianchi identity of the GR-like part of the split (footnote at
line 328).

**Three properties that matter for the implementation:**

1. Both expressions are functions of `(ϖ, ϖ̇, H, a)` — **not closed-form in `a`**. Producing
   a table means first integrating `eq:TorCF2` + `eq:varpiEOM` as a coupled ODE system
   (line 344).
2. `ρ_Λ^eff` carries an overall `1/ϖ²`, so the standard density parameters cannot be
   separated from the torsion contribution — the whole dark sector moves together.
3. Both are singular at `ϖ = 1` (the `(ϖ² − 1)` denominators), which is exactly the ΛCDM
   limit. This is why the prior excludes `ϖ_r = 1` (line 419 footnote).

### 1.4 The pole — why forks exist at all

CAMB's stock dark-energy interface accepts an equation of state `w(a)`. The paper cannot
use it (line 360):

> attempting to incorporate the TorC extension through the equation of state
> `w_Λ^eff(a)` introduces poles when the dark energy density changes sign, as when `ϖ`
> exceeds unity in `eq:effrho` … This obstructs exploration of the full parameter space of
> `Ω_Λ` and `ϖ_r`. However, the issue arises only when using `w_Λ^eff`, and is resolved
> when pressure and density are implemented separately.

Restated at line 372 (the effective EOS "exhibits poles") and line 724. **This is the sole
reason both forks exist**, and it fixes the signature of our optional hook: the feature
takes `(a, ρ, P)` and must never be expressed as `w`. §8 quantifies where the poles
actually appear.

### 1.5 Priors, likelihoods, sampler

Uniform priors (`tab:prior`, line 395):

| Parameter | Range |
| --- | --- |
| `Ω_Λ` | `[0.1, 1.5]` |
| `ϖ_r` | `[0.1, 1.5]` |
| `h` | `[0.2, 1.0]` |
| `Ω_b h²` | `[0.005, 0.1]` |
| `Ω_c h²` | `[0.005, 0.99]` |
| `τ_reio` | `[0.01, 0.4]` |
| `n_s` | `[0.885, 1.04]` |
| `log A_s` | `[2.5, 3.7]` |

The ΛCDM ranges are Cobaya defaults. Two extra constraints (line 419 footnote): `ϖ_r`
restricted away from the `ϖ_r = 1` singularity, and `Ω_m = Ω_b + Ω_c < 1` imposed as a
prior "to avoid segmentation faults".

Likelihoods (line 421): Planck 2018 high-`ℓ` and low-`ℓ` TT, TE, EE, and lensing; plus
SH0ES 2020 `H_0 = 73.04 ± 1.04 km s⁻¹ Mpc⁻¹` as a Gaussian on `H_0`. The exact Cobaya
identifiers are not in the paper but are recoverable from the chains — see §5.3.

Sampler: PolyChord with **1000 live points** (line 469), driven through Cobaya. The paper
does not state `num_repeats`, the precision criterion, or the seed; the archive does not
contain a run configuration either (§5.4).

### 1.6 Results — O1's acceptance targets

Bayes factors, `Δ log Z = log Z_TorC − log Z_ΛCDM`, decomposed as
`Δ⟨log L⟩_P − Δ D_KL`:

| Dataset | `Δ log Z` | `Δ⟨log L⟩_P` | `Δ D_KL` | Reference |
| --- | --- | --- | --- | --- |
| Planck 2018 | `−5.627 ± 0.269` | `0.417 ± 0.114` | `6.044 ± 0.265` | `eq:LogZPlanck`, line 537 |
| Planck + SH0ES | `+0.725 ± 0.280` | `6.792 ± 0.113` | `6.071 ± 0.263` | `eq:LogZPlanckSH0ES`, line 544 |

Tension `R`-statistic between Planck and SH0ES (line 556): ΛCDM `log R = −5.60 ± 0.28`
(strong tension), TorC `log R = 0.70 ± 0.27` (consistent). Computed with `anesthetic`.

Qualitative posterior results: `H_0` anticorrelates with `ϖ_r` and correlates with `Ω_Λ`
(line 514); TorC broadens the CMB `H_0` posterior enough that SH0ES ceases to be an outlier
(line 512); `S_8` shifts low but is flagged as unreliable without perturbations (line 519).

**There is no published table of marginal parameter values** — only contour figures. The
usable numeric targets are the four evidences (recovered exactly in §5.3) and the two
`log R` values.

### 1.7 Two statements that bear directly on the program

**The perturbation caveat** (line 383, repeated at 211 and 260) — the limitation the whole
program exists to lift:

> Note that this analysis modifies the background expansion in `CAMB` using the TorC
> effective dark energy density and pressure, while the perturbation equations remain those
> of standard ΛCDM. […] As a initial study, only the background evolution of TorC is being
> investigated, cosmological perturbation theory for TorC will be developed in future work.

**`ΔN_eff`** (line 716) — relevant to the program's spectator-validity flags:

```text
ΔN_{ϖ_r,eff} = (ϖ_r⁻² − 1) · ( (8/7)(11/4)^{4/3} + N_eff )
```

Torsion behaves as early dark radiation. The main runs **exclude** this from the
cosmological evolution to avoid double-counting the same effect already carried by the
modified dark energy (line 720); a separate validation run routed it through PArthENoPE's
`Y_P` only, and did not change conclusions.

---

## 2. `slegner/CAMB` — what was patched and how

Fork of `cmbant/CAMB`, default branch `master`, last pushed 2026-06-10.

### 2.1 Isolating the paper-era patch

The fork has continued past the paper. The publication state is commit **`2fb908af`**
("merged to master", 2025-07-10) — one day before the Zenodo deposit of 2025-07-11. Use it
as the reference point:

```bash
gh api "repos/cmbant/CAMB/compare/master...slegner:CAMB:2fb908af"   # paper-era: 12 commits, 21 files
gh api "repos/cmbant/CAMB/compare/master...slegner:CAMB:master"     # current:   18 commits, 24 files
```

| | Paper-era (`2fb908af`) | Current `master` |
| --- | --- | --- |
| Commits ahead | 12 | 18 |
| Files changed | 21 | 24 |
| Behind upstream | 180 | 98 |

**Post-paper only — exclude from any baseline** (commits `08d9675a` 2026-04-20 "adding
massive neutrinos and hmcode using tabulated rho and P", `a9458be7` 2026-06-09):
`camb/massive_neutrinos.py` (+405, new), `fortran/massive_neutrinos.f90` (+19/−1), the
enlarged `fortran/halofit.f90` (+67/−3 vs +1 at paper-era), the `CAMB_GetMassiveNuRhoP`
C-callable in `fortran/camb_python.f90` (+36 vs +1), `.gitignore`, and a stray
`.claude/settings.local.json`. Conversely `camb/camb.py` (a one-line `read_ini` fix)
appears at paper-era but not in the current diff — upstream absorbed it.

### 2.2 The mechanism, end to end

The design adds a dark-energy model that reports **density and pressure directly**, instead
of reporting `w` and having CAMB reconstruct `P = w·ρ`.

**New Fortran module `fortran/DarkEnergyPressure.f90` (+169).**
`TDarkEnergyDensityAndPressure` extends `TDarkEnergyModel` and holds two `TCubicSpline`s
(`density`, `pressure`) sampled in `log a`. `SetrhoTable` / `SetPTable` both `error stop`
unless the table ends at `a = 1`. Inside the table, `grho_de(a)` returns `ρ·a⁴` and
`P_de(a)` returns the spline value; outside it, both extrapolate by freezing
`w = P/ρ` at the nearest endpoint and continuing `grho_de ∝ a^{1−3w}` from there
(equivalently `ρ ∝ a^{−3(1+w)}`, the usual constant-`w` scaling). `Init` sets
`is_cosmological_constant = .false.` and `is_no_mod_P = .false.`.
`Effective_w_wa` synthesizes a `(w, wa)` pair at `a = 1` purely so that HMcode/halofit —
which insist on a `w`,`wa` parameterization — have something to consume.

**New Fortran module `fortran/DarkEnergyPressurePPF.f90` (+140).**
`TDarkEnergyPressurePPF` extends the above with the PPF perturbation closure
(`c_Gamma_ppf = 0.4`, `PerturbedStressEnergy`, `diff_rhopi_Add_Term`), supports `cs2 = 1`
only, and registers Python name `DarkEnergyPressurePPF`. **This is the class actually
used.**

**`fortran/DarkEnergyInterface.f90` (+36/−8) — the load-bearing change.** Three parts:

1. Two new flags on `TDarkEnergyModel`: `is_no_mod_w`, `is_no_mod_P`, both defaulting
   `.true.`, plus a new `P_de` type-bound procedure (default `−1`).
2. `BackgroundDensityAndPressure` gains an **optional `P` output**:
   ```fortran
   subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w, P)
   ```
   with the branch that matters — if `is_no_mod_w` is false, `w = w_de(a)` and `P` is
   derived from it (the old behavior); if `is_no_mod_P` is false, **`P = P_de(a)` is taken
   directly from the spline and `w` is derived from `P`**. The pole is sidestepped here:
   `w` may be arbitrarily large or undefined, and nothing downstream needs it.
3. `PerturbedStressEnergy` and `diff_rhopi_Add_Term` gain a `Pgrhova2` argument.

**`fortran/equations.f90` (+28/−7) — consumption in the solver.** In `derivs`, a new
`Pgrhova2 = P_dark_energy_t · State%grhov · a²` is computed alongside `grhov_t`, and the
two total-pressure sites change from

```fortran
gpres = gpres_noDE + w_dark_energy_t*grhov_t
```
to
```fortran
gpres = gpres_noDE + Pgrhova2
```

`Pgrhova2` is then threaded into `PerturbedStressEnergy` and `diff_rhopi_Add_Term`. Also
adds an early-abort guard when `tauend` exceeds every switch time — with a bare
`print *, "tauend larger than all switch error encountered"` left in.

**Signature churn, not new physics.** Because the two abstract signatures changed,
`DarkEnergyFluid.f90` (+6/−5), `DarkEnergyPPF.f90` (+4/−4) and
`DarkEnergyQuintessence.f90` (+4/−4) all change — they only add the unused `Pgrhova2`
argument. (`DarkEnergyFluid.f90` additionally sets `is_no_mod_w` for
`TAxionEffectiveFluid`.)

**Python surface.**

| File | Change |
| --- | --- |
| `camb/dark_energy.py` (+96/−1) | Adds `DarkEnergyPressure` and `@fortran_class DarkEnergyPressurePPF`, with `set_rho_a_table(a, rho)` / `set_P_a_table(a, P)` validating equal lengths, `a[-1] ≈ 1`, and `a > 0`. Registers short name `"pressureppf"`. Adds the two new flags to `DarkEnergyModel._fields_` |
| `camb/model.py` (+17/−1) | Adds `CAMBparams.set_dark_energy_rho_p_a(a_rho, rho, a_p, p, dark_energy_model="pressureppf")` — the public entry point |
| `camb/results.py` (+5/−4) | Widens `DarkEnergyStressEnergy` to return `P` as well, turning `get_dark_energy_rho_w` from a 2-tuple into a 3-tuple. **A breaking change to a public upstream API** |

**Build plumbing** (mechanical): `fortran/Makefile_main` (adds `DarkEnergyPressure` to
`SOURCEFILES` and `DarkEnergyPressurePPF` to `DARKENERGY_FILES` — module order matters),
`fortran/Makefile`, `fortran/camb.cbp`, `setup.py`, `fortran/camb.f90` (accepts
`DarkEnergyModel = 'PRESSUREPPF'` from an ini file), and one-line `use` additions in
`fortran/camb_python.f90`, `fortran/cmbmain.f90`, `fortran/model.f90`,
`fortran/halofit.f90`.

**All 21 paper-era files accounted for:** 2 new modules; 1 interface change; 1 solver
change; 3 signature-churn DE modules; 3 Python files; 10 build/plumbing/one-liners
(`Makefile`, `Makefile_main`, `camb.cbp`, `setup.py`, `camb.f90`, `camb_python.f90`,
`cmbmain.f90`, `model.f90`, `halofit.f90`, `camb/camb.py`) plus `fortran/results.f90`
(+3/−3, the `DarkEnergyStressEnergy` Fortran side).

### 2.3 Reading

There is also an ini-file route: `TDarkEnergyDensityAndPressure_ReadParams` loads two
two-column text files named by `rhoafile` and `Pafile`. The Cobaya path does not use it,
but it is the simplest way to drive the feature from a standalone CAMB run — useful for
building a reference binary in O1.

---

## 3. `slegner/cobaya` — three-way attribution

Fork chain: `CobayaSampler/cobaya` (source) → `AdamOrmondroyd/cobaya` (parent) →
`slegner/cobaya`. Last pushed 2025-07-10. Against upstream it is `ahead_by 39`,
`behind_by 95`, 8 files — but **only 2 of those files are TorC's own work**.

```bash
gh api "repos/AdamOrmondroyd/cobaya/compare/master...slegner:cobaya:master"   # TorC delta: 7 commits, 2 files
gh api "repos/CobayaSampler/cobaya/compare/master...slegner:cobaya:master"    # everything: 39 commits, 8 files
```

### 3.1 TorC's own delta — 2 files

**`cobaya/theories/camb/camb.py` (+72/−2)** adds two boolean options, `external_wa` and
`external_rhopa` (declared in `camb.yaml`, +2, both defaulting `False`). When
`external_rhopa` is set:

- `must_provide` returns `{"dark_energy": None}`, so Cobaya requires **some other
  component** to provide a `dark_energy` product;
- in `set()`, before any cosmology is configured,
  ```python
  de = self.provider.get_dark_energy()
  a_rho, a_p, rho, p = de["a_rho"], de["a_p"], de["rho"], de["p"]
  ...
  params.set_dark_energy_rho_p_a(**darkenergypressure(a_rho, rho, a_p, p, **self.extra_args))
  ```
  The ordering is deliberate and commented — *"DarkEnergy has to be set before cosmology if
  theta is used instead of H0"* — because `cosmomc_theta` solving depends on the background;
- a post-condition asserts `type(params.DarkEnergy).__name__` ends in `"PressurePPF"`,
  raising `LoggedError` otherwise;
- two module-level helpers, `darkenergy(a, w, dark_energy_model, **kwargs)` and
  `darkenergypressure(a_rho, rho, a_p, p, dark_energy_model, **kwargs)`, assemble the
  kwargs and pop `dark_energy_model` out of `extra_args`.

Note what this fork does **not** contain: the component that computes `ρ(a), P(a)` from
`(Ω_Λ, ϖ_r)`. The fork is pure plumbing; the physics lives in the (unpublished) provider
component, whose reference implementation is the Zenodo `TorC_rhopa.py` (§5.2).

> **Template value for WS5.** This provider→consumer pattern — a `Theory` declaring
> `{"product": None}` in `must_provide`, pulling it via `self.provider.get_*()` in
> `calculate`/`set`, and asserting the downstream state — is exactly the shape our own
> `Theory` class needs when it chains off `CAMBdata`. Worth copying as wiring even though
> the physics direction is opposite.

### 3.2 Inherited from `AdamOrmondroyd`, not TorC's

**`cobaya/samplers/polychord/polychord.py` (+14/−27)** is the substantive inherited change.
Stock Cobaya passes PolyChord the *posterior* on a linearly rescaled unit hypercube, adding
back `logvolume = log Π(scales)` to compensate. The fork deletes that and instead:

- transforms the cube through the true inverse CDF,
  `theta[i] = self.model.prior.pdf[i].ppf(x)`;
- passes the *likelihood* (`max(loglikes.sum(), logzero)`), not the posterior;
- dumps `.paramnames` up front.

This makes `log Z` a genuine Bayesian evidence under non-uniform priors — a prerequisite
for the `anesthetic` `R`-statistic the paper reports, and hence load-bearing for the
results even though it is not TorC-specific. Also `polychord.yaml` `nlives: {}` → `nlives:`
(+1/−1) and setuptools-scm versioning (`cobaya/__init__.py`, `pyproject.toml`). The
eighth changed file, `.gitignore` (+1), is trivial.

### 3.3 A robustness patch worth stealing

**`cobaya/likelihoods/base_classes/planck_clik.py` (+5)** returns `−inf` when any `C_ℓ`
array contains NaN, instead of letting `clik` crash:

```python
if np.any(np.isnan(cl_array)):
    self.log.error("nans in cl['%s']: returning logzero and carrying on." % cl_key)
    return -np.inf
```

That is a direct symptom of exotic backgrounds reaching the likelihood in a nested-sampling
prior box. Our own component should fail the same way — as a flagged rejection, not a
crash.

---

## 4. Upstream CAMB today

Latest release **2.0.3**, published 2026-08-24. Checked directly:

- `camb/dark_energy.py` exposes only `DarkEnergyEqnOfState.set_w_a_table(a, w)` — a `w(a)`
  table, exactly the pole-prone interface the paper had to abandon. Model classes are
  `DarkEnergyFluid`, `DarkEnergyPPF`, `AxionEffectiveFluid`, `Quintessence`,
  `EarlyQuintessence`.
- `fortran/DarkEnergyInterface.f90` still declares
  `BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)` — **no `P` argument** — and has
  no `P_de` procedure.
- A code search for `set_dark_energy_rho` across `cmbant/CAMB` returns **0 hits**.

**Upstream has not gained `(ρ, P)` table support.**

---

## 5. The Zenodo record

[doi:10.5281/zenodo.15866507](https://doi.org/10.5281/zenodo.15866507) — "Chains and
supplementary files for *Alleviating the Hubble tension with Torsion Condensation (TorC)*",
Legner, Handley & Barker, deposited 2025-07-11.

One file: `run_chains.zip`, **623,630,974 bytes**, md5
`ba9a10f698465d5535caaabf60b05586`.

### 5.1 Reading it without downloading it

Zenodo's file endpoint honors HTTP Range requests (verified: `Range: bytes=-100` → `206`),
even though it does not advertise `Accept-Ranges`. The archive is a classic (non-Zip64)
zip, so its index is cheap to read:

- **116 entries**; central directory at offset **623,614,488**, size **16,464** bytes.

```bash
URL="https://zenodo.org/records/15866507/files/run_chains.zip?download=1"
curl -sL -H "Range: bytes=623614488-623630951" "$URL" > cd.bin   # 16 KB: the whole index
```

To extract one member, range-fetch from its local-header offset for
`compressed_size + len(name) + 200` bytes, skip the 30-byte header plus name and extra
fields, and `zlib.decompress(data, -15)`. Every file quoted in this report was obtained
this way; total transferred, well under 100 KB.

### 5.2 Supplementary material

| Path | Size | What it is |
| --- | --- | --- |
| `run_chains/Readme.md` | 3.0 KB | Contact `sl2091@cam.ac.uk`. Names dependencies: modified CAMB, and **the `Tension` branch of `handley-lab/anesthetic`** (not the release) |
| `run_chains/TorC_equation_derivations.nb` | 630 KB | xAct/xTras derivation of the field equations and of `ρ_Λ^eff`, `P_Λ^eff` |
| `run_chains/Plotting/TorC_rhopa.py` | 6.4 KB | **The table generator** — see below |
| `run_chains/Plotting/Plot_chains.py` | 15.8 KB | Posteriors, model comparison, tension analysis |
| `run_chains/Plotting/Plot_cpr_OmL_effects.py` | 13.1 KB | `C_ℓ` and comoving-horizon effect figures |

**`TorC_rhopa.py` is the reference implementation O1 needs.** It is pure numpy/scipy — no
CAMB, no xAct, no TorC framework. Its structure:

1. `omega_r(h)` computes `Ω_r` from CMB temperature 2.7255 K with `N_eff = 3.046`.
2. `aseries`/`cpseries` (and their derivatives) are the paper's power-series ICs, evaluated
   at `tauini = small(cpr, omega_m, h)` — the paper's `τ_ini`.
3. `coupdiff` is the coupled `(a, ϖ, a′, ϖ′)` conformal-time system, integrated with
   `solve_ivp` at `rtol=1e-10`.
4. **Two-pass normalization**: integrate with `a0 = 1` until the event `y[2] − 100h = 0`
   fires (i.e. `H = H_0`), read off `a0` there, then re-integrate with that `a0` until
   `a = 1`. This is the paper's `α` rescaling, done numerically.
5. `rho_deriveda0` / `P_deriveda0` evaluate `ρ_Λ^eff`, `P_Λ^eff` in the code's conformal
   variables, on a **10,000-point `logspace` grid in `a`**.
6. **The tables handed to CAMB are normalized:**
   ```python
   rho_renorm = rho_derived / rho_derived[-1]
   P_renorm   = P_deriveda0(...) / rho_derived[-1]
   ```
   i.e. both are divided by `ρ(a=1)`, so the table satisfies `ρ(1) = 1`. **Our optional
   hook must adopt the same convention** or the CAMB normalization will be wrong.

### 5.3 Chains — seven runs

| Run | Model | Data |
| --- | --- | --- |
| `TorC_Planck_Lense` | TorC | Planck 2018 (+ lensing) |
| `TorC_SH0ES_Planck_Lense` | TorC | Planck 2018 + lensing + SH0ES 2020 |
| `TorC_SH0ES` | TorC | SH0ES 2020 only |
| `LCDM_Planck_Lense` | ΛCDM | Planck 2018 (+ lensing) |
| `LCDM_SH0ES_Planck_Lense` | ΛCDM | Planck 2018 + lensing + SH0ES 2020 |
| `LCDM_SH0ES` | ΛCDM | SH0ES 2020 only |
| `LCDM_Freedman` | ΛCDM | Freedman 2020 TRGB |

(`Lense` in the directory name marks inclusion of the Planck 2018 lensing likelihood.)

Each run directory holds `<name>.txt`, `_dead.txt`, `_dead-birth.txt`,
`_equal_weights.txt`, `_phys_live.txt`, `_phys_live-birth.txt`, `_prior.txt`, `.stats`,
`.paramnames`, `.prior_info`, plus a `clusters/` subdirectory duplicating the single
cluster found in each run.

**Evidences, read directly from the `.stats` files:**

| Run | `log Z` | `n_posterior` |
| --- | --- | --- |
| `TorC_Planck_Lense` | `−1447.31778 ± 0.22115` | 72,207 |
| `TorC_SH0ES_Planck_Lense` | `−1449.30915 ± 0.22053` | 72,540 |
| `LCDM_Planck_Lense` | `−1441.71360 ± 0.20726` | 66,123 |
| `LCDM_SH0ES_Planck_Lense` | `−1450.08732 ± 0.20671` | 66,130 |

Differencing gives `Δ log Z^Planck = −5.604` and `Δ log Z^joint = +0.778`, against the
paper's `−5.627 ± 0.269` and `+0.725 ± 0.280` — consistent, the small offsets being the
paper's `anesthetic` resampling over 1000 simulated data points (figure caption, line 530).
**These four numbers are the acceptance targets** if a full reproduction is ever attempted.

**The `.paramnames` files give the exact run configuration**, which the paper does not.
From `TorC_Planck_Lense.paramnames`:

- sampled: `OmL`, `cpr` (`ϖ_r`), `logA`, `ns`, `h`, `ombh2`, `omch2`, `tau`, plus **21
  Planck nuisance parameters** (`A_planck`, `calib_100T`, `calib_217T`, `A_cib_217`,
  `xi_sz_cib`, `A_sz`, `ksz_norm`, four `gal545_A_*`, four `ps_A_*`, six `galf_TE_A_*`);
- likelihoods: `planck_2018_lowl.TT`, `planck_2018_lowl.EE`,
  `planck_2018_highl_plik.TTTEEE`, `planck_2018_highl_plik.SZ`,
  `planck_2018_lensing.clik`;
- priors: `0` (the base prior), `omegam_prior`, `cpr_prior`, `SZ`;
- derived include `H0`, `omegam`, `omega_de`, `sigma8`, `S_8` variants, `rdrag`, `age`,
  `YHe`, `DHBBN`.

`TorC_Planck_Lense.prior_info` records `nprior = 10000`, `ndiscarded = 159228`.

### 5.4 Two gaps

1. **No Cobaya run configuration is archived** — no `.yaml`, no `.updated.yaml`, no `.ini`.
   The configuration must be reconstructed from `.paramnames` plus the paper's prior table.
   In particular `num_repeats`, the precision criterion and the seed are unrecoverable.
2. **No published marginal-parameter table**, in paper or archive — posteriors exist only as
   samples and contour figures.

### 5.5 What a reproduction would actually need

Per run: `.paramnames` (~2 KB), `.stats` (~4 KB), `_equal_weights.txt` (~10 MB for the
Planck runs). For the four headline runs that is roughly **40 MB**, versus 623 MB for the
archive.

Skippable — about 600 MB of the total: `_dead.txt` and `_dead-birth.txt` (~100–137 MB
each), `_phys_live*`, `_prior.txt`, the full `.txt` weighted chains (needed only for
`anesthetic` evidence *recomputation*, not for posteriors), and all of `clusters/`, which
duplicates the single cluster of each run.

The `.stats` values in §5.3 are already extracted, so **for the O1a route chosen below
nothing further needs downloading at all.**

---

## 6. Decisions

### R1 — Scope of O1: fixed-table pass-through

O1 is a **plumbing gate**, not a TorC reproduction. Two rungs were considered.

**O1a — fixed-table pass-through. CHOSEN.** Generate *one* `(a, ρ, P)` table offline by
running `TorC_rhopa.py` at a chosen `(Ω_Λ, ϖ_r)`, feed it through our optional
tabulated-background hook, and check we drive CAMB to the same `H(a)` and `C_ℓ` as a
reference build. This exercises the hook, the Cobaya wiring and the CAMB seam — everything
O1 exists to test — with no TorC physics inside our package and no ODE in any hot path.
Table generation is a one-off data-preparation step, not a package feature.

**O1b — full posterior reproduction. REJECTED for now.** Regenerating the table per
likelihood call and sampling all 8 + 21 parameters with PolyChord against the Planck
likelihood is the only route that reproduces the published numbers. It costs a full
nested-sampling campaign for a plumbing check, and it puts the nonlinear TorC background —
a regime this program is explicitly not investigating — inside the fast block. The four
evidences in §5.3 stay on record should this ever be wanted.

> **For the orchestrator.** This narrows `docs/COSMOLOGY_PROGRAM.md`'s stated O1 validation
> ("against published chains, Zenodo 10.5281/zenodo.15866507") to a `C_ℓ`- and
> background-level check. The program document is not edited here; see §9.

---

## 7. Decision R2 — how the CAMB patch gets made

Three options, stated so the rejected ones are on record.

**(a) Inherit / rebase `slegner/CAMB`** — depend on their fork, rebased onto current
upstream. **Rejected.** The paper-era state is 180 commits behind upstream (current
`master`, 98). It carries unrelated post-paper work (massive neutrinos, HMcode), a breaking
change to the public `get_dark_energy_rho_w` signature, a leftover debug `print *`, and a
stray `.claude/settings.local.json`. Pinning a stale Fortran fork also contradicts the
program's latest-upstream-CAMB policy, for a feature that is meant to be optional and off by
default.

**(c) Use upstream CAMB's own dark-energy interface** — no patch at all. **Rejected** by
§4: only `set_w_a_table` exists, and that is precisely the pole-prone route the paper had
to abandon. §8 shows the poles are real and unbounded.

**(b) Re-apply the design cleanly against current upstream. CHOSEN.** Study the fork
carefully to understand the change, then implement our own minimal `(ρ, P)`-table
dark-energy model against 2.0.3 — keeping the design, not the fork. The genuinely new code
is two self-contained Fortran modules (~309 lines) plus ~115 lines of Python surface; the
rest of the 21-file footprint is signature churn from one optional argument, plus build
plumbing. Re-applying lets us drop the breaking API change, drop the post-paper additions,
and stay rebaseable against upstream releases.

### 7.1 `2fb908af` is a specification, not a patch

Isolating the paper-era change is just a matter of choosing the revision — the compare
command in §2.1 does it. What choosing the revision does **not** give us is something
`git apply` can consume: the snapshot is 180 commits behind 2.0.3 and the change touches
shared code (`DarkEnergyInterface.f90`, `equations.f90`, three DE-subclass signatures) that
upstream has moved under. So `2fb908af` defines *what to build*; the re-apply is real work.

Four steps, each independently checkable:

1. Add the optional `P` output to `BackgroundDensityAndPressure`, the `P_de` procedure and
   the two `is_no_mod_*` flags, propagating the signature to the existing DE subclasses.
   Mechanical; mirrors the fork's churn.
2. Port `DarkEnergyPressure.f90` and `DarkEnergyPressurePPF.f90` largely as-is — they are
   self-contained and barely touch upstream surface.
3. Switch the two `gpres` sites in `equations.f90` to `Pgrhova2` and thread it through
   `PerturbedStressEnergy` / `diff_rhopi_Add_Term`.
4. Python surface and build plumbing — but add a **new** accessor rather than widening
   `get_dark_energy_rho_w`, so no public upstream API breaks.

Verification: with a constant table `P = −ρ` the patched build must reproduce stock ΛCDM
`C_ℓ` to machine precision; and a `TorC_rhopa.py` table must reproduce a reference build's
`H(a)`. §8 adds a third, independent check.

### 7.2 Distribution

Fork `cmbant/CAMB` under our own account and branch **from the `2.0.3` release tag**, not
from `master`, so the base is a released, reproducible state. Land the change as a small
readable commit series whose messages and a `README` note state plainly that it re-applies
the design of `slegner/CAMB@2fb908af` (Legner, Handley & Barker, arXiv:2507.09228) to
current upstream, with attribution.

Our package depends on that branch **only** when the optional tabulated-background feature
is enabled; the default install takes stock CAMB from PyPI with no Fortran build.

Rejected alternative: carrying `.patch` files in our own repo and applying them to a pinned
upstream tag at build time. It is reviewable in-tree, but it makes installation fragile and
gives the change no stable identity to cite.

The residual cost is that users of the optional path need a Fortran toolchain. That is
acceptable precisely because the path is optional — and it is the strongest argument for
eventually offering the optional `P` output upstream, which would retire the fork entirely.
Flag that as a follow-up, not a dependency.

### 7.3 Feature shape

The hook takes `(a, ρ, P)` with `ρ` and `P` normalized so that `ρ(a=1) = 1` (matching
`TorC_rhopa.py`, §5.2), requires the grid to end at `a = 1`, is declared per run, and is
inert unless explicitly enabled. It is never expressed as `w`.

---

## 8. Where the poles actually are — a free cross-check

The paper asserts that `w_Λ^eff` develops poles; it does not say where. Running Zenodo's
`TorC_rhopa.py` unmodified, at Planck-like `h = 0.674`, `Ω_m = 0.314`, over the 10,000-point
grid it produces:

| `ϖ_r` | `Ω_Λ` | sign changes in `ρ_Λ^eff` | `max \|w\|` |
| --- | --- | --- | --- |
| 0.8 | 0.685 | **0** | 0.998 |
| 0.95 | 0.685 | **0** | 0.9995 |
| 1.05 | 0.685 | 1 | 2.7e4 |
| 1.2 | 0.3 | 1 | 4.2e2 |
| 0.5 | 1.0 | 2 | 2.0e3 |

`(ϖ_r, Ω_Λ) = (0.8, 0.685)` is the paper's own fiducial point (figure caption, line 369).
The poles are real and unbounded, confirming §1.4 — and they appear only once `ϖ_r > 1` or
`Ω_Λ` is pushed to the prior edge, exactly as the paper describes.

**This does not argue for skipping the fork.** The `(ρ, P)` interface is what makes the
feature general; without it the prior box is unreachable. What the result buys is a
**second, independent route to the same physics at the fiducial point** — stock upstream
CAMB fed the same table through `set_w_a_table` — and therefore the cleanest possible oracle
for our re-applied patch, at zero cost:

- O1a runs on the patched `(ρ, P)` path, as the feature is meant to be used;
- the stock-CAMB `w(a)` run at `ϖ_r = 0.8` is a **cross-check**: the two must agree on
  `H(a)` and `C_ℓ` to solver tolerance. Disagreement localizes a bug in our re-apply rather
  than in the physics;
- moving to `ϖ_r > 1`, where the `w(a)` route provably breaks, then demonstrates the
  feature's whole reason for existing, and is the natural regression test.

Caveat: the probe used the paper's fiducial point. Any other point intended as a `w(a)`
cross-check must be re-verified for sign changes and bounded `|w|` first.

---

## 9. Open items

**For the orchestrator (`docs/COSMOLOGY_PROGRAM.md`, not edited here):**

- O1's validation criterion currently reads "against published chains (Zenodo …)". Under R1
  it becomes a `C_ℓ`/background-level check against a reference CAMB build at a fixed
  `(Ω_Λ, ϖ_r)`. The four evidences (§5.3) remain the targets for O1b if it is ever revived.
- The CAMB policy line ("`slegner/CAMB` fork is an *optional hook*; H1 decides
  inherit-vs-reapply") is settled: re-apply, own fork off the `2.0.3` tag (§7).

**For WS5 (#494):**

- The provider→consumer pattern in §3.1 is the wiring template for our `Theory` class.
- Adopt the `planck_clik` NaN guard (§3.3) as a flagged rejection path.
- Cobaya's stock PolyChord sampler does *not* give correct evidences under non-uniform
  priors without the Ormondroyd patch (§3.2). If we ever need evidences rather than
  posteriors, check whether that fix has reached upstream Cobaya.

**For H2 (observable ladder):**

- O1 as scoped here tests plumbing only and produces no new physics; H2 should weigh it
  purely as a gate.
- `ΔN_eff` (§1.7) is already parameterized for the torsion sector and is a ready-made input
  to the program's spectator-validity flags.

**Unresolved:**

- PolyChord `num_repeats`, precision criterion and seed are not recorded anywhere (§5.4).
- No marginal-parameter *table* was published (§1.6, §5.4) — but the posteriors themselves
  are recoverable: `_equal_weights.txt` is archived for all seven runs and `.paramnames`
  carries the derived parameters, so marginals are obtainable via `anesthetic`. What the
  missing run configuration rules out is *exact* reproduction, not *statistical*
  comparison.

## Sources

- `literature/2507.09228/paper_Qtorsion.tex` (line references throughout)
- `slegner/CAMB` @ `2fb908af` and `master`; `cmbant/CAMB` @ 2.0.3
- `slegner/cobaya`, `AdamOrmondroyd/cobaya`, `CobayaSampler/cobaya`
- Zenodo [10.5281/zenodo.15866507](https://doi.org/10.5281/zenodo.15866507),
  `run_chains.zip` md5 `ba9a10f698465d5535caaabf60b05586`
