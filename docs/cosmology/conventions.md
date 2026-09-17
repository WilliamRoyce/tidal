# Conventions — signature, ε, CAMB variables, PSALTer input, legacy comparisons

> **Status: CANONICAL since 2026-09-15** — written by R-1 (#566), adopted by the user, recorded as
> the register row *Conventions* in `docs/COSMOLOGY_PROGRAM.md`. Every entry is pinned to a source
> line that was re-read for this document. Where this file and `spectrum_design.md` §4.3
> disagree about CAMB, this file is the corrected reading (§2). The reasoning trail is in
> `r1_planning_record.md`; the decision memo is `interfaces_decision.md`.

<!-- cspell:words etak dgrho dgpi gpres Scal adotoa clxcdot hdot epsilonG uₐuᵃ Challinor Lasenby dxⁱdxʲ ayprime hexic ipynb DEVEL -->

**Pinned versions.** PSALTer v2.0.2 at `bb45adb0` (installed under Wolfram's
`$UserBaseDirectory/Applications/xAct/PSALTer`; paths below are relative to it). CAMB 2.0.4:
the installed `camb` Python package, which is byte-identical to tag `2.0.4` =
`cmbant/CAMB@a6de8cc59124c8bbe3924f1c546f96cef73dcabb`; Fortran and docs paths below are at that
commit. Legacy TIDAL at this repository's `feat/cosmology-program`.

## 1. The rule: each quantity takes the convention of the tool that defines it

| quantity | owner | why |
| --- | --- | --- |
| metric signature | PSALTer — and CAMB agrees (§2) | the ghost verdict's parity factor is defined by it (`spectrum_design.md:320-327`) |
| ε orientation | PSALTer — CAMB defines none (§2) | parity-odd operators in our Lagrangians flip with it |
| perturbation variables, gauge/frame names, conformal time, units | CAMB (§3) | the solver hands values to and from CAMB |
| sampled parameter names | Cobaya (`params:`, `renames`) | users already know them |
| Lagrangian input form and validity | PSALTer's requirements, enforced by us (§4) | PSALTer does not check most of them itself |

Nothing is converted at run time (#513). A formula taken from a source with a different
convention is rewritten **once, where it is imported**, with the source cited and a test.
Every derived artifact records `conventions` and every consumer asserts it (§5).

## 2. Signature and ε — no competing standard among the tools we call

| source | signature | ε | evidence |
| --- | --- | --- | --- |
| **PSALTer (code)** | (+,−,−,−) | ε₀₁₂₃ = +1 | `PSALTer.m:101` — "the West Coast signature (1,-1,-1,-1)"; `Sources/ReloadPackage/DefGeometry.m:69-72` sets `G` to `DiagonalMatrix@{1,-1,-1,-1}` and `epsilonG` lower components to `$epsilonSign*LeviCivitaTensor`, with `$epsilonSign` defaulting to 1 (`xTensor.m:1834`) |
| PSALTer papers | (+,−,−,−) | ε₀₁₂₃ = +1 | arXiv 2406.09500 `Manuscript.tex:102`; arXiv 2506.02111 `Manuscript.tex:100` ("'mostly minus' signature … ε₀₁₂₃ = 1") |
| **CAMB (declared)** | (+,−,−,−) | none | `camb/symbolic.py:460` (Newtonian: ds² = a²((1+2Ψ_N)dt² − (1−2Φ_N)δᵢⱼdxⁱdxʲ)); `docs/source/variables_guide.rst:223` (synchronous: ds² = a²(τ)[dτ² − (δᵢⱼ+hᵢⱼ)dxⁱdxʲ]); `docs/ScalEqs.ipynb` cell at lines 278-285; the CAMB notes (cosmologist.info/notes/CAMB.pdf) "Using the uₐuᵃ = 1 signature" |
| **CAMB (Fortran)** | signature-agnostic | none | `fortran/equations.f90:2365, 2427, 2430-2431, 2442` are Ma–Bertschinger's synchronous-gauge equations in variables that do not change under g → −g (`:38-39` points to the notes). A grep of every `.f90`/`.py` at the tag finds no Levi-Civita convention; the notes assume a parity-symmetric ensemble (C_ℓ^TB = C_ℓ^EB = 0) |
| Cobaya | none | none | a sampler |

**Corrections to `spectrum_design.md` §4.3 (reported, not edited here):** `:301` labels CAMB
"(−,+,+,+), Ma–Bertschinger" but cites Ma & Bertschinger, not CAMB; `:354` repeats it for the
solver branch; `:312-313` says ε·ε identities "carry the sign of det g" — det g = −1 in both
signatures, so that is not the mechanism (see §6). One convention, PSALTer's, serves both the
spectrum and the solver branch.

**Sources that do differ, and where they enter this programme:**

| source | signature | ε | where it enters |
| --- | --- | --- | --- |
| Ma & Bertschinger (astro-ph/9506072) | (−,+,+,+) | — | R-C validates FRW equations against it (`handoffs/R-C.md:69, 114`); `9506072.tex:295` (synchronous), `:361-363` (conformal Newtonian) |
| legacy TIDAL | (−,+,+,+) | ε₀₁₂… = −1 | M3's comparison with the frozen legacy oracle (§6); `tidal/wolfram/CommonUtilities.wl:30-34` |
| SymBoltz.jl v1.7.0 | (−,+,+,+) | — | a survey item for R-C only (`docs/src/conventions.md:7`) |
| Challinor & Lasenby (astro-ph/9804301) | (+,−,−,−) | η₀₁₂₃ = −√(−g), i.e. **opposite** to PSALTer | CAMB's formalism lineage; matters only if a parity-odd CMB formula (O4) is ever taken from it — none is planned (`observable_ladder.md:336-339` rotates CAMB's C_ℓ instead) |

### 2.1 Curvature signs, the connection's slots, and the contortion family

Added 2026-09-16 from R-C (#567), which found that **nothing in the project stated a curvature
convention**: this file fixed the signature and `ε` only, PSALTer defines no curvature and sets
none of xTensor's sign globals, and CAMB carries the Einstein equations in a fixed form with no
Riemann tensor. Adopted, in R-C's words:

> Curvature signs are xTensor's defaults, `$RiemannSign = $RicciSign = $TorsionSign =
> $epsilonSign = +1` (`xTensor.m:287-289`, defaults at `:1837-1843`), with the derivative
> index in the middle slot of the connection and of the contortion. They are asserted at the
> start of every kernel, re-asserted after every `Needs` and after any induced decomposition,
> and printed into every artifact header. xPand needs no adjustment; **xMAG sets
> `$RiemannSign = −1` when it loads** (`xMAG.m:110`, undocumented) and
> `$ExtrinsicKSign = $AccelerationSign = −1` inside `StartInducedDecomposition`
> (`xMAG.m:1792`), so a kernel that uses xMAG resets them.

**The contortion, in xTensor's slot meaning** (`perturbation_tooling.md` §7.1, sentinels
`RC_T1_A6_…` in run `t1/`): `K^a{}_{bc} = ½(T^a{}_{bc} + T_b{}^a{}_c + T_c{}^a{}_b)`, which is
`identical` to xMAG's `Contorsion`, has `K^a{}_{bc} − K^a{}_{cb} = T^a{}_{bc}` and
`K_abc + K_cba = 0`, and enters the rewrite as `ChristoffelCDCDT → −K` because
`CDT_b v^a − CD_b v^a = −ChristoffelCDCDT^a{}_{bs} v^s`. Einstein–Cartan check:
`R̃ = R − 6v² + 6∇·v` for `T^a{}_{bc} = δ^a_c v_b − δ^a_b v_c`. **Legacy's identity fails both
tests (#582)** — see `tests_cosmo/data/oracles/README.md` for which frozen specs carry it.

**Import map for "family A" sources** — papers whose connection carries the derivative index in
a different slot (Aoki et al. 2310.16007, Nikiforova–Damour 1804.09215, Shapiro
hep-th/0103093). Applied **once, at import, with the source cited and a test** (§1): their
torsion tensor is **minus** xTensor's `TorsionCDT`, their Riemann slot order is remapped, and
their contortion is re-expressed by the identity above. This applies to the SVT parametrization
R-C transcribed from Aoki et al.: the map flips the sign of its 24 potentials, which leaves the
counting and representability untouched but must travel into M3's rule set
(`perturbation_tooling.md` §7.2).

## 3. CAMB perturbation variables — adopt verbatim at the CAMB seam

CAMB exposes two naming families: **covariant** (`eta`, `hdot`, `sigma`, `z`) and
**Ma–Bertschinger-style synchronous** (`eta_s`, `hdot_s`), with stated conversions. Our FRW
derivation and all seam code are written in `camb.symbolic`'s named variables, CAMB's native
"CDM" frame. Sign hazards are in the definitions, not the signature:

| quantity | definition (as CAMB states it) | evidence |
| --- | --- | --- |
| frame "CDM" | zero acceleration, comoving with CDM: "Equivalent to the synchronous gauge but using the covariant variable names" | `symbolic.py:480-491` |
| `etak` / η | `etak = kη_s = −kη/2` in the CDM frame; `camb.symbolic` maps η via `camb_sub="-2*etak/k"` | `variables_guide.rst:133`; `symbolic.py:183-188` |
| η_s, ḣ_s | `eta = −2·K_fac·eta_s`, `hdot = hdot_s/6`; `ḣ_s = 6ḣ = 2kz` | `symbolic.py:505`; `variables_guide.rst:252` |
| covariant vs Ma–Bertschinger | "The covariant hdot and eta variables differ from Ma & Bertschinger by factors (the latter have _s after the name)" | `symbolic.py:481` |
| `z`, σ | `z = (0.5·dgrho/k + etak)/adotoa`; flat: `sigma = z + 1.5·dgq/k²`; `ayprime(ix_etak) = 0.5·dgq`; CDM: `clxcdot = −k·z` | `equations.f90:2427, 2430-2431, 2442` |
| φ | the **Weyl potential**; Φ_N = φ + ½a²κΠ/k², Ψ_N = φ − ½a²κΠ/k²; Newtonian gauge φ = (Φ_N+Ψ_N)/2, A = −Ψ_N, η = −2Φ_N·K_fac | `symbolic.py:181, 439-440, 443-449` |
| Φ_N, Ψ_N signs | "(default, as defined by Ma and Bertschinger …)"; Hu et al. use (1+2Φ_N) in the spatial part, "corresponding to a sign change in Φ_N" | `ScalEqs.ipynb:278-285` |
| `grho`, `gpres`, `dgrho`, `dgq`, `dgpi` | a²κρ, a²κp, a²κδρ, a²κq (heat flux), a²κπ (anisotropic stress) | `equations.f90:1119-1125` |
| κ | 8πG | `constants.f90:49` |
| `adotoa` | conformal Hubble rate, `sqrt(grho/3)` | `equations.f90:2365` |
| time | conformal time: `camb.symbolic`'s `t` = Fortran `tau` = TorC's τ = our η | `repo_reshape.md:443-445` |

**Competing definitions, named so none is chosen by accident:** Ma–Bertschinger / CLASS use
`h`, `η` (= CAMB's `hdot_s`, `eta_s` family) and ψ, φ in (−,+,+,+); Hu et al. flip Φ_N;
CAMB's φ is the Weyl potential, not a Newtonian potential. **A test that our variables match
`camb.symbolic`'s `camb_var`/`camb_sub` mapping** is the CAMB-side check (§5); it belongs to
I-532 / M1a.

## 4. PSALTer input requirements — what the Stage-1 package enforces

| requirement | source | enforced by |
| --- | --- | --- |
| quadratic in the perturbed fields, **linear in the couplings** | `PSALTer.m:80` (`ParticleSpectrum::usage`); arXiv 2406.09500 `Manuscript.tex:130` | one operator per coupling key (`lagrangian: {coupling: operator}`); an operator containing a coupling is refused |
| linear **parameterization** — not one coupling per monomial | same | an operator may be a composite expression (e.g. a torsionful curvature invariant expanded by the package); all its monomials share the coupling |
| couplings are symbolic constants; no function of a coupling as a coefficient (`1/kappa^2` is refused; Einstein–Hilbert gets its own generic coupling) | `ValidateLagrangian.m:32`; `stage1_engineering_plan.md:492-500` | name and whitelist checks |
| exact coefficients | PSALTer works symbolically | real-number atoms are refused (`0.5` → use `1/2`) |
| coupling and field names: letters and digits, starting with a letter, **no underscores** (`m_phi` parses as `Pattern[m, Blank[phi]]`) | Wolfram *OperatorInputForms* | name check before parsing |
| names must not be symbols PSALTer/xAct already own (`G`, `CD`, …) | `UpdateTheoryAssociation.m:8`; `DefGeometry.m:5-7` | context check after a held parse |
| field and coupling symbols live in `Global` | `DefField.m:37, 52`; `ValidateLagrangian.m:37` | the package creates them there |
| abstract indices a–z are PSALTer's public symbols | `DefGeometry.m:5` (read in the public context) | whitelist |
| `MaxLaurentDepth` ∈ {1,2,3}: 1 finds quadratic 1/k² poles; 2, 3 also quartic/hexic poles (always sick, expensive) | `PSALTer.m:82`; `Manuscript.tex:132` | derivation option, part of the fingerprint |
| PSALTer never checks linearity (`NonLinearCouplings` is declared, never thrown) and throws uncaught (`ParticleSpectrum.m:23`); it can `Quit[]` silently (`ConjectureInverse.m:37-39`) | source | refuse before PSALTer; `Catch` in the driver; verify the artifact, never the exit status |

## 5. Carriage: record, assert, and the checks that catch a slip

- Every derived spectrum's manifest records `"conventions": {"signature": [1,-1,-1,-1],
  "epsilon0123": 1}`; the Cobaya Theory refuses any other value (prototype case 4 in
  `scripts/research/interfaces/`).
- **Spectrum side:** a reference theory with a published healthy/ghost verdict runs through the
  whole pipeline in CI; a flipped sign flips the verdict and fails the test.
- **CAMB side:** the variable-mapping test of §3.

## 6. Comparing against legacy TIDAL or Ma–Bertschinger: physics, not bytes

The new package changes conventions deliberately, so its equations will not match legacy
text. What must match is the physics, and **every difference must be attributed** — to a
convention change or to a real change. M3's written mapping (`repo_reshape.md:787-790`)
therefore records, per term, a sign factor with three sources:

1. **metric contractions** — each inverse metric flips under g → −g, e.g. T_{λμν}T^{λμν} picks
   up (−1)³ (`docs/tex/pgt_stability_priors.tex:136-139`);
2. **field definitions** — how each field is defined relative to the metric (a perturbation
   defined as δg_ab itself flips);
3. **ε orientation** — legacy ε₀₁₂… = −1 against PSALTer's +1; the parity-odd oracle specs
   (`examples/torsion_gertsenshtein/theory_parity_odd.toml`) are where this appears.

The per-term table is derived and checked against the oracle **by M3**; this file states the
rule and its sources only. `tests_cosmo/data/oracles/README.md` carries the drift classes,
including the one added 2026-09-16: legacy's contortion identity is wrong, so the frozen torsion
specs encode wrong `R̃` terms and M3 attributes that difference rather than reporting a
regression (#582).

### 6.1 Running xPand in `(+,−,−,−)` — the measured rule (#586)

R-C ran xPand natively in our signature; the rule, from Probe C (`c/20260916T114552Z`):

- `SetSlicing[g, n, +1, h, cd, {…}, "FLFlat"]` — the normal's norm is `+1`; accepted with no
  messages, `n·n = 1`.
- The map from the mostly-plus literature has **exactly three sites**: the lapse flip, the
  Laplacian term, and the determinant sign. Dropping the Laplacian site leaves the residual
  `4εD²φ + 8εD²ψ`; dropping the lapse flip is `proved-different` (both controls run).
- **Two package sites hard-code `−1`** and are done by hand: `ExtractComponents`' Time
  projection (it returns `−V0`; by hand `+V0`) and the built-in `ToxPand` fluid path (use
  `SplitMatter[…, +1, …]`).
- Confirmation that this is CAMB-compatible: xPand's `00` constraint equals `camb.symbolic`'s
  own text in both gauges, `proved-equal` with `c = 1`.

A kernel that loads xMAG must re-assert the curvature globals afterwards (§2.1).

## 7. Open items

- The CAMB notes were read as extracted text; their equations were not machine-checked.
- No convention yet for the sign of the birefringence angle β or polarization handedness (O4).
- The per-term legacy flip table (§6) is M3's.
