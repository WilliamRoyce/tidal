# The cosmological magnetic field background for O3

**Status:** literature survey, 2026-08-30, produced for H2's O3 rung. O3 (Gertsenshtein
graviton↔photon mixing on FRW, `docs/COSMOLOGY_PROGRAM.md`) cannot be posed without an
assumed magnetic field: the mixing term is *linear* in an external coherent field, so with
no `B̄` there is no `h ↔ a` coupling at all. This document records what the
expanding-universe Gertsenshtein literature actually assumes, how it justifies the
assumption, and what O3 should therefore adopt.

> **Nothing here is a result.** It is a survey with the conventions and traps marked, so
> that O3's magnetic-field model is a recorded decision rather than an inherited habit.

---

## 1. Why a background field is needed, and why it is still a spectator

Two distinct objects must not be conflated:

- **The CMB** is the photon *gas* — the quanta being converted, and the distribution whose
  distortion is the observable.
- **The magnetic field** is a separate, coherent classical background that mediates the
  conversion. Cosmologically it is a **primordial magnetic field (PMF)**.

Adding a PMF does not take us outside the spectator limit. Flux freezing in the highly
conducting cosmic plasma keeps the comoving field `a²B` constant, so `B ∝ a⁻²` and
`ρ_B ∝ a⁻⁴`: the field is a **radiation-like** component whose relative size is

```text
r_B ≡ ρ_B(t)/ρ_γ(t) = B₀²/(8πρ_γ0) ≈ 10⁻⁷ B₋₉²     (Subramanian, arXiv:1504.02311, L649)
```

where `B₋₉` is the present-day field in nG; equivalently `ρ_B = ρ_γ` at `B₀ ~ 3.2 μG`
(L639–641). Paoletti et al. state the same relation with a
different rounding (arXiv:2204.06302 line 207):
`ρ_B(z) = ⟨B²⟩(1+z)⁴/8π ≈ 9.5×10⁻⁸ (⟨B²⟩/nG²) ρ_γ(z)`, i.e. equipartition at 3.24 μG.

At any observationally allowed amplitude (§4) this gives `r_B ≲ 10⁻¹⁰–10⁻⁶`, hence
`ΔN_eff ≲ 10⁻⁵` — utterly negligible in the Friedmann equation. **The PMF is a legitimate
spectator background in exactly the program's sense**, and its own validity is a cheap,
checkable number rather than an assertion (§6).

## 2. What each paper assumes — the survey

Line references are to the local TeX sources under `literature/<id>/`.

| Paper | Field model | `B₀` used | Bound cited | `λ_B` | `a`-scaling | Transverse factor |
| --- | --- | --- | --- | --- | --- | --- |
| Dolgov & Ejlli, `1211.0500` | primordial, uniform within `λ_B` | 3 nG (→ 3×10⁻³ G at recombination); figures also use 0.12 G and 1.2 G | 3 nG for `λ_B < 1 Mpc` (Paoletti & Finelli); COBE 5×10⁻⁹ G (Barrow 1997); Faraday 6×10⁻⁸–2×10⁻⁶ G | assumed `≫ λ_photon`, never specified (L266, L330) | `B(t_rec) = B(t₀)(1+z)²` (L631) | none; "we omit index `T`" (L518) |
| Cembranos et al., `2302.08186` | cosmic, uniform, static or slowly varying, global | 5×10⁻⁹ G (= the bound itself, L379); figure at recombination quotes `B = 3×10³ G` (L529) | COBE / Barrow 1997 (L366) | **none** — field treated as globally uniform, path length `η₀` | `B² = B_T²/a⁴` (L200) | `B_T` symbol only |
| Domcke & Garcia-Cely, `2006.01161` | primordial, stochastic `P_B(k)`, reduced to a domain-averaged rate × `𝓕` | 47 pG–nG scan | **Jedamzik & Saveliev 47 pG**; blazar lower limit `10⁻¹⁶ G`; MHD damping (L194) | `λ_B = ∫dλ B_λ²/⟨B²⟩`; `Δz₀ = min[λ_EQ, λ_B⁰]`, `λ_EQ/2π = 95 Mpc` (L204) | `ρ_B = ρ_B0(1+z)⁴`, `λ_B = λ_B⁰/(1+z)` (L192–193) | **explicit `B² → 2B²/3`** (L192 footnote) |
| He et al., `2312.17636` | primordial (cosmological) and astrophysical (cluster) | 0.1 nG at `λ_B ~ 1 Mpc`; 0.1 μG at 10 kpc for clusters | IGM `B < 4 nG` (O'Sullivan), `< 40 nG` (Vernstrom); filaments 30–60 nG; blazar `> 10⁻¹⁴ G` in voids (L228–235) | 1 Mpc / 10 kpc | **never written**; inherited through the borrowed rate formula | none explicit |
| Addazi et al., `2401.15965` | primordial, **stochastic Gaussian** with a broken-power-law `P_B(k)` | 10⁻¹² G fiducial; scan to 10⁻⁹ G | Durrer & Neronov exclusion plane (CMB + MHD + blazar) | `λ_B = 2πρ_B⁻¹∫ρ_B k⁻¹dk`; scanned; patches `Δl ~ 1–10 Mpc` | `B = B₀/a²`, `P_B(z,k) = (1+z)⁴P_B0(k)`; "diluted by expansion, neglecting its dynamical evolution" (L642–648) | full `k̂` projection (L403–406) |
| Kushwaha & Jain, `2502.12517` | primordial, domain-uniform | 47 pG – 1 nG at `ℓ₀ = 1 Mpc` (L279) | attributed to Planck 2015 — **mis-attribution**, 47 pG is Jedamzik & Saveliev | 1 Mpc | `B̄(z) = B̄₀(1+z)²`, `ℓ = ℓ₀/(1+z)` (L158) | folded into the symbol `B̄` |
| Cembranos et al., `1806.11020` | external, static, homogeneous, **Minkowski** | none given | none | none | none | **`B_T = B_e sin φ`** (L86) — the clearest statement of the projection |
| Paoletti et al., `2204.06302` | primordial, stochastic non-helical `P_B = A_B k^{n_B}` | derives `√⟨B²⟩ < 0.69 nG` | own result (Planck 2018 / SROLL2) | via damping scale `k_D(⟨B²⟩, n_B)`; **1 Mpc normalization** | `B(x,τ) = B(x)/a²(τ)`, `ρ_B ∝ a⁻⁴` (L162) | n/a |
| Lella et al., `2406.17853` | **Galactic** (Jansson–Farrar regular + turbulent) | μG scale, `l_corr ~ 10–100 pc` | — | — | — | `s_φ, c_φ` (L223) |
| Domcke et al., `2507.16609` | **static fields, flat spacetime** (neutron-star magnetospheres) | — | — | — | — | — |

The last two are listed to record that they are *not* cosmological: they solve the same
mixing problem in the Milky Way and in flat space respectively, and carry no `a`-scaling,
no `λ_B` cosmology, and no CMB distortion calculation.

### Conversion-probability scalings

Every paper gets `P ∝ B₀²`; they differ in how the coherence length enters.

| Paper | Result |
| --- | --- |
| `1211.0500` | `P ≈ 10⁻¹⁵ (ω_i/10⁵ eV)² (B_i/3×10⁻³ G)² (40.5 cm⁻³/n_e)²` (L664) |
| `2302.08186` | `P = Δ_M² e^{−Δ_g^I η}/(α²+β²)·[sinh²(βη) + sin²(αη)]`, `Δ_M = κ_eff B_T/2` (L352); GR value `8.6×10⁻¹²` |
| `2006.01161` | `𝓕 = 6.3×10⁻¹⁹ (B₀/nG)²(ω₀/T₀)²(Mpc/Δz₀)(I(z_ini)/10⁶)` (L206) — `∝ 1/λ_B` |
| `2312.17636` | `P ≈ 3.78×10⁻²⁰ (B/0.1 nG)²(f/f_eq)²(1 Mpc/Δl₀)(I/6×10⁶)` (L202) |
| `2502.12517` | `P⁺ = 4.04×10⁻¹⁹ h⁻¹ (Mpc/ℓ₀)(B̄₀/nG)²(ω₀/T₀)²(I/10⁶)`, `I(1100) ≈ 6.31×10⁶` |
| `2401.15965` | three regimes: `κ²B²Δl²` (short), `κ²B²l_osc²` (saturated, `λ_B > l_osc`), and **`κ²B²λ_B Δl`** — linear growth in `λ_B`, maximal at `λ_B = π l_osc` (the paper's resonance result) |

The first three cosmological results are the same formula in different conventions
(factors of `h`, the transverse `2/3`, and the definition of the patch length).

## 3. How the backreaction question is actually handled — a gap

**Not one of the Gertsenshtein papers surveyed justifies neglecting the field's effect on
the expansion history.** Each takes a fixed classical background and imports an
observational upper bound. Specifically:

- `1211.0500` is the only one that engages with a backreaction bound at all, and it
  *declines* it: Caprini & Durrer's anisotropic-stress/GW-generation limits
  (`B ≲ 10⁻²⁷ G` for an electroweak-phase-transition field, `10⁻³⁹ G` for inflationary)
  are cited at L533–539 and set aside at L541 — "we assume validity of the CMB bounds on
  the large scale magnetic fields, quoted above" — on the strength of a published
  criticism of those limits.
- `2006.01161` comes closest to a self-consistency check, but it checks a *different*
  quantity: the field's effect on the **gravitational-wave propagation speed**,
  `1 − v_g/c = c²κ²B²/2ω_pl² = 10⁻⁴⁶ (B/nG · Hz/ω_pl)²`, against the neutron-star-merger
  bound of `10⁻¹⁶` (Supplementary L136–139).
- `2401.15965` argues qualitatively from isotropy that the coupling "can be treated as a
  perturbative interaction" (L331). Its only explicit backreaction discussion concerns the
  *inflaton*, not the expansion (L760).
- `2302.08186`, `2312.17636`, `2502.12517` make no statement.
- `2204.06302` — a PMF-constraint paper, not a mixing paper — is the only source in the set
  that writes `ρ_B/ρ_γ` down, and it does so as a modeling ingredient for plasma heating.
  Its "backreaction" caveat (L161) is the *fluid's* effect on the field, not the field's on
  the expansion.

So the justification exists and is easy (§1), but it lives in the PMF reviews rather than
in the conversion literature. **This is the same "asserted, never enforced" pattern the
program already documents for the spectator limit itself**
(`docs/cosmology/spectator_route.md` §2), and it is handled the same way: by computing the
number per run (§6).

## 4. The standard PMF facts, with sources

From Durrer & Neronov (arXiv:1303.7121) and Subramanian (arXiv:1504.02311), which are the
sources every paper above cites:

- **Scaling.** Maxwell theory is conformally invariant and the plasma is highly
  conducting, so flux freezing gives `B ∝ a⁻²` and `ρ_B ∝ a⁻⁴` (Durrer & Neronov
  `section2.tex` L160, L167: `B_i ∝ 1/a`, `B^i ∝ a⁻²` in the orthonormal basis, `B² ∝ a⁻⁴`). The comoving field `a²B`
  is the constant, and a PMF therefore does not change the *functional form* of the
  expansion history — it adds a radiation-like component of relative size `r_B`.
- **Energy budget.** `r_B ≈ 10⁻⁷ B₋₉²`; equipartition with the CMB at `B₀ ~ 3.2 μG`
  (Subramanian L649, L639). At generation, Durrer & Neronov write
  `ρ_B/ρ_rad ≈ (2/g_eff)^{1/3}(B_*/3×10⁻⁶ G)²`.
- **Smoothing convention.** The field is quoted smoothed on a scale `λ`, with
  `B_λ² ≡ 8π P_B(2π/λ)/λ³` and coherence length `λ_B = ∫dλ B_λ²/⟨B²⟩`. The community
  normalization scale is **1 Mpc**; "`B_{1 Mpc} < few nG`" always means this.
- **Observational window.** Lower bound `B ≳ 10⁻¹⁶ G` on Mpc coherence from the
  non-observation of TeV-blazar electromagnetic cascades; upper bounds of a few nG on
  1 Mpc from the CMB, tightened to **47 pG (95% C.L.) for a scale-invariant spectrum** by
  Jedamzik & Saveliev (arXiv:1804.06115, PRL 123, 021301) via magnetically induced baryon
  clumping at recombination; plus an MHD-turbulence damping bound that cuts a diagonal
  through the `(B₀, λ_B)` plane at small scales.
- **Planck 2015 XIX** (arXiv:1502.01594) gives `B_{1Mpc} < 4.4 nG` (zero helicity, CMB
  spectra), `< 2.0 nG` near-scale-invariant (`< 0.9 nG` including ionization-history
  effects), and a much weaker `< 1380 nG` from Faraday rotation.

### Two structure-dependent refinements

- **Small-scale power (Domcke & Garcia-Cely, Supplementary L294–295).** For
  `P_B(k) ∝ k^{−α}`, the rate acquires a boost `𝓕 ~ (l_osc/λ_B)^{α−3}`. A scale-invariant
  spectrum (`α = 3`) gives `𝓕 ~ 1`; MHD damping makes `α ≈ 4` "the more realistic (and
  conservative) scenario", giving `𝓕 ~ l_osc/λ_B` and recovering the simple domain formula
  with `Δz = λ_B`.
- **Resonance (Addazi et al.).** When `λ_B ≲ l_osc` the conversion grows *linearly* in
  `λ_B` rather than saturating, peaking at `λ_B = π l_osc`. Under the cosmological plasma
  conditions of §5 this regime is not reached, but it is the reason a stochastic treatment
  can differ from a domain model by orders of magnitude.

## 5. Why the plasma, not the field, sets the oscillation length

At allowed amplitudes the magnetic term is negligible in the oscillation length compared
with the plasma detuning:

```text
l_osc⁻¹ = ½√( ω²(1−μ)² + κ²B² ) ,    μ = √(1 − ω_pl²/ω²)
        ≈ (1+z)² X_e(z) ω_pl,0²/(4 ω₀)          for B₀ ≲ nG
```

(`2006.01161` L212; `2502.12517` L156; `2401.15965` L367). The consequence is structural,
and it is what forces O3's solver design: `l_osc ≪ 1 pc ≪ λ_B ~ Mpc`, so the conversion
probability oscillates many times within one coherence patch and the literature
**averages over patches and accumulates a rate along the line of sight** rather than
resolving the oscillation:

```text
⟨Γ_{g↔γ}⟩ = c|K₁₂|² l_osc²/(2Δz) ,    𝓕 = ∫_{l.o.s.} ⟨Γ⟩ dt
```

with `Δz = min[λ_EQ, λ_B⁰]/(1+z)`. Electron fraction and density enter through
`n_e(z) = n_b0(1+z)³X_e(z)`, `n_b0 = 0.251 m⁻³`, with `X_e = 1, 0.68, 0.0002, 0.15` at
`z = 0, 10, 20, 1100`; the redshift integral is
`I(z_ini) = ∫₀^{z_ini} dz (1+z)^{−3/2} X_e^{−2}(z)`, dominated by `z ~ 10`.

## 6. What O3 adopts

**Baseline model** (all of it an *assumed input*, declared per run, never derived):

| Quantity | Choice | Source |
| --- | --- | --- |
| Amplitude | `B₀ = 47 pG` baseline; `1 nG` as an optimistic column | Jedamzik & Saveliev; the range used by `2006.01161` and `2502.12517` |
| Smoothing / coherence | `λ_B⁰ = 1 Mpc` | community convention |
| Patch length | `Δz(η) = min[λ_EQ, λ_B⁰]/(1+z)`, `λ_EQ/2π = 95 Mpc` | `2006.01161` L204 |
| Scaling | comoving-constant `F̄_ij`; `B_phys ∝ a⁻²`, `λ_B ∝ a` | Durrer & Neronov |
| Geometry | transverse projection `B² → 2B²/3` for an isotropic average | `2006.01161` L192 |
| Structure | uniform-within-patch; stochastic `P_B(k)` with `α ≈ 4` as a later refinement | `2006.01161`, `2401.15965` |

Because `P ∝ B₀²` and `∝ 1/λ_B`, the choice is not a detail: the literature's spread from
47 pG to 5 nG is a factor `~10⁴` in conversion probability. Any O3 bound must be quoted
*with* its assumed field, and the report's cost/target tables do so.

**Enforced flags** (the program's honest-flags style, `gauge_certificate` lineage) — the
point of this document is that these are computed, not asserted:

- `r_B = ρ_B/ρ_γ ≈ 10⁻⁷ (B₀/nG)²` reported per run, against the `ΔN_eff ≲ 0.1` budget
  already used for the new sector;
- `l_osc ≪ Δz` — the precondition for patch averaging, checked rather than assumed;
- `P ≪ 1` and the number of patches `N = D/Δz` with `P·N ≪ 1`, the condition the domain
  sum requires (`2502.12517` L130).

## 7. Traps and discrepancies found in the survey

- **`2502.12517` mis-attributes the 47 pG bound** to Planck 2015; it is Jedamzik &
  Saveliev 2019. Cite the latter.
- **`2302.08186`'s figure value `B = 3×10³ G` at recombination** (L529) is `10⁶` times the
  recombination-era field implied by Dolgov & Ejlli's 3 nG today. Either a different
  normalization or a typo; do not adopt the number without re-deriving it from `B₀`.
- **`1211.0500`'s larger benchmarks (0.12 G, 1.2 G at recombination) exceed the 3 nG
  present-day bound the same paper cites.** They are illustrative, not allowed.
- **Two `κ` conventions** circulate in this literature; the local project convention and
  the known normalization errors are documented in `docs/tex/gertsenshtein_formula.tex`
  (§`sec:kappa-conventions`). The papers surveyed here all agree with `κ² = 16πG`.
- **Galactic and flat-space results do not transfer.** `2406.17853` (Jansson–Farrar
  Galactic field) and `2507.16609` (static fields, flat spacetime) solve the same mixing
  algebra under entirely different backgrounds.

## 8. Gaps

- No paper in this set quantifies the error of the test-field treatment of the magnetic
  field — the same gap the program records for the spectator approximation generally.
- The stochastic-field treatment (`P_B(k)`, the `𝓕` boost, the `λ_B ≈ π l_osc` resonance)
  is a refinement over the domain model, not the baseline; adopting it changes the answer
  by orders of magnitude in some corners of the `(B₀, λ_B)` plane and should be revisited
  before any published O3 bound.
- Magnetic-field evolution across epochs (MHD processing between generation and
  recombination) is explicitly flagged as future work by He et al. (L236–239) and is not
  modeled by anyone in the set.

---

## See also

- `docs/cosmology/observable_ladder.md` — H2's feasibility study; O3's rung consumes this
  document's §6.
- `docs/cosmology/spectator_route.md` — the spectator limit and its validity criteria.
- `docs/tex/gertsenshtein_formula.tex` — the project's verified flat-space formula and the
  `κ` conventions.
