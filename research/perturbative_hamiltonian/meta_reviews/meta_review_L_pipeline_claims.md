# Meta-Review L — Independent Verification of Review 3's TIDAL-Pipeline Claims

**Auditor**: independent meta-reviewer, posture = "verify, don't trust"
**Date**: 2026-04-27
**Subject**: `research/perturbative_hamiltonian/reviews/review3_project_relevance.md`

---

## TL;DR

Review 3's high-level conclusion — *"Path A and Path B do not produce a TIDAL-pipeline-usable Hamiltonian for the metric h₄/h₇/h₉ subspace, so for the actual project the academic results are interesting but not operational"* — is broadly **correct**, but several of its detailed pipeline claims are **wrong or misleading**, and at least one of its "showstoppers" is materially weaker than presented. The most important reliability issues:

| Claim | Verdict | Notes |
|------|---------|-------|
| L1 — `_resolve_time_derivative` returns None for time_order ≥ 3 | **CORRECT** (literally) | But irrelevant to the actual JSON, which only carries `mixed_2_*` |
| L2 — operator schema only supports `mixed_T_*` with T ≤ 2 in practice | **PARTIALLY CORRECT** | Wolfram emits `mixed_T2`/`mixed_T3`/`mixed_T4` and modal solver registers up to T4; only the validator's *warning gate* and `_resolve_time_derivative` cap at T=2 |
| L3 — TIDAL is "polynomial-coefficient", no rational denominators | **WRONG** | Existing `coefficient_symbolic` strings include `1/(2*kappa^2)`, `B0^2/8`, `1/8*B0^2` — i.e. rational in *parameters*. `eval()` based resolution accepts arbitrary Python expressions |
| L4 — h_4/h_7/h_9 are metric (not torsion) and are the constraint-promoted PU fields | **CORRECT** | JSON: `time_order=4, kc=2*b5`; Hamiltonian terms 75/76/77 are `(∂_t² h_X)²` |
| L5 — Path B's metric exclusion is the Gertsenshtein bottleneck | **WRONG (or at minimum unproven)** | The actual conversion observable in every sweep is `h_5 → a_1`. h_5 is `time_order=2, kc=-kappa^(-2)` (standard graviton), NOT one of the constraint-promoted fields. Both the h_5 EOM and the a_1 EOM only reference `{a_1, h_5}` |
| L6 — `_conversion.py` is amplitude-based, not Hamiltonian-based | **WRONG** | `_conversion.py` calls `compute_energy_timeseries`, which directly evaluates `canonical.hamiltonian_terms` (energy ratio) |
| L7 — cross-sector couplings unmodelled by toys | **CORRECT** but with caveat | The active observable's two equations are sector-isolated; cross-sector content lives in the *gauge / non-propagating* h_0…h_3 algebraic equations and the b5-promoted h_4/h_7/h_9 self-terms |
| L8 — t_end-independence is a documented workaround | **CORRECT** | `CLAUDE.md:67` codifies it for distinguishing genuine amplification from tachyonic artefact |
| L9 — campaigns show NULL amplification | **CORRECT, with major caveats Review 3 omits** | The 276-run dark-photon null was retracted as solver bug #257 (rank-deficient eigenbasis), and the propagating-PGT NULL holds only inside a `|δ₁|<0.005` window |

The net effect: Review 3 **correctly identifies that the Round 1-3 papers don't deliver an operational Hamiltonian for h₄/h₇/h₉**, but it **mis-identifies which fields actually carry the project's headline observable**, and it **mis-reports what `_conversion.py` actually computes**. Both errors point in the same direction — overstating how blocked the project is.

---

## Per-claim audit (file:line evidence)

### Claim L1 — `_resolve_time_derivative` and the time_order ≥ 3 gap

**Verdict: literally CORRECT, but its consequence as Review 3 frames it overstates the problem.**

`tidal/measurement/_energy.py:430-464`:

```python
def _resolve_time_derivative(
    data: SimulationData,
    field: str,
    t_idx: int,
    time_order: int,
    params: dict[str, float],
) -> NDArray[np.float64] | None:
    """Resolve the Nth time derivative of a field from stored data + EOM.
    ...
    ≥ 3        Not supported (returns None)
    """
    if time_order == 0: return _resolve_term_target(...)
    if time_order == 1: return data.velocities[field][t_idx] or _differentiate_constraint(...)
    if time_order == 2:
        accel = _compute_acceleration_from_eom(...)
        if accel is not None: return accel
        return _differentiate_constraint_twice(...)
    # time_order >= 3: not currently supported
    return None
```

So the literal claim — "returns None for time_order ≥ 3" — is true (line 463-464). But Review 3 implicitly suggests this rules out energy evaluation of the b5·R̃² metric Hamiltonian. The actual Hamiltonian terms for h_4/h_7/h_9 use `mixed_2_0_0_0` (i.e. `(∂_t² h)²`, time_order = 2), which is fully supported via `_compute_acceleration_from_eom`. The "EOM-defines-Hamiltonian" critique is real, but it is not blocked by L1.

A jet-5 polynomial like Path A's L_VT *would* require time_order ≥ 3 evaluation. So Review 3's L1-as-applied-to-VT is correct; L1-as-applied-to-the-current-JSON is a non-issue.

### Claim L2 — operator schema's max T

**Verdict: PARTIALLY CORRECT.** The schema *parser* and the *Wolfram emitter* support T up to 4, and the modal solver explicitly registers `mixed_T4_*`. Only the validator's LPS warning gate and `_resolve_time_derivative` cap at T ≤ 2.

Evidence:

- `tidal/wolfram/ExportJSON.wl:1239`:
  ```mathematica
  parts = {"mixed_T" <> If[timeOrder > 1, ToString[timeOrder], ""]};
  ```
  No upper cap on `timeOrder`.
- `tidal/solver/modal.py:180-185` registers `mixed_T4_S1x`, `mixed_T4_S1y`, `mixed_T4_S1z`, `mixed_T4_S2x/y/z` as eigendecomposable operators. So the modal solver explicitly handles T=4.
- `tidal/symbolic/json_loader.py:96-99` accepts both legacy numeric (`mixed_1_0_0_1`) and Wolfram-style (`mixed_T2_S2x`) names with no T cap.
- `tidal/cli/_validate.py:189-207` issues a *warning* (not error) for `mixed_T_* with T ≥ 2`, calling it "the irreducible LPS residue". This is documented as architectural barrier #321.

So Review 3's claim that "TIDAL's JSON schema operator names support only `mixed_T_S1_..._SD`" is true; its claim that T is effectively capped at 2 is *partially* true (the energy evaluator caps at 2, the validator warns at ≥ 2), but the broader schema and the modal solver handle up to T=4.

### Claim L3 — Routhian denominators / polynomial-coefficient pipeline

**Verdict: WRONG.** TIDAL accepts arbitrary Python-`eval`-able coefficient strings.

Evidence:

- `tidal/symbolic/json_loader.py:807-837`:
  ```python
  def _resolve_symbolic_coeff(sym: str, parameters: Mapping[str, float]) -> float | None:
      ...
      try:
          py_expr = sym.replace("^", "**")
          result = eval(py_expr, {"__builtins__": {}}, dict(parameters))  # noqa: S307
          value = float(result)
      except (NameError, SyntaxError, TypeError, ValueError, ZeroDivisionError):
          ...
          return None
  ```
  The `eval` accepts any expression Python can parse, including rationals.
- The actual `torsion_gertsenshtein.json` already contains `coefficient_symbolic` strings such as
  `"-1/2*1/kappa^2"`, `"3/kappa^2"`, `"1/(2*kappa^2)"`, `"-1/8*B0^2"`. These are field-rational in the *parameters* — exactly the form Review 3 says TIDAL "does not support".

What TIDAL does *not* support directly is *spatially-rational* coefficients (rational in the field values themselves), which is what a Routhian projector denominator like `2·M_1²·M_2²·M_3²` would actually be in cases where M_i depend on dynamical fields. If M_i are masses (parameters), the denominator is a parameter-rational, fully supported. So Review 3 conflates parameter-rational (supported) with field-rational (not directly supported) — the L_VT obstruction may still be real, but as stated L3 is incorrect.

### Claim L4 — h_4, h_7, h_9 are metric components, not torsion

**Verdict: CORRECT.** Confirmed by direct JSON inspection.

Evidence (jq queries on `examples/data/torsion_gertsenshtein.json`):

```
field h_4: time_order=4, kinetic_coefficient_symbolic="2*b5"
field h_7: time_order=4, kinetic_coefficient_symbolic="2*b5"
field h_9: time_order=4, kinetic_coefficient_symbolic="2*b5"
```

and (Hamiltonian terms 75/76/77):
```
{coefficient: 1, factor_a: {field: h_4, operator: mixed_2_0_0_0},
 factor_b:                   {field: h_4, operator: mixed_2_0_0_0},
 term_class: "self", coefficient_symbolic: "-b5"}
```
(same for h_7 and h_9.)

The JSON has 38 fields total: `a_0…a_3` (4 photon), `h_0…h_9` (10 metric), `t_0…t_23` (24 torsion). The b5-promoted `time_order=4` fields are 3 of the 10 metric components. Per the `[torsion]` TOML section, torsion lives in `t_*` not `h_*`. So h_4/h_7/h_9 are metric.

### Claim L5 — Path B's metric exclusion is THE bottleneck

**Verdict: WRONG (or at minimum undefended).** The active conversion observable does *not* go through h_4/h_7/h_9.

Evidence:

- All sweep scripts in `examples/torsion_gertsenshtein/` use `--source h_5 --target a_1`:
  ```
  examples/torsion_gertsenshtein/sweep_B0_scaling.sh:22:  --measure peak_conversion --source h_5 --target a_1
  examples/torsion_gertsenshtein/sweep_alpha2_hires.sh:18:  --measure peak_conversion --source h_5 --target a_1
  examples/torsion_gertsenshtein/sweep_delta1.sh:63:  --source h_5 --target a_1
  examples/torsion_gertsenshtein/sweep_phase1_1d.sh:23:    --source h_5 --target a_1
  examples/torsion_gertsenshtein/sweep_phase2_pairs.sh:18:    --source h_5 --target a_1
  ```
- `h_5` is **not** a constraint-promoted field. Per JSON: `time_order=2, kinetic_coefficient_symbolic="-kappa^(-2)"`. It is a standard graviton component — the h_+ TT-like polarisation in this gauge.
- Both the `d2_t(h_5)` and `d2_t(a_1)` equations couple ONLY to `{h_5, a_1}` (verified by jq listing `rhs.terms[].field | unique`). The b5-promoted h_4/h_7/h_9 do not enter these equations.
- The h_4/h_7/h_9 self-energy Hamiltonian terms (75/76/77) are *self-self*: `(∂_t²h_4)·(∂_t²h_4)`. They do not enter the h_5 or a_1 self-energy or the cross-coupling interaction energy that drives conversion.

Consequence: the conversion observable `P(h_5 → a_1)` can be computed by `_conversion.py` **without ever touching the constraint-promoted subspace**. The total energy denominator (interaction + all per-field) does include h_4/h_7/h_9 contributions, but `compute_conversion_probability` (line 155) divides by *source-field* energy at t=0, not total — so even the denominator is sector-clean.

This is the most important reliability issue in Review 3. Its showstopper #1 ("Path B doesn't cover the metric h₄/h₇/h₉ subspace") only matters if those fields gate the headline observable. They don't — the headline observable runs through the standard-kinetic h_5 ↔ a_1 channel.

(Caveat: the *full* graviton energy budget includes h_4/h_7/h_9 self-energies, which use `mixed_2_*` operators and therefore do depend on `_compute_acceleration_from_eom`. So if the user wanted *total system energy conservation*, h_4/h_7/h_9 do enter and are evaluated via the EOM-from-RHS pathway — the very pattern §pr-cb-conclusion criticises. This is real but it is not a blocker for `P(h_5 → a_1)`.)

### Claim L6 — `_conversion.py` is amplitude-based

**Verdict: WRONG.** It is energy-based and Hamiltonian-based.

Evidence: `tidal/measurement/_conversion.py:1-78`:

- Module docstring (line 4-7): `P(t) = E_target(t) / E_source(0)`.
- `_per_field_energy_timeseries` (line 66-78) calls `compute_energy_timeseries` from `_energy.py`, with this comment:
  > "Uses `compute_energy_timeseries` which correctly handles volume weights (`sqrt|g_spatial|`), operator-aware gradient axes, BC types, and position-dependent masses."
- Line 92-95 of `compute_conversion_probability`:
  > "Energy is computed via the canonical Hamiltonian (kinetic + gradient + mass), including the spatial volume element `sqrt|g_spatial|` for curved coordinates."
- Line 127-128:
  > "Use Hamiltonian energy (volume-weighted, operator-aware gradient axes)."

Review 3 even cites `_conversion.py` as "the operationally usable pathway" for amplitude-based conversion (line 653-654, 632-637). That recommendation is doubly wrong: `_conversion.py` is **already** Hamiltonian-based, and its energy ratio is what the entire amplification campaign uses (and reports as 6.7e-6 precision against `sin²(κB₀t/2)`).

If Review 3's *intent* was "use a non-Hamiltonian observable like raw amplitude ratios", then the recommendation is to *write* such code — not to flag the existing `_conversion.py` as "already amplitude-based". The recommendation conflates "energy ratio of dynamical fields" with "amplitude ratio".

### Claim L7 — Cross-sector couplings unmodelled by toys

**Verdict: CORRECT for the toys, but the conversion observable's two equations are sector-isolated.**

Evidence: cross-sector coupling does exist in this theory (the JSON's a_2, h_0, h_3 algebraic constraints reference each other and the b5-promoted block), but the `h_5 ↔ a_1` channel that the campaign actually uses is a clean 2-equation closed sub-system. This was already noted in `dark_photon_amplification_campaign_v0.31.md` (the sub-system was bit-exact to `sin²(κB₀t/2)` to 1e-6 across 276 runs).

So the toy isolation criticism is academically valid but operationally moot for the headline observable. For total-energy conservation tests it would matter.

### Claim L8 — t_end independence test is documented

**Verdict: CORRECT.** `CLAUDE.md:67`:

> "**t_end independence test for conversion amplification**: After measuring P_torsion/P_GR, ALWAYS verify at two different t_end values (e.g., t and 2t). If A(2t)/A(t) ≈ 1 → genuine amplification. If A(2t)/A(t) >> 1 → tachyonic instability artifact (see #238). B₀ scaling does NOT distinguish amplification from instability (growth rate is B₀-independent). IC amplitude must be ≪ B₀ for valid linearization."

Plus `CLAUDE.md:68` adds the perturbative-P-regime requirement (P_max ≪ 1) for trust in A. This is exactly the kind of pragmatic workaround Review 3 advocates and accurately describes what the project already practises.

### Claim L9 — Existing campaigns show NULL amplification

**Verdict: CORRECT for the headline result, but Review 3 omits major plot twists that materially weaken its citation.**

Evidence:

- `dark_photon_amplification_campaign_v0.31.md:9-14`: "the 276-run null result is not a bug — it is the exact prediction of the Lagrangian, a manifestation of Holdom's 1986 theorem". Then in lines 428-449 the same memory file says "the apparent null is an artefact of issue #257 [rank-deficient eigenbasis decoupling]" — i.e. an internal retraction. Then again at lines 9-14 a *re*-confirmation as physical Holdom triviality. The file is structurally a multi-stage investigation log, not a single coherent verdict.
- `propagating_model_finding.md:10-13`: "stability window is `|δ₁| < 0.005`. Within this window, A = P/P_EM = 1.0 exactly (zero amplification)". So "A=1.0" only holds inside an extremely narrow stability window; the propagating model is generically unstable and the null is not a robust finding across the full parameter space.

Both files do support "null amplification at the chosen-stable points", but neither cleanly says "the physics answer is established and the Hamiltonian gap doesn't matter". The campaign findings are themselves contested by their own authors.

Review 3 cites these to argue urgency-reduction ("the science answer is no amplification, so the Hamiltonian question is less pressing") — but this is over-reading the citations. The campaigns are partial null observations that *would* benefit from a clean Hamiltonian observable to confirm.

---

## Overall assessment

**Reliability of Review 3: medium.** Its high-level conclusion is right (the academic results don't operationalise into a TIDAL pipeline artefact for the most stringent definition of "Hamiltonian observable"), but its detailed pipeline claims contain at least three significant errors:

1. **`_conversion.py` is Hamiltonian-based, not amplitude-based** (Review 3 line 631-637 and 653-654 are wrong).
2. **The metric h₄/h₇/h₉ subspace is NOT the bottleneck for the actual conversion observable** (Review 3 lines 117-172 and showstopper #2 are misframed). The active observable is h_5↔a_1, which is sector-clean.
3. **TIDAL is not "polynomial-coefficient"** — the existing pipeline already supports rational-in-parameters coefficients (Review 3 line 79-81 is wrong).

These errors all push in the same direction: making the project look more blocked than it is. A reader would conclude from Review 3 that the project genuinely cannot proceed; in fact, the existing `_conversion.py` + `h_5→a_1` measurement pathway is already operational, has run a 276-run HPC campaign, and produced quantitatively trustworthy results to ~1e-6 precision against the analytical Boccaletti formula. The b5-promoted h_4/h_7/h_9 sector is a real architectural concern for *total system energy conservation* checks but is not the gate to the headline physics.

**Where Review 3 is right and useful**:
- The toy-to-PGT inferential gap is real (gauge structure, cross-sector couplings, parity-odd content).
- The Path A L_VT inheriting PU structure is a genuine theoretical limitation: at b5≠0 the phase space jumps from 6 → 30, and no current Path B sector covers that.
- The recommendation to update the TeX writeup with the strengthened no-go (Round 2 Agent G + quantitative phase-space count from Agent I) is sensible.
- Tier 1-2 pragmatic engineering recommendations (document non-Hamiltonian fallback observable, build a TT-projected energy variant) are actionable. (Although the fallback is already the de facto practice.)

**What Review 3 missed**:
- The `h_5 → a_1` pathway is sector-clean and already provides quantitatively trustworthy P(t) — so the user's claim "without a working Hamiltonian we cannot do anything with the simulations" is itself slightly overstated; the project HAS been doing things with simulations, with real measurement-grade rigor, for a year.
- The `compute_conversion_probability` energy-ratio observable is sound for the headline channel; the constraint-promotion residue affects only the *total system energy* and the h_4/h_7/h_9 self-energy contributions to it.
- The `mixed_T4_*` modal-solver registrations and `mixed_T2_*`-supported energy evaluation already exist — the schema is more capable than Review 3 implies.

**Net recommendation to the user**:
1. Treat Review 3's *direction* as correct — the constraint-promotion barrier for h_4/h_7/h_9 is genuine theoretical territory and Path A/B don't bridge it for the metric subspace.
2. Treat Review 3's *details* with caution — the conversion observable is already operational on the standard-kinetic graviton channel and does not require the Hamiltonian fix to produce trustworthy P(t) timeseries for `h_5 → a_1`.
3. The remaining real architectural concern is *total system energy conservation* in the presence of h_4/h_7/h_9 self-energy terms — which is a separate, narrower problem than Review 3 frames it as.
4. The TeX writeup update Review 3 recommends (strengthened no-go + phase-space count) is appropriate. The "redirect investigation to amplitude-based observables" recommendation is unnecessary because the existing energy-based observable is already trustworthy on the channel that matters.

---

## File evidence cited

- `tidal/measurement/_energy.py:430-464` — `_resolve_time_derivative`, T≥3 returns None
- `tidal/measurement/_energy.py:1080-1126` — `_prepare_hamiltonian_context`
- `tidal/measurement/_energy.py:1195-1280` — `_evaluate_single_hamiltonian_term`
- `tidal/measurement/_energy.py:1284-1352` — `_compute_hamiltonian_per_field`
- `tidal/measurement/_energy.py:1425-1479` — `compute_energy_timeseries`
- `tidal/measurement/_conversion.py:1-172` — Hamiltonian-based conversion probability
- `tidal/wolfram/ExportJSON.wl:1227-1252` — `BuildMixedOperatorName` (no T cap)
- `tidal/solver/modal.py:158-185` — modal-solver `mixed_T1`/`T2`/`T3`/`T4` registrations
- `tidal/symbolic/json_loader.py:80-101` — operator name validation (T-uncapped)
- `tidal/symbolic/json_loader.py:807-837` — `_resolve_symbolic_coeff` (Python `eval` based)
- `tidal/cli/_validate.py:145-228` — `_check_perturbative_consistency` (warns at T≥2)
- `examples/torsion_gertsenshtein/theory.toml:117-148` — Lagrangian, gauge, perturbation declaration
- `examples/torsion_gertsenshtein/sweep_*.sh` — all sweeps use `--source h_5 --target a_1`
- `examples/data/torsion_gertsenshtein.json` — h_4/h_7/h_9 `time_order=4, kc=2*b5`; h_5 `time_order=2, kc=-kappa^(-2)`; 78 hamiltonian_terms total, 3 `mixed_2_0_0_0` self-self on h_4/h_7/h_9
- `CLAUDE.md:67-68` — t_end independence + perturbative-P regime guidance
- `~/.claude/projects/-workspaces-torsion-gertsenshtein/memory/dark_photon_amplification_campaign_v0.31.md` — multi-stage investigation log; null result alternately attributed to Holdom triviality and to solver bug #257
- `~/.claude/projects/-workspaces-torsion-gertsenshtein/memory/propagating_model_finding.md` — A=1.0 holds only inside `|δ₁|<0.005` window
