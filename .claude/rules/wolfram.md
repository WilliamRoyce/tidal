---
paths:
  - "tidal/wolfram/**"
  - "tidal/symbolic/_derive.py"
  - "examples/**/theory.toml"
  - "examples/**/*.wls"
---
> **SCOPE (2026-09-04): two program-level policies apply on top of everything below.**
>
> **1. Wolfram runs at derivation time only.** Nothing symbolic happens at sampling time.
> This is what the two-stage spectrum architecture (#495) and the symbolic eikonal reduction
> (#504) both exist to guarantee: run Wolfram once, export a numeric contract, evaluate that
> per sample. A design that needs a kernel inside a likelihood call is the wrong design.
>
> **2. #513 — the new package emits CAMB and PSALTer conventions NATIVELY.** When porting
> `ExportJSON.wl` or any of this pipeline into `tidalcosmo/`, do **not** faithfully preserve
> legacy TIDAL notation: discarding it is the point. Legacy has no backward compatibility to
> protect. Conformal time agrees with `camb.symbolic` already (their `t` = Fortran `tau` =
> our `eta`) — do not build a conversion layer at that seam. Gauge becomes an explicit named
> input (covariant / Newtonian / synchronous), applied symbolically in Wolfram, recorded in
> the spec as first-class metadata, and asserted at the CAMB seam.
>
> Consequence: the port gate **cannot be a byte diff** — `tidal inspect OLD --diff NEW` would
> report every intended change as a failure. Equivalence is semantic, against the frozen
> oracle (#525), recorded as a written mapping.

# Wolfram/xAct Rules

## Symbol Management
- **No underscores in constant names** — `X_Y` = `Pattern[X, Blank[Y]]` in Mathematica. Use `Bpeak` not `B_peak`.
- Check existence before defining: `If[!xTensorQ[M2], DefManifold[...]]`
- Use `DefConstantSymbol` for mass/coupling (not bare Symbol)
- Standard shared manifolds: M2 (1+1D), M3 (2+1D), M4 (3+1D)

## Critical Patterns
- `SeparateMetric` BEFORE `ToBasis` — ensures correct signs on temporal components
- `MetricOfCovD[CD]` (SINGULAR) returns metric for a CovD
- Pattern variable pitfall: never reassign function parameter `eom_` — use `Module` local variable
- Parenthesize multiline Lagrangians: `L = (term1 + term2);`

## xCoba Rules vs DownValues
- `ComponentValue`/`MetricInBasis` use xCoba internal rules — do NOT auto-evaluate
- Direct `Set` (e.g., `g[{0,-chart},{1,-chart}] = 0`) creates auto-evaluating DownValues
- `SetMetricDownValues` and `SetBackgroundFieldDownValues` in CommonUtilities.wl use direct `Set`
- Critical for performance: eliminates off-diagonal zeros during TraceBasisDummy

## Background Fields
- `EvaluatePDBackgroundField` (3 rules) for derivative evaluation
- Applied after metric eval, after ConvertCDToDerivatives, and in batch loop
- Threaded via `"BackgroundFieldRules"` option in DecomposeToComponents

## Matter Perturbation (xPert)
- `SetupMetricPerturbation` MUST always be called (xPert global requirement)
- Matter-only: call it, then zero h terms: `l2Raw /. hpert[LI[_], idx___] :> 0`
- Naming collision: field "A" → `geA`, perturbation "a" → `gea` (uncapitalized)

## Cross-Field Decomposition
- Pass ALL coupled fields in `additionalFields` list to `DecomposeToComponents`
- Without this, cross-field coupling terms are silently dropped

## Equation Generation: Measure and LHS Normalization
These two produce **silently wrong equations** — no error, no warning. Both are
invisible on flat Minkowski, so any change here must be validated on at least one
curved-metric AND one multi-component gauge-field theory before landing.

- **`ComponentEulerLagrange` varies a DENSITY, not a Lagrangian** (#394). A non-constant
  `√|g|` must be passed via the `"Measure"` option. Omitting it drops Christoffel
  first-derivative terms for spatially curvilinear metrics, and *sign-flips* them for
  time-dependent ones (de Sitter came out `+2H` anti-damping instead of `−H`).
- **Build the measure as `PowerExpand[Sqrt[metricDetSign * Det[g]]]` — never `Abs`.**
  `Abs` is not symbolically differentiable: `D[Sqrt[Abs[u]], x]` leaves an unevaluated
  `Derivative[1][Abs][u]` in the coefficient, which operator identification cannot match,
  so the term is generated and then **silently dropped on export**. Symptom: the fix
  appears to do nothing at all.
- **LHS normalization must dispatch on `NumericQ[lhsCoeff]`, never `Abs[lhsCoeff] =!= 1`**
  (#381). The latter sends both `+1` and `−1` down a "no normalization needed" path, so
  the `lhsCoeff = -1` that xAct emits for a temporal component under (−,+,+,+) is never
  divided through — yielding a temporal-only tachyon (`∂²ₜa₀ = −∇²a₀`).
- **Do NOT "fix" these two asymmetries — both are correct physics:**
  - `canonical.hamiltonian_terms` gives the temporal photon `a_0 = -0.5` against
    `a_3 = +0.5`. That is the negative-norm temporal photon of Lorenz-gauge-fixed
    Maxwell. The negative norm belongs in the *energy*; the *equation* is `□A_μ = 0` for
    every μ.
  - `h_5` carries `laplacian_x = -kappa^(-2)` matched by
    `lhs.kinetic_coefficient_symbolic = -kappa^(-2)`, i.e. correctly normalized.
- **When comparing coefficients across components, compare EFFECTIVE signs** —
  `sign(coefficient) × sign(kinetic_coefficient)`, summed over *all* matching terms.
  Raw comparison misreads EH-class theories (two self-laplacian terms per component) and
  reports false differences when a re-derivation merely introduces an explicit kinetic
  coefficient. Torsion-family sign non-uniformity is **physical** (rank-3 irreducible
  components carry different normalizations) — do not flag it.

## Common Errors
- All-zero components → field strength not expanded before decomposition
- Cross terms missing → other fields not in additionalFields
- Epsilon not evaluating → chart name mismatch or mixed index signs
- Package function unevaluated → missing `::usage` declaration in public section
- Christoffel first-derivative term missing on a curved metric → `"Measure"` not passed,
  or built with `Abs` (see Equation Generation above)

## Perturbative Reduction
- **`order_in_eps`** tagging lives in `ExportJSON.wl` — `Max[Total[Exponent]]` over `small_parameters`. Always `∈ {0, 1}` by architecture (quadratic Lagrangians × linear couplings). `order=2` is permanently gated by `NotImplementedError` (#273).
- **Power-of-contraction rewrite** (#271): user Lagrangians containing `Power[tensor_contraction, n]` (e.g. `(F·F)^2` for Euler–Heisenberg) must be wrapped in `Hold[]` at assignment, then ReplaceRepeated `Power[X, n] → Scalar[X]^n` applied before `ReleaseHold`. Mathematica's built-in `Power[Times[a,b,...], n] → a^n·b^n·...` fires at parse time and produces repeated abstract indices → xAct `Validate::repeated + Throw[Null]`. See `_derive.py::_wls_lagrangian` and `docs/tex/perturbative_reduction.tex` §Power-of-Contraction (sec:pr-eh-cd).
- **`RenameDummies` is deterministic** — never use `Product[RenameDummies[X], n]` to "freshen" copies; they come out byte-identical and reproduce the clash. Use `Scalar[X]^n` instead (xAct treats it as opaque).
- **CD ComponentValue precompute gate** (#271): `_wls_precompute_cd_component_values` must run whenever a single-dyn-field Lagrangian has derivative-only dependence (matter-only Maxwell-type). Gate: runs if `len(dyn_fields) ≥ 2` OR any `derived_field` with `"CD["` in its definition is used in `L` OR any dyn-field name appears as `CD[...][name[...]]` in the user expression. The legacy `len < 2` skip was a performance heuristic that silently broke EH-style theories.
- **`$CDShorthandReverseRules`** is generated at `_derive.py:1407` but only applied as a FreeQ-gated safety net before Component-E-L field detection. Never apply it unconditionally: a naïve `lagComp //. rules` regresses gertsenshtein by ~11s. The FreeQ guard keeps it O(LeafCount[lagComp]).
- **Validity monitor** (`tidal/solver/perturbative_driver.py`): checks ε·|λ|·t (secular error) and max Re(λ)·t (base-theory stability). Missing small parameter in runtime `params` dict emits a warning (#284–#285).
