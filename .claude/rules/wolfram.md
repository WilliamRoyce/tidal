---
paths:
  - "tidal/wolfram/**"
  - "tidal/symbolic/_derive.py"
  - "examples/**/theory.toml"
  - "examples/**/*.wls"
---

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

## Common Errors
- All-zero components → field strength not expanded before decomposition
- Cross terms missing → other fields not in additionalFields
- Epsilon not evaluating → chart name mismatch or mixed index signs
- Package function unevaluated → missing `::usage` declaration in public section

## Perturbative Reduction
- **`order_in_eps`** tagging lives in `ExportJSON.wl` — `Max[Total[Exponent]]` over `small_parameters`. Always `∈ {0, 1}` by architecture (quadratic Lagrangians × linear couplings). `order=2` is permanently gated by `NotImplementedError` (#273).
- **Power-of-contraction rewrite** (#271): user Lagrangians containing `Power[tensor_contraction, n]` (e.g. `(F·F)^2` for Euler–Heisenberg) must be wrapped in `Hold[]` at assignment, then ReplaceRepeated `Power[X, n] → Scalar[X]^n` applied before `ReleaseHold`. Mathematica's built-in `Power[Times[a,b,...], n] → a^n·b^n·...` fires at parse time and produces repeated abstract indices → xAct `Validate::repeated + Throw[Null]`. See `_derive.py::_wls_lagrangian` and `docs/tex/perturbative_reduction.tex` §Power-of-Contraction (sec:pr-eh-cd).
- **`RenameDummies` is deterministic** — never use `Product[RenameDummies[X], n]` to "freshen" copies; they come out byte-identical and reproduce the clash. Use `Scalar[X]^n` instead (xAct treats it as opaque).
- **CD ComponentValue precompute gate** (#271): `_wls_precompute_cd_component_values` must run whenever a single-dyn-field Lagrangian has derivative-only dependence (matter-only Maxwell-type). Gate: runs if `len(dyn_fields) ≥ 2` OR any `derived_field` with `"CD["` in its definition is used in `L` OR any dyn-field name appears as `CD[...][name[...]]` in the user expression. The legacy `len < 2` skip was a performance heuristic that silently broke EH-style theories.
- **`$CDShorthandReverseRules`** is generated at `_derive.py:1407` but only applied as a FreeQ-gated safety net before Component-E-L field detection. Never apply it unconditionally: a naïve `lagComp //. rules` regresses gertsenshtein by ~11s. The FreeQ guard keeps it O(LeafCount[lagComp]).
- **Validity monitor** (`tidal/solver/perturbative_driver.py`): checks ε·|λ|·t (secular error) and max Re(λ)·t (base-theory stability). Missing small parameter in runtime `params` dict emits a warning (#284–#285).
