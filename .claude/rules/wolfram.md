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
