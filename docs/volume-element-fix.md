# Volume Element Fix: √|g| in Lagrangian Linearization

## Problem (2026-03-04)

Massive gravity (Fierz-Pauli) simulation was unstable — all fields grew exponentially. Root cause: the linearization pipeline computed `Perturbation[L, 2]` instead of `Perturbation[√|g|·L, 2]`.

The action integral is `S = ∫ √|g| L d^n x`. The second-order perturbation of the integrand is `δ²(√|g| L)`, not just `δ²(L)`. The missing volume element contribution means VarD gives `R^(1)_{ab}` (linearized Ricci) instead of the full Einstein tensor `G^(1)_{ab} = R^(1)_{ab} - ½η_{ab}R^(1)`.

**Why gravitational_waves was unaffected**: TT gauge imposes `h = 0` (traceless), so the missing trace term `h·R^(1)` vanishes.

**Only 3 theories use linearization path**: gravitational_waves, gravitational_waves_1d, massive_gravity.

## Solution: Native xPert Detg[] Perturbation

### Research Findings

- xPert natively supports `Perturbation[Sqrt[-Detg[]], n]` — the metric determinant symbol created by `DefMetric`
- xPert paper (arXiv:0807.0824) and example notebook `Lagrangian-variation-xPert-VarD.nb` confirm the canonical approach:

  ```wolfram
  LGR = Sqrt[-Detg[]] * RicciScalarCD[];
  l2Raw = Perturbation[LGR, 2] // ExpandPerturbation;
  eom = 1/Sqrt[-Detg[]] * VarD[h[LI[1], a, b], CD] @ l2Raw;
  ```

- `∇_i(√|g|) = 0` (covariantly constant for Levi-Civita connection), so `√|g₀|` passes through VarD as a constant
- xAct naming convention: `DefMetric[-1, mgEta[-a,-b], mgCD]` creates `DetmgEta[]`. Pattern: `Det` + metric symbol name.

### Implementation (commit ca256cb)

In `_wls_linearize_from_lagrangian()` in `tidal/cli/_derive.py`:

1. **Multiply by √(-g)**: `lDensity = Sqrt[-DetmgEta[]] * lOriginal`
2. **Perturb density**: `l2Raw = Perturbation[lDensity, 2]; l2Raw = ExpandPerturbation[l2Raw]`
3. **Divide out background**: `l2Raw = l2Raw / Sqrt[-DetmgEta[]]`
4. **Evaluate determinant**: `l2Raw = l2Raw /. DetmgEta[] -> Det[mgMetricMatrix]; l2Raw = Simplify[l2Raw]`

The `det_sym = f"Det{ctx.metric}"` gives the correct symbol for any prefix.

### Why This Works for All Backgrounds

- **Minkowski**: `Det[MetricMatrix] = -1`, `Sqrt[-(-1)] = 1` — factor vanishes
- **Curved** (polar, spherical, etc.): `√|g₀|` factor correctly captured by xPert and divided out
- **TT gauge**: `h = 0` → determinant perturbation vanishes → GW results unchanged
- **Non-curvature Lagrangians**: If `L₀ = 0` (all terms O(h)), correction vanishes naturally

### Previous Approaches (Superseded)

1. **Handmade formula** (commit 4d7c317): `δ²L + h·δL` where `h = η^{ab}h_{ab}`. Only correct for flat Minkowski with `√|g₀| = 1` and `L₀ = 0`. Replaced by native approach.
2. **Phantom d²_t detection** (commits a8ccf29, 45fb7df): Elaborate code in ComponentDecompose.wl to cancel "phantom" d²_t terms. Wrong hypothesis — the terms were genuinely non-zero because of the missing volume element.

### Related Files Changed

- `tidal/cli/_derive.py` — Core fix: native Detg[] perturbation
- `tidal/wolfram/ComponentDecompose.wl` — Reverted phantom d²_t detection to `Expand[Simplify[componentEq]]` (commit 4d7c317)
- `tidal/wolfram/ExportJSON.wl` — Cleaned up DIAG prints, kept equation reassignment logic (commit 4d7c317)

### Key References

- xPert paper: arXiv:0807.0824 (Brizuela, Martin-Garcia, Mena Marugan 2009)
- xPert example: `Lagrangian-variation-xPert-VarD.nb` (xAct-contrib/examples on GitHub)

- xPert docs: `xact.es/Documentation/HTML/HTMLLinks/xPertDoc.nb_5.html`
- VarL (xTras): `xact.es/xTras/documentation/ref/VarL.html` — alternative approach (not used)
