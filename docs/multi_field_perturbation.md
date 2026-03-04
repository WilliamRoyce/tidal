# Multi-Field Perturbation Theory: Architecture Design

## Motivation

The current `[linearization]` pipeline handles a single perturbation field — the metric:

```
g_{μν} = η_{μν} + ε h_{μν}
```

via `SetupMetricPerturbation` + `Perturbation[L, 2]` + `VarD[h]`.

Many physically interesting theories require **multiple perturbation fields** expanded simultaneously around their backgrounds. The Gertsenshtein effect (Einstein-Maxwell) is the immediate use case, but the architecture is general-purpose.

## Use Cases

| Theory | Perturbation Fields | Background |
| ------ | ------------------- | ---------- |
| **Einstein-Maxwell** (Gertsenshtein) | metric h_{μν} + EM vector a_μ | η_{μν} + Ā_μ (uniform B₀) |
| **Einstein-Proca** | metric h_{μν} + massive vector a_μ | η_{μν} |
| **Einstein-Yang-Mills** | metric h_{μν} + gauge field a_μ^I | η_{μν} + Ā_μ^I |
| **Einstein-Cartan** | metric h_{μν} + torsion T^λ_{μν} | η_{μν} + T̄ |
| **Scalar-tensor gravity** | metric h_{μν} + scalar φ | η_{μν} + φ̄ (VEV) |
| **Matter-only** (fixed geometry) | vector a_μ (or scalar, tensor) | η_{μν} (fixed) + Ā_μ |

## TOML Configuration

### Full multi-field perturbation (e.g., Einstein-Maxwell):

```toml
[linearization]
perturbation_field = "h"   # Metric perturbation (existing, unchanged)

[[linearization.matter_perturbations]]
field = "A"                # Must match a [[fields]] entry
perturbation_name = "a"    # Name for the 1st-order perturbation tensor
background = "Abar"        # Resolved via matching [[background_fields]] entry
```

### Matter-only perturbation (fixed metric):

```toml
[linearization]
# perturbation_field omitted → metric is NOT perturbed
# Internally uses: Perturbation[Bg[__], ___] := 0

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"
```

### Multiple matter fields:

```toml
[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"

[[linearization.matter_perturbations]]
field = "phi"
perturbation_name = "dphi"
background = "phi0"
```

### Backward compatibility

When `matter_perturbations` is absent, behavior is identical to the current pipeline. No existing TOML configs are affected.

## Pipeline Architecture

### Derived Fields and Perturbation Order

**Critical design detail**: The Lagrangian uses derived fields (e.g., F), but perturbation acts on fundamental fields (A).

**Pipeline sequence:**

1. **User writes Lagrangian** with derived field F:
   ```
   L = -1/4 F[-a,-b] η[a,c] η[b,d] F[-c,-d]
   ```

2. **MakeRule expansion** (existing `_wls_derived_fields()` → `_wls_lagrangian()`):
   F → `CD[-a][A[-b]] - CD[-b][A[-a]]`. The Lagrangian now has `CD[-a][A[-b]]` terms.

3. **xPert perturbation** of the expanded Lagrangian:
   `Perturbation[√|g| L_expanded, 2]` + `ExpandPerturbation[]`.
   xPert correctly produces:
   - 0th order: `CD[-a][A[-b]]` → background field strength F̄
   - 1st order: `CD[-a][aPert[LI[1], -b]]` → perturbation field strength f

4. **Background resolution** via xCoba during `ToBasis` + `TraceBasisDummy`:
   xCoba evaluates `CD[-a][A[-b]]` using `ComponentValue`-defined background.
   For flat space: `∂_z(A_y) = ∂_z(-B₀z) = -B₀`.

**Key point**: The existing `_validate_lagrangian` check rejects `CD[-a][G[]]` in the user's Lagrangian expression. This does NOT trigger because the user writes `F[-a,-b]` (not `CD[-a][A[-b]]`). The derivatives only appear after MakeRule expansion.

### Perturbation Order Truncation

TIDAL operates exclusively in the linearized regime. After `Perturbation[L, 2]` + `ExpandPerturbation`:

| Term Type | Example | Keep? | Reason |
| --------- | ------- | ----- | ------ |
| Products of 1st-order (same field) | `h[LI[1]] × h[LI[1]]` | **YES** | Quadratic action |
| Products of 1st-order (cross) | `h[LI[1]] × a[LI[1]]` | **YES** | Coupling terms |
| Pure 2nd-order | `h[LI[2]]` or `a[LI[2]]` | **NO** | Linearized: truncate at 1st order |

Dropping `LI[2]` is equivalent to setting field = background + ε·perturbation^(1) with no ε² terms. This is correct for all linearized theories.

## Wolfram Code Generation

### New function: `_wls_matter_perturbation_setup()`

Generates `DefTensorPerturbation` + `ComponentValue` for each matter perturbation.

```mathematica
(* For each [[linearization.matter_perturbations]] entry: *)

(* 1. Define xPert perturbation for field A *)
DefTensorPerturbation[aPert[LI[order], -a], A[-a], M4];

(* 2. Set background values for A (0th order) *)
(* From matching [[background_fields]] entry *)
ComponentValue[A[{0, -Cart}], 0];
ComponentValue[A[{1, -Cart}], 0];
ComponentValue[A[{2, -Cart}], -B0 * z[]];
ComponentValue[A[{3, -Cart}], 0];
```

For tensor matter perturbations (rank-2):
```mathematica
DefTensorPerturbation[TPert[LI[order], -a, -b], T[-a, -b], M4];
```

For scalar matter perturbations:
```mathematica
DefTensorPerturbation[dphiPert[LI[order]], phi[], M4];
```

### Modified: `_wls_linearize_from_lagrangian()`

After `ExpandPerturbation`, drop LI[2] and replace LI[1] for ALL perturbation fields:

```mathematica
(* Existing: metric perturbation *)
SetupMetricPerturbation[g, hpert, epsilon];

(* NEW: each matter perturbation *)
DefTensorPerturbation[aPert[LI[order], -a], A[-a], M4];

(* 2nd-order perturbation of action density *)
l2Raw = Perturbation[Sqrt[-Detg[]] * Lagrangian, 2];
l2Raw = ExpandPerturbation[l2Raw];

(* Truncate ALL fields at 1st order *)
l2Raw = l2Raw /. hpert[LI[2], __] :> 0;
l2Raw = l2Raw /. aPert[LI[2], __] :> 0;

(* Replace xPert notation with declared field tensors *)
l2Raw = l2Raw /. hpert[LI[1], idx__] :> H[idx];
l2Raw = l2Raw /. aPert[LI[1], idx__] :> a[idx];

(* VarD for EACH dynamical field *)
eomH = VarD[H[-a,-b], CD][l2ForVarD];
eomA = VarD[a[-a], CD][l2ForVarD];
```

### Multi-field VarD and Component Decomposition

Each dynamical field gets its own `VarD` call. The results feed into `DecomposeToComponents` with all other dynamical fields listed as `additionalFields`:

```mathematica
(* Decompose graviton equation — a components are additional fields *)
DecomposeToComponents[eomH, H, {a}, ...]

(* Decompose photon equation — H components are additional fields *)
DecomposeToComponents[eomA, a, {H}, ...]
```

The existing `additionalFields` mechanism in `DecomposeToComponents` already handles cross-field terms correctly (verified in coupled_scalars, scalar_vector_coupling examples).

## Background Field Strength: How F̄ = dĀ Gets Evaluated

For the Gertsenshtein effect, the coupling coefficients contain the background field strength F̄. Here is the evaluation chain:

1. **TOML**: `[[background_fields]]` with components `["0", "0", "-B0 * z[]", "0"]`
2. **Wolfram**: `ComponentValue[A[{2, -Cart}], -B0 * z[]]`
3. **After ExpandPerturbation**: 0th-order terms contain `CD[-a][A[-b]]`
4. **ToBasis**: xCoba evaluates `CD[-a][A[-b]]` → `∂_a(Ā_b)` using known component values
5. **Result**: For Ā_y = -B₀z, xCoba computes ∂_z(Ā_y) = -B₀

For uniform B₀, the coupling coefficients are constants (proportional to B₀). For localized B₀(z), they become position-dependent coefficients — already supported by the solver's coordinate-dependent coefficient system (`IsCoordinateDependentCoefficient` → `_mathematica_to_python()` → grid evaluation → L2 cache).

## Files Modified

| File | Change |
| ---- | ------ |
| `tidal/cli/_derive.py` | Extend `_validate_linearization` for `matter_perturbations`; new `_wls_matter_perturbation_setup()` helper; extend `_wls_linearize_from_lagrangian` for multi-field LI[2] drop + LI[1] replacement + multi-field VarD |
| `tidal/cli/_derive.py` | Modify `_wls_metadata_and_export` to handle multi-field EOM |
| `tidal/wolfram/ComponentDecompose.wl` | May need extension for multi-primary-field decomposition (verify `additionalFields` suffices) |
| `tests/test_cli.py` | New tests for TOML validation and WLS generation |

## Validation

1. **Regression**: All 1343+ existing tests pass; all 25 examples derive/simulate correctly
2. **TOML validation**: `matter_perturbations` accepted/rejected correctly
3. **WLS generation**: `--dry-run` produces correct `DefTensorPerturbation` + multi-field `VarD`
4. **End-to-end**: Einstein-Maxwell derives to correct JSON with h-a coupling terms
5. **Physics**: Conversion probability matches analytical formula to within 1%

## Design Principles

1. **General-purpose**: Not Gertsenshtein-specific — any theory with multiple perturbation fields
2. **Lagrangian-first**: All equations derived from the action, never hardcoded
3. **xPert-native**: Uses `DefTensorPerturbation`, not ad-hoc splitting
4. **Backward-compatible**: No changes when `matter_perturbations` is absent
5. **Composable**: Works with existing gauge fixing, plane-wave reduction, background fields
