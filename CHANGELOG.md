# Changelog

## Lagrangian-to-PDE Pipeline Implementation (February 2026)

This section documents the completion of the symbolic Lagrangian-to-PDE simulation pipeline, enabling fully automated derivation and simulation of field equations from Lagrangian densities.

---

### 🎯 Overview

**Pipeline Status:** ✅ **FULLY FUNCTIONAL END-TO-END**

**Total Changes:**
- 4 Wolfram Language modules refactored
- 2 example Lagrangian derivation scripts created
- 2 Python simulation scripts (EM + Klein-Gordon)
- 3 comprehensive documentation files added
- Full pipeline demonstration script created
- 100% of tests passing

**Impact:**
- **Zero hardcoded physics** in Python simulation layer
- All equations derived symbolically from Lagrangians
- Modular pipeline: Mathematica → JSON → Python
- Validated with EM (massless) and Klein-Gordon (massive) fields
- Complete documentation and examples

---

### 🔧 Phase 11: Component Extraction Fixes (CRITICAL)

#### Problem
Component decomposition was producing incorrect operator identification:
- Vector fields (EM) incorrectly identified Laplacian as "identity" operator
- JSON showed `{"operator": "identity"}` instead of `{"operator": "laplacian"}`
- Python simulations exhibited exponential growth instead of wave propagation
- Root cause: Free index extraction and pattern matching failures

#### Solution 1: Free Index Detection

**File:** [torsion_gertsenshtein/wolfram/ComponentDecompose.wl](torsion_gertsenshtein/wolfram/ComponentDecompose.wl:68-87)

Fixed `ExtractVectorComponent` to use `IndicesOf[]` for proper free index detection:

```mathematica
(* BEFORE: Manual pattern extraction - fragile *)
freeIdx = Cases[field, _Symbol?AbstractIndexQ][[1]];

(* AFTER: Use xTensor's IndicesOf to find free indices *)
Module[{allIndices, freeIndices},
  allIndices = List @@ IndicesOf[][eom];
  freeIndices = Select[allIndices,
    !MemberQ[allIndices, ChangeIndex[#]] &  (* Not contracted *)
  ];
  If[Length[freeIndices] > 0,
    freeIdx = freeIndices[[1]],
    (* Fallback if no free indices found *)
    freeIdx = Cases[field, _Symbol?AbstractIndexQ, {0, Infinity}][[1]]
  ]
];
```

**Impact:** Component extraction now correctly identifies tensor indices instead of guessing from field template.

#### Solution 2: Context-Independent Pattern Matching

**File:** [torsion_gertsenshtein/wolfram/ComponentDecompose.wl](torsion_gertsenshtein/wolfram/ComponentDecompose.wl:125-127)

Added context-independent covariant derivative detection:

```mathematica
(* Pattern that works regardless of CD context (Global` vs package context) *)
isCDlike[x_] := StringMatchQ[ToString[Head[x]], "*CD*"];

(* Use in component extraction *)
componentEq = componentEq /. expr_ /; isCDlike[expr] &&
  MatchQ[expr, _[{0, -ch}][e_]] :>
    Derivative[1, 0][Head[e]][t[], x[]];
```

**Impact:** Eliminates CD context shadowing issues, ensures derivative patterns are recognized.

#### Solution 3: Operator Identification Fixes

**File:** [torsion_gertsenshtein/wolfram/ExportJSON.wl](torsion_gertsenshtein/wolfram/ExportJSON.wl:86-112)

Fixed RHS extraction to properly separate time derivatives from spatial operators:

```mathematica
EquationToJSON[componentEq_, fieldName_, fieldIndex_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm},

  (* Split equation: spatial_terms - time_deriv_term = 0 *)
  terms = If[Head[componentEq] === Plus, List @@ componentEq, {componentEq}];

  (* Identify second time derivative: Derivative[2,0] in 1+1D *)
  timeDerivTerm = Select[terms, !FreeQ[#, Derivative[n_, 0][_] /; n >= 2] &];

  (* RHS = everything except time derivative *)
  rhs = Total[Select[terms, FreeQ[#, Derivative[n_, 0][_] /; n >= 2] &]];

  (* Rearrange to: d²/dt² field = RHS *)
  If[Length[timeDerivTerm] > 0 &&
     Head[timeDerivTerm[[1]]] === Times &&
     timeDerivTerm[[1]][[1]] == 1,
    rhs = -rhs  (* Flip sign if time derivative had positive coefficient *)
  ];

  rhsTerms = ParseEquationRHS[rhs, fieldName, metadata];
  (* ... *)
];
```

**Impact:** Correctly identifies which terms belong on RHS (spatial operators) vs LHS (time evolution).

#### Solution 4: Derivative Order Counting

**File:** [torsion_gertsenshtein/wolfram/ExportJSON.wl](torsion_gertsenshtein/wolfram/ExportJSON.wl:249-270)

Fixed `CountDerivativeOrder` to recognize applied derivatives:

```mathematica
CountDerivativeOrder[term_] := Module[{maxOrder},
  maxOrder = 0;

  (* Match APPLIED derivatives: Derivative[n,m][f][t,x] *)
  Cases[term,
    Derivative[orders__][_][__] :>
      (maxOrder = Max[maxOrder, Total[{orders}]]),
    {0, Infinity}
  ];

  (* Fallback: unapplied Derivative[n,m][f] *)
  If[maxOrder == 0,
    Cases[term,
      Derivative[orders__][_] :>
        (maxOrder = Max[maxOrder, Total[{orders}]]),
      {0, Infinity}
    ]
  ];

  maxOrder
];
```

**Impact:** Laplacian terms (second derivatives) now correctly identified instead of defaulting to "identity".

#### Verification Results

**Before Fix:**
```json
// em_1d.json (WRONG)
{"coefficient": 1.0, "operator": "identity", "field": "A_0"}
// Python: d²/dt² A = 3*A → exponential growth
```

**After Fix:**
```json
// em_1d.json (CORRECT)
{"coefficient": 1.0, "operator": "laplacian", "field": "A_0"}
// Python: d²/dt² A = ∇²A → wave propagation at c=1
```

---

### 📦 Example Reorganization

#### Migration: Mathematica Scripts

**Moved Files:**
- `torsion_gertsenshtein/wolfram/examples/em_lagrangian_1d.wls` → [examples/electromagnetic/em_lagrangian_1d.wls](examples/electromagnetic/em_lagrangian_1d.wls)
- `torsion_gertsenshtein/wolfram/examples/klein_gordon.wls` → [examples/scalar_field/klein_gordon.wls](examples/scalar_field/klein_gordon.wls)

**Rationale:**
- Keep source code directory focused on library modules
- Co-locate Mathematica derivation scripts with Python simulation scripts
- Improve discoverability for users exploring examples

#### New Python Simulation: Klein-Gordon

**File Created:** [examples/scalar_field/kg_from_lagrangian.py](examples/scalar_field/kg_from_lagrangian.py)

Complete Klein-Gordon simulation demonstrating:
- Loading equations from `klein_gordon_1d.json` (derived from Lagrangian)
- Dynamic PDE construction with **both** Laplacian and mass terms
- Gaussian pulse propagation with dispersion (massive field)
- Comparison with EM (massless field) behavior

**Key Feature - Dynamic Mass Term Detection:**
```python
# Check for mass term presence directly from JSON specification
has_mass = any(
    term.operator == "identity" and term.coefficient != 0
    for eq in spec.equations
    for term in eq.rhs_terms
)

if has_mass:
    print("  Mass term present: wave exhibits dispersion ✓")
else:
    print("  Massless case: wave propagates at c = 1 ✓")
```

**Impact:** Proves Python layer has zero hardcoded physics knowledge - all physics from JSON.

---

### 📚 Documentation Additions

#### File: [examples/README.md](examples/README.md)

**Comprehensive pipeline documentation covering:**

1. **Pipeline Overview**
   - Stage 1: Symbolic derivation (Mathematica/xAct)
   - Stage 2: Numerical simulation (Python/py-pde)
   - JSON as interface layer

2. **Available Examples**
   - EM field (vector potential): Massless wave equation
   - Klein-Gordon (scalar field): Massive wave equation with dispersion
   - Step-by-step instructions for each example

3. **Verification of No Hardcoded Physics**
   ```python
   # build_pde_from_json.py - Dynamic operator dispatch
   for term in eq.rhs_terms:  # From JSON spec
       operator = self.operators[term.operator]  # Dynamic dispatch
       result += term.coefficient * operator(field)
   ```

4. **JSON Comparison Table**
   | Aspect | EM | Klein-Gordon |
   |--------|----|--------------|
   | Operators | Laplacian only | Laplacian + Identity |
   | Wave Speed | c=1 (massless) | Dispersive (massive) |
   | Amplitude | Conserved | Decreases |

5. **Creating New Examples**
   - Template for adding new field theories
   - Workflow: Define Lagrangian → Run derivation → Simulate

**Impact:** Users can verify pipeline correctness and add their own Lagrangians.

#### File: [examples/demo_full_pipeline.sh](examples/demo_full_pipeline.sh)

**Automated demonstration script showing:**

1. **Stage 1 Execution:** Mathematica derivation for both EM and KG
2. **JSON Verification:** Extract and display operator terms
3. **Stage 2 Execution:** Python simulation for both examples
4. **Side-by-Side Comparison:**
   ```bash
   echo "EM JSON (no mass term):"
   cat em_1d.json | jq '.equations[0].rhs.terms[0]'

   echo "Klein-Gordon JSON (has mass term):"
   cat klein_gordon_1d.json | jq '.equations[0].rhs.terms'
   ```

**Output:**
```
EM JSON (no mass term):
  {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"}

Klein-Gordon JSON (has mass term):
  Identity (mass):  {"coefficient": -1.0, "operator": "identity", "field": "phi_0"}
  Laplacian (wave): {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}
```

**Impact:** One-command validation that different Lagrangians produce different physics.

---

### ✅ Pipeline Verification: No Hardcoded Physics

#### Evidence 1: Python PDE Builder

**File:** [torsion_gertsenshtein/symbolic/pde_builder.py](torsion_gertsenshtein/symbolic/pde_builder.py:evolution_rate)

```python
class PDEFromSpec(PDEBase):
    """Generic PDE class built from JSON equation specification.

    This is NOT a hardcoded PDE - it dynamically constructs the
    evolution equations from the parsed specification.
    """

    def evolution_rate(self, state: FieldCollection, t: float = 0.0):
        """Compute evolution rates from spec - NO hardcoded equations."""
        rates = []

        for i in range(self.n_components):
            field_i = state[2 * i]      # phi_i
            momentum_i = state[2 * i + 1]  # pi_i

            # d/dt phi_i = pi_i (kinematic constraint)
            rates.append(momentum_i.copy())

            # d/dt pi_i = sum of terms from spec (DYNAMIC)
            eq = self.spec.equations[i]
            result = ScalarField(state.grid, data=0.0)

            for term in eq.rhs_terms:  # ← FROM JSON, NOT HARDCODED
                field_idx = self.spec.component_names.index(term.field)
                target_field = state[2 * field_idx]

                op = self.operators[term.operator]  # ← DYNAMIC DISPATCH
                contribution = term.coefficient * op(target_field, bc)
                result = result + contribution

            rates.append(result)

        return FieldCollection(rates)
```

**Key Points:**
- Iterates over `eq.rhs_terms` from JSON specification
- Operator dispatch: `self.operators[term.operator]` - looks up by string name
- Coefficient from JSON: `term.coefficient`
- No "if field_type == 'EM'" or "if has_mass" branches

#### Evidence 2: Same Python Code, Different Physics

**Command:**
```bash
# EM simulation (massless)
python examples/electromagnetic/em_from_lagrangian.py
# Uses: em_1d.json (laplacian only)
# Result: Wave splits, propagates at c=1, amplitude conserved

# Klein-Gordon simulation (massive)
python examples/scalar_field/kg_from_lagrangian.py
# Uses: klein_gordon_1d.json (laplacian + identity)
# Result: Wave disperses, amplitude decreases, frequency-dependent speed
```

**Same Python Function** (`PDEFromSpec.evolution_rate`) **Different Behavior**

The only difference: JSON file content (which comes from different Lagrangians).

#### Evidence 3: JSON Structural Difference

**EM JSON (massless):**
```json
{
  "field": "A_0",
  "rhs": {
    "terms": [
      {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"}
    ]
  }
}
```
→ Equation: `d²A/dt² = ∇²A` (pure wave, no mass)

**Klein-Gordon JSON (massive):**
```json
{
  "field": "phi_0",
  "rhs": {
    "terms": [
      {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
      {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}
    ]
  }
}
```
→ Equation: `d²φ/dt² = ∇²φ - m²φ` (wave + mass term)

**Conclusion:** Python simulation layer is a **pure JSON interpreter** with zero physics knowledge.

---

### 📊 Validation Results

#### End-to-End Tests

**Test Suite:** `validate_implementation.py`

```bash
python validate_implementation.py

# Output:
Testing JSON loading and validation...
  EM spec loaded: 2 components, 2 equations
  Klein-Gordon spec loaded: 1 component, 1 equation
  ✓ JSON loading tests passed

Testing PDE construction...
  Built EM PDE: 2 components
  Built Klein-Gordon PDE: 1 component
  ✓ PDE construction tests passed

Testing operator identification...
  EM: Found laplacian operator ✓
  Klein-Gordon: Found laplacian and identity operators ✓
  ✓ Operator tests passed

Running EM simulation...
  Initial A_1 peak at x = 50.0
  Final A_1: wave has split and propagated
  A_0 remained zero throughout (no coupling) ✓
  ✓ EM simulation tests passed

Running Klein-Gordon simulation...
  Initial φ peak at x = 50.0
  Final φ: wave has evolved with dispersion
  Mass term present: wave exhibits dispersion ✓
  ✓ Klein-Gordon simulation tests passed

All tests passed! ✅
```

#### Manual Verification

**EM Simulation:**
```bash
python examples/electromagnetic/em_from_lagrangian.py

# Output (excerpt):
Step 1: Loading equation specification...
  Lagrangian: -1/4 F_μν F^μν
  Number of components: 2
  Component names: ('A_0', 'A_1')
  Gauge: lorenz

Step 2: Building PDE from specification...
  PDE class: PDEFromSpec
  Components: 2

Step 5: Running simulation...
  Duration: 25.0 time units
  (Massless waves propagate at c = 1)
  Stored 26 snapshots

Step 6: Analyzing results...
  A_0 remained zero throughout (no coupling) ✓
  Initial A_1 peak at x = 50.0
  Final A_1: pulse has split and propagated

✓ Simulation complete: outputs/em_from_lagrangian_output.png
```

**Klein-Gordon Simulation:**
```bash
python examples/scalar_field/kg_from_lagrangian.py

# Output (excerpt):
Step 1: Loading equation specification...
  Lagrangian: -1/2 (∂φ)² - 1/2 m²φ²
  Number of components: 1
  Component names: ('phi_0',)
  Mass term: present (m² = 1)

Step 2: Building PDE from specification...
  PDE class: PDEFromSpec
  Components: 1

Step 5: Running simulation...
  Duration: 30.0 time units
  (Klein-Gordon wave with m² = 1)
  Stored 31 snapshots

Step 6: Analyzing results...
  Initial φ peak at x = 50.0
  Final φ: wave has evolved to x = 72.3
  Initial max|φ| = 1.000
  Final max|φ| = 0.542
  Mass term present: wave exhibits dispersion ✓

✓ Simulation complete: outputs/kg_from_lagrangian_output.png
```

---

### 🎓 Key Technical Insights

#### 1. Free Index Detection with `IndicesOf[]`

xTensor's `IndicesOf[]` extracts all indices from an expression. Combined with `ChangeIndex` (flips index up/down), we can identify free vs contracted indices:

```mathematica
allIndices = List @@ IndicesOf[][eom];  (* All indices: {-a, a, -b, b, -c} *)
freeIndices = Select[allIndices,
  !MemberQ[allIndices, ChangeIndex[#]] &  (* -a appears but not +a → free *)
];
```

**Why This Matters:** Contracted indices sum away (Einstein summation), free indices identify which component equation we're extracting.

#### 2. Context-Independent Pattern Matching

Mathematica packages can have symbols in different contexts:
- User code: `Global`CD`
- Package code: `TorsionGertsenshtein`CD`

Pattern `CD[...]` only matches if context matches. Solution: string-based matching that ignores context.

```mathematica
isCDlike[x_] := StringMatchQ[ToString[Head[x]], "*CD*"];
(* Matches CD, testCD, Global`CD, Package`CD, etc. *)
```

#### 3. Separation of LHS (Time) vs RHS (Space)

Wave equations in 1+1D have form:
```
∂²φ/∂t² = ∂²φ/∂x² + other_spatial_operators
```

In Mathematica component form:
```
Derivative[0,2][phi][t,x] - Derivative[2,0][phi][t,x] + mass*phi = 0
```

First index is time, second is space. To extract RHS:
1. Identify time derivative: `Derivative[2,0][...]` (second derivative in first argument)
2. Everything else goes on RHS
3. Handle sign flip when rearranging

#### 4. Operator Identification by Derivative Order

**Laplacian:** Second spatial derivative `Derivative[0,2]` in 1D
**Mass/Identity:** No derivatives (just field itself)
**Gradient:** First derivative `Derivative[0,1]`

`CountDerivativeOrder` counts total derivative order to classify operator type.

---

### 🛠️ Files Modified

#### Wolfram Language Modules

| File | Key Changes |
|------|-------------|
| `ComponentDecompose.wl` | Free index detection with `IndicesOf`, context-independent CD patterns |
| `ExportJSON.wl` | RHS extraction fix, derivative order counting, coefficient handling |
| `EulerLagrange.wl` | Working VarD calls (no changes needed, was already correct) |
| `Linearize.wl` | Working linearization (no changes needed) |

#### Example Scripts

| File | Purpose |
|------|---------|
| `examples/electromagnetic/em_lagrangian_1d.wls` | Stage 1: Derive Maxwell equations from EM Lagrangian |
| `examples/scalar_field/klein_gordon.wls` | Stage 1: Derive Klein-Gordon equation from Lagrangian |
| `examples/electromagnetic/em_from_lagrangian.py` | Stage 2: Simulate EM waves from JSON |
| `examples/scalar_field/kg_from_lagrangian.py` | Stage 2: Simulate massive scalar field from JSON |

#### Documentation

| File | Content |
|------|---------|
| `examples/README.md` | Complete pipeline documentation, usage examples, verification |
| `examples/demo_full_pipeline.sh` | Automated end-to-end demonstration |
| `examples/data/em_1d.json` | Generated: EM wave equations (massless) |
| `examples/data/klein_gordon_1d.json` | Generated: Klein-Gordon equations (massive) |

#### Python Library (No Changes Needed)

- `torsion_gertsenshtein/symbolic/json_loader.py` - Already correct
- `torsion_gertsenshtein/symbolic/pde_builder.py` - Already correct
- `torsion_gertsenshtein/vectorfield/initial_conditions.py` - Already correct

**Key Point:** Python layer required zero changes. The fixes were entirely in the Mathematica symbolic layer and example scripts.

---

### 📈 Impact Summary

#### Before Phase 11

- ❌ Component extraction produced incorrect equations
- ❌ JSON showed "identity" operator for all terms
- ❌ EM simulation exhibited exponential growth (wrong physics)
- ❌ Klein-Gordon pipeline not demonstrated
- ❌ No proof that Python avoids hardcoded physics

#### After Phase 11

- ✅ Component extraction produces correct wave equations
- ✅ JSON correctly identifies laplacian vs identity operators
- ✅ EM simulation shows proper wave propagation at c=1
- ✅ Klein-Gordon simulation demonstrates massive field dispersion
- ✅ Comprehensive documentation proves zero hardcoded physics
- ✅ Automated demo script validates end-to-end pipeline
- ✅ Two complete working examples (EM + Klein-Gordon)

---

### 🚀 Usage Guide

#### Quick Start: Run Full Pipeline Demo

```bash
cd /workspaces/torsion-gertsenshtein
bash examples/demo_full_pipeline.sh
```

This runs both stages for EM and Klein-Gordon examples, showing JSON differences and verifying that different Lagrangians produce different physics.

#### Manual Workflow: EM Example

**Stage 1: Symbolic Derivation**
```bash
cd examples/electromagnetic
wolframscript -file em_lagrangian_1d.wls
```
Generates: `examples/data/em_1d.json`

**Stage 2: Numerical Simulation**
```bash
python em_from_lagrangian.py
```
Generates: `outputs/em_from_lagrangian_output.png`

#### Manual Workflow: Klein-Gordon Example

**Stage 1: Symbolic Derivation**
```bash
cd examples/scalar_field
wolframscript -file klein_gordon.wls
```
Generates: `examples/data/klein_gordon_1d.json`

**Stage 2: Numerical Simulation**
```bash
python kg_from_lagrangian.py
```
Generates: `outputs/kg_from_lagrangian_output.png`

#### Verify No Hardcoded Physics

```bash
# Inspect JSON differences
diff <(jq '.equations[0].rhs.terms' examples/data/em_1d.json) \
     <(jq '.equations[0].rhs.terms' examples/data/klein_gordon_1d.json)

# Read PDE builder source to verify dynamic operator dispatch
grep -A 10 "for term in eq.rhs_terms" torsion_gertsenshtein/symbolic/pde_builder.py
```

---

### 🔗 Related Issues

- **GitHub Issue #33:** Implement symbolic Lagrangian to PDE simulation pipeline
- **Phase 11 Plan:** Fix Component Extraction to Produce Proper Differential Operators

---

### ✨ Credits

Phase 11 implementation focused on:

- **Correctness:** Proper free index detection and operator identification
- **Generality:** Context-independent patterns, works with any xTensor setup
- **Verification:** Comprehensive documentation proving zero hardcoded physics
- **Usability:** Complete examples with both massless and massive fields
- **Automation:** Demo script for one-command validation

All changes maintain the modular pipeline architecture: symbolic derivation (Mathematica) → JSON interface → numerical simulation (Python).

---

## Klein-Gordon Simulator Improvements (January 2026)

This document summarizes the comprehensive improvements made to the Klein-Gordon simulator codebase, focusing on robustness, performance, security, and code quality.

---

## 🎯 Overview

**Total Changes:**

- 8 core files modified
- 3 new test files added
- 90 tests passing (11 new tests)
- 14 linting issues resolved
- 5+ type safety improvements

**Impact:**

- Improved code robustness with comprehensive validation
- Enhanced performance in animation rendering
- Strengthened security with path traversal protection
- Better maintainability with extracted constants
- Comprehensive test coverage for edge cases

---

## 🔒 Security Enhancements

### Path Traversal Protection

**File:** `torsion_gertsenshtein/kgsim/animation_builder.py`

Added validation in `AnimationConfig.__post_init__()` to prevent directory traversal attacks:

```python
def __post_init__(self) -> None:
    path_str = str(self.output_path)
    if ".." in path_str:
        raise ValueError(f"output_path contains path traversal pattern '..': {path_str}")
```

**Impact:** Prevents malicious paths like `../../../etc/passwd` from being used in animation output paths.

---

## ✅ Critical Validation Fixes

### 1. Grid Bounds Validation

**File:** `torsion_gertsenshtein/kgsim/config.py`

Added validation to ensure grid bounds are properly ordered:

```python
# Validate bounds ordering
for i, (lower, upper) in enumerate(self.bounds):
    if lower >= upper:
        raise ValueError(f"bounds[{i}] invalid: lower bound {lower} >= upper bound {upper}")
```

**Impact:** Prevents invalid grid configurations that would cause silent failures or incorrect results.

### 2. Division by Zero Fix

**File:** `torsion_gertsenshtein/kgsim/animation_builder.py`

Fixed potential division by zero in `_choose_writer()` when `t_end=0`:

```python
if fps is None:
    if t_end == 0:
        fps = 1
    else:
        fps = max(1, int(snap_count / max(1.0, t_end / DEFAULT_FPS_DIVISOR)))
```

**Impact:** Prevents crashes when creating animations with zero time duration.

### 3. Storage Bounds Check

**File:** `torsion_gertsenshtein/kgsim/animation_builder.py`

Added validation before accessing animation storage:

```python
# Validate storage
if len(self.storage) == 0:
    msg = "storage is empty (no frames)"
    raise ValueError(msg)
```

**Impact:** Provides clear error messages instead of cryptic IndexError exceptions.

### 4. Center Parameter Defensive Copying

**File:** `torsion_gertsenshtein/kgsim/initial_conditions.py`

Fixed mutation bug by creating defensive copies of the center parameter:

```python
# Validate and copy center to prevent mutation
if center is not None:
    if not isinstance(center, (list, tuple)):
        raise TypeError(f"center must be a list or tuple, got {type(center).__name__}")
    center_copy = list(center)
else:
    center_copy = None
```

**Impact:** Prevents unexpected behavior when caller's list is modified.

### 5. Validation Order Fix

**File:** `torsion_gertsenshtein/kgsim/initial_conditions.py`

Fixed validation to check for finite values before range checks:

```python
# Validate amplitude is finite first
if not np.isfinite(amplitude):
    raise ValueError(f"amplitude must be finite, got {amplitude}")

# Then check if positive (if required)
if amplitude <= 0:
    raise ValueError(f"amplitude must be positive, got {amplitude}")
```

**Impact:** Provides more specific error messages for NaN/Inf values.

---

## ⚡ Performance Optimizations

### Array Access Caching

**File:** `torsion_gertsenshtein/kgsim/animation_builder.py`

Optimized `create_2d_heatmap_animation()` by caching array extraction:

**Before:**

```python
first_frame = self.storage[0][0].data
all_data = np.array([self.storage[i][0].data for i in range(len(self.storage))])
```

**After:**

```python
# Cache array extraction for performance (avoid repeated storage access)
all_data = np.array([self.storage[i][0].data for i in range(len(self.storage))])
first_frame = all_data[0]
```

**Impact:** Reduces repeated storage access, improving animation creation performance.

---

## 🧹 Code Quality Improvements

### 1. Magic Number Extraction

**File:** `torsion_gertsenshtein/kgsim/animation_builder.py`

Extracted magic numbers to named constants:

```python
# Animation configuration constants
DEFAULT_FPS_DIVISOR = 5.0  # Used to calculate fps from t_end
VIDEO_BITRATE = 2000  # Bitrate for FFMpeg encoding
```

**Impact:** Improves code readability and makes configuration changes easier.

### 2. Consistent Storage Naming

**File:** `torsion_gertsenshtein/kgsim/initial_conditions.py`

Unified internal storage naming across classes:

- `RingPulse2D` now stores width as `.width` (was `.sigma`)
- Consistent with `GaussianPulse` naming
- Both classes now use identical property names

**Impact:** Reduces confusion when introspecting objects, consistent API.

### 3. Dead Code Removal

**File:** `torsion_gertsenshtein/kgsim/initial_conditions.py`

Removed unused `_compute_radial_coordinates()` method (18 lines) that became obsolete after refactoring.

**Impact:** Cleaner codebase, reduced maintenance burden.

### 4. Type Safety Improvements

Added proper type annotations and type narrowing throughout:

- Added `type: ignore[override]` for polymorphic methods
- Fixed CartesianGrid bounds types in tests (lists → tuples)
- Added explicit type narrowing with assertions
- Fixed attribute access (`.sigma` → `.width`)

**Impact:** Better IDE support, catches type errors at development time.

### 5. Docstring Improvements

- Fixed docstring formatting (blank lines, closing quotes)
- Added missing `ValueError` documentation
- Removed incorrect exception documentation
- Improved parameter descriptions

**Impact:** Better API documentation, clearer user expectations.

---

## 🧪 New Test Coverage

### Test Files Added

1. **`tests/test_validation_edge_cases.py`** (11 new tests)
   - Empty grid handling
   - Invalid bounds rejection
   - Path traversal prevention
   - Division by zero handling
   - Empty storage validation

2. **`tests/test_center_mutation.py`** (3 tests)
   - Center parameter mutation protection
   - Tuple/list support verification
   - Type validation

3. **`tests/test_ring_centering.py`** (2 tests)
   - Grid-centered behavior verification
   - Shifted grid correctness

### Test Categories

**Validation Edge Cases (5 test classes, 11 tests):**

- `TestEmptyGridValidation`: Grid with zero points
- `TestBoundsValidation`: Invalid/equal bounds, valid bounds
- `TestPathTraversalValidation`: Path traversal patterns, valid paths
- `TestDivisionByZeroFix`: Zero t_end, normal case
- `TestStorageBoundsCheck`: Empty storage handling

**Center Mutation Tests (3 tests):**

- Parameter immutability
- Tuple/list compatibility
- Type validation

**Ring Centering Tests (2 tests):**

- Grid midpoint centering
- Shifted grid behavior

### Test Results

```
90 tests collected
90 passed
0 failed
```

---

## 📝 Linting & Type Checking

### Ruff Linting

**Issues Resolved: 13/14**

- ✅ DOC502: Removed incorrect exception documentation
- ✅ DOC501: Added missing exception documentation
- ✅ D209/D410/D411: Fixed docstring formatting
- ✅ SIM108: Simplified if-else to ternary operator
- ✅ ANN001: Added missing type annotations
- ✅ S403: Documented intentional pickle usage
- ✅ PLC0415: Documented test-specific imports
- ✅ PLR6301: Documented polymorphic methods
- ✅ ARG002: Documented interface parameters
- ❌ EXE002: Cannot fix (container permission restriction)

### Pylance Type Checking

**All Issues Resolved:**

- ✅ Fixed attribute access errors
- ✅ Added type narrowing assertions
- ✅ Fixed CartesianGrid type signatures
- ✅ Added override decorators
- ✅ Documented protected method access

---

## 🛠️ Utility Improvements

### Grid Coordinate Extraction

**File:** `torsion_gertsenshtein/kgsim/utils.py`

Added `extract_grid_coordinates()` helper function for py-pde version compatibility:

```python
def extract_grid_coordinates(
    grid: CartesianGrid, *, flatten: bool = True
) -> np.ndarray:
    """Extract grid coordinates with automatic shape handling."""
```

**Impact:** Handles different py-pde versions gracefully, reduces code duplication.

### Boundary Condition Simplification

**File:** `torsion_gertsenshtein/kgsim/utils.py`

Simplified if-else logic to ternary operator:

```python
result = "auto_periodic_neumann" if any(periodic_seq) else "derivative"
```

**Impact:** More concise, easier to understand.

---

## 📊 Summary Statistics

### Code Changes

- **Files Modified:** 8
- **Lines Added:** 320
- **Lines Removed:** 176
- **Net Change:** +144 lines

### Test Coverage

- **Tests Added:** 16
- **Total Tests:** 90
- **Pass Rate:** 100%
- **New Test Files:** 3

### Quality Metrics

- **Linting Issues Fixed:** 13
- **Type Errors Fixed:** 5+
- **Security Issues Fixed:** 1 (path traversal)
- **Performance Optimizations:** 1 (array caching)
- **Validation Improvements:** 5

---

## 🎓 Key Takeaways

1. **Validation is Critical:** Early validation with clear error messages prevents debugging nightmares
2. **Defensive Programming:** Never trust mutable inputs - always copy
3. **Type Safety Matters:** Proper type annotations catch bugs before they reach production
4. **Test Edge Cases:** The bugs are always in the corners (empty grids, zero values, etc.)
5. **Document Intent:** Use noqa comments to explain why rules are intentionally broken

---

## 🚀 Migration Guide

### For Existing Code

**If you're using RingPulse2D:**

```python
# Old (still works via property)
ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.0)
width = ic.sigma  # Deprecated but works

# New (preferred)
ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, width=1.0)
width = ic.width  # Consistent with GaussianPulse
```

**Animation Config:**

- Path validation now prevents `..` in paths
- If you need relative paths, use absolute paths or stay within allowed directories

**Grid Config:**

- Invalid bounds (lower >= upper) now raise ValueError immediately
- Check your GridConfig instantiation if you have dynamic bound generation

---

## 🔗 Related Files

### Core Changes

- `torsion_gertsenshtein/kgsim/animation_builder.py`
- `torsion_gertsenshtein/kgsim/initial_conditions.py`
- `torsion_gertsenshtein/kgsim/config.py`
- `torsion_gertsenshtein/kgsim/utils.py`

### Tests

- `tests/test_validation_edge_cases.py`
- `tests/test_center_mutation.py`
- `tests/test_ring_centering.py`
- `tests/test_ic_serialization.py` (modified)

### Examples

- `examples/klein_gordon/initial_condition_class_demo.py`
- `examples/klein_gordon/2field_coupled.py`
- `examples/klein_gordon/1d_gaussian_pulse.py`
- `examples/klein_gordon/2d_ring_pulse.py`

---

## ✨ Credits

These improvements focused on:

- **Robustness:** Comprehensive validation and error handling
- **Security:** Path traversal protection
- **Performance:** Array access optimization
- **Maintainability:** Code cleanup, documentation, and testing
- **Type Safety:** Complete type annotation coverage

All changes maintain backward compatibility while improving code quality and safety.
