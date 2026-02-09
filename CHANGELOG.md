# Changelog

## Scalar-Vector Coupling Stress Test (February 2026)

**Status:** ✅ COMPLETE

**Summary:** End-to-end stress test of the pipeline with mixed-rank cross-field coupling: scalar phi + vector A_mu in 2+1D, combining Klein-Gordon, Proca, Chern-Simons, and divergence coupling terms in a single Lagrangian with 4 symbolic constants.

**Files Created:**
- `examples/scalar_vector_coupling/theory.toml` — TOML config with `[[derived_fields]]` for F_ab
- `examples/scalar_vector_coupling/scalar_vector_coupling.wls` — Wolfram derivation script
- `examples/scalar_vector_coupling/scalar_vector_coupling_simulation.py` — Python simulation with 3x2 plot layout
- `examples/scalar_vector_coupling/run.sh` — CLI-equivalent workflow script
- `examples/data/scalar_vector_coupling.json` — JSON equation spec (4 fields, 4x4 matrices)
- `tests/test_scalar_vector_physics.py` — 22 physics tests

**Key Features:**
- Mixed-rank cross-field: scalar phi_0 + vector (A_0, A_1, A_2) with cross-field first_derivative_t and gradient operators
- 4 symbolic constants (phim2, Am2, kCS, gSV) all preserved symbolically for runtime parameter sweeps
- 4x4 mass/coupling matrices auto-computed from identity terms
- A_0 as constraint equation (time_derivative_order=0)
- Demonstrates scalar→vector energy transfer via gradient coupling
- Uses `[[derived_fields]]` TOML feature for field strength tensor F_ab

**Tests:** 743 Python tests passing (22 new physics tests + 147 CLI tests)

---

## CLI (`tg` Command) Implementation (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Unified command-line interface for the Lagrangian-to-PDE pipeline with 5 subcommands, zero new dependencies (stdlib argparse + tomllib), and 147 CLI-specific tests.

**Subcommands:**
- `tg derive theory.toml` — Generate .wls from TOML config, run wolframscript to produce JSON
- `tg simulate spec.json` — Full simulation with smart defaults, plotting, and parameter override
- `tg inspect spec.json` — Display equation system info (fields, operators, parameters)
- `tg list` — Discover available JSON specs in examples/data/
- `tg validate spec.json` — Validate JSON equation specification structure

**Key Features:**
- `theory.toml` configuration with `[[derived_fields]]` for intermediate tensor definitions
- IC presets: `gaussian`, `plane-wave`, `zero`, `formula` via `--ic` flag
- Per-axis boundary conditions via `--bc neumann,periodic`
- `--mode constraint` for single constraint solve (no time evolution)
- `--scheme scipy` (adaptive) or `--scheme runge-kutta` (explicit)
- `--ic-formula` hardened with AST validation
- Version via `importlib.metadata.version()` (single source from pyproject.toml)
- `py.typed` marker included in sdist/wheel

**Architecture:**
- `_WlsContext` dataclass bundles WLS generation state
- `PlotContext` dataclass bundles plot arguments
- Plotting extracted to `_plot.py` (separate from `_simulate.py`)
- `_validate.py` for JSON spec validation

**Files Created:**
- `torsion_gertsenshtein/cli/` — 8 modules (__init__, __main__, _derive, _simulate, _inspect, _list, _plot, _validate)
- `tests/test_cli.py` + `tests/test_cli_parsing.py` — 147 CLI tests
- 14 `theory.toml` files across examples

**Tests:** 743 Python tests passing (147 CLI tests), 0 ruff violations, 0 pyright errors

---

## Derived Fields TOML Feature (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Added `[[derived_fields]]` support to `theory.toml` for defining intermediate tensor fields (e.g., field strength F_ab) that are automatically expanded during Wolfram script generation.

**Usage:**
```toml
[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[lagrangian]
expression = "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d]"
```

**Implementation:** `_wls_derived_fields()` generates `DefTensor` + `MakeRule`; `_wls_lagrangian()` expands definitions via `/. rules`.

---

## Critical Review Pass 1 (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Comprehensive review addressing fail-fast enforcement, constraint solver ordering, dynamic axis naming, operator plugin API, all-JSON parametrized tests, physics validation, and epsilon tensor fixes.

**Key Changes:**
- Fail-fast symbolic coefficient validation
- Constraint solver ordering fixes
- Dynamic axis naming (supports up to 6D)
- Operator plugin API
- All-JSON parametrised tests
- Physics validation (energy conservation + analytical solutions)
- Epsilon `sqrt(|det(g)|)` fix

---

## Critical Review Pass 2 (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Second review pass focusing on code quality, stability, and extensibility.

**Key Changes:**
- `evolution_rate` refactored into 3 helper methods
- Configurable `constraint_eps` parameter
- CFL `check_stability()` method
- `MinkowskiMetricFactor` signature parameter
- Tensor symmetry reduction for rank 3+

---

## Massive 3-Form Example (February 2026)

**Status:** ✅ COMPLETE

**Summary:** End-to-end example of rank-3 antisymmetric tensor field decomposition, demonstrating symmetry reduction from 64 to 4 independent components in 3+1D.

**Key Fixes:**
- xAct's `Cycles` context: `xAct`xPerm`Cycles` vs `System`Cycles` — all pattern matching updated
- Explicit metric Lagrangian with `DefConstantSymbol` for mass parameter
- `EnumerateComponentTuples` symmetry reduction for antisymmetric rank-3

**Files Created:**
- `examples/massive_3form/massive_3form.wls` — Wolfram derivation
- `examples/massive_3form/theory.toml` — TOML config
- `examples/data/massive_3form.json` — 4-component KG system

**Tests:** 3 physics tests in dedicated test file

---

## Issue #71: 3+1D Klein-Gordon Working Example (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Created complete end-to-end 3+1D Klein-Gordon example demonstrating full 4D spacetime support with Wolfram derivation and Python simulation.

**Files Created:**
- `examples/scalar_field_3d/klein_gordon_3d.wls` — Wolfram derivation script (4D manifold, Minkowski metric)
- `examples/scalar_field_3d/kg_3d_simulation.py` — Python simulation (32³ grid, 4-panel visualization)

**Files Updated:**
- `examples/data/klein_gordon_3d.json` — Replaced hand-written version with pipeline-derived format
- `tests/test_3d_validation.py` — Updated field name from `"phi"` to `"phi_0"`

**Key Features:**
- Full Lagrangian → Euler-Lagrange → Component decomposition → JSON → Simulation pipeline in 3+1D
- Uses KG-prefixed xAct symbols (`kgM4`, `kgEta`, `kgCD`, `kgCart`) to avoid kernel conflicts
- Symbolic mass parameter `m2 = Symbol["m2"]` resolved at simulation time
- 32³ = 32,768 cell grid with periodic boundary conditions
- 4-panel visualization: z-profile (initial/final), x-y slices (initial/final), amplitude decay over time
- Demonstrates expected 3D physics: Gaussian pulse amplitude decays from 0.964 to 0.362 (ratio ~0.375)

**Impact:** Proves 3+1D capabilities are fully functional end-to-end, not just claimed.

**Tests:** 496 Python tests passing (including updated test_3d_validation.py)

---

## Issue #75: Operator Dimension Validation at Construction Time (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Added two-stage validation to fail fast when operators require dimensions not supported by the equation specification, catching errors at PDE construction instead of runtime.

**Implementation:**

**New Methods:**
- `PDEFromSpec._operator_min_dim(operator_name: str) -> int` — Static method returning minimum spatial dimension required by an operator
- `PDEFromSpec._validate_operator_dimensions()` — Instance method validating all operators in spec against spatial dimension

**Validation Logic:**
- Explicit registry: `gradient_z`, `laplacian_z`, `cross_derivative_xz` require 3D
- Dynamic patterns: `derivative_N_z` requires 3D, `derivative_Nx_My` requires 2D if y present
- Fail-fast: Raises `ValueError` at `__init__` with clear message indicating which operator and field failed

**Test Coverage:**
- 9 new tests in `TestOperatorDimensionValidation`
- Updated `test_unknown_operator_raises` to expect error at construction

**Example Error:**
```
ValueError: Operator 'gradient_z' in equation for 'phi_0' requires at least 3D spatial grid,
but the spec has spatial_dimension=2 (from 3D spacetime).
```

**Impact:** Developer experience improvement — errors caught immediately with clear diagnostic messages instead of cryptic runtime failures.

**Tests:** 496 Python tests passing (9 new)

---

## Issue #67: Replace assert isinstance with Explicit TypeErrors (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Replaced 15 `assert isinstance(...)` statements in `pde_builder.py` with explicit `if not isinstance(...): raise TypeError(...)` checks to prevent silent failures when Python runs with `-O` optimization flag.

**Changes:**
- 15 isinstance checks converted to explicit TypeError raises
- 9 `# pyright: ignore[reportUnnecessaryIsInstance]` comments added (pyright doesn't understand the pattern)
- 4 docstrings updated to document raised TypeErrors
- 2 tests updated to expect TypeError instead of ValueError for complex numbers
- 1 pyright fix: `float(cast("SupportsFloat", result))` for line 835

**Linting Fixes:**
- ruff EM102: Extracted 15 inline f-strings to `msg = f"..."` variables
- ruff DOC501: Added TypeError to 5 method docstrings
- ruff TRY004: Changed ValueError to TypeError for complex number check (autofix side effect)
- Removed stale `# noqa: C901`, added `PLR0912` to `evolution_rate`

**Rationale:** Python's `-O` flag strips all assert statements for performance. Using assert for validation means code silently breaks in production when optimizations are enabled.

**Tests:** 487 Python tests passing after changes

---

## Issue #85: Unified Derivative Classification (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Refactored 5 dimension-specific derivative wrappers into a single dimension-agnostic `ExtractDerivativeProfile` function, saving ~108 lines of code and improving maintainability.

**Replaced Functions:**
- `Extract1DSpatialDerivativeForm` (1+1D specific)
- `Extract2DSpatialDerivativeForm` (2+1D specific)
- `Extract3DSpatialDerivativeForm` (3+1D specific)
- `ExtractTimeDerivativeOrder` (dimension-agnostic but redundant)
- `ExtractSpatialDerivativeOrders` (dimension-agnostic but redundant)

**New Unified Function:**
```mathematica
ExtractDerivativeProfile[term_, dim_] := Module[{...},
  (* Returns: <|"time" -> n, "space" -> {m_x, m_y, m_z, ...}|> *)
]
```

**Key Improvements:**
- Single source of truth for derivative classification across all dimensions
- Handles 1+1D (2-arg), 2+1D (3-arg), 3+1D (4-arg) Derivative forms uniformly
- Eliminates code duplication and dimension-specific logic branches
- ~108 lines of code removed

**Test Coverage:**
- 10 new Wolfram tests in `test_common_utilities.wls` covering 1D/2D/3D cases
- All 496 Python tests passing after refactor

**Impact:** Cleaner codebase, easier to extend to 4+1D or higher dimensions in the future.

---

## Issue #79: Mixed Time-Spatial Derivatives (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Added support for mixed time-space derivative terms like `∂_t ∂_x A` (common in Hubble friction, curved spacetime, and momentum gradient terms).

**Wolfram Side:**
- `ClassifySpatialProfile[spatialOrders]` — Identifies gradient direction from spatial derivative pattern
- `ExtractSpatialOperatorFromMixed[term, dim]` — Extracts spatial operator name from mixed derivative term

**Python Side:**
- Multi-axis handler in `_apply_operator` for `gradient_x(pi_i)` terms
- Validates that referenced momentum field exists
- Applies spatial gradient to time derivative (momentum) field

**Use Cases:**
- Hubble friction: `H * ∂_t φ` appears as `gradient_t(phi)` → handled as `first_derivative_t`
- Momentum gradients: `∂_x (∂_t φ)` appears as `gradient_x(pi_i)` in curved spacetime

**Test Coverage:**
- Wolfram tests for mixed derivative classification
- Python tests for momentum field references

**Impact:** Enables realistic curved spacetime simulations (de Sitter, conformal, S²) with proper friction terms.

**Tests:** All 496 Python tests passing

---

## Issue #68: Eval Validation for Coefficient Expressions (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Added `_validate_eval_result` guard to catch NaN, Inf, and complex values in coefficient evaluation, preventing silent physics errors from invalid mathematical expressions.

**Implementation:**
```python
def _validate_eval_result(self, result: complex, expression: str, context: dict) -> None:
    """Validate that evaluated coefficient is finite and real."""
    if np.isnan(result):
        msg = f"Coefficient expression '{expression}' evaluated to NaN with {context}"
        raise ValueError(msg)
    if np.isinf(result):
        msg = f"Coefficient expression '{expression}' evaluated to Inf with {context}"
        raise ValueError(msg)
    if np.iscomplex(result) or result.imag != 0:
        msg = f"Coefficient must be real, got {result} from '{expression}' with {context}"
        raise TypeError(msg)
```

**Called From:**
- `_resolve_coefficient_at_point` (position-dependent coefficients)
- After every `eval()` of Mathematica expressions

**Test Coverage:**
- 13 new tests in `TestEvalValidation`:
  - `test_validate_eval_result_nan` / `test_nan_in_position_dependent_raises`
  - `test_validate_eval_result_inf` / `test_inf_in_position_dependent_raises` / `test_overflow_to_inf_raises`
  - `test_validate_eval_result_complex` / `test_sqrt_negative_gives_nan_raises`
  - `test_validate_eval_result_valid_float` / `test_validate_eval_result_valid_int` / `test_validate_eval_result_zero`

**Impact:** Physics errors (like `1/0` or `sqrt(-1)` in coefficients) now fail loudly with clear diagnostics instead of producing silent NaN propagation.

**Tests:** 496 Python tests passing (13 new)

---

## Phase 13: Rank 3+ Tensor Support (February 2026)

**Status:** ✅ COMPLETE (Issue #70)

**Summary:** Extended component extraction to support rank-3 and higher tensors (epsilon tensors, Riemann curvature, field strength tensors, etc.).

**New Function:**
```mathematica
ReplaceHigherRankFieldComponents[componentEq_, fieldTemplate_, chart_] := Module[{...},
  (* Detects rank-3+ tensors and replaces with component functions *)
  (* Example: epsilon[a,b,c] -> epsilon012[t,x,y] *)
]
```

**Supported Ranks:**
- Rank 0: Scalars (phi[])
- Rank 1: Vectors (A[-a])
- Rank 2: 2-tensors (h[-a,-b])
- **Rank 3+: NEW** — Epsilon tensors (ε[-a,-b,-c]), Riemann tensors (R[-a,-b,-c,-d])

**Integration:**
- Called automatically in `DecomposeToComponents` after standard rank-1/2 replacement
- Uses xAct introspection to detect tensor rank
- Constructs component names from index values

**Test Coverage:**
- 27 new Wolfram tests in `test_component_decompose.wls`
- 4 new Python tests for higher-rank tensor JSON loading

**Impact:** Enables topological field theories (Chern-Simons with ε_{abc}), gravitational theories (Riemann curvature), and non-Abelian gauge theories.

**Tests:** ~100 Wolfram tests + 496 Python tests passing (31 new across both)

---

## Phase 12: Auto-Computed Mass/Coupling Matrices (February 2026)

**Status:** ✅ COMPLETE

**Summary:** Automated computation of mass and coupling matrices from equation terms with symbolic coefficient preservation, eliminating manual matrix construction and enabling runtime parameter sweeps.

**Convention:**
```
mass_matrix[i][j] = -(coefficient of identity(field_j) in equation_i)
coupling_matrix[i][j] = -(coefficient of identity(field_j) in equation_i, where i≠j)
```

**Implementation:**

**Wolfram Side (`ExportJSON.wl`):**
- `ExtractMassCouplingFromEquations[fieldEquations]` — Parses RHS terms to extract matrix coefficients
- Preserves symbolic expressions: `"coefficient": -1.0, "coefficient_symbolic": "-m2"`
- Defense-in-depth: Both numeric and symbolic coefficients exported

**Python Side (`json_loader.py`):**
- `EquationSystem._compute_matrices_from_terms` — Computes matrices from equation terms at load time
- Returns 4-tuple: `(mass_numeric, coupling_numeric, mass_symbolic, coupling_symbolic)`
- `__post_init__` guard: UserWarning if constructor-provided matrices inconsistent with terms

**Symbolic Preservation:**
- `mass_matrix_symbolic` / `coupling_matrix_symbolic` preserve exact Mathematica expressions
- Evaluated at runtime using `_mathematica_to_python` expression evaluator
- Enables parameter sweeps without re-deriving equations

**Test Coverage:**
- 453 tests passing after Phase 12 implementation
- Tests for matrix auto-computation, symbolic preservation, consistency validation

**Impact:** Users can sweep parameters (`m2`, `g`, `H`) in simulations without regenerating JSON — symbolic coefficients evaluated dynamically.

**Example:**
```json
{
  "coefficient": 1.0,
  "coefficient_symbolic": "m2",
  "operator": "identity",
  "field": "phi_0"
}
```

**Tests:** 453 Python + Wolfram tests passing after implementation

---

## Phase 4: Pipeline Robustness & Testing (February 2026)

### 🎯 Overview

**Status:** ✅ **ALL CRITICAL IMPLEMENTATION COMPLETE**

**Total Changes:**
- 4 Wolfram test files fixed (81 tests passing)
- 3 utility functions in CommonUtilities.wl enhanced
- 5 development utility scripts created
- 3 private helper function usage strings added
- All modules now have comprehensive header documentation
- Test runner script with kernel caching support

**Impact:**
- **100% test pass rate** (81 Wolfram tests + 186 Python tests)
- Robust Wolfram test infrastructure with proper xAct symbol management
- Complete pipeline validation scripts for development workflows
- Professional documentation coverage

---

### 🔧 Issue 18: Wolfram Test Symbol Conflicts ✅ RESOLVED

#### Problem
Running `./scripts/run_wolfram_tests.sh` failed with xAct kernel caching issues:
- `ValidateSymbol::used: Symbol TestM2 is already used as a manifold`
- Multiple test files defining same symbols
- Bash arithmetic exit code issues

#### Root Causes
1. **DefMetric incorrect syntax**: Passing manifold instead of covariant derivative name as 3rd argument
2. **Pattern matching gaps**: `ExtractNumericCoefficient` missing patterns for bare tensors (`f[]`, `f[_]`)
3. **xAct introspection issues**: `IsCovDOperator` not handling applied CD forms like `CD[-a][phi[]]`
4. **Test runner bug**: `((PASSED++))` returns exit code 1 when PASSED=0, triggering `set -e`

#### Solutions

**1. Fixed DefMetric Syntax**

**File:** [tests/wolfram/test_euler_lagrange.wls](tests/wolfram/test_euler_lagrange.wls:40)

```mathematica
(* BEFORE: Incorrect - passes manifold as 3rd arg *)
DefMetric[-1, elTestEta[-a, -b], ELTestM2, SymbolOfCovD -> {";", "D"}];

(* AFTER: Correct - passes covariant derivative name *)
DefMetric[-1, elTestEta[-a, -b], elTestCD, SymbolOfCovD -> {";", "D"}, PrintAs -> "g"];
```

**Impact:** Eliminates `ValidateSymbol::used` errors from incorrect DefMetric usage.

**2. Enhanced ExtractNumericCoefficient**

**File:** [torsion_gertsenshtein/wolfram/CommonUtilities.wl](torsion_gertsenshtein/wolfram/CommonUtilities.wl:380-383)

```mathematica
(* Added patterns for xAct tensor forms *)
f_[] /; StringContainsQ[ToString[f], ToString[fieldName]] :> 1,  (* Scalar: phi[] *)
f_[_] /; StringContainsQ[ToString[f], ToString[fieldName]] :> 1,  (* Vector: A[-a] *)
```

**Impact:** Correctly extracts coefficients from bare xAct tensors, not just coordinate-form functions.

**3. Fixed IsCovDOperator for Applied Forms**

**File:** [torsion_gertsenshtein/wolfram/CommonUtilities.wl](torsion_gertsenshtein/wolfram/CommonUtilities.wl:139-150)

```mathematica
IsCovDOperator[expr_] := Quiet[
  Module[{h = Head[expr], baseSymbol},
    If[Head[h] === Symbol,
      TrueQ[xAct`xTensor`CovDQ[h]],
      (* Extract base: CD from CD[-a][phi[]] *)
      baseSymbol = Head[h];
      TrueQ[xAct`xTensor`CovDQ[baseSymbol]]
    ]
  ],
  {xAct`xTensor`CovDQ::argx, General::stop}
];
```

**Impact:** Properly detects covariant derivatives in applied form, not just bare CD symbols.

**4. Fixed Test Runner Bash Arithmetic**

**File:** [scripts/run_wolfram_tests.sh](scripts/run_wolfram_tests.sh:38)

```bash
# BEFORE: ((PASSED++)) returns exit 1 when PASSED=0, triggers set -e
if wolframscript -file "$TEST_PATH"; then
    ((PASSED++))

# AFTER: PASSED=$((PASSED + 1)) always returns exit 0
if wolframscript -file "$TEST_PATH"; then
    PASSED=$((PASSED + 1))
```

**Impact:** Test runner no longer exits prematurely on first passing test.

#### Verification Results

**Before Fixes:**
```bash
$ ./scripts/run_wolfram_tests.sh
ValidateSymbol::used: Symbol TestM2 is already used as a manifold.
TorsionQ::unknown: Unknown covariant derivative TestM2.
Throw::nocatch: Uncaught Throw[Null] returned to top level.
[Exit 1]
```

**After Fixes:**
```bash
$ ./scripts/run_wolfram_tests.sh
=== Running Wolfram Tests ===
Running: tests/wolfram/test_euler_lagrange.wls
*** ALL TESTS PASSED *** (7 tests)
Running: tests/wolfram/test_common_utilities.wls
*** ALL TESTS PASSED *** (24 tests)
Running: tests/wolfram/test_export_json.wls
*** ALL TESTS PASSED *** (50 tests)
=== Wolfram Test Summary ===
Passed: 3
Failed: 0
*** ALL WOLFRAM TESTS PASSED ***
```

**Kernel Caching Verification:**
```bash
$ ./scripts/run_wolfram_tests.sh && ./scripts/run_wolfram_tests.sh
# Both runs pass - kernel caching properly handled
```

---

### 🛠️ Issue 17: Development Utility Scripts ✅ COMPLETE

Created 5 helper scripts to streamline local development and testing workflows.

#### Scripts Created

**1. `scripts/run_wolfram_tests.sh`**
- Runs all Wolfram unit tests with summary output
- Tracks pass/fail counts
- Lists failed tests for debugging
- Exit codes: 0 (all pass), 1 (any fail)

**2. `scripts/run_examples.sh`**
- Regenerates all JSON files from Lagrangian derivations
- Runs: Klein-Gordon, EM, Coupled Scalars, Chern-Simons, Navier-Cauchy
- Validates symbolic derivation pipeline

**3. `scripts/full_test.sh`**
- Complete test suite: Python (186 tests) + Wolfram (81 tests)
- One-command verification before commits
- Sequential execution with clear output separation

**4. `scripts/validate_pipeline.sh`**
- End-to-end pipeline validation
- Derives Klein-Gordon equations → JSON → Python simulation
- Validates JSON file creation and structure

**5. `scripts/lint_wolfram.sh`**
- Basic syntax checking for Wolfram modules
- Loads each module and reports errors
- Quick verification after Wolfram changes

#### Usage

```bash
# Run just Wolfram tests
./scripts/run_wolfram_tests.sh

# Regenerate all example JSON files
./scripts/run_examples.sh

# Full test suite before committing
./scripts/full_test.sh

# Validate entire pipeline works end-to-end
./scripts/validate_pipeline.sh

# Check Wolfram module syntax
./scripts/lint_wolfram.sh
```

**Impact:** Developers can quickly validate changes without memorizing commands or paths.

---

### 📚 Issue 15: Module Header Documentation ✅ COMPLETE

All Wolfram modules now have comprehensive header comments with MODULE, PURPOSE, DEPENDENCIES, and DATA FLOW sections.

#### Documentation Added

**ExportJSON.wl** - Complete module documentation including:
- Supported operators (identity, laplacian, laplacian_x/y/z, gradient_x/y/z, cross_derivative_xy/xz/yz)
- PDE types (elliptic, parabolic, hyperbolic) via LHS structure
- Momentum gradient handling for mixed time-space derivatives
- Data flow from ComponentDecompose → ExportEquationSystem

**ComponentDecompose.wl** - Complete module documentation including:
- Tensor to scalar component decomposition process
- additionalFields parameter for cross-field coupling
- Dimension-agnostic design using GetChartDimension
- Automatic epsilon tensor evaluation for topological terms

**Linearize.wl** - Complete module documentation including:
- Zero background vs custom background linearization
- Polynomial degree selection vs xPert integration
- Gauge fixing considerations
- Usage patterns with examples

**CommonUtilities.wl** - Already had comprehensive documentation:
- CD → Derivative conversion process
- Dimension validation and max supported dimensions
- Sign conventions for Minkowski spacetime
- xAct introspection helpers

**Impact:** Developers can understand module responsibilities and data flow without reading implementation details.

---

### 📝 Issue 14: Private Helper Function Documentation ✅ COMPLETE

Added `::usage` strings for 3 private helper functions in Linearize.wl for completeness.

**File:** [torsion_gertsenshtein/wolfram/Linearize.wl](torsion_gertsenshtein/wolfram/Linearize.wl:48-60)

```mathematica
SelectLinearTerms::usage =
  "SelectLinearTerms[expr, field] extracts terms that are linear (degree-1) in the field. \
For polynomial expressions, returns Coefficient[expr, field, 1] * field. For complex \
expressions, selects terms with exactly one field factor.";

HasExactlyOneFieldFactor::usage =
  "HasExactlyOneFieldFactor[term, field] returns True if term contains exactly one \
occurrence of the field or its derivatives. Used to identify linear terms in non-polynomial \
expressions.";

ExpandAroundBackground::usage =
  "ExpandAroundBackground[eom, field, background] performs manual epsilon expansion \
around a non-zero background: field -> background + epsilon*field, then extracts O(epsilon) \
terms. Fallback method when xPert systematic perturbation is not applicable.";
```

**Impact:** Complete API documentation coverage, including internal helpers for maintainers.

---

### 📊 Phase 4 Summary Statistics

#### Test Coverage
- **Wolfram Tests:** 81 tests across 3 files (100% pass rate)
- **Python Tests:** 186 tests (100% pass rate)
- **Total:** 267 tests passing
- **Kernel Caching:** Verified with consecutive runs

#### Documentation
- **Module Headers:** 4/4 complete (ExportJSON, ComponentDecompose, Linearize, CommonUtilities)
- **Public API Functions:** 100% documented with `::usage` strings
- **Private Helpers:** 3/3 documented (Linearize.wl)
- **Scripts:** 5/5 with comprehensive README documentation

#### Code Quality
- **Fixed Functions:** 4 (DefMetric syntax, ExtractNumericCoefficient, ExtractCoefficientWithSymbolic, IsCovDOperator)
- **Fixed Scripts:** 1 (run_wolfram_tests.sh bash arithmetic)
- **New Scripts:** 5 (run_wolfram_tests, run_examples, full_test, validate_pipeline, lint_wolfram)

---

### 🎓 Key Technical Insights

#### DefMetric Signature
```mathematica
DefMetric[signdet, metric[-a, -b], covd, options]
(*         ^        ^                ^     ^
           |        |                |     Optional: SymbolOfCovD, PrintAs
           |        |                Covariant derivative NAME (not manifold!)
           |        Metric tensor with abstract indices
           Signature determinant: -1 for Minkowski, +1 for Euclidean
*)
```

**Common Error:** Passing manifold as 3rd argument instead of covariant derivative name.

#### xAct Applied Forms
When xAct applies a covariant derivative:
```mathematica
CD[-a][phi[]]  (* Head = CD[-a], not CD *)
```

To detect: `Head[Head[expr]]` gives the base symbol `CD`.

#### Bash Arithmetic Exit Codes
```bash
((PASSED++))  # Returns exit 1 when PASSED=0 (falsy value)
PASSED=$((PASSED + 1))  # Always returns exit 0
```

With `set -e`, the first form causes premature exit.

---

### 🛠️ Files Modified

| File | Changes |
|------|---------|
| `torsion_gertsenshtein/wolfram/CommonUtilities.wl` | Enhanced coefficient extraction, fixed IsCovDOperator |
| `torsion_gertsenshtein/wolfram/Linearize.wl` | Added 3 private helper usage strings |
| `tests/wolfram/test_euler_lagrange.wls` | Fixed DefMetric syntax |
| `tests/wolfram/test_common_utilities.wls` | Verified with fixed utilities |
| `tests/wolfram/test_export_json.wls` | Verified with fixed patterns |
| `scripts/run_wolfram_tests.sh` | Fixed bash arithmetic, enhanced output |
| `scripts/run_examples.sh` | NEW: Regenerate all JSON files |
| `scripts/full_test.sh` | NEW: Complete test suite |
| `scripts/validate_pipeline.sh` | NEW: End-to-end validation |
| `scripts/lint_wolfram.sh` | NEW: Wolfram syntax checking |

---

### ✨ Impact Summary

#### Before Phase 4
- ❌ Wolfram tests failing due to kernel caching issues
- ❌ No development utility scripts
- ❌ Incomplete documentation for private helpers
- ❌ Manual test execution required memorizing commands

#### After Phase 4
- ✅ 81 Wolfram tests + 186 Python tests passing (100%)
- ✅ 5 utility scripts for streamlined development
- ✅ Complete documentation coverage (modules + APIs + helpers)
- ✅ One-command validation for all workflows
- ✅ Robust kernel caching handling
- ✅ Professional test infrastructure

---

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

This update adds complete symbolic tensor algebra capabilities to the project through automated Wolfram Engine and xAct/xCoba package installation, enabling derivation of linearized field equations for the Gertsenshtein effect.

---

## 🚀 Symbolic Computing Infrastructure

**New Components:**

- 5 installation/verification scripts
- Wolfram Engine 14.3 headless integration
- xAct 1.2.1 tensor algebra suite
- xCoba coordinate-based tensor computations
- GLIBC compatibility fixes for xPerm binary
- Comprehensive verification and smoke testing

**Impact:**

- Enables symbolic derivation of coupled EM/gravity/torsion field equations
- Automated container setup for symbolic computations
- Reproducible installation across development environments
- Ready for xAct-based coordinate calculations and Christoffel symbol derivations

### Wolfram Engine Integration

**Files:** `scripts/install-wolfram-engine.sh`, `scripts/activate-wolfram.sh`

Automated installation and activation of Wolfram Engine for symbolic computations:

- Downloads and installs Wolfram Engine 14.3
- Handles system dependencies (libglib2.0, X11 libraries, etc.)
- Interactive license activation with free Wolfram ID
- Verification of installation and activation status

**Container Integration:** Dev container automatically installs if installer placed in `third_party/`

### xAct/xCoba Tensor Algebra Suite

**Files:** `scripts/install-xact-xcoba.sh`, `scripts/xact_smoke.wl`

Complete tensor algebra capabilities for General Relativity calculations:

- Official xAct 1.2.1 package suite (xCore, xPerm, xTensor, xCoba)
- GLIBC compatibility fix: recompiles xPerm binary using Wolfram's MathLink compiler
- Installs to user Applications directory (`~/.WolframEngine/Applications/`)
- Smoke test verifies manifold/chart definitions, metric tensors, Riemann tensor operations

```wolfram
DefManifold[M, 4, IndexRange[a, z]];
DefChart[spherical, M, {0, 1, 2, 3}, {t[], r[], θ[], φ[]}];
DefMetric[-1, g[-a, -b], CD];
```

### Comprehensive Verification

**File:** `scripts/verify-wolfram-setup.sh`

Multi-stage verification process:

1. **Wolfram Engine**: binary availability, activation, version/license check
2. **xAct Installation**: package directory structure, core components present
3. **xPerm Compatibility**: binary dependencies, GLIBC version compatibility
4. **Package Loading**: successful loading of xCore, xPerm, xTensor, xCoba
5. **Smoke Test**: complete tensor operation workflow verification

**Container Automation:** `postAttachCommand` auto-installs xAct if missing and prompts for verification

---

## Previous Improvements (January 2026)

This section covers the earlier Klein-Gordon simulator improvements focusing on robustness, performance, security, and code quality.

---

## 🎯 Overview

**Total Changes:**

- 8 core files modified
- 3 new test files added
- 90 tests passing (11 new tests)
- 14 linting issues resolved
- 5+ type safety improvements
- **NEW:** 5 symbolic computing scripts added
- **NEW:** Wolfram Engine & xAct integration

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
