#!/usr/bin/env python3
"""Script to create GitHub issues from the issue tracking plan.

This script creates all 25 issues identified in the codebase review.
Requires GitHub CLI (gh) to be installed and authenticated.

Usage:
    python scripts/create_github_issues.py [--dry-run]
"""

import subprocess  # noqa: S404
import sys

ISSUES = [
    # CRITICAL BUGS & SECURITY
    {
        "title": "Replace Assertions with Explicit Error Handling in PDE Operators",
        "body": """## Problem
13 critical type checks use `assert isinstance(...)` statements in operator functions. Assertions are stripped with Python's `-O` optimization flag, causing silent failures in production.

## Location
`torsion_gertsenshtein/symbolic/pde_builder.py` (lines 61, 74, 76, 89, 91, 100, 102, 120, 621, 647, 715, 781, 822)

## Example
```python
# Current (problematic):
assert isinstance(component, ScalarField)

# Should be:
if not isinstance(component, ScalarField):
    raise TypeError(f"Expected ScalarField, got {type(component).__name__}")
```

## Impact
In optimized Python, operator functions silently accept wrong types, causing cryptic downstream errors.

## Affected Functions
- `_op_gradient()`, `_op_directional_laplacian()`, `_op_cross_derivative()`
- `_op_biharmonic()`, `_op_nth_derivative()`
- `_get_field_from_state()`, `_compute_rhs_for_component()`, `evolution_rate()`

## Fix
Replace all assertions with explicit `if not ... raise TypeError(...)` checks.
""",
        "labels": ["bug", "priority: critical", "type-safety"],
    },
    {
        "title": "Strengthen Mathematica Expression Evaluation Security",
        "body": """## Problem
Uses `eval()` to evaluate Mathematica-generated expressions in Python with potential security and validation issues:
- No validation that evaluated coefficients are finite (NaN/Inf)
- No protection against complex number results when float expected
- Edge cases in namespace restrictions

## Location
`torsion_gertsenshtein/symbolic/pde_builder.py` (line 492)

## Current Code
```python
result = eval(py_expr, {"__builtins__": {}}, namespace)  # noqa: S307
```

## Edge Cases
1. **Malformed expressions:** `sqrt(-1)` returns complex, not float
2. **Floating-point errors:** `1.0/0.0` returns inf
3. **Syntax injection:** While namespace is restricted, edge cases may exist

## Recommendations
1. Add post-evaluation validation:
   ```python
   if np.isnan(result).any() or np.isinf(result).any():
       raise ValueError(f"Coefficient evaluation produced NaN/Inf: {result}")
   ```
2. Type-check result before returning
3. Consider replacing `eval()` with `ast.literal_eval()` where possible or sympy evaluation
4. Add unit tests for edge cases
""",
        "labels": ["security", "priority: critical", "validation"],
    },
    {
        "title": "Add Wolfram Tests to GitHub Actions CI",
        "body": """## Problem
Wolfram tests only run locally (81 tests across 4 test files). Changes to Wolfram code can break silently on PR merge without automated validation of Wolfram pipeline correctness.

## Location
`.github/workflows/test.yml`, `scripts/run_wolfram_tests.sh`

## Blockers
- Wolfram Engine license activation (requires interactive input or pre-provisioned license)
- xAct installation (240MB download)
- Build time concerns

## Recommendations
1. Create optional workflow triggered by label (`run-wolfram-tests`) or manual dispatch
2. Cache xAct installation in GitHub Actions
3. Use self-hosted runner with pre-activated Wolfram license (if available)
4. Alternative: Weekly scheduled workflow instead of per-PR

## Benefits
Catch Wolfram pipeline breakage before merging, ensure JSON generation correctness.
""",
        "labels": ["ci-cd", "priority: critical", "testing"],
    },

    # CRITICAL FEATURES
    {
        "title": "Support Rank-3+ Tensor Decomposition",
        "body": """## Current Limitation
Explicit error thrown for rank ≥ 3 tensors. Only handles rank 0 (scalar), rank 1 (vector), rank 2 (tensor).

## Location
`torsion_gertsenshtein/wolfram/ComponentDecompose.wl` (lines 377-383)
`torsion_gertsenshtein/wolfram/ExportJSON.wl` (line 339)

## Error Message
```mathematica
Throw["ReplaceTensorFieldComponents: Rank-3 field replacement is not yet implemented..."]
```

## Impact
- Cannot decompose stress-energy tensor (rank 2) with derivatives (effectively rank 3)
- Riemann curvature tensor (rank 4) unsupported
- Higher-rank gauge fields blocked

## Use Cases
- Elastic continua with stress tensor derivatives
- Full general relativity (Riemann, Ricci tensors)
- Yang-Mills field strength (rank 2 with gauge indices)

## Implementation Notes
- Need pattern-matching rules for `ReplaceTensorFieldComponents` at rank 3+
- Symmetry reduction logic in `GenerateIndexTuples` needs generalization
- May require significant xAct introspection work
""",
        "labels": ["enhancement", "priority: critical", "wolfram", "tensors"],
    },
    {
        "title": "Add 3+1D Spacetime Examples",
        "body": """## Current State
- Documentation mentions 3+1D support
- `klein_gordon_3d.json` exists but no working example
- Only 1+1D and 2+1D examples implemented
- `$MaxSupportedDimension = 4` in code but not demonstrated

## Missing Examples
1. `examples/scalar_field_3d/` - Basic Klein-Gordon in full 4D spacetime
2. `examples/em_3d/` - Electromagnetic waves in 3D space
3. `examples/gravity_3d/` - Gravitational wave propagation

## Recommendations
1. Create 3+1D Klein-Gordon example with:
   - Wolfram derivation script
   - Python simulation with 3D visualization
   - Validation against analytical solutions
2. Add to CI/CD pipeline
3. Document performance considerations (3D grids scale as N³)

## Benefits
Validates dimension-generalization claims, provides template for users.
""",
        "labels": ["documentation", "priority: critical", "examples", "3d"],
    },
    {
        "title": "Implement Automatic Gauge Fixing",
        "body": """## Current State
`DeDonderGaugeRules` infrastructure exists but returns empty rules `{}`. Users must manually apply gauge conditions.

## Location
`torsion_gertsenshtein/wolfram/Linearize.wl` (lines 170-191)

## Impact
- Gauge theory examples cannot derive gauge-fixed equations of motion automatically
- Lorenz gauge for electromagnetism requires manual constraint
- Coulomb gauge for quantum field theory blocked

## Implementation Path
1. Integrate with xPert's gauge transformation utilities
2. Add Lorenz gauge (∂^μ A_μ = 0) as first implementation
3. Support Coulomb gauge (∇·A = 0) for non-relativistic
4. Document pattern for users to add custom gauges

## Priority Justification
Listed in MEMORY.md as "Future Priority #3".
""",
        "labels": ["enhancement", "priority: critical", "gauge-theory", "wolfram"],
    },

    # HIGH PRIORITY
    {
        "title": "Add Animation Module Test Coverage",
        "body": """## Problem
Animation features completely untested (250+ lines, 0 tests). FFMpeg vs Pillow writer selection could fail silently. Frame rate calculations and output path validation unchecked.

## Location
`torsion_gertsenshtein/kgsim/animation_builder.py`

## Test Requirements
1. AnimationBuilder initialization with various configs
2. Writer selection (FFMpeg available vs fallback to GIF)
3. 1D and 2D animation rendering (at least smoke tests)
4. Frame rate edge cases (very fast, very slow)
5. Output path validation (invalid paths, permissions)
6. Integration with simulation results

## Recommendation
Create `tests/test_animation_builder.py` with ~10-15 test cases.
""",
        "labels": ["testing", "priority: high", "animation"],
    },
    {
        "title": "Add Code Coverage Reporting to CI",
        "body": """## Current State
`pytest-cov` installed but not used in CI. No coverage metrics visible in PR checks. Unknown which code paths are exercised.

## Location
`.github/workflows/test.yml`

## Implementation
1. Add to pytest step:
   ```yaml
   - run: uv run pytest --cov=torsion_gertsenshtein --cov-report=xml --cov-report=term
   ```
2. Upload coverage to Codecov or similar service
3. Set minimum coverage threshold (suggest 80%)
4. Add coverage badge to README.md

## Benefits
- Visibility into untested code paths
- Prevent coverage regressions in PRs
- Identify high-risk areas
""",
        "labels": ["ci-cd", "priority: high", "testing"],
    },
    {
        "title": "Validate Grid Dimensions During PDE Construction (Not Runtime)",
        "body": """## Problem
Grid dimension validated in `evolution_rate()` (runtime) not `__init__` (construction). Mismatch causes cryptic py-pde operator errors before clear validation message. Users waste time debugging wrong issues.

## Location
`torsion_gertsenshtein/symbolic/pde_builder.py` (lines 695, 789-797)

## Current Flow
```python
# Line 695 - zero field created without dimension check
result = ScalarField(grid, data=0.0)

# Lines 789-797 - dimension checked later
if grid_dim != expected_dim:
    raise ValueError(...)  # Too late!
```

## Fix
1. Add `_validate_grid_dimension()` method in `PDEFromSpec.__init__()`
2. Check grid.dim against `spec.spacetime.dimension`
3. Validate operator requirements (e.g., `gradient_z` needs 3D) at construction time
4. Pre-flight check prevents wasted computation
""",
        "labels": ["bug", "priority: high", "validation"],
    },
    {
        "title": "Improve Coefficient Resolution Performance",
        "body": """## Problem
Coefficient resolution happens per-term per-timestep. Repeated evaluation of same expressions (e.g., `m2(t)` evaluated 100+ times per step). Spatial coordinates recomputed for every term.

## Location
`torsion_gertsenshtein/symbolic/pde_builder.py` (lines 698-742)

## Current Bottleneck
```python
for term in eq.rhs_terms:
    coefficient = self._resolve_coefficient_at_point(term, t, grid)
    # This evaluates Mathematica expression every time
```

## Optimization Strategies
1. **Cache time-dependent coefficients:** Evaluate `m2(t)` once per timestep, reuse for all terms
2. **Pre-evaluate spatial coordinates:** Compute meshgrid once, pass to all terms
3. **Vectorize multiple terms:** If multiple terms share coordinate dependence, batch evaluation
4. **Memoization:** Use `functools.lru_cache` for expensive expressions

## Expected Impact
10-20% speedup for position-dependent systems (e.g., `sphere_kg/`).

## Recommendation
Profile `sphere_kg` example, implement caching layer for coefficients.
""",
        "labels": ["performance", "priority: high", "optimization"],
    },
    {
        "title": "Document JSON Schema with Detailed Guide",
        "body": """## Problem
JSON schema defined in code (`torsion_gertsenshtein/symbolic/json_loader.py`) but not documented. No explanation of field constraints, valid ranges, operator format. Users unclear how to extend schema or migrate old JSON.

## Needed Documentation
1. **Field-by-field reference:** Each JSON key explained with examples
2. **Operator format guide:** How to specify `laplacian`, `gradient_x`, custom operators
3. **Coefficient types:** Constant vs time-dependent vs position-dependent
4. **Cross-field references:** Syntax for momentum fields, coupled terms
5. **Validation rules:** What makes a valid JSON (required fields, constraints)
6. **Versioning strategy:** How schema evolves, migration path
7. **Extension guide:** Adding custom operators

## Target File
`docs/JSON_SCHEMA_GUIDE.md` (NEW)

## Benefits
Users can generate JSON programmatically, troubleshoot schema errors, extend for new physics.
""",
        "labels": ["documentation", "priority: high", "json-schema"],
    },

    # MEDIUM PRIORITY
    {
        "title": "Add Tests for Observers, Profiling, Runners Modules",
        "body": """## Modules Without Tests
1. `kgsim/observers.py` - Energy conservation calculations, custom observers (0 tests)
2. `kgsim/profiling.py` - Timer context managers, profiling callbacks (0 tests)
3. `kgsim/runners.py` - Solver execution pipeline (0 tests)

## Test Requirements
- Energy computation accuracy (compare against analytical)
- Observer state handling (serialization, callbacks)
- Timer correctness (mock time.time())
- Runner integration with different solver configs

## Estimated Effort
~20-30 test cases total across 3 modules.
""",
        "labels": ["testing", "priority: medium"],
    },
    {
        "title": "Handle Mixed Time-Space Cross-Derivatives Properly",
        "body": """## Problem
Terms like `∂_t ∂_x ∂_y` (time + 2 spatial) converted to momentum gradient with warning. `ExtractSpatialGradientFromMixed` defaults to `gradient_x(pi)` but loses direction info. In 2+1D, `Derivative[1,1,1]` is ambiguous.

## Location
`torsion_gertsenshtein/wolfram/ExportJSON.wl` (lines 533-542)

## Current Handling
```mathematica
Print["WARNING: Multiple spatial derivatives..."];
(* Falls back to first spatial direction *)
```

## Better Approach
1. Detect pattern `Derivative[1, dx, dy, ...]` where dx+dy+... > 1
2. Generate nested operators: `gradient_y(gradient_x(pi_0))`
3. Or create new operator type: `mixed_derivative_xy_t`
4. Validate in Python that nested operators are supported

## Impact
Currently minor (no examples hit this case), but blocks complex elasticity/fluid dynamics.
""",
        "labels": ["bug", "priority: medium", "wolfram", "derivatives"],
    },
    {
        "title": "Expand _mathematica_to_python Function Set",
        "body": """## Problem
Only ~7 math functions supported: `Sqrt`, `Sin`, `Cos`, `Tan`, `Exp`, `Log`. Missing inverse trig, hyperbolic, power functions. Unknown operators silently become Python identifiers (causes NameError at runtime).

## Location
`torsion_gertsenshtein/symbolic/pde_builder.py` (lines 390-397)

## Missing Functions
- Inverse trig: `ArcSin`, `ArcCos`, `ArcTan`, `ArcTan2`
- Hyperbolic: `Sinh`, `Cosh`, `Tanh`, `ArcSinh`, `ArcCosh`, `ArcTanh`
- Power: `Power[x, n]` → `x**n`
- Special: `Erf`, `BesselJ`, `BesselY` (for advanced physics)

## Implementation Example
```python
conversions = {
    "Sqrt": "np.sqrt",
    "ArcSin": "np.arcsin",
    "Sinh": "np.sinh",
    "Power": lambda a, b: f"({a}**{b})",
    # ... etc
}
```

## Benefits
Supports more complex coefficient expressions (e.g., damping with `tanh`, special functions).
""",
        "labels": ["enhancement", "priority: medium", "mathematica"],
    },
    {
        "title": "Add Non-Cartesian Coordinate System Examples",
        "body": """## Current State
Position-dependent coefficients documented as working with non-Cartesian. Only one example: `sphere_kg/` (stereographic coordinates, 2+1D). Unclear if full 3+1D spherical/cylindrical works.

## Missing Examples
1. **Spherical coordinates (3+1D):** Full 4D spacetime in spherical
2. **Cylindrical coordinates (3+1D):** Axisymmetric systems
3. **Polar coordinates (2+1D):** 2D polar grid

## Implementation Notes
- Requires py-pde grid support (already available: PolarSymGrid, SphericalSymGrid)
- Wolfram derivation uses standard xCoba charts
- Main challenge: Visualization in non-Cartesian grids

## Benefits
Validates coordinate-independence, enables realistic physics (e.g., black hole spacetimes).
""",
        "labels": ["documentation", "priority: medium", "examples", "coordinates"],
    },
    {
        "title": "Add Convergence and Stability Stress Tests",
        "body": """## Missing Test Categories
1. **Convergence studies:** Verify solution converges as grid/dt refined
2. **Stability bounds:** Test extreme parameters (m→0, m≫1, strong coupling)
3. **CFL violation:** Ensure solver fails gracefully with too-large dt
4. **Grid extremes:** Very large (10k+ points), very small (8-16 points), anisotropic (1000x10)

## Example Test
```python
def test_convergence_klein_gordon():
    \"\"\"Verify KG solution converges to analytical solution as grid refined.\"\"\"
    for N in [32, 64, 128, 256]:
        grid = CartesianGrid([[0, 10]], N)
        # Run simulation, compare to analytical
        # Assert error decreases as O(N^-2)
```

## Benefits
Numerical reliability validation, user confidence in results.
""",
        "labels": ["testing", "priority: medium", "validation", "convergence"],
    },
    {
        "title": "Create Architecture Diagrams",
        "body": """## Current State
Architecture described in text (PIPELINE_README.md). No visual diagrams for module dependencies, data flow, state structure.

## Needed Diagrams
1. **Data flow:** Lagrangian → xAct → JSON → Python → Simulation
2. **Module dependencies:** Which Python modules import which
3. **State structure:** How multi-field state is organized in `FieldCollection`
4. **Execution flow:** From user script to solver to output

## Format
Use diagrams.net (.drawio) or mermaid.js (inline in markdown)

## Target Location
`docs/architecture/` (NEW)

## Benefits
New contributors understand architecture faster, easier onboarding.
""",
        "labels": ["documentation", "priority: medium", "architecture"],
    },
    {
        "title": "Add Full Pipeline Validation to CI",
        "body": """## Problem
Mathematica derivations run locally; JSON committed to repo. No CI check that re-deriving produces identical JSON. Risk: Wolfram code changes make JSON stale.

## Proposed Workflow
1. Optional workflow (manual trigger or label-based)
2. Run `scripts/run_examples.sh` to regenerate all JSON
3. Git diff to compare against committed versions
4. Fail if differences detected (or warn with diff in PR comment)

## Alternative
Weekly scheduled job instead of per-PR (lighter weight).

## Benefits
Ensures JSON stays synchronized with Wolfram code.
""",
        "labels": ["ci-cd", "priority: medium", "validation"],
    },

    # LOW PRIORITY
    {
        "title": "Refactor ContainsTimeDerivative and IsMixedTimeSpaceDerivative",
        "body": """## Opportunity
Unify into single `ExtractDerivativeProfile` function that returns structured data: `{dt, dx, dy, dz}` instead of boolean. Reduces duplication, easier to extend to 4D.

## Location
`torsion_gertsenshtein/wolfram/ExportJSON.wl`
""",
        "labels": ["refactoring", "priority: low", "code-quality"],
    },
    {
        "title": "Add Parameter Sweep Examples",
        "body": """## Feature
Demonstrate batch simulation with varying parameters. Example: Vary mass m from 0.1 to 10.0, plot amplitude decay. Use joblib or multiprocessing for parallelization.

## Target Location
`examples/parameter_sweep/` (NEW)
""",
        "labels": ["documentation", "priority: low", "examples"],
    },
    {
        "title": "Add Python 3.12+ Testing to CI Matrix",
        "body": """## Current State
Python 3.11 pinned (`requires-python = ">=3.11,<3.12"`).

## Recommendation
Add Python 3.12 test job (allow_failure: true) to detect future issues.

## Location
`.github/workflows/test.yml`
""",
        "labels": ["ci-cd", "priority: low", "future-proofing"],
    },
    {
        "title": "Add PGF Export Module Tests",
        "body": """## Problem
PGF export functionality untested (150+ lines, 0 tests). May produce invalid LaTeX.

## Location
`torsion_gertsenshtein/plot_pgf.py`

## Test Requirements
- Validate PGF output is valid LaTeX
- Test figure metadata extraction
- Edge cases: empty figures, very large datasets
""",
        "labels": ["testing", "priority: low", "visualization"],
    },
    {
        "title": "Cache Coordinate Symbols and Chart Dimension",
        "body": """## Optimization
`GetCoordinateSymbols` and `GetChartDimension` called repeatedly. Memoize results per chart to avoid repeated xAct introspection. Estimated speedup: 10-20% for large systems.

## Location
`torsion_gertsenshtein/wolfram/CommonUtilities.wl`
""",
        "labels": ["performance", "priority: low", "optimization", "wolfram"],
    },
    {
        "title": "Add Debugging and Performance Tuning Guides",
        "body": """## Content Needed
- Step-by-step Wolfram derivation debugging
- JSON visualization before simulation
- Performance profiling workflows
- Solver selection guide (RK vs implicit)
- Grid sizing recommendations

## Target Files
- `docs/DEBUGGING.md` (NEW)
- `docs/PERFORMANCE_TUNING.md` (NEW)
""",
        "labels": ["documentation", "priority: low"],
    },
    {
        "title": "Support Elliptic PDE Solving (Constraint Equations)",
        "body": """## Current Limitation
Constraint equations (time_order=0) just return zero evolution. No mechanism to solve elliptic problems (Poisson equation).

## Implementation Path
- Integrate py-pde implicit solvers
- Add example: Electrostatic potential from charge distribution

## Location
`torsion_gertsenshtein/symbolic/pde_builder.py`
""",
        "labels": ["enhancement", "priority: low", "elliptic-pde"],
    },
]


def create_issue(issue_data: dict, *, dry_run: bool = False) -> None:
    """Create a GitHub issue using gh CLI."""
    title = issue_data["title"]
    body = issue_data["body"]
    labels = ",".join(issue_data["labels"])

    if dry_run:
        print(f"\n{'=' * 80}")
        print(f"Title: {title}")
        print(f"Labels: {labels}")
        print(f"\nBody:\n{body}")
        return

    # Create issue using gh CLI
    cmd = [
        "gh", "issue", "create",
        "--title", title,
        "--body", body,
        "--label", labels
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✓ Created: {title}")
        print(f"  URL: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to create: {title}")
        print(f"  Error: {e.stderr}")
    except FileNotFoundError:
        print("Error: gh CLI not found. Install from https://cli.github.com/")
        sys.exit(1)


def main() -> None:
    dry_run = "--dry-run" in sys.argv

    if dry_run:
        print("DRY RUN MODE - No issues will be created")
        print(f"Total issues to create: {len(ISSUES)}")

    print(f"\n{'=' * 80}")
    print("Creating GitHub Issues for Torsion-Gertsenshtein Pipeline")
    print(f"{'=' * 80}\n")

    for i, issue in enumerate(ISSUES, 1):
        print(f"[{i}/{len(ISSUES)}] {issue['title']}")
        create_issue(issue, dry_run=dry_run)

    if not dry_run:
        print(f"\n✓ Successfully created {len(ISSUES)} issues!")
    else:
        print("\nDry run complete. Run without --dry-run to create issues.")


if __name__ == "__main__":
    main()
