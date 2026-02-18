# Gauge Fixing: Optional Per-Field Toolkit

Gauge fixing simplifies equation structure for theories with gauge symmetry
(massless vectors, linearized gravity). It is **never required** — TIDAL's
existing pipeline handles gauge-invariant theories correctly, and all
measurement quantities (energy, conversion, mixing) are gauge-invariant.

## Quick Start

Add a `[[gauge]]` section to your `theory.toml`:

```toml
[[fields]]
name = "A"
type = "vector"

[[gauge]]
field = "A"
type = "lorenz"
xi = 1.0              # optional, default 1.0 (Feynman gauge)
```

This adds the Lorenz gauge-fixing term `-(1/2xi)(d_mu A^mu)^2` to the
Lagrangian before Euler-Lagrange derivation. For Maxwell equations, this
reduces coupled equations to uncoupled wave equations for each component.

No `[[gauge]]` section means no gauge fixing (the default). Existing
examples continue to work unchanged.

## Named Presets

| Preset | Mechanism | Fields | Expression | Effect |
|--------|-----------|--------|------------|--------|
| `lorenz` | lagrangian_term | vector | `-(1/2xi)(d_mu A^mu)^2` | Maxwell -> uncoupled wave equations |
| `de_donder` | lagrangian_term | sym. rank-2 | `-(1/2xi)(d_a h^a_b - 1/2 d_b h)^2` | Lin. Einstein -> uncoupled waves |
| `temporal` | constraint | vector | `A_0 = 0` | Eliminates temporal component |
| `coulomb` | constraint | vector | `div A = 0` | Transversality constraint |
| `axial` | constraint | vector | `A_n = 0` | Eliminates one spatial component |

Presets are syntactic sugar over the expression-based mechanism described
below. Each maps to a builder function in `GaugeFix.wl`.

## Custom Gauge Terms

For full flexibility, use `type = "custom"` with an arbitrary Wolfram
expression:

```toml
# Type A: add a term to the Lagrangian (before EL derivation)
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
expression = "-(1/(2*xi)) * eta[a,b] CD[-a][A[-c]] eta[c,d] CD[-b][A[-d]]"

# Type B: impose a constraint (after EOM derivation)
[[gauge]]
field = "A"
type = "custom"
mechanism = "constraint"
expression = "eta[a,b] CD[-a][A[-b]]"   # set to zero
```

The `expression` field uses the same syntax as `[lagrangian].expression` and
`[[derived_fields]].definition` — field names, `eta`, `CD`, and chart
placeholders are all substituted automatically.

### Worked Example: Lorenz Gauge Two Ways

**Using the preset:**
```toml
[[gauge]]
field = "A"
type = "lorenz"
xi = 1.0
```

**Equivalent custom expression:**
```toml
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
expression = "-(1/2) * eta[a,b] CD[-a][A[-c]] eta[c,d] CD[-b][A[-d]]"
```

Both produce the same gauge-fixed Lagrangian. The preset is more readable;
the custom form shows the underlying mechanism and serves as a template for
non-standard gauge choices.

## Two Mechanisms

**Type A — Lagrangian term** (`mechanism = "lagrangian_term"`):
- An expression is added to L *before* Euler-Lagrange derivation
- Changes the structure of the equations of motion
- Example: Lorenz gauge adds `-(1/2xi)(div A)^2`, turning Maxwell into wave equations

**Type B — Constraint** (`mechanism = "constraint"`):
- A constraint is imposed on the EOM *after* derivation
- Eliminates degrees of freedom without modifying the Lagrangian
- Example: temporal gauge sets `A_0 = 0`, eliminating the temporal component

For named presets, the mechanism is inferred automatically (e.g., `lorenz`
is always `lagrangian_term`, `temporal` is always `constraint`).

## Multi-Field Example

Different fields can have different gauge choices, or no gauge at all:

```toml
[[fields]]
name = "A"
type = "vector"

[[fields]]
name = "B"
type = "vector"

# Lorenz gauge on A only
[[gauge]]
field = "A"
type = "lorenz"

# B has no gauge fixing (default)
```

## Adding New Presets (Developer Guide)

Adding a new named gauge preset is a 3-step process:

### Step 1: Write a builder function in `tidal/wolfram/GaugeFix.wl`

```mathematica
BuildMyGaugeTerm[field_, metric_, covd_, xi_:1] := Module[
  {gaugeTerm},

  If[!NumericQ[xi] || xi <= 0,
    Throw["BuildMyGaugeTerm: xi must be positive, got " <> ToString[xi]]
  ];

  (* Your gauge-fixing expression here *)
  gaugeTerm = (* ... *);

  gaugeTerm
];
```

The function must return an *expression* (a Wolfram symbolic expression),
not perform an action. The pipeline calls `AddGaugeFixingTerm[L, yourTerm]`
to inject it into the Lagrangian.

### Step 2: Add a registry entry in `tidal/cli/_derive.py`

Add to the `_GAUGE_PRESETS` dict:

```python
_GAUGE_PRESETS["my_gauge"] = {
    "mechanism": "lagrangian_term",   # or "constraint"
    "builder": "BuildMyGaugeTerm",    # Wolfram function name
    "requires": "vector",             # or "tensor", "scalar"
}
```

### Step 3: (Optional) Add validation rules

If your gauge has specific requirements beyond field type (e.g., requires
a specific symmetry, extra parameters), add checks in `_validate_gauge()`.

That's it. Users can now write `type = "my_gauge"` in their TOML.

## FAQ

**Do I need gauge fixing?**
No. All existing TIDAL examples work without it. Gauge fixing is a
computational convenience — it simplifies equation structure but doesn't
change the physics. Massless vectors already get implicit constraint
structure (A_0 becomes time_order=0). Massive vectors have no gauge
freedom at all.

**Are measurements affected by gauge choice?**
No. Energy density, conversion probability, mixing length, and spectral
quantities are all gauge-invariant physical observables.

**Can I mix presets and custom gauges?**
Yes. Each `[[gauge]]` entry is independent. You can use `type = "lorenz"`
on one field and `type = "custom"` on another.

**What happens if I gauge-fix a field that has no gauge symmetry?**
For massive fields (Proca), the mass term already breaks gauge symmetry.
Adding a gauge-fixing term is mathematically valid but physically
unnecessary and may produce unexpected equations. Validation prevents
gauge-fixing scalar fields (which never have gauge symmetry).
