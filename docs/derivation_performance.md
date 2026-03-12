# Derivation Performance Notes

## Spherical 4D Einstein-Maxwell + TT Gauge (theory_radial.toml)

Total wall time: ~40 minutes. Peak memory: ~840 MB.

### Phase breakdown

| Phase | Time (est.) | Memory delta | Notes |
|-------|-------------|-------------|-------|
| xPert L^(2) expansion | ~3 min | +5 MB | Fast |
| h decomp (6 components) | **~25 min** | +690 MB | **DOMINANT** |
| — BatchedTBD per component | ~2 min | +3-5 MB | 1000-2300 terms each |
| — Christoffel (first) | ~30s | +285 MB | Cache miss; subsequent ~5s |
| a decomp (4 components) | ~10 min | +50 MB | Rank-1, smaller |
| ExportJSON | ~1 min | minimal | |

### Bottlenecks

1. **BatchedTraceBasisDummyWithMetric** (~20 min total for 10 components)
   - Fused TraceBasisDummy + Expand + MetricEval in adaptive batch groups
   - Generates 1000-2300 intermediate terms per component
   - Single largest time sink — see "Adaptive batch sizing" below

2. **Memory growth** (154->840 MB during h decomp)
   - `Share[]` between components helps but doesn't fully reclaim

3. **Christoffel first-call cost** (+285 MB, ~30s)
   - Already cached via `GetCachedChristoffels` -- subsequent calls O(1)
   - The 285 MB is the Christoffel array itself (4x4x4 with symbolic entries)

### Flat 4D Einstein-Maxwell (theory.toml)

Total: ~40 minutes. Peak: ~170 MB (Wolfram heap) / ~1.3 GB (OS RSS).

| Phase | Time | Notes |
|-------|------|-------|
| EOM h decomposition (6 comp) | ~14 min | 332 & 226 terms per comp |
| EOM a decomposition (4 comp) | ~1.5 min | 99 terms per comp |
| Linearization total | ~16 min | Includes xPert + gauge |
| **Canonical Lagrangian decomp** | **~24 min** | **34 additive terms, DOMINANT** |
| Legendre + parse | <1s | Trivial |

The canonical Lagrangian decomposition is the dominant cost (24 of 40 min).
This is intrinsic to the full Legendre transform — each of the 34 Lagrangian
terms must be independently decomposed via `DecomposeScalarExpression`.
Previously (~10 min) the canonical pipeline was skipped; now it is always
computed for correct coordinate-invariant energy measurement.

### Wolfram License Constraint

Only ONE wolframscript session at a time (engine license). Running two concurrent
derivations exhausts the license and causes "license error" exit. Always serialize
derivation commands. The license error can also occur mid-derivation if the WLS runs
for too long; the Python post-processing (plane-wave reduction, constraint elimination)
then fails because it checks `ret == 0`.

**Workaround**: If the WLS completes but Python post-processing is skipped due to
license error, run the post-processing manually:

```python
from tidal.symbolic.reduction import reduce_spec, eliminate_degenerate_constraints
import json, tomllib

spec = json.loads(Path('examples/data/file.json').read_text())
with open('examples/.../theory.toml', 'rb') as f:
    config = tomllib.load(f)
reduced = reduce_spec(spec, config['reduction'])
reduced = eliminate_degenerate_constraints(reduced)
Path('examples/data/file.json').write_text(json.dumps(reduced, indent='\t'))
```

## Optimizations

### Timing instrumentation

All major derivation phases print `[TIMING]` lines to stderr with elapsed
seconds. Phases instrumented:

- `Linearization (xPert L^(2) + EOM decomposition)` — overall linearization
- `EOM decomposition (<field>)` — per-field `DecomposeToComponents`
- `Canonical Lagrangian decomposition` — `DecomposeScalarExpression` loop
- `Legendre transform (momenta + H)` — canonical momentum + H = Σπv - L
- `ParseHamiltonianExpression` — structured quadratic term parsing

Each canonical Lagrangian term also reports individual time and memory.

### Derivation caching

`tidal derive` hashes the generated WLS script (SHA-256) and stores it in
the output JSON's `metadata.derivation_hash`. On subsequent runs, if the
TOML config produces an identical WLS script and the output JSON exists with
a matching hash, wolframscript is skipped entirely.

```bash
# First run: full derivation (~40 min for theory_radial.toml)
uv run tidal derive examples/gertsenshtein/theory_radial.toml

# Second run: cache hit, instant
uv run tidal derive examples/gertsenshtein/theory_radial.toml
# → "Derivation cache hit: gertsenshtein_radial.json"

# Force re-derivation
uv run tidal derive examples/gertsenshtein/theory_radial.toml --force-derive
```

The hash captures ALL inputs: TOML fields, parameters, Wolfram pipeline code
paths, gauge config, linearization options, reduction settings. Any change
to the TOML or to the Python code that generates the WLS script invalidates
the cache automatically.

### Adaptive batch sizing in BatchedTraceBasisDummyWithMetric

`BatchedTraceBasisDummyWithMetric` (ComponentDecompose.wl) fuses
`TraceBasisDummy` + `Expand` + metric evaluation in batches to bound peak
memory. Previously used a fixed batch size of 50.

Now uses **adaptive sizing**: the first batch runs at the default size (50),
then the expansion factor (peak intermediate terms / input terms) is measured.
If the expansion is high enough, subsequent batches are **decreased** in size
to keep peak intermediate terms below ~2000:

```
batch_size = min(default, max(5, floor(2000 / expansion_factor)))
```

**Critical**: Batch size is NEVER increased above the default. `Expand[]`
scales super-linearly with batch size due to cross-term interactions during
simplification — empirically, doubling batch size can increase processing
time by 10-20x, not 2x. The default batch size of 50 is near-optimal for
most theories.

- Simple scalar theories (expansion ~1-5x): batch stays at 50 (default)
- Complex tensor theories (expansion ~40-50x): batches decrease to 20-30
- Per-batch diagnostics: `TraceBD`, `Expand`, `Eval` term counts + timing

Ref: ComponentDecompose.wl, `BatchedTraceBasisDummyWithMetric` function.

### EOM vs canonical decomposition (no caching possible)

Investigation confirmed that the EOM pass decomposes the **equations of
motion** (`VarD[L]`) via `DecomposeToComponents`, while the canonical path
decomposes the **Lagrangian** itself via `DecomposeScalarExpression`. These
are fundamentally different algebraic objects — the Legendre transform
requires L_comp, not EOM_comp. The Wolfram kernel state (Christoffels,
metric DownValues, background field DownValues) IS already reused between
passes via `Share[]`.
