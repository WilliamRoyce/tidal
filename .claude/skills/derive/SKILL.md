---
name: derive
description: Run tidal derive with safety checks. Blocks parallel wolframscript. Use when deriving PDEs from a theory TOML.
---

# Safe Wolfram Derivation

## Running wolframscript processes
!`pgrep -a wolframscript 2>/dev/null || echo "none running"`

## Instructions

### Pre-flight checks
1. **Check wolframscript**: If any process is shown above, STOP immediately. Tell the user a wolframscript is already running (single Wolfram Engine license — only ONE at a time).
2. **Validate TOML** at the path given in $ARGUMENTS:
   - Constant names must NOT contain underscores (`B0_peak` → corrupts equations silently)
   - Check `[[background_fields]]` and `[[gauge]]` sections are well-formed
   - Verify the TOML file exists

### Run derivation
```bash
uv run tidal derive $ARGUMENTS
```
The default 10-minute timeout (--timeout 600) applies automatically. Set the Bash tool timeout to 600000ms to match.

If derivation times out: **do NOT increase the timeout or change the physics.** Instead, investigate how to optimize the Wolfram pipeline code itself — identify which stage is the bottleneck (decomposition, basis transformation, canonical pipeline) and optimize the algorithm, caching, or code path.

### Post-flight
1. Verify JSON output was created in `examples/data/`
2. **Validate JSON structure** — check operators, field references, parameters, and stability:
```bash
uv run tidal validate <json_path> --stability
```
3. **Check canonical.hamiltonian_terms exists** — this is CRITICAL. Without hamiltonian_terms, ALL energy conservation measurements and energy-based analyses silently fail. Check with:
```bash
python3 -c "import json; d=json.load(open('<json_path>')); h=d.get('canonical',{}).get('hamiltonian_terms',[]); print(f'{len(h)} hamiltonian terms'); assert len(h)>0, 'MISSING hamiltonian_terms — canonical pipeline failed or was skipped!'"
```
If missing, the derivation's canonical pipeline likely failed — do NOT proceed to simulation.
4. Run a quick smoke simulation to confirm the JSON runs:
```bash
uv run tidal simulate <json_path> --grid-shape 64 --bounds 0:10 --periodic \
  --ic gaussian --t-end 1.0 --output /tmp/tidal_derive_smoke
```
5. Report: derivation time, JSON file path, number of fields, hamiltonian terms count, solver auto-selected

### Multiple TOMLs
If the user asks to derive multiple examples, run them SEQUENTIALLY in a loop. Never use parallel execution or background processes for wolframscript.