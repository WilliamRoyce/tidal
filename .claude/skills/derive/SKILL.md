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

### Post-flight
1. Verify JSON output was created in `examples/data/`
2. Run a quick smoke simulation to confirm the JSON is valid:
```bash
uv run tidal simulate <json_path> --grid-shape 64 --bounds 0:10 --periodic \
  --ic gaussian --t-end 1.0 --output /tmp/tidal_derive_smoke
```
3. Report: derivation time, JSON file path, number of fields found, solver auto-selected

### Multiple TOMLs
If the user asks to derive multiple examples, run them SEQUENTIALLY in a loop. Never use parallel execution or background processes for wolframscript.