# Scripts

Utility scripts for project setup and maintenance.

## Wolfram Engine Setup

The project uses Wolfram Engine for symbolic computation (xCoba tensor algebra).

### Prerequisites

1. **Download the Wolfram Engine installer** from [wolfram.com/engine](https://www.wolfram.com/engine/)
2. Place the downloaded `.sh` file in the `third_party/` directory
3. The installer should be named like: `WolframEngine_14.3.0_LIN.sh`

### Installation

**No lifecycle hook installs Wolfram.** Container creation wires up an engine that is already
on the mount; installing it is a manual, one-time step. The full path — six steps, ending at
the certified configuration — is in
[`.devcontainer/docs/WOLFRAM_GUIDE.md`](../.devcontainer/docs/WOLFRAM_GUIDE.md), which is the
single source. This file documents the individual scripts.

```bash
# Installs into $HOME/.local/wolfram/engine/14.3 -- the dev container's bind mount.
# No sudo: only the system-library step escalates, and the container has those already.
./scripts/install-wolfram-engine.sh

# Or if the installer is already in third_party/:
./scripts/install-wolfram-engine.sh --skip-download
```

### Activation

After installation, you must activate Wolfram Engine with a free Wolfram ID:

```bash
# Interactive activation (recommended)
./scripts/activate-wolfram.sh

# Or directly:
wolframscript -activate
```

You will be prompted for your Wolfram ID credentials. If you don't have one, create a free account at [account.wolfram.com](https://account.wolfram.com/).

## xAct & xCoba Setup

After Wolfram Engine is installed and activated, install xAct for tensor algebra:

```bash
# Install xAct and xCoba packages
./scripts/install-xact-xcoba.sh
```

This script:

- Downloads the official xAct package suite
- Recompiles xPerm binary for GLIBC compatibility
- Installs to the Wolfram user Applications directory
- Verifies installation with a test

### xAct Usage Examples

```wolfram
Needs["xAct`xCoba`"];
DefManifold[M, 4, IndexRange[a, z]];
DefChart[spherical, M, {0, 1, 2, 3}, {t[], r[], θ[], φ[]}];
```

### Verification

Check that everything is working:

```bash
# Check activation status
./scripts/activate-wolfram.sh --check

# Or test directly:
wolframscript -code "1+1"
# Should output: 2

# Test symbolic computation:
wolframscript -code "Integrate[x^2, {x, 0, 1}]"
# Should output: 1/3
```

### Troubleshooting

**"Wolfram Engine is not activated"**

- Run `wolframscript -activate` and enter your Wolfram ID credentials

**"wolframscript not found"**

- Ensure installation completed successfully
- Check that `/usr/local/bin` is in your PATH

**License limitations**

- Free Wolfram Engine license has 2GB memory limit
- For larger computations, consider a commercial license

**GLIBC compatibility errors (xPerm)**

- If you see a `GLIBC_… not found` error, the xPerm binary shipped in the tarball is newer
  than this image's GLIBC (2.36) and must be rebuilt from source
- Run `bash .devcontainer/scripts/build-xperm.sh` — the route that produced the certified
  binary. It runs `mprep` on xAct's own `xperm.tm`, compiles, and installs a wrapper that puts
  MathLink's shared libraries on `LD_LIBRARY_PATH`

## Verification

Run the comprehensive verification script to check all components:

```bash
./scripts/verify-wolfram-setup.sh
```

Eleven checks: the binary on `PATH`; activation, version, license and — critically — that
`$InstallationDirectory` is the **mounted** engine rather than a cloud fallback; the userbase;
the xAct packages and their four-version fingerprint; the xPerm binary; loading each package;
the PSALTer install and its pinned commit; PSALTer's own functions; a headless PDF smoke test;
the certified identities of the two registered resource functions on master and subkernel; and
an xAct smoke test.

```bash
./scripts/verify-wolfram-setup.sh --require-psalter   # the certification gate: exit 0
```

Exit **0** = certified, **2** = works but degraded (not certifiable), **1** = broken.

### Smoke Test

Run the xAct/xCoba smoke test directly:

```bash
wolframscript -file scripts/xact_smoke.wl
```

Expected output includes:

- Package loading messages
- Manifold and chart definitions
- Metric tensor definition
- Riemann tensor antisymmetry verification
- "SMOKE TEST PASSED" message

## Development Utility Scripts

Scripts for local development and testing workflows.

### Version Management

Update version numbers across all project files atomically:

```bash
# Interactive mode (prompts for new version)
python scripts/bump_version.py

# Direct mode - specify version
python scripts/bump_version.py 0.3.0

# Preview changes without modifying files
python scripts/bump_version.py 0.3.0 --dry-run

# Automatically create git commit after updating
python scripts/bump_version.py 0.3.0 --commit

# Allow running with uncommitted changes
python scripts/bump_version.py 0.3.0 --allow-dirty
```

The script updates version numbers in:
- `pyproject.toml` - Package version metadata
- `tidal/__init__.py` - Python module `__version__`
- `CITATION.cff` - Citation metadata (version + release date)
- `uv.lock` - Regenerated via `uv lock` command

**Features:**
- ✅ Validates semantic versioning format (X.Y.Z or X.Y.Z-suffix)
- ✅ Creates `.bak` backups before modification
- ✅ Automatic rollback on failure
- ✅ Detects and reports version inconsistencies
- ✅ Updates CITATION.cff release date automatically
- ✅ Optional git commit with conventional format

**Example workflow:**
```bash
# Check current version and preview changes
python scripts/bump_version.py 0.3.0 --dry-run

# Apply changes
python scripts/bump_version.py 0.3.0

# Review and commit
git diff
git add -A && git commit -m "chore: bump version to 0.3.0"
git tag v0.3.0
git push && git push --tags
```

### Testing

```bash
# Run all Wolfram unit tests
./scripts/run_wolfram_tests.sh

# Run full test suite (Python + Wolfram)
./scripts/full_test.sh
```

### Regenerating Equations

```bash
# Regenerate all JSON equation files from Lagrangians
./scripts/run_examples.sh
```

### Validation

```bash
# End-to-end pipeline validation (derive → JSON → simulate)
./scripts/validate_pipeline.sh

# Check Wolfram module syntax (no tests, just load verification)
./scripts/lint_wolfram.sh
```

## Files

| Script                      | Purpose                                           |
| --------------------------- | ------------------------------------------------- |
| `install-wolfram-engine.sh` | Downloads and installs Wolfram Engine             |
| `activate-wolfram.sh`       | Helps with license activation                     |
| `install-xact-xcoba.sh`     | Installs xAct/xCoba with GLIBC compatibility      |
| `install-psalter.sh`        | Installs PSALTer at a pinned commit               |
| `install-xpand.sh`          | Installs xPand 0.4.4, the FRW derivation engine (D-C) |
| `verify-wolfram-setup.sh`   | Comprehensive verification of all components      |
| `xact_smoke.wl`             | Wolfram Language smoke test for xAct/xCoba        |
| `psalter_smoke.wl`          | PSALTer smoke test, incl. headless PDF export     |
| `psalter/`                  | PSALTer Tier-1 install gate, probes and diff tooling — **see `psalter/README.md`** |
| `oracles/`                  | Frozen legacy oracle (M0.5, #525) — **see `oracles/README.md`** |
| `bump_version.py`           | Atomic version updates across project files       |
| `run_wolfram_tests.sh`      | Run all Wolfram unit tests                        |
| `run_examples.sh`           | Regenerate JSON files from example derivations    |
| `full_test.sh`              | Run complete test suite (Python + Wolfram)        |
| `validate_pipeline.sh`      | End-to-end pipeline validation                    |
| `lint_wolfram.sh`           | Check Wolfram module syntax                       |

Two subdirectories are listed by directory rather than by file, because both carry their own
README that stays current as their contents change:

- **`scripts/psalter/`** — `run_tier1_gate.sh` (the gate), `tier1_diff.wls`,
  `summarize_diff.py` (**read-only** readback of a recorded verdict — the gate's `--diff-only`
  path recomputes and would overwrite it), `extract_checkpoints.py`, `stamp_lines.py`,
  `vector_smoke.wls`, and the `probe_52*.wls` scripts behind #521–#523.
- **`scripts/oracles/`** — `freeze_legacy_oracle.py`, the **one place allowed to touch
  legacy**. Regenerates or verifies the 185 committed fixtures under
  `tests_cosmo/data/oracles/`. `make oracle-check` runs the verify path; CI runs it
  path-filtered via `.github/workflows/oracle.yml`.

**Standing rule:** if `tidal/` or `examples/data/` changes, re-run `scripts/oracles/` in the
same commit, so the frozen oracle keeps describing what legacy actually produces.

## Environment Variables

| Variable              | Default                            | Description                       |
| --------------------- | ---------------------------------- | --------------------------------- |
| `WOLFRAM_VERSION`     | `14.3.0`                           | Wolfram Engine version to install |
| `WOLFRAM_INSTALL_DIR` | `$HOME/.local/wolfram/engine`      | Install root (the bind mount)     |
| `XACT_VERSION`        | `1.3.0`                            | xAct version to install           |
| `PSALTER_COMMIT`      | `bb45adb0…` (v2.0.2)               | PSALTer revision to install       |
| `QT_QPA_PLATFORM`     | unset                              | Set to `offscreen` before running PSALTer |

`QT_QPA_PLATFORM` is not a preference. PSALTer exports a PDF through the Wolfram
front end every time a field is declared, and the call is neither guarded nor
time-limited. When no Qt platform plugin can be initialized, the front end aborts
and that call **blocks indefinitely** rather than failing — which, inside a long
run, is indistinguishable from PSALTer merely being slow, and holds the
single-license Wolfram lane open forever. `offscreen` is the plugin most likely to
be satisfiable on a bare container. `scripts/install-psalter.sh` and the Tier-1
gate set it themselves; anything else you write must too.

## Container Rebuild Behavior

On container rebuild:

1. **initializeCommand** (on the host): creates the bind-mount source directories
2. **postCreateCommand**: installs system dependencies, then runs
   `.devcontainer/scripts/setup-wolfram-links.sh` to wire the engine, license and
   `WolframScript.conf`, and finally `scripts/psalter/ensure_registered.sh` to re-register
   PSALTer's two resource functions
3. **Manual**: `./scripts/verify-wolfram-setup.sh --require-psalter` for full verification

There is no post-attach hook, and nothing installs Wolfram, xAct or PSALTer for you.

**Note**: activation is **not** per-container. `mathpass` is written into the bind-mounted
userbase, so you activate once per machine, not once per rebuild. What the resource registry
under `~/.Wolfram` needs after a rebuild is re-registration, which step 2 above does.
