# Scripts

Utility scripts for project setup and maintenance.

## Wolfram Engine Setup

The project uses Wolfram Engine for symbolic computation (xCoba tensor algebra).

### Prerequisites

1. **Download the Wolfram Engine installer** from [wolfram.com/engine](https://www.wolfram.com/engine/)
2. Place the downloaded `.sh` file in the `third_party/` directory
3. The installer should be named like: `WolframEngine_14.3.0_LIN.sh`

### Installation

#### Automatic (Dev Container)

If you place the installer in `third_party/` before building the dev container, it will be installed automatically during container creation.

#### Manual Installation

```bash
# Run the installation script (requires sudo)
sudo ./scripts/install-wolfram-engine.sh

# Or if installer is already downloaded:
sudo ./scripts/install-wolfram-engine.sh --skip-download
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

- If you see `GLIBC_2.38 not found`, the xPerm binary needs recompilation
- Run `./scripts/install-xact-xcoba.sh` which handles this automatically

## Verification

Run the comprehensive verification script to check all components:

```bash
./scripts/verify-wolfram-setup.sh
```

This checks:

- Wolfram Engine installation and activation
- xAct package installation (xCore, xPerm, xTensor, xCoba)
- xPerm binary compatibility
- Full smoke test with tensor operations

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
| `verify-wolfram-setup.sh`   | Comprehensive verification of all components      |
| `xact_smoke.wl`             | Wolfram Language smoke test for xAct/xCoba        |
| `run_wolfram_tests.sh`      | Run all Wolfram unit tests                        |
| `run_examples.sh`           | Regenerate JSON files from example derivations    |
| `full_test.sh`              | Run complete test suite (Python + Wolfram)        |
| `validate_pipeline.sh`      | End-to-end pipeline validation                    |
| `lint_wolfram.sh`           | Check Wolfram module syntax                       |

## Environment Variables

| Variable              | Default                            | Description                       |
| --------------------- | ---------------------------------- | --------------------------------- |
| `WOLFRAM_VERSION`     | `14.3.0`                           | Wolfram Engine version to install |
| `WOLFRAM_INSTALL_DIR` | `/usr/local/Wolfram/WolframEngine` | Installation directory            |
| `XACT_VERSION`        | `1.2.1`                            | xAct version to install           |

## Container Rebuild Behavior

On container rebuild:

1. **postCreateCommand**: Installs system dependencies including build tools for xPerm recompilation
2. **postAttachCommand**: Checks Wolfram activation, auto-installs xAct/xCoba if missing
3. **Manual**: Run `./scripts/verify-wolfram-setup.sh` for full verification

**Note**: Wolfram Engine activation is per-container and needs to be redone after each rebuild.
