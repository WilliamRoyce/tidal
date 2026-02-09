# DevContainer Configuration

This directory contains the complete development container setup for the TIDAL project, with integrated Wolfram Engine 14.3 and xAct tensor computation framework.

## Directory Structure

```
.devcontainer/
├── README.md                 # This file
├── devcontainer.json         # Main VS Code dev container configuration
├── test-noble.json          # Test configuration (legacy)
├── docs/                     # Documentation files
│   ├── WOLFRAM_GUIDE.md     # Complete Wolfram Engine setup and usage guide
│   └── XACT_TESTS.md        # xAct test suite documentation
├── scripts/                  # Utility scripts for setup and maintenance
│   ├── setup_wolfram_engine.sh        # Interactive Wolfram Engine 14.3 setup (NEW)
│   ├── setup_xact.sh                  # Interactive xAct package installation (NEW)
│   ├── validate-setup.sh              # Comprehensive setup validation (NEW)
│   ├── build-xperm.sh                 # Build xPerm from source for current GLIBC
│   ├── check-wolfram.sh               # Health check for Wolfram Engine
│   ├── fix-xperm.sh                   # Fix xPerm installation issues
│   ├── install-extensions-final.sh    # Install VS Code extensions
│   ├── notify-install-extensions.sh   # Extension installation notification
│   └── wolfram-activation-manager.sh  # Manage Wolfram licensing and activation
└── tests/                    # Comprehensive test suite
    ├── test-all-xact.sh              # Master test runner (runs all tests below)
    ├── test-xtensor.wls              # Test xTensor (abstract tensor algebra)
    ├── test-xcoba.wls                # Test xCoba (coordinate-based computations)
    ├── test-xperm.wls                # Test xPerm (permutation algorithms + MathLink)
    ├── test-xpert.wls                # Test xPert (perturbation theory)
    ├── test-integration.wls          # Integration test (Schwarzschild computation)
    ├── test-xact.wls                 # Combined xAct test (legacy)
    └── comprehensive-xact-test.wls   # Comprehensive test suite (legacy)
```

## First-Time Setup

If this is your first time using this devcontainer and you don't have Wolfram Engine installed:

### 1. Install Wolfram Engine 14.3

```bash
bash .devcontainer/scripts/setup_wolfram_engine.sh
```

This interactive wizard will:
- Guide you through downloading the Wolfram Engine installer (~4 GB)
- Run the installer and activate with your Wolfram ID (free)
- Verify installation and create backups

**Time:** ~15-20 minutes (mostly download)

### 2. Install xAct Packages

```bash
bash .devcontainer/scripts/setup_xact.sh
```

This will:
- Download xAct 1.3.0 (~6 MB)
- Install all tensor computation packages
- Optionally compile xPerm for performance

**Time:** ~5 minutes

### 3. Validate Setup

```bash
bash .devcontainer/scripts/validate-setup.sh
```

Runs 9 comprehensive checks and provides a health report.

**Time:** ~3-5 minutes

See [WOLFRAM_GUIDE.md](docs/WOLFRAM_GUIDE.md) for detailed setup instructions.

---

## Quick Start

### Running Tests

Run the complete xAct test suite:

```bash
cd .devcontainer/tests
./test-all-xact.sh
```

Individual test files can be run directly with `wolframscript`:

```bash
cd .devcontainer/tests
wolframscript test-xtensor.wls
wolframscript test-xcoba.wls
wolframscript test-xperm.wls
wolframscript test-xpert.wls
wolframscript test-integration.wls
```

### Wolfram Engine Management

Check Wolfram Engine health:

```bash
bash .devcontainer/scripts/check-wolfram.sh
```

Manage activation and licensing:

```bash
bash .devcontainer/scripts/wolfram-activation-manager.sh status
bash .devcontainer/scripts/wolfram-activation-manager.sh backup
bash .devcontainer/scripts/wolfram-activation-manager.sh restore
```

Install or update VS Code extensions:

```bash
bash .devcontainer/scripts/install-extensions-final.sh
```

## Key Features

### Wolfram Engine Integration

- **Version**: 14.3 (July 31, 2025)
- **Licensing**: Dual system with offline mathpass and cloud activation
- **Persistence**: All licensing and activation data persists across rebuilds
- **MathLink**: Enabled for high-performance computations with xPerm

### xAct Tensor Framework

- **xTensor**: Abstract tensor algebra
- **xCoba**: Coordinate-based computations
- **xPerm**: Advanced permutation algorithms (MathLink enabled)
- **xPert**: Systematic perturbation theory

### VS Code Configuration

- **Extensions**: 12 pre-configured extensions (Python, Ruff, Prettier, GitHub Copilot, etc.)
- **Python Support**: Full debugging, linting, and testing setup
- **Settings**: Format on save, auto-import organization, pytest configuration

## Development Workflow

### First-Time Setup (New Users)
1. **Wolfram Engine Install**: Run `setup_wolfram_engine.sh` (one-time, ~15-20 min)
2. **xAct Install**: Run `setup_xact.sh` (one-time, ~5 min)
3. **Validation**: Run `validate-setup.sh` to verify (~3-5 min)

### Every Container Build
1. **Container Creation**: devcontainer.json runs `onCreateCommand` then `postCreateCommand`
2. **Mount Activation**: Wolfram Engine and xAct automatically available via mounts
3. **Activation Restore**: License and activation automatically restored from backup
4. **Extension Installation**: Manual via `install-extensions-final.sh` (runs once per rebuild)
5. **Ready to Use**: Everything works immediately - no reinstallation needed!

## Important Notes

- **Mounts**: 4 critical mounts ensure Wolfram Engine, licensing, and caches persist:
  - Engine: `~/.local/wolfram/engine/14.3`
  - User base: `~/.local/wolfram/userbase`
  - Machine ID: `/etc/machine-id` (read-only)
  - Cache: `~/.cache/Wolfram`

- **GLIBC Compatibility**: xPerm is compiled from source to match container's GLIBC version (2.36)

- **File Paths**: All scripts and documentation reference files using relative paths (e.g., `.devcontainer/docs/WOLFRAM_GUIDE.md`)

## Documentation

See the comprehensive guides in the `docs/` folder:

- **WOLFRAM_GUIDE.md**: Complete Wolfram Engine setup, licensing, activation, and troubleshooting
- **XACT_TESTS.md**: xAct test suite documentation and usage

## Troubleshooting

If tests fail or Wolfram Engine isn't working:

1. Run health check: `bash .devcontainer/scripts/check-wolfram.sh`
2. Check activation: `bash .devcontainer/scripts/wolfram-activation-manager.sh status`
3. Review logs: `cat .devcontainer/docs/WOLFRAM_GUIDE.md`
4. Restore activation: `bash .devcontainer/scripts/wolfram-activation-manager.sh restore`

## Contact & Updates

For issues or updates related to xAct framework, visit:

- xAct homepage: http://xact.es
- xPerm MathLink: Part of xAct 1.3.0+
