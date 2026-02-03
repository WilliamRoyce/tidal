# DevContainer Configuration

This directory contains the complete development container setup for the torsion-gertsenshtein project, with integrated Wolfram Engine 14.3 and xAct tensor computation framework.

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
│   ├── build-xperm.sh       # Build xPerm from source for current GLIBC
│   ├── check-wolfram.sh     # Health check for Wolfram Engine
│   ├── fix-xperm.sh         # Fix xPerm installation issues
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

1. **Setup Phase**: devcontainer.json runs `onCreateCommand` then `postCreateCommand`
2. **Extension Installation**: Manual via `install-extensions-final.sh` (runs once per rebuild)
3. **Wolfram Activation**: Automatically restored from persistent backup if available
4. **xPerm Compilation**: Automatically compiled from source on first use
5. **Testing**: Run `test-all-xact.sh` to verify all functionality

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
