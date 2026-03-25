# Contributing to TIDAL

Thank you for your interest in contributing to the TIDAL project! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Development Workflow](#development-workflow)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

This project follows a standard code of conduct:

- Be respectful and inclusive
- Welcome newcomers and help them learn
- Focus on constructive feedback
- Maintain professional communication

## Getting Started

### Prerequisites

- Python 3.11 (required for `tomllib` stdlib + tested configuration)
- [uv](https://github.com/astral-sh/uv) package manager
- Git
- A GitHub account

### Finding Issues to Work On

- Check the [issue tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues)
- Look for issues labeled `good first issue` or `help wanted`
- Comment on an issue to let others know you're working on it

## Development Setup

### 1. Fork and Clone

```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/torsion-gertsenshtein.git
cd torsion-gertsenshtein
```

### 2. Set Up Development Environment

```bash
# Pin Python version
uv python pin 3.11

# Install dependencies including dev tools
uv sync --all-extras

# Install pre-commit hooks (optional but recommended)
uv run pre-commit install
```

### 3. Verify Setup

```bash
# Run tests (1,700 Python tests)
uv run pytest

# Run linter
uv run ruff check .

# Run type checker
uv run pyright

# Verify CLI works
tidal list
tidal validate examples/data/klein_gordon_1d.json

# Run a pipeline example
cd examples/scalar_field && bash run.sh
```

### Dev Container (Alternative)

If you use VS Code, you can use the dev container:

1. Install "Dev Containers" extension
2. Open project in VS Code
3. Press F1 → "Dev Containers: Reopen in Container"

## Development Workflow

### 1. Create a Branch

```bash
# Always create a new branch for your changes
git checkout -b feature/your-feature-name
# or
git checkout -b fix/issue-description
```

Branch naming conventions:

- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation improvements
- `refactor/` - Code refactoring
- `test/` - Test additions/improvements

### 2. Make Changes

- Write clear, focused commits
- Follow the coding standards (see below)
- Add tests for new functionality
- Update documentation as needed

### 3. Test Your Changes

```bash
# Run all tests
uv run pytest

# Run with coverage
uv run pytest --cov=tidal --cov-report=html

# Run specific test file
uv run pytest tests/test_your_feature.py

# Run linter and formatter
uv run ruff check . --fix
uv run ruff format .

# Type check
uv run pyright
```

### 4. Commit Changes

```bash
git add .
git commit -m "feat: add new feature description"
```

Commit message format:

- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `test:` - Test additions/changes
- `refactor:` - Code refactoring
- `perf:` - Performance improvements
- `chore:` - Build/tooling changes

### 5. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a pull request on GitHub.

## Coding Standards

### Style Guide

This project uses strict code quality tools:

- **Ruff**: Linting and formatting (configured in `pyproject.toml`)
- **Pyright**: Type checking in strict mode
- **NumPy docstring convention**: For all public APIs

### Type Hints

All code must include type hints:

```python
from __future__ import annotations  # Required at top of every file

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

def my_function(x: float, items: Sequence[int]) -> bool:
    """Docstring here."""
    ...
```

### Docstrings

Use NumPy docstring format:

```python
def example_function(param1: float, param2: str) -> int:
    """
    Brief one-line description.

    Longer description if needed, explaining the function's purpose,
    behavior, and any important details.

    Parameters
    ----------
    param1 : float
        Description of param1.
    param2 : str
        Description of param2.

    Returns
    -------
    int
        Description of return value.

    Raises
    ------
    ValueError
        When parameter validation fails.

    Examples
    --------
    >>> example_function(1.0, "test")
    42
    """
    ...
```

### Dataclasses for Configuration

Use frozen dataclasses with validation:

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class MyConfig:
    """Configuration for X."""

    value: float = 1.0

    def __post_init__(self) -> None:
        """Validate configuration."""
        if self.value <= 0:
            msg = "value must be positive"
            raise ValueError(msg)
```

### Import Organization

```python
from __future__ import annotations

# Standard library
import logging
from dataclasses import dataclass

# Third-party
import numpy as np

# Local
from tidal.symbolic.json_loader import load_equation_system

# Type-only imports
if TYPE_CHECKING:
    from collections.abc import Sequence
```

### Error Messages

Always use descriptive error messages:

```python
# Good
if x <= 0:
    msg = f"x must be positive, got {x}"
    raise ValueError(msg)

# Bad
if x <= 0:
    raise ValueError("invalid")
```

## Testing

### Writing Tests

- Place tests in `tests/` directory
- Use `test_*.py` naming convention
- Use fixtures from `tests/conftest.py` when possible

Example test:

```python
from __future__ import annotations

import numpy as np


def test_my_feature() -> None:
    """Test that my feature works correctly."""
    # Arrange
    expected_value = 42

    # Act
    result = my_feature()

    # Assert
    assert result == expected_value
```

### Test Coverage

- Aim for >80% code coverage
- All new features must include tests
- All bug fixes should include regression tests

```bash
# Generate coverage report
uv run pytest --cov=tidal --cov-report=html
# View in browser: open htmlcov/index.html
```

### Test Categories

- **Unit tests**: Test individual functions/classes in isolation
- **Integration tests**: Test component interactions
- **Smoke tests**: Basic "does it run" tests
- **Edge case tests**: Boundary conditions and error handling

## Documentation

### Code Documentation

- All public APIs must have docstrings
- Include examples in docstrings when helpful
- Document parameters, return values, and exceptions

### User Documentation

If adding new features, update:

- `README.md` - If it affects setup or basic usage
- `docs/source/` - Add tutorial or example if appropriate
- `CHANGELOG.md` - Document changes under "Unreleased"

### Building Docs Locally

```bash
cd docs
uv run sphinx-apidoc -f -o source/ ../tidal/
uv run make html
# Open build/html/index.html in browser
```

## Submitting Changes

### Pull Request Checklist

Before submitting a PR, ensure:

- [ ] All tests pass (`uv run pytest`)
- [ ] Code passes linting (`uv run ruff check .`)
- [ ] Code passes type checking (`uv run pyright`)
- [ ] New features include tests
- [ ] Documentation is updated
- [ ] CHANGELOG.md is updated
- [ ] Commits follow conventional commit format
- [ ] Branch is up to date with main

### Pull Request Description

Include in your PR description:

1. **Summary**: What does this PR do?
2. **Motivation**: Why is this change needed?
3. **Changes**: List of key changes
4. **Testing**: How was this tested?
5. **Breaking Changes**: Any backwards-incompatible changes?
6. **Related Issues**: Links to related issues

Example:

```markdown
## Summary

Add support for 3D Klein-Gordon simulations

## Motivation

Users have requested 3D simulation capabilities for modeling more complex systems.

## Changes

- Extended GridConfig to support dim=3
- Added 3D visualization utilities
- Updated documentation with 3D examples

## Testing

- Added 15 new tests for 3D functionality
- All existing tests pass
- Ran performance benchmarks

## Breaking Changes

None

## Related Issues

Closes #42
```

### Review Process

1. Automated checks run on your PR (tests, linting, type checking)
2. Maintainers review your code
3. Address any feedback
4. Once approved, maintainers will merge

### After Your PR is Merged

- Your contribution will be acknowledged in the changelog
- Delete your feature branch
- Pull latest main: `git checkout main && git pull upstream main`

## Getting Help

- **Questions**: Open a [discussion](https://github.com/WilliamRoyce/torsion-gertsenshtein/discussions)
- **Bugs**: Open an [issue](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues)
- **Feature Requests**: Open an issue with the `enhancement` label

## Project Structure

Understanding the codebase:

```
tidal/                   # TIDAL project root
├── tidal/               # Main package
│   ├── symbolic/             # Lagrangian-to-PDE pipeline (Python side)
│   │   ├── json_loader.py   # Load equations from JSON → EquationSystem
│   │   ├── _eval_utils.py   # Mathematica→Python expression conversion
│   │   └── __init__.py
│   ├── solver/               # PDE time-integration backends
│   │   ├── ida.py           # SUNDIALS IDA (DAE solver)
│   │   ├── cvode.py         # SUNDIALS CVODE (adaptive BDF ODE)
│   │   ├── leapfrog.py      # Störmer-Verlet symplectic integrator
│   │   ├── modal.py         # Fourier modal solver (eigendecomposition)
│   │   ├── fields.py        # FieldSet typed container
│   │   ├── coefficients.py  # CoefficientEvaluator (4-level cache)
│   │   ├── rhs.py           # RHSEvaluator (operator+coefficient application)
│   │   ├── operators.py     # Pure numpy spatial operators (FD stencils)
│   │   ├── grid.py          # GridInfo minimal grid dataclass
│   │   ├── state.py         # StateLayout (field→slice mapping)
│   │   ├── validation.py    # SpecValidator (CFL, mass sign, dimensions)
│   │   ├── constraint_solve.py  # Three-tier constraint pre-solve
│   │   └── progress.py     # Simulation progress bar (tqdm)
│   ├── cli/                  # CLI (`tidal` command, 9 subcommands)
│   │   ├── __init__.py      # Entry point + argument parsing
│   │   ├── _derive.py       # tidal derive: TOML → .wls → wolframscript
│   │   ├── _simulate.py     # tidal simulate: JSON → PDE → solve → plot
│   │   ├── _measure.py      # tidal measure: post-hoc analysis
│   │   ├── _inspect.py      # tidal inspect: display equation system info
│   │   ├── _list.py         # tidal list: discover available JSON specs
│   │   ├── _validate.py     # tidal validate: JSON spec validation
│   │   ├── _plot.py         # Plotting utilities for simulate
│   │   ├── _plot_command.py # tidal plot: standalone plotting
│   │   ├── _sweep.py        # tidal sweep: parameter sweeps
│   │   └── _analyze.py      # tidal analyze: sensitivity analysis
│   ├── wolfram/              # Mathematica/xAct pipeline modules
│   │   ├── EulerLagrange.wl
│   │   ├── ComponentDecompose.wl
│   │   ├── ExportJSON.wl
│   │   ├── GaugeFix.wl
│   │   ├── CommonUtilities.wl
│   │   └── ...
│   └── measurement/          # Post-hoc analysis (energy, conversion, spectra)
├── tests/                    # Test suite (1,701 Python tests)
│   ├── conftest.py          # Shared fixtures
│   ├── test_cli.py          # CLI integration tests
│   ├── test_solver_ida.py   # IDA solver tests
│   ├── test_solver_leapfrog.py  # Leapfrog tests
│   └── test_*.py            # Other test modules
├── examples/                 # 20 working pipeline examples
│   ├── data/                # Generated JSON specifications
│   └── {example}/           # Each has theory.toml, run.sh
├── docs/                     # Documentation
│   ├── tex/                 # LaTeX technical docs (25 files)
│   ├── figures/             # TikZ diagrams (17 files)
│   └── source/              # Sphinx API source files
└── pyproject.toml           # Project configuration
```

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to TIDAL! 🎉
