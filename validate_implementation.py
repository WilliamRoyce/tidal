#!/usr/bin/env python3
"""Validation script to test the Lagrangian-to-PDE pipeline implementation.

This script performs basic sanity checks to ensure all modules import correctly
and the pipeline works end-to-end.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
    PDEFromSpec,
    build_pde_from_json,
    create_initial_state,
    load_equation_system,
)
from torsion_gertsenshtein.vectorfield import (
    ComponentFieldParams,
    ComponentGaussianPulse,
    ComponentPlaneWave,
    create_gaussian_pulse_state,
)

LOGGER = logging.getLogger(__name__)
LOG_SEPARATOR = "=" * 60
EM_COMPONENTS = 2
KG_COMPONENTS = 1
EM_STATE_FIELDS = 4
EM_GRID_DOMAIN = (0, 100)
EM_GRID_RESOLUTION = 64
DEFAULT_DT = 0.01
SHORT_SIM_TIME = 1.0
A0_TOL = 1e-10


def test_imports() -> bool:
    """Test that all modules import correctly."""
    LOGGER.info(LOG_SEPARATOR)
    LOGGER.info("Testing imports...")
    LOGGER.info(LOG_SEPARATOR)

    try:
        _ = (
            ComponentEquation,
            EquationSystem,
            OperatorTerm,
            PDEFromSpec,
            build_pde_from_json,
            create_initial_state,
            load_equation_system,
            ComponentFieldParams,
            ComponentGaussianPulse,
            ComponentPlaneWave,
            create_gaussian_pulse_state,
            CartesianGrid,
            FieldCollection,
            ScalarField,
        )
    except ImportError:
        LOGGER.exception("Import error")
        return False
    else:
        LOGGER.info("✓ symbolic package imports OK")
        LOGGER.info("✓ vectorfield package imports OK")
        LOGGER.info("✓ py-pde imports OK")
        return True


def test_json_loading() -> bool:
    """Test loading JSON equation specifications."""
    LOGGER.info("\n%s", LOG_SEPARATOR)
    LOGGER.info("Testing JSON loading...")
    LOGGER.info(LOG_SEPARATOR)

    try:
        # Test loading EM spec
        em_path = Path("examples/data/em_1d.json")
        if not em_path.exists():
            LOGGER.error("✗ EM JSON file not found: %s", em_path)
            return False

        em_spec = load_equation_system(em_path)
        LOGGER.info("✓ Loaded EM spec: %s components", em_spec.n_components)
        assert em_spec.n_components == EM_COMPONENTS
        assert em_spec.component_names == ("A_0", "A_1")

        # Test loading KG spec
        kg_path = Path("examples/data/klein_gordon_1d.json")
        if not kg_path.exists():
            LOGGER.error("✗ KG JSON file not found: %s", kg_path)
            return False

        kg_spec = load_equation_system(kg_path)
        LOGGER.info("✓ Loaded KG spec: %s component(s)", kg_spec.n_components)
        assert kg_spec.n_components == KG_COMPONENTS
        assert kg_spec.component_names == ("phi",)
    except (OSError, ValueError, TypeError):
        LOGGER.exception("✗ Error loading JSON")
        return False
    else:
        return True


def test_pde_construction() -> bool:
    """Test building PDE from specification."""
    LOGGER.info("\n%s", LOG_SEPARATOR)
    LOGGER.info("Testing PDE construction...")
    LOGGER.info(LOG_SEPARATOR)

    try:
        # Build EM PDE
        em_path = Path("examples/data/em_1d.json")
        em_pde = build_pde_from_json(em_path)
        LOGGER.info("✓ Built EM PDE: %s", type(em_pde).__name__)
        assert em_pde.n_components == EM_COMPONENTS

        # Create grid and initial state
        grid = CartesianGrid([EM_GRID_DOMAIN], EM_GRID_RESOLUTION, periodic=True)
        spec = load_equation_system(em_path)
        state = create_initial_state(grid, spec)
        LOGGER.info("✓ Created initial state: %s fields", len(state))
        assert len(state) == EM_STATE_FIELDS  # A_0, Pi_0, A_1, Pi_1

        # Test evolution rate computation
        rates = em_pde.evolution_rate(state)
        LOGGER.info("✓ Computed evolution rates: %s fields", len(rates))
        assert len(rates) == EM_STATE_FIELDS
    except (OSError, ValueError, TypeError, RuntimeError):
        LOGGER.exception("✗ Error constructing PDE")
        return False
    else:
        return True


def test_vectorfield() -> bool:
    """Test vectorfield utilities."""
    LOGGER.info("\n%s", LOG_SEPARATOR)
    LOGGER.info("Testing vectorfield utilities...")
    LOGGER.info(LOG_SEPARATOR)

    try:
        # Load spec
        em_path = Path("examples/data/em_1d.json")
        spec = load_equation_system(em_path)

        # Create params
        params = ComponentFieldParams.from_equation_system(spec)
        LOGGER.info(
            "✓ Created ComponentFieldParams: %s components",
            params.n_components,
        )
        assert params.is_massless

        # Create Gaussian pulse
        grid = CartesianGrid([EM_GRID_DOMAIN], EM_GRID_RESOLUTION, periodic=True)
        pulse = ComponentGaussianPulse(
            center=(50.0,), width=5.0, active_components={"A_1": 1.0}
        )
        state = pulse.create(grid, spec)
        LOGGER.info("✓ Created Gaussian pulse state: %s fields", len(state))

        # Verify pulse is only in A_1
        assert np.allclose(state[0].data, 0.0)  # A_0 = 0
        assert np.max(state[2].data) > 0  # A_1 has pulse

    except (OSError, ValueError, TypeError, RuntimeError):
        LOGGER.exception("✗ Error in vectorfield")
        return False
    else:
        return True


def test_simulation() -> bool:
    """Test running a short simulation."""
    LOGGER.info("\n%s", LOG_SEPARATOR)
    LOGGER.info("Testing simulation...")
    LOGGER.info(LOG_SEPARATOR)

    try:
        # Load EM PDE
        em_path = Path("examples/data/em_1d.json")
        pde = build_pde_from_json(em_path)

        # Load spec for initial condition
        spec = load_equation_system(em_path)

        # Create grid and initial state
        grid = CartesianGrid([EM_GRID_DOMAIN], EM_GRID_RESOLUTION, periodic=True)
        pulse = ComponentGaussianPulse(
            center=(50.0,), width=5.0, active_components={"A_1": 1.0}
        )
        initial = pulse.create(grid, spec)

        # Run short simulation
        result = pde.solve(initial, t_range=SHORT_SIM_TIME, dt=DEFAULT_DT)
        field_result = result[0] if isinstance(result, tuple) else result
        if field_result is None:
            LOGGER.error("✗ Simulation returned no field data")
            return False

        assert isinstance(field_result, FieldCollection)
        LOGGER.info("✓ Simulation completed: %s fields", len(field_result))

        # Verify result is finite
        assert np.all(np.isfinite(field_result.data))

        # Verify A_0 stayed zero (no coupling)
        assert np.allclose(field_result.data[0], 0.0, atol=A0_TOL)

        LOGGER.info("✓ Simulation produced valid results")
    except (OSError, ValueError, TypeError, RuntimeError):
        LOGGER.exception("✗ Error in simulation")
        return False
    else:
        return True


def main() -> int:
    """Run all validation tests."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    LOGGER.info("\nLagrangian-to-PDE Pipeline Validation")
    LOGGER.info(LOG_SEPARATOR)
    LOGGER.info("")

    tests = [
        ("Imports", test_imports),
        ("JSON Loading", test_json_loading),
        ("PDE Construction", test_pde_construction),
        ("VectorField Utilities", test_vectorfield),
        ("Simulation", test_simulation),
    ]

    results: list[tuple[str, bool]] = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except (OSError, ValueError, TypeError, RuntimeError):
            LOGGER.exception("\n✗ Unexpected error in %s", name)
            results.append((name, False))

    # Summary
    LOGGER.info("\n%s", LOG_SEPARATOR)
    LOGGER.info("Validation Summary")
    LOGGER.info(LOG_SEPARATOR)

    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        symbol = "✓" if passed else "✗"
        LOGGER.info("%s %s: %s", symbol, name, status)

    all_passed = all(passed for _, passed in results)
    LOGGER.info("\n%s", LOG_SEPARATOR)
    if all_passed:
        LOGGER.info("✓ All validation tests PASSED")
        return 0
    LOGGER.info("✗ Some validation tests FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
