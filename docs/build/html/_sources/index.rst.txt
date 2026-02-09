TIDAL Documentation
====================

**TIDAL** (Tensor Integration and Derivation for Any Lagrangian) is a symbolic physics
pipeline that derives PDEs from Lagrangians using xAct/Mathematica and simulates them
numerically with py-pde. All equation structure comes from symbolic computation —
no physics is hardcoded.

Key Features
------------

* **Symbolic Pipeline**: Lagrangian → Euler-Lagrange → Component Decomposition → JSON → PDE Simulation
* **CLI Tool**: ``tidal`` command with derive, simulate, inspect, list, and validate subcommands
* **18 Working Examples**: Spanning 1+1D through 3+1D spacetimes
* **Multi-Field Support**: Scalars, vectors, and rank-3+ tensors with cross-field coupling
* **Curvilinear Coordinates**: Automatic Christoffel symbol computation for non-Cartesian grids
* **Parameter Sweeps**: Override symbolic coefficients at runtime without re-deriving
* **743 Python Tests + ~108 Wolfram Tests**: Comprehensive validation of both pipeline stages

Quick Start
-----------

Install and run:

.. code-block:: bash

   # Install
   uv sync --all-extras

   # List available equation specifications
   tidal list

   # Simulate a Klein-Gordon field
   tidal simulate examples/data/klein_gordon_1d.json --t-end 20 --ic gaussian

   # Derive equations from a Lagrangian (requires wolframscript)
   tidal derive examples/scalar_field/theory.toml --run

   # Inspect an equation system
   tidal inspect examples/data/chern_simons.json

   # Run a pipeline example directly
   uv run python examples/scalar_field/kg_from_lagrangian.py

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   pipeline
   cli
   examples
   TEX_SUPPORT

.. toctree::
   :maxdepth: 2
   :caption: Legacy Implementation Notes

   2D_COUPLED_IMPLEMENTATION
   ADVANCED_KLEIN_GORDON_IMPLEMENTATION
   RESTRUCTURING_PHASE2

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   modules

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
