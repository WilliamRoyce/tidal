TIDAL Documentation
====================

**TIDAL** (Tensor Integration and Derivation for Any Lagrangian) is a symbolic physics
pipeline that derives PDEs from Lagrangians using xAct/Mathematica and integrates them
numerically with SUNDIALS (IDA/CVODE), an exact spectral modal solver, symplectic
leapfrog, or scipy. All equation structure comes from symbolic computation —
no physics is hardcoded.

Key Features
------------

* **Symbolic Pipeline**: Lagrangian → Euler-Lagrange → Component Decomposition → JSON → PDE Simulation
* **CLI Tool**: ``tidal`` command with derive, simulate, measure, inspect, list, and validate subcommands
* **Worked Examples**: Spanning 1+1D through 3+1D spacetimes
* **Multi-Field Support**: Scalars, vectors, and rank-3+ tensors with cross-field coupling
* **Curvilinear Coordinates**: Automatic Christoffel symbol computation for non-Cartesian grids
* **Parameter Sweeps**: Override symbolic coefficients at runtime without re-deriving

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

   # Run a pipeline example
   cd examples/scalar_field && bash run.sh

User Guide
----------

.. note::

   **This section is not written yet.** The narrative documentation — pipeline
   internals, CLI reference, JSON schema, solver backends, gauge fixing,
   background fields, inference and troubleshooting — currently lives as LaTeX
   sources under ``docs/tex/`` in the repository and has not yet been migrated
   into this site. Until it is, the API Reference below is the only
   documentation published here.

   Progress is tracked in `issue #416
   <https://github.com/WilliamRoyce/tidal/issues/416>`_.

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   modules

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
