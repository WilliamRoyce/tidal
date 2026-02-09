TIDAL Documentation
====================

Klein-Gordon PDE simulations for electromagnetic-gravitational wave conversion research.

Overview
--------

The **tidal** package provides a comprehensive framework for simulating 
Klein-Gordon equations in various configurations, including:

* Single-field and multi-field coupled systems
* Homogeneous and inhomogeneous mass distributions
* Advanced PDE variants (anisotropic, higher-order dispersion)
* Comprehensive visualization and animation tools
* Type-safe configuration system with runtime validation

This package is built on the `py-pde <https://py-pde.readthedocs.io/>`_ framework and 
leverages Numba for high-performance numerical computations.

Quick Start
-----------

Install dependencies and run an example:

.. code-block:: bash

   uv sync --all-extras
   uv run python examples/klein_gordon/1d_gaussian_pulse.py

Features
--------

* **Type-Safe Configuration**: Frozen dataclasses with runtime validation
* **Comprehensive Testing**: 90+ tests with strict type checking
* **Multiple PDE Variants**: Klein-Gordon, coupled fields, inhomogeneous, anisotropic
* **Visualization**: Built-in animation builder with spacetime heatmaps
* **Performance**: Numba JIT compilation for CPU-intensive operations
* **Flexible API**: Both function-based and class-based initial condition APIs

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   examples
   tutorials

.. toctree::
   :maxdepth: 2
   :caption: Implementation Details

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