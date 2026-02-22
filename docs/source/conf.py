"""Sphinx configuration for TIDAL documentation."""

from importlib.metadata import version as _get_version

project = "TIDAL"
author = "William Royce"
version = _get_version("tidal")
release = version

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",  # Auto-generate summary tables
    "sphinx.ext.doctest",      # Test docstring examples
    "myst_parser",
]

# MyST Markdown settings
myst_enable_extensions = [
    "colon_fence",
    "fieldlist",
]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# Napoleon settings (NumPy docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Intersphinx mappings
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pde": ("https://py-pde.readthedocs.io/en/stable/", None),
}

# Autosummary settings
autosummary_generate = True

# Doctest settings
doctest_global_setup = """
import numpy as np
from tidal.symbolic import load_equation_system
"""

templates_path = ["_templates"]
exclude_patterns = []

# HTML theme
html_theme = "pydata_sphinx_theme"
html_logo = "../TIDAL_Logo_TikZ_Figure.svg"
html_static_path = []

html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/WilliamRoyce/torsion-gertsenshtein",
            "icon": "fa-brands fa-github",
        },
    ],
    "show_toc_level": 2,
}
