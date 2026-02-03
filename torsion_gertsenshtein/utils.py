"""Utility functions."""

from __future__ import annotations

# Type alias for boundary condition data compatible with py-pde
# Can be: string for auto BC, dict for explicit BC mapping, or None
BCDescriptor = str | dict[str, dict[str, str | float]] | None
