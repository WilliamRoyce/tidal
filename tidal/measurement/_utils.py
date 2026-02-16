"""Shared private utilities for the measurement package.

Functions here are used by sibling modules (_conversion, _spectral_conversion)
within the measurement package.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence


def _normalize_group(  # pyright: ignore[reportUnusedFunction]
    fields: str | Sequence[str],
) -> tuple[str, ...]:
    """Normalize a string or sequence to a tuple of field names."""
    if isinstance(fields, str):
        return (fields,)
    return tuple(fields)
