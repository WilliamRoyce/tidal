"""Profiling and timing utilities for simulation performance analysis."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from time import perf_counter
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from collections.abc import Callable

    from pde import FieldCollection


@dataclass
class Timer:
    """Simple timer utility for profiling code sections."""

    t0: float = field(default_factory=perf_counter)
    marks: dict[str, float] = field(
        default_factory=lambda: cast("dict[str, float]", {})
    )

    def mark(self, label: str) -> None:
        """Record the current time under the given label.

        Parameters
        ----------
        label : str
            Identifier for this timing mark.
        """
        self.marks[label] = perf_counter()

    def since(self, label: str) -> float:
        """Return the time elapsed since the given label was marked.

        Parameters
        ----------
        label : str
            The label to measure time from.

        Returns
        -------
        float
            Elapsed time in seconds since the label was marked.
        """
        return perf_counter() - self.marks[label]

    def summary(self) -> dict[str, float]:
        """Return a summary of all marked times.

        Returns
        -------
        dict[str, float]
            Dictionary with keys showing interval labels and deltas from
            previous mark, plus a 'total (since start)' entry.
        """
        out: dict[str, float] = {}
        last = self.t0
        for k, t in self.marks.items():
            out[f"{k} (Δ from prev)"] = t - last
            last = t
        out["total (since start)"] = perf_counter() - self.t0
        return out


def print_summary(d: dict[str, float], prefix: str = "[prof] ") -> None:
    """Log a summary of timing information using the logging module.

    Parameters
    ----------
    d : dict[str, float]
        Dictionary of label-to-duration mappings to log.
    prefix : str, optional
        Prefix for each log line (default: "[prof] ").
    """
    logger = logging.getLogger(__name__)
    for k, v in d.items():
        # use logging formatting (deferred interpolation) instead of f-strings
        logger.info("%s%s: %.3fs", prefix, k, v)


def first_tick_tracker(
    report: dict[str, float], label: str
) -> Callable[[Any, float], dict[str, float]]:
    """Create a tracker callback that records the time of first invocation.

    This is useful for measuring the initialization delay between calling
    pde.solve() and when the first time step actually executes.

    Parameters
    ----------
    report : dict[str, float]
        Dictionary where the timestamp will be stored.
    label : str
        Key under which to store the timestamp in `report`.

    Returns
    -------
    Callable[[Any, float], dict[str, float]]
        A callback function compatible with py-pde's CallbackTracker.
        On first call, stores perf_counter() in report[label].

    Examples
    --------
    >>> prof = {}
    >>> tracker = first_tick_tracker(prof, "t_first_tick")
    >>> # Pass tracker to CallbackTracker, then after solve:
    >>> init_delay = prof["t_first_tick"] - prof["t_call_solve"]
    """
    seen = {"done": False}

    def _cb(_state: FieldCollection, _t: float) -> dict[str, float]:
        if not seen["done"]:
            report[label] = perf_counter()
            seen["done"] = True
        return {}

    return _cb
