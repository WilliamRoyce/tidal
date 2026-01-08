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
        """Record the current time under the given label."""
        self.marks[label] = perf_counter()

    def since(self, label: str) -> float:
        """Return the time elapsed since the given label was marked."""
        return perf_counter() - self.marks[label]

    def summary(self) -> dict[str, float]:
        """Return a summary of all marked times."""
        out: dict[str, float] = {}
        last = self.t0
        for k, t in self.marks.items():
            out[f"{k} (Δ from prev)"] = t - last
            last = t
        out["total (since start)"] = perf_counter() - self.t0
        return out


def print_summary(d: dict[str, float], prefix: str = "[prof] ") -> None:
    """Log a summary of timing information (use logging, not print)."""
    logger = logging.getLogger(__name__)
    for k, v in d.items():
        # use logging formatting (deferred interpolation) instead of f-strings
        logger.info("%s%s: %.3fs", prefix, k, v)


def first_tick_tracker(
    report: dict[str, float], label: str
) -> Callable[[Any, float], dict[str, float]]:
    """
    Tracker that stores the time of the first callback under `report[label]`.

    You'll set report['t_call_solve'] right before pde.solve(...).
    """
    seen = {"done": False}

    def _cb(_state: FieldCollection, _t: float) -> dict[str, float]:
        if not seen["done"]:
            report[label] = perf_counter()
            seen["done"] = True
        return {}

    return _cb
