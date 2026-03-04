"""Simulation progress bar for TIDAL solvers (tqdm-based).

Provides a thin wrapper around tqdm that tracks simulation time ``t`` and
displays elapsed time, ETA, and percentage completion.  Designed to be
called from solver step loops or RHS functions with negligible overhead.

Key features:

- ``disable=True`` makes tqdm a complete no-op (zero overhead when quiet).
- ``mininterval=0.5`` throttles terminal writes to ~2 Hz.
- Monotonic ``t`` tracking prevents backwards progress from IDA Jacobian
  finite-difference evaluations.
- Output goes to ``sys.stderr`` so it doesn't interfere with piped stdout.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from tqdm import tqdm  # pyright: ignore[reportMissingTypeStubs]

if TYPE_CHECKING:
    from tqdm import tqdm as TqdmType  # noqa: N812


class SimulationProgress:
    """tqdm-based progress bar driven by simulation time.

    Parameters
    ----------
    t_start : float
        Simulation start time.
    t_end : float
        Simulation end time.
    solver_name : str
        Label for the progress bar (e.g. ``"IDA"``, ``"CVODE"``).
    disable : bool
        If ``True``, completely disables output (zero overhead).
    """

    def __init__(
        self,
        t_start: float,
        t_end: float,
        *,
        solver_name: str = "",
        disable: bool = False,
    ) -> None:
        self._t_start = t_start
        self._max_t = t_start
        self._pbar: TqdmType[None] = tqdm(  # type: ignore[assignment]  # tqdm stub overload
            total=t_end - t_start,
            desc=f"  {solver_name}" if solver_name else "  Simulating",
            unit="t",
            bar_format=(
                "{desc}: {percentage:3.0f}%|{bar}| "
                "t={n:.3g}/{total:.3g} [{elapsed}<{remaining}]"
            ),
            disable=disable,
            file=sys.stderr,
            mininterval=0.5,
            leave=True,
        )

    def set_phase(self, phase: str) -> None:
        """Update the progress bar description to show current phase.

        Useful for indicating solver initialization before the step loop
        begins, so users know the solver hasn't hung on large systems.
        """
        self._pbar.set_description(f"  {phase}")
        self._pbar.refresh()

    def update(self, t: float) -> None:
        """Report current simulation time.

        Monotonic: ignores ``t`` values less than the previous maximum,
        which can occur during IDA Jacobian finite-difference evaluations.
        Clamps to total to prevent tqdm warnings when leapfrog overshoots.
        """
        if t > self._max_t:
            self._max_t = t
            n = t - self._t_start
            if self._pbar.total is not None:
                n = min(n, self._pbar.total)
            self._pbar.n = n
            self._pbar.refresh()

    def finish(self) -> None:
        """Set progress to 100% and close the bar."""
        if self._pbar.total is not None:
            self._pbar.n = self._pbar.total
        self._pbar.refresh()
        self._pbar.close()
