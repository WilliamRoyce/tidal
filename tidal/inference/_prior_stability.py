"""Post-hoc prior stability sweep for nested-sampling visualisation.

Nested-sampling chains (PolyChord) only contain accepted samples — points
where the pre-flight stability guard returned ``-inf`` never enter the
chain.  To visualise the unstable parameter region on a corner plot we
sample the prior independently and run only the cheap eigenvalue check
on each draw, saving the rejected ones as a CSV side file.

The stability check is a Fourier-block eigenvalue analysis with no
simulation — empirically ~1 ms per sample on the dark_photon_plasma
4-parameter model.  N=5000 prior draws cost ~5 s on a single thread,
< 0.3% of a typical 30-minute INTR inference run.

Following anesthetic conventions (Handley 2019, JOSS 4(37) 1414) for
visualising inaccessible regions of parameter space.
"""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Callable
    from pathlib import Path


def run_prior_stability_sweep(
    *,
    base_args: Namespace,
    spec_path: Path,
    param_names: list[str],
    prior_transform: Callable[[Any], Any],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    output_path: Path,
    n_samples: int = 5000,
    seed: int = 0,
    quiet: bool = False,
) -> int:
    """Sample the prior and run the stability check on each draw.

    Writes a CSV with columns ``[*param_names, run_status,
    tachyonic_excess]`` to ``output_path``.  Only stability-rejected
    rows have ``run_status='tachyonic'`` and a finite excess; the rest
    are tagged ``'success'``.

    Parameters
    ----------
    base_args : Namespace
        CLI args carrying ``grid_shape``, ``bounds``, ``param`` etc.
        — same setup as the inference pre-flight guard.
    spec_path : Path
        Path to the equation system JSON.
    param_names : list[str]
        Parameter names matching the prior transform output.
    prior_transform : callable
        Same prior transform as nested sampling: takes a uniform
        hypercube vector and returns a parameter vector.
    source, target : tuple[str, ...] | None
        Source / target field names for the conversion channel.
    output_path : Path
        Destination CSV.  Parent must exist.
    n_samples : int
        Number of prior draws.  Defaults to 5000 — dense enough for a
        scatter-overlay visualisation, ~5 s wall time on one thread.
    seed : int
        Reproducible sampling seed.
    quiet : bool
        If False, print the rejection rate at the end.

    Returns
    -------
    int
        Number of rejected (tachyonic) draws written.
    """
    if not source or not target:
        # Stability check is only meaningful for a conversion channel.
        # Write an empty file so callers can still rely on it existing.
        output_path.write_text("", encoding="utf-8")
        return 0

    from tidal.cli._simulate import _parse_params  # pyright: ignore[reportPrivateUsage]
    from tidal.measurement._stability import check_conversion_stability
    from tidal.solver.grid import GridInfo
    from tidal.symbolic._spec_cache import load_spec_cached

    raw_spec = load_spec_cached(spec_path)
    base_p = _parse_params(list(getattr(base_args, "param", []) or []), raw_spec)
    grid_n = int(getattr(base_args, "grid_shape", 256))
    bounds_str = getattr(base_args, "bounds", "0:100")
    if isinstance(bounds_str, str):
        bparts = bounds_str.split(":")
        bounds_tuple: tuple[tuple[float, float], ...] = (
            (float(bparts[0]), float(bparts[1])),
        )
    else:
        bounds_tuple = tuple(bounds_str)
    grid = GridInfo(shape=(grid_n,), bounds=bounds_tuple, periodic=(True,))
    ic_wavevector_str = getattr(base_args, "ic_wavevector", None)
    ic_k: float | None = None
    if ic_wavevector_str:
        try:
            ic_k = float(str(ic_wavevector_str).split(",")[0])
        except (ValueError, IndexError):
            ic_k = None

    rng = np.random.default_rng(seed)
    n_rejected = 0
    n_evaluated = 0

    # The stability check warns on every ill-conditioned eigendecomposition
    # (cond(V) > 1e12).  Across 5000 prior draws on PGT/CDT models this is
    # tens of thousands of identical warnings — suppress them; the
    # conservative=True path correctly handles the ill-conditioned case
    # without us needing the per-call diagnostic.
    import warnings as _warnings

    with (
        output_path.open("w", encoding="utf-8", newline="") as f,
        _warnings.catch_warnings(),
    ):
        _warnings.filterwarnings(
            "ignore",
            message="Stability analysis: cond.V.",
        )
        writer = csv.writer(f)
        writer.writerow(
            [
                *param_names,
                "run_status",
                "tachyonic_excess",
                "stability_profile",
                "borderline_stability",
            ],
        )
        for _ in range(n_samples):
            cube = list(rng.random(len(param_names)))
            try:
                theta = list(prior_transform(cube))
            except Exception:  # noqa: BLE001, S112
                continue
            params = {**base_p, **dict(zip(param_names, theta, strict=False))}
            try:
                stability = check_conversion_stability(
                    raw_spec,
                    grid,
                    params,
                    source=source[0],
                    target=target[0],
                    ic_wavevector=ic_k,
                )
            except Exception:  # noqa: BLE001, S112
                continue
            n_evaluated += 1
            if stability.stable:
                writer.writerow(
                    [
                        *theta,
                        "success",
                        float(stability.max_excess),
                        stability.profile_name,
                        int(bool(stability.borderline)),
                    ],
                )
            else:
                n_rejected += 1
                writer.writerow(
                    [
                        *theta,
                        "tachyonic",
                        float(stability.max_excess),
                        stability.profile_name,
                        0,
                    ],
                )

    if not quiet:
        if n_evaluated:
            pct = 100.0 * n_rejected / n_evaluated
            print(
                f"  Prior stability sweep: {n_rejected}/{n_evaluated} rejected "
                f"({pct:.1f}%); saved to {output_path.name}"
            )
        else:
            print("  Prior stability sweep: no valid prior draws.")

    return n_rejected
