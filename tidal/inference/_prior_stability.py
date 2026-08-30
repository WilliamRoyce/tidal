"""Post-hoc prior stability sweep for nested-sampling visualization.

Nested-sampling chains (PolyChord) only contain accepted samples — points
where the pre-flight stability guard returned ``-inf`` never enter the
chain.  To visualize the unstable parameter region on a corner plot we
sample the prior independently and run only the cheap eigenvalue check
on each draw, saving the rejected ones as a CSV side file.

The stability check is a Fourier-block eigenvalue analysis with no
simulation — empirically ~1 ms per sample on the dark_photon_plasma
4-parameter model.  N=5000 prior draws cost ~5 s on a single thread,
< 0.3% of a typical 30-minute INTR inference run.

Following anesthetic conventions (Handley 2019, JOSS 4(37) 1414) for
visualizing inaccessible regions of parameter space.
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
        scatter-overlay visualization, ~5 s wall time on one thread.
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

    from tidal.measurement._run_stages import (
        parse_bounds,
        parse_grid_shape,
        parse_params,
    )
    from tidal.measurement._stability import check_conversion_stability
    from tidal.solver.grid import GridInfo
    from tidal.symbolic._spec_cache import load_spec_cached

    raw_spec = load_spec_cached(spec_path)
    base_p = parse_params(list(getattr(base_args, "param", []) or []), raw_spec)
    # Resolve the grid the way the simulation does (GH #479).  This used
    # to read ``int(getattr(base_args, "grid_shape", 256))``, and since
    # ``--grid-shape`` defaults to None in every subparser the getattr
    # default never fired: it was ``int(None)``, a TypeError, caught by
    # the caller at tidal/cli/_sample.py and reported as "Prior stability
    # sweep skipped".  Every `tidal sample` run without an explicit
    # --grid-shape therefore produced no _rejected_prior.csv and no
    # prior-stability overlay on its corner plot.
    spatial_dim = raw_spec.spatial_dimension
    raw_shape = getattr(base_args, "grid_shape", None)
    shape = parse_grid_shape(
        str(raw_shape) if raw_shape is not None else None, spatial_dim
    )
    raw_bounds = getattr(base_args, "bounds", None)
    if raw_bounds is None or isinstance(raw_bounds, str):
        bounds_list = parse_bounds(raw_bounds, spatial_dim)
    else:
        bounds_list = [tuple(b) for b in raw_bounds]
    grid = GridInfo(shape=(shape[0],), bounds=(bounds_list[0],), periodic=(True,))
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
