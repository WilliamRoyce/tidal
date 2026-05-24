# ruff: noqa: N803, N806
"""Phase E — Wavepacket-transit stability diagnostics.

Naming exception: the dimensional quantities ``sigB``, ``L``,
``A_plateau_*``, ``A_t_check_*`` follow the physics conventions in
``scripts/hpc_submit_drafts/v3e_localised/_geometry.env`` and the
underlying literature on Gertsenshtein conversion. Per-file noqa lifts
the lowercase rule so identifiers stay one-to-one with the env file.


Post-simulation analysis of a localised wavepacket passing through a
localised B-field region. Distinguishes between

* a *finite-interaction-time* growth (tachyonic-during-transit; the
  conversion factor saturates once the wavepacket exits the field) and

* a *catastrophic* instability (growth in regions the wavepacket never
  reached; theory ill-posed in vacuum).

Four metrics, each O(grid · snapshots) and sub-second:

1. Sup-norm of the dynamical state at t_check_2 (overflow flag).
2. Wavepacket-norm conservation: ratio of post-transit ‖h‖² in the
   exit window to the IC ‖h‖² in the entry window. Catastrophic
   instabilities push this ratio ≫ 1.
3. Pre-arrival vacuum check: ‖h‖² in a window the wavepacket NEVER
   crosses (between the wrap-around tail and the first B-field peak).
   A clean run keeps this at numerical floor; a catastrophic instability
   pumps energy into vacuum.
4. A-plateau ratio between two checkpoint times: in a finite-interaction
   regime, the conversion amplitude saturates so the ratio approaches 1.
   A persistent monotonic growth flags background-source pumping.

This module does NOT decide PASS / SOFT-PENALIZED / CATASTROPHIC by
itself — it returns the raw numbers and a struct of pass/fail flags
relative to user thresholds. The likelihood pipeline maps the flags
onto soft penalties.

The geometry parameters here (``zc1``, ``sigB``, ``x_c``, ``sigma_w``,
``t_check_1``, ``t_check_2``) live in
``scripts/hpc_submit_drafts/v3e_localised/_geometry.env`` and must be
forwarded by the simulate / sample CLI as runtime params.
"""

from __future__ import annotations

import dataclasses
import json
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclasses.dataclass(frozen=True)
class PhaseETransitResult:
    """Outputs of the four Phase E wavepacket-transit diagnostics.

    All four `*_pass` booleans must be True for a sample to be considered
    fully clean. Otherwise the likelihood pipeline applies a soft penalty
    proportional to the worst violation.
    """

    # Diagnostic 1 — sup-norm overflow guard
    sup_norm_t_check_2: float
    sup_norm_pass: bool

    # Diagnostic 2 — wavepacket norm conservation
    norm_ic: float
    norm_post_transit: float
    norm_ratio: float
    norm_ratio_pass: bool

    # Diagnostic 3 — pre-arrival vacuum check
    vacuum_norm: float
    vacuum_floor: float  # threshold below which vacuum is "clean"
    vacuum_pass: bool

    # Diagnostic 4 — A plateau test
    A_t_check_1: float
    A_t_check_2: float
    A_plateau_ratio: float
    A_plateau_pass: bool

    # Aggregate verdict (string for human consumption)
    verdict: str  # "PASS", "SOFT-PENALIZED", "CATASTROPHIC"

    def to_dict(self) -> dict[str, Any]:
        raw = dataclasses.asdict(self)
        # Coerce numpy scalar types to native Python so json.dumps works.
        out: dict[str, Any] = {}
        for key, value in raw.items():
            coerced = value.item() if isinstance(value, np.generic) else value
            if isinstance(coerced, bool):
                out[key] = bool(coerced)
            elif isinstance(coerced, (int, float)):
                out[key] = float(coerced)
            else:
                out[key] = coerced
        return out


def _integrate_norm_window(
    field_snapshot: NDArray[np.float64],
    z_coords: NDArray[np.float64],
    z_lo: float,
    z_hi: float,
) -> float:
    """Trapezoidal ∫|h(z)|² dz over [z_lo, z_hi].

    `field_snapshot` is 1D along z (multi-D arrays must be reduced to a
    z-axis cut by the caller).
    """
    mask = (z_coords >= z_lo) & (z_coords <= z_hi)
    if not np.any(mask):
        return 0.0
    z_sub = z_coords[mask]
    f_sub = field_snapshot[mask]
    return float(np.trapezoid(f_sub * f_sub, z_sub))


def _z_coords(data: SimulationData) -> NDArray[np.float64]:
    """Return the 1D z-coordinate array for a Phase E simulation.

    Phase E uses 1+1D (single spatial axis); for multi-D runs the first
    axis is taken to be z by convention.
    """
    (zmin, zmax) = data.grid_bounds[0]
    n_pts = data.fields[next(iter(data.fields))].shape[1]
    return np.linspace(zmin, zmax, n_pts, dtype=np.float64)


def _project_to_z(snapshot_slice: NDArray[np.float64]) -> NDArray[np.float64]:
    """Reduce a snapshot field-slice to a 1D z-cut.

    For 1D grids this is a no-op. For higher-D grids we take the sup over
    non-z axes — robust to which transverse mode is excited.
    """
    if snapshot_slice.ndim == 1:
        return snapshot_slice
    axes = tuple(range(1, snapshot_slice.ndim))
    return np.max(np.abs(snapshot_slice), axis=axes)


def _find_snapshot_index(times: NDArray[np.float64], target: float) -> int:
    """Return the snapshot index closest to `target` time."""
    return int(np.argmin(np.abs(times - target)))


def compute_transit_diagnostics(  # noqa: PLR0913, PLR0914
    data: SimulationData,
    *,
    source_field: str,
    target_field: str,
    zc1: float,
    sigB: float,
    x_c: float,
    sigma_w: float,
    t_check_1: float,
    t_check_2: float,
    h0: float,
    sup_norm_limit: float = 1e6,
    norm_ratio_limit: float = 10.0,
    vacuum_floor_factor: float = 1e-3,
    A_plateau_tol: float = 0.05,
) -> PhaseETransitResult:
    """Compute the four Phase E wavepacket-transit diagnostics.

    Parameters
    ----------
    data
        SimulationData from a Phase E run (must have snapshots at or
        near t_check_1 and t_check_2 — `--snapshots 2` from the
        ``feedback_snapshots_mandatory`` rule guarantees this).
    source_field
        Field name carrying the incoming wavepacket (typically a TT
        metric slot, e.g. ``"h_5"``).
    target_field
        Field name carrying the conversion product (typically a photon
        slot, e.g. ``"a_1"``).
    zc1, sigB
        Position and width of the first B-field peak. The post-transit
        window starts at ``zc1 + 3·sigB`` (right edge of the field).
    x_c, sigma_w
        Wavepacket centre and envelope width at t=0.
    t_check_1, t_check_2
        Check times in code units.
    h0
        IC amplitude (for normalising the vacuum-floor threshold).
    sup_norm_limit
        Above this, flag catastrophic overflow.
    norm_ratio_limit
        Above this, wavepacket norm has blown up.
    vacuum_floor_factor
        Threshold multiplier on ``h0²``; vacuum-region norm above
        ``vacuum_floor_factor · h0² · (vacuum window width)`` is flagged.
    A_plateau_tol
        Required closeness of A(t_check_2)/A(t_check_1) to 1.
    """
    times = data.times
    z = _z_coords(data)
    L = float(z[-1] - z[0])

    # ------------------------------------------------------------------
    # Snapshots at the two check times.
    # ------------------------------------------------------------------
    idx_1 = _find_snapshot_index(times, t_check_1)
    idx_2 = _find_snapshot_index(times, t_check_2)
    src_field = data.fields[source_field]  # (n_snap, *grid)
    tgt_field = data.fields[target_field]
    src_t2 = _project_to_z(src_field[idx_2])
    src_ic = _project_to_z(src_field[0])
    tgt_t1 = _project_to_z(tgt_field[idx_1])
    tgt_t2 = _project_to_z(tgt_field[idx_2])

    # ------------------------------------------------------------------
    # 1. Sup-norm overflow check.
    # ------------------------------------------------------------------
    sup_norm_t2 = float(np.max(np.abs(src_t2)))
    sup_norm_pass = sup_norm_t2 < sup_norm_limit and np.isfinite(sup_norm_t2)

    # ------------------------------------------------------------------
    # 2. Wavepacket-norm conservation.
    # ------------------------------------------------------------------
    # IC window: ±3·sigma_w around x_c.
    ic_lo, ic_hi = x_c - 3.0 * sigma_w, x_c + 3.0 * sigma_w
    # Post-transit window: from right edge of first B-field to L/2 minus
    # the envelope tail, so we capture the wavepacket but not the second
    # B-field region (~3·sigB before zc2).
    post_lo = zc1 + 3.0 * sigB
    post_hi = max(post_lo + sigma_w, 0.5 * L - sigma_w)
    norm_ic = _integrate_norm_window(src_ic, z, ic_lo, ic_hi)
    norm_post = _integrate_norm_window(src_t2, z, post_lo, post_hi)
    norm_ratio = norm_post / norm_ic if norm_ic > 0 else float("nan")
    norm_ratio_pass = np.isfinite(norm_ratio) and norm_ratio < norm_ratio_limit

    # ------------------------------------------------------------------
    # 3. Pre-arrival vacuum check.
    # ------------------------------------------------------------------
    # Region the wavepacket NEVER crosses: between the wrap-around tail
    # (extends from x_c+3*sigma_w by ~v·t_check_2 forward) and zc1-3*sigB
    # on the entry side, AND the region that the back of the wavepacket
    # has already vacated. We pick a simple, robust band: 2*sigma_w wide,
    # immediately to the left of the first B-field, but to the LEFT of
    # where the wavepacket back lives at t_check_2.
    # Assume group velocity c=1; back of packet at t=t_check_2 is at
    # roughly x_c - sigma_w + t_check_2 (envelope back edge advanced).
    back_at_t2 = x_c - sigma_w + t_check_2
    vac_lo = back_at_t2 + sigma_w
    vac_hi = zc1 - 3.0 * sigB
    if vac_hi - vac_lo < sigma_w:
        # Degenerate band — pick the smallest legal slab.
        vac_lo = max(0.0, zc1 - 3.0 * sigB - 2.0 * sigma_w)
        vac_hi = zc1 - 3.0 * sigB
    vacuum_norm = _integrate_norm_window(src_t2, z, vac_lo, vac_hi)
    vacuum_window_width = max(vac_hi - vac_lo, 1e-12)
    vacuum_floor = vacuum_floor_factor * h0 * h0 * vacuum_window_width
    vacuum_pass = vacuum_norm < vacuum_floor

    # ------------------------------------------------------------------
    # 4. A-plateau test.
    # ------------------------------------------------------------------
    # A is proportional to the target-field peak amplitude. Without an
    # independent P_GR denominator here, we use the *ratio* of target
    # peaks at the two check times — which equals A(t2)/A(t1) modulo a
    # constant that cancels.
    peak_t1 = float(np.max(np.abs(tgt_t1)))
    peak_t2 = float(np.max(np.abs(tgt_t2)))
    A_ratio = peak_t2 / peak_t1 if peak_t1 > 0 else float("inf") if peak_t2 > 0 else 1.0
    A_plateau_pass = abs(A_ratio - 1.0) < A_plateau_tol

    # ------------------------------------------------------------------
    # Aggregate verdict.
    # ------------------------------------------------------------------
    if not (sup_norm_pass and np.isfinite(norm_ratio)) or not (
        vacuum_pass and norm_ratio_pass
    ):
        verdict = "CATASTROPHIC"
    elif not A_plateau_pass:
        verdict = "SOFT-PENALIZED"
    else:
        verdict = "PASS"

    return PhaseETransitResult(
        sup_norm_t_check_2=sup_norm_t2,
        sup_norm_pass=sup_norm_pass,
        norm_ic=norm_ic,
        norm_post_transit=norm_post,
        norm_ratio=norm_ratio,
        norm_ratio_pass=norm_ratio_pass,
        vacuum_norm=vacuum_norm,
        vacuum_floor=vacuum_floor,
        vacuum_pass=vacuum_pass,
        A_t_check_1=peak_t1,
        A_t_check_2=peak_t2,
        A_plateau_ratio=A_ratio,
        A_plateau_pass=A_plateau_pass,
        verdict=verdict,
    )


def write_stability_json(result: PhaseETransitResult, path: Path) -> None:
    """Serialize a transit result to stability.json (Phase E convention)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(result.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
    )
