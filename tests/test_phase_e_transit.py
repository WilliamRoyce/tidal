"""Tests for tidal.measurement._phase_e_transit.

Constructs synthetic SimulationData representing the three regimes we
want the diagnostics to distinguish:

* CLEAN — Gaussian wavepacket cleanly transits and exits; vacuum stays
  at zero; conversion saturates.
* CATASTROPHIC — wavepacket plus growing exponential everywhere; the
  vacuum-region check trips.
* PERSISTENT-GROWTH — A(t_check_2) is much larger than A(t_check_1);
  the plateau check trips but the vacuum check does not.
"""

from __future__ import annotations

import dataclasses
import json
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._phase_e_transit import (
    PhaseETransitResult,
    compute_transit_diagnostics,
    write_stability_json,
)

if TYPE_CHECKING:
    from pathlib import Path


@dataclasses.dataclass
class _FakeSim:
    """Minimal SimulationData stand-in for testing.

    Real SimulationData has many more attributes, but compute_transit_diagnostics
    only reads .times, .fields, .grid_bounds — duck-typed access.
    """

    times: np.ndarray
    fields: dict[str, np.ndarray]
    grid_bounds: tuple[tuple[float, float], ...]


# Phase E canonical geometry (matches _geometry.env).
GEOM = {
    "zc1": 50.0,
    "sigB": 5.0,
    "x_c": 15.0,
    "sigma_w": 5.0,
    "t_check_1": 60.0,
    "t_check_2": 80.0,
    "h0": 0.01,
}


def _make_sim(
    src_snapshots: list[np.ndarray],
    tgt_snapshots: list[np.ndarray],
    times: list[float],
    L: float = 200.0,
) -> _FakeSim:
    src = np.stack(src_snapshots, axis=0)
    tgt = np.stack(tgt_snapshots, axis=0)
    return _FakeSim(
        times=np.array(times, dtype=np.float64),
        fields={"h_5": src, "a_1": tgt},
        grid_bounds=((0.0, L),),
    )


def _gaussian(z: np.ndarray, *, centre: float, width: float, amp: float) -> np.ndarray:
    return amp * np.exp(-((z - centre) ** 2) / (2.0 * width**2))


def test_clean_wavepacket_passes() -> None:
    """A clean Gaussian transiting through the B-field gives PASS verdict."""
    z = np.linspace(0.0, 200.0, 512)
    # IC at x_c=15 with h0=0.01
    src_ic = _gaussian(z, centre=GEOM["x_c"], width=GEOM["sigma_w"], amp=GEOM["h0"])
    # Post-transit at t_check_1=60 (group velocity c=1): centre at 15+60=75, well past B-field
    src_t1 = _gaussian(z, centre=75.0, width=GEOM["sigma_w"], amp=GEOM["h0"])
    # Late post-transit at t_check_2=80: centre at 15+80=95
    src_t2 = _gaussian(z, centre=95.0, width=GEOM["sigma_w"], amp=GEOM["h0"])
    # Target conversion saturated post-transit (same A at both checkpoints)
    tgt_t0 = np.zeros_like(z)
    tgt_t1 = _gaussian(z, centre=75.0, width=GEOM["sigma_w"], amp=GEOM["h0"] * 0.01)
    tgt_t2 = _gaussian(z, centre=95.0, width=GEOM["sigma_w"], amp=GEOM["h0"] * 0.01)
    sim = _make_sim(
        [src_ic, src_t1, src_t2], [tgt_t0, tgt_t1, tgt_t2], [0.0, 60.0, 80.0]
    )

    result = compute_transit_diagnostics(
        sim, source_field="h_5", target_field="a_1", **GEOM
    )

    assert result.sup_norm_pass
    assert result.vacuum_pass, (
        f"vacuum_norm={result.vacuum_norm}, floor={result.vacuum_floor}"
    )
    assert result.norm_ratio_pass
    assert result.A_plateau_pass, f"A_ratio={result.A_plateau_ratio}"
    assert result.verdict == "PASS"


def test_catastrophic_vacuum_growth_flagged() -> None:
    """A growing background everywhere trips the vacuum-region check."""
    z = np.linspace(0.0, 200.0, 512)
    src_ic = _gaussian(z, centre=GEOM["x_c"], width=GEOM["sigma_w"], amp=GEOM["h0"])
    # Add a uniform-amplitude growing exponential to t_check_2
    src_t1 = _gaussian(z, centre=55.0, width=GEOM["sigma_w"], amp=GEOM["h0"]) + 0.1
    src_t2 = _gaussian(z, centre=75.0, width=GEOM["sigma_w"], amp=GEOM["h0"]) + 0.5
    tgt = np.zeros_like(z)
    sim = _make_sim([src_ic, src_t1, src_t2], [tgt, tgt, tgt], [0.0, 60.0, 80.0])

    result = compute_transit_diagnostics(
        sim, source_field="h_5", target_field="a_1", **GEOM
    )

    assert not result.vacuum_pass
    assert result.verdict == "CATASTROPHIC"


def test_persistent_growth_soft_penalty() -> None:
    """Wavepacket norm OK but A(t_check_2)/A(t_check_1) ≫ 1 → soft penalty."""
    z = np.linspace(0.0, 200.0, 512)
    src_ic = _gaussian(z, centre=GEOM["x_c"], width=GEOM["sigma_w"], amp=GEOM["h0"])
    src_t1 = _gaussian(z, centre=55.0, width=GEOM["sigma_w"], amp=GEOM["h0"])
    src_t2 = _gaussian(z, centre=75.0, width=GEOM["sigma_w"], amp=GEOM["h0"])
    tgt_t0 = np.zeros_like(z)
    # Target keeps growing in the post-transit window — A_plateau fails
    tgt_t1 = _gaussian(z, centre=55.0, width=GEOM["sigma_w"], amp=GEOM["h0"] * 0.01)
    tgt_t2 = _gaussian(
        z, centre=75.0, width=GEOM["sigma_w"], amp=GEOM["h0"] * 0.05
    )  # 5x growth
    sim = _make_sim(
        [src_ic, src_t1, src_t2], [tgt_t0, tgt_t1, tgt_t2], [0.0, 60.0, 80.0]
    )

    result = compute_transit_diagnostics(
        sim, source_field="h_5", target_field="a_1", **GEOM
    )

    assert result.sup_norm_pass
    assert result.vacuum_pass
    assert result.norm_ratio_pass
    assert not result.A_plateau_pass, f"A_ratio={result.A_plateau_ratio}"
    assert result.verdict == "SOFT-PENALIZED"


def test_sup_norm_overflow_flagged() -> None:
    """A blown-up field (NaN/inf or > sup_norm_limit) flags CATASTROPHIC."""
    z = np.linspace(0.0, 200.0, 512)
    src_ic = _gaussian(z, centre=GEOM["x_c"], width=GEOM["sigma_w"], amp=GEOM["h0"])
    src_t1 = _gaussian(z, centre=55.0, width=GEOM["sigma_w"], amp=GEOM["h0"])
    src_t2 = np.ones_like(z) * 1e12  # overflowed
    tgt = np.zeros_like(z)
    sim = _make_sim([src_ic, src_t1, src_t2], [tgt, tgt, tgt], [0.0, 60.0, 80.0])

    result = compute_transit_diagnostics(
        sim, source_field="h_5", target_field="a_1", **GEOM
    )

    assert not result.sup_norm_pass
    assert result.verdict == "CATASTROPHIC"


def test_write_stability_json_roundtrip(tmp_path: Path) -> None:
    """write_stability_json produces a readable JSON of all fields."""
    z = np.linspace(0.0, 200.0, 512)
    src_ic = _gaussian(z, centre=GEOM["x_c"], width=GEOM["sigma_w"], amp=GEOM["h0"])
    src_t = _gaussian(z, centre=75.0, width=GEOM["sigma_w"], amp=GEOM["h0"])
    tgt = np.zeros_like(z)
    sim = _make_sim([src_ic, src_t, src_t], [tgt, tgt, tgt], [0.0, 60.0, 80.0])
    result = compute_transit_diagnostics(
        sim, source_field="h_5", target_field="a_1", **GEOM
    )

    out = tmp_path / "stability.json"
    write_stability_json(result, out)

    parsed = json.loads(out.read_text())
    expected_keys = {f.name for f in dataclasses.fields(PhaseETransitResult)}
    assert set(parsed) == expected_keys
