"""Tests for ``tidal simulate --resume`` checkpoint restart functionality."""

from __future__ import annotations

import json
import math
from typing import TYPE_CHECKING

import numpy as np
import pytest

from tidal.cli._simulate import (
    ResumeState,
    _load_resume_state,
    _validate_resume_grid,
)

if TYPE_CHECKING:
    from pathlib import Path

# ==================== Helper: create a synthetic snapshot directory ====================


def _make_snapshot_dir(
    tmp_path: Path,
    *,
    fields: list[str] | None = None,
    velocities: list[str] | None = None,
    grid_shape: tuple[int, ...] = (32,),
    grid_bounds: list[list[float]] | None = None,
    periodic: list[bool] | None = None,
    n_snapshots: int = 11,
    t_end: float = 5.0,
    parameters: dict[str, float] | None = None,
    bc_types: list[str] | None = None,
    dt: float | None = None,
) -> Path:
    """Create a synthetic snapshot directory with metadata and .npy files."""
    if fields is None:
        fields = ["phi_0"]
    if velocities is None:
        velocities = ["phi_0"]
    if grid_bounds is None:
        grid_bounds = [[0.0, 10.0]] * len(grid_shape)
    if periodic is None:
        periodic = [True] * len(grid_shape)
    if parameters is None:
        parameters = {"m2": 1.0}
    if bc_types is None:
        bc_types = ["periodic"] * len(grid_shape)

    out = tmp_path / "snapshot_dir"
    out.mkdir()

    # Write metadata
    meta = {
        "version": 1,
        "n_snapshots": n_snapshots,
        "grid_shape": list(grid_shape),
        "grid_bounds": grid_bounds,
        "grid_spacing": [
            (b[1] - b[0]) / s for b, s in zip(grid_bounds, grid_shape, strict=True)
        ],
        "periodic": periodic,
        "parameters": parameters,
        "fields": fields,
        "velocities": velocities,
        "dtype": "float64",
        "bc_types": bc_types,
        "dt": dt,
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2))

    # Write times
    times = np.linspace(0.0, t_end, n_snapshots)
    np.save(out / "times.npy", times)

    # Write field data (ramp so snapshots are distinguishable)
    rng = np.random.default_rng(42)
    for fname in fields:
        data = rng.standard_normal((n_snapshots, *grid_shape))
        np.save(out / f"{fname}.npy", data)

    for vname in velocities:
        data = rng.standard_normal((n_snapshots, *grid_shape))
        np.save(out / f"v_{vname}.npy", data)

    return out


# ==================== Inline spec fixture ====================


@pytest.fixture
def kg_spec(tmp_path: Path) -> Path:
    """Minimal inline 1D Klein-Gordon spec for resume tests."""
    spec = {
        "metadata": {
            "source": "inline-test",
            "lagrangian_expr": "KG",
            "derived_from": "Euler-Lagrange",
            "gauge": "none",
            "linearized": False,
            "parameters": {"m2": 1.0},
        },
        "spacetime": {
            "dimension": 2,
            "signature": [-1, 1],
            "coordinates": ["t", "x"],
        },
        "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {
                    "expression": "d2_t(phi_0)",
                    "order": {"time": 2, "space": 0},
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "identity",
                            "field": "phi_0",
                            "coefficient_symbolic": "-m2",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi_0",
                        },
                    ],
                },
            },
        ],
    }
    p = tmp_path / "kg_test.json"
    p.write_text(json.dumps(spec, indent=2))
    return p


# ==================== Unit tests: _load_resume_state ====================


class TestLoadResumeState:
    """Unit tests for _load_resume_state()."""

    def test_load_last_snapshot(self, tmp_path: Path, kg_spec: Path) -> None:
        """Loading without snapshot_index returns the last snapshot."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, n_snapshots=11, t_end=5.0)
        spec = load_equation_system(kg_spec)

        state = _load_resume_state(snap_dir, spec)

        assert state.snapshot_index == 10
        assert state.t_start == pytest.approx(5.0)
        assert state.grid_shape == (32,)
        assert state.parameters == {"m2": 1.0}
        assert state.y0.shape[0] > 0

    def test_load_specific_snapshot(self, tmp_path: Path, kg_spec: Path) -> None:
        """Loading a specific snapshot index returns the correct time and data."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, n_snapshots=11, t_end=5.0)
        spec = load_equation_system(kg_spec)

        state = _load_resume_state(snap_dir, spec, snapshot_index=5)

        assert state.snapshot_index == 5
        assert state.t_start == pytest.approx(2.5)

        # Verify the loaded field matches snapshot 5
        raw = np.load(snap_dir / "phi_0.npy")
        expected_field = raw[5]
        # Field should be in the first grid_size elements of y0
        actual_field = state.y0[: np.prod(state.grid_shape)]
        np.testing.assert_array_equal(actual_field, expected_field.ravel())

    def test_snapshot_index_out_of_range(self, tmp_path: Path, kg_spec: Path) -> None:
        """Snapshot index beyond range raises ValueError."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, n_snapshots=5)
        spec = load_equation_system(kg_spec)

        with pytest.raises(ValueError, match="out of range"):
            _load_resume_state(snap_dir, spec, snapshot_index=10)

    def test_missing_metadata_raises(self, tmp_path: Path, kg_spec: Path) -> None:
        """Missing metadata.json raises FileNotFoundError."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path)
        (snap_dir / "metadata.json").unlink()
        spec = load_equation_system(kg_spec)

        with pytest.raises(FileNotFoundError, match=r"metadata\.json"):
            _load_resume_state(snap_dir, spec)

    def test_field_mismatch_raises(self, tmp_path: Path) -> None:
        """Fields missing from checkpoint that spec requires raise ValueError."""
        from tidal.symbolic import load_equation_system

        # Create spec with phi_0, but snapshot only has chi_0
        spec_data: dict[str, object] = {
            "metadata": {
                "source": "test",
                "lagrangian_expr": "KG",
                "derived_from": "Euler-Lagrange",
                "gauge": "none",
                "linearized": False,
                "parameters": {},
            },
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec_path = tmp_path / "mismatch_spec.json"
        spec_path.write_text(json.dumps(spec_data, indent=2))
        spec = load_equation_system(spec_path)

        snap_dir = _make_snapshot_dir(tmp_path, fields=["chi_0"], velocities=["chi_0"])

        with pytest.raises(ValueError, match="missing fields"):
            _load_resume_state(snap_dir, spec)

    def test_parameter_inheritance(self, tmp_path: Path, kg_spec: Path) -> None:
        """Parameters from metadata are preserved in ResumeState."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(
            tmp_path, parameters={"m2": math.pi, "coupling": 0.5},
        )
        spec = load_equation_system(kg_spec)

        state = _load_resume_state(snap_dir, spec)
        assert state.parameters == {"m2": math.pi, "coupling": 0.5}

    def test_bc_types_preserved(self, tmp_path: Path, kg_spec: Path) -> None:
        """BC types from metadata are preserved."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, bc_types=["neumann"])
        spec = load_equation_system(kg_spec)

        state = _load_resume_state(snap_dir, spec)
        assert state.bc_types == ("neumann",)

    def test_load_snapshot_zero(self, tmp_path: Path, kg_spec: Path) -> None:
        """Loading snapshot index 0 returns t=0 and first snapshot data."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, n_snapshots=5, t_end=4.0)
        spec = load_equation_system(kg_spec)

        state = _load_resume_state(snap_dir, spec, snapshot_index=0)

        assert state.snapshot_index == 0
        assert state.t_start == pytest.approx(0.0)

        raw = np.load(snap_dir / "phi_0.npy")
        expected = raw[0]
        actual = state.y0[: np.prod(state.grid_shape)]
        np.testing.assert_array_equal(actual, expected.ravel())

    def test_negative_snapshot_index_raises(
        self, tmp_path: Path, kg_spec: Path,
    ) -> None:
        """Negative snapshot indices raise ValueError."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, n_snapshots=5)
        spec = load_equation_system(kg_spec)

        with pytest.raises(ValueError, match="out of range"):
            _load_resume_state(snap_dir, spec, snapshot_index=-1)

    def test_single_snapshot_directory(self, tmp_path: Path, kg_spec: Path) -> None:
        """Resume from a directory with exactly 1 snapshot."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path, n_snapshots=1, t_end=0.0)
        spec = load_equation_system(kg_spec)

        state = _load_resume_state(snap_dir, spec)

        assert state.snapshot_index == 0
        assert state.t_start == pytest.approx(0.0)
        assert state.y0.shape[0] > 0

    def test_missing_npy_warns(self, tmp_path: Path, kg_spec: Path) -> None:
        """Missing field .npy files produce a warning."""
        from tidal.symbolic import load_equation_system

        snap_dir = _make_snapshot_dir(tmp_path)
        # Delete the velocity file
        (snap_dir / "v_phi_0.npy").unlink()
        spec = load_equation_system(kg_spec)

        with pytest.warns(UserWarning, match="v_phi_0"):
            state = _load_resume_state(snap_dir, spec)

        # Velocity should be zero (default from FieldSet.from_dict)
        n_pts = np.prod(state.grid_shape)
        velocity_part = state.y0[n_pts : 2 * n_pts]
        np.testing.assert_array_equal(velocity_part, 0.0)


# ==================== Unit tests: _validate_resume_grid ====================


class TestValidateResumeGrid:
    """Unit tests for _validate_resume_grid()."""

    def test_matching_grid_passes(self) -> None:
        """Identical grid shape and bounds passes validation."""
        from tidal.solver.grid import GridInfo

        state = ResumeState(
            y0=np.zeros(64),
            t_start=5.0,
            parameters={},
            grid_shape=(32,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            bc_types=None,
            dt=None,
            snapshot_index=0,
        )
        grid = GridInfo(bounds=((0.0, 10.0),), shape=(32,), periodic=(True,))
        _validate_resume_grid(state, grid)  # should not raise

    def test_shape_mismatch_raises(self) -> None:
        """Different grid shapes raise ValueError."""
        from tidal.solver.grid import GridInfo

        state = ResumeState(
            y0=np.zeros(64),
            t_start=5.0,
            parameters={},
            grid_shape=(32,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            bc_types=None,
            dt=None,
            snapshot_index=0,
        )
        grid = GridInfo(bounds=((0.0, 10.0),), shape=(64,), periodic=(True,))

        with pytest.raises(ValueError, match="Grid shape mismatch"):
            _validate_resume_grid(state, grid)

    def test_bounds_mismatch_raises(self) -> None:
        """Different grid bounds raise ValueError."""
        from tidal.solver.grid import GridInfo

        state = ResumeState(
            y0=np.zeros(64),
            t_start=5.0,
            parameters={},
            grid_shape=(32,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            bc_types=None,
            dt=None,
            snapshot_index=0,
        )
        grid = GridInfo(bounds=((0.0, 20.0),), shape=(32,), periodic=(True,))

        with pytest.raises(ValueError, match="Grid bounds mismatch"):
            _validate_resume_grid(state, grid)


# ==================== Integration tests: CLI roundtrip ====================


class TestSimulateResume:
    """Integration tests for --resume flag through the CLI."""

    def test_resume_roundtrip(self, inline_kg_1d_json: Path, tmp_path: Path) -> None:
        """Run a simulation, then resume from checkpoint and verify continuity."""
        from tidal.cli import main

        # Step 1: Run initial simulation to t=5
        ckpt_dir = tmp_path / "checkpoint"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--t-end",
                "5.0",
                "--grid-shape",
                "32",
                "--output",
                str(ckpt_dir),
                "--quiet",
            ],
        )
        assert ret == 0
        assert ckpt_dir.is_dir()

        # Step 2: Resume from checkpoint to t=10
        resume_dir = tmp_path / "resumed"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--resume",
                str(ckpt_dir),
                "--t-end",
                "10.0",
                "--output",
                str(resume_dir),
                "--quiet",
            ],
        )
        assert ret == 0
        assert resume_dir.is_dir()

        # Step 3: Verify the resumed simulation starts at t=5
        resumed_times = np.load(resume_dir / "times.npy")
        assert resumed_times[0] == pytest.approx(5.0, abs=0.01)
        assert resumed_times[-1] == pytest.approx(10.0, abs=0.01)

        # Step 4: Verify field + velocity continuity — last snapshot of
        # checkpoint matches first snapshot of resumed run (exact disk I/O)
        ckpt_field = np.load(ckpt_dir / "phi_0.npy")[-1]
        resumed_field = np.load(resume_dir / "phi_0.npy")[0]
        np.testing.assert_allclose(ckpt_field, resumed_field, atol=1e-10)

        ckpt_vel = np.load(ckpt_dir / "v_phi_0.npy")[-1]
        resumed_vel = np.load(resume_dir / "v_phi_0.npy")[0]
        np.testing.assert_allclose(ckpt_vel, resumed_vel, atol=1e-10)

    def test_resume_specific_snapshot(
        self, inline_kg_1d_json: Path, tmp_path: Path,
    ) -> None:
        """Resume from a specific mid-run snapshot."""
        from tidal.cli import main

        # Run initial simulation
        ckpt_dir = tmp_path / "checkpoint"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--t-end",
                "5.0",
                "--grid-shape",
                "32",
                "--output",
                str(ckpt_dir),
                "--quiet",
            ],
        )
        assert ret == 0

        # Resume from snapshot 5 (mid-run)
        resume_dir = tmp_path / "resumed"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--resume",
                str(ckpt_dir),
                "--snapshot",
                "5",
                "--t-end",
                "10.0",
                "--output",
                str(resume_dir),
                "--quiet",
            ],
        )
        assert ret == 0

        # Verify start time is at snapshot 5
        ckpt_times = np.load(ckpt_dir / "times.npy")
        resumed_times = np.load(resume_dir / "times.npy")
        assert resumed_times[0] == pytest.approx(ckpt_times[5], abs=0.1)

    def test_resume_t_additional(self, inline_kg_1d_json: Path, tmp_path: Path) -> None:
        """--t-additional computes correct absolute t_end."""
        from tidal.cli import main

        # Run to t=5
        ckpt_dir = tmp_path / "checkpoint"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--t-end",
                "5.0",
                "--grid-shape",
                "32",
                "--output",
                str(ckpt_dir),
                "--quiet",
            ],
        )
        assert ret == 0

        # Resume with --t-additional 3 (should run to t=8)
        resume_dir = tmp_path / "resumed"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--resume",
                str(ckpt_dir),
                "--t-additional",
                "3.0",
                "--output",
                str(resume_dir),
                "--quiet",
            ],
        )
        assert ret == 0

        resumed_times = np.load(resume_dir / "times.npy")
        assert resumed_times[-1] == pytest.approx(8.0, abs=0.1)

    def test_resume_nonexistent_dir(self, inline_kg_1d_json: Path) -> None:
        """--resume with nonexistent directory returns error."""
        from tidal.cli import main

        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--resume",
                "/nonexistent/dir",
                "--t-end",
                "10.0",
                "--quiet",
            ],
        )
        assert ret == 1

    def test_resume_parameter_override(
        self, inline_kg_1d_json: Path, tmp_path: Path,
    ) -> None:
        """CLI --param overrides checkpoint metadata parameters."""
        from tidal.cli import main

        # Run with m2=1.0
        ckpt_dir = tmp_path / "checkpoint"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "2.0",
                "--grid-shape",
                "32",
                "--output",
                str(ckpt_dir),
                "--quiet",
            ],
        )
        assert ret == 0

        # Resume with different m2
        resume_dir = tmp_path / "resumed"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--resume",
                str(ckpt_dir),
                "--param",
                "m2=4.0",
                "--t-end",
                "4.0",
                "--output",
                str(resume_dir),
                "--quiet",
            ],
        )
        assert ret == 0

        # Verify the resumed metadata has the overridden parameter
        with (resume_dir / "metadata.json").open() as f:
            meta = json.load(f)
        assert meta["parameters"]["m2"] == pytest.approx(4.0)
