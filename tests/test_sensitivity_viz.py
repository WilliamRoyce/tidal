"""Tests for F4 (sensitivity analysis) and F8 (advanced visualization).

F4: Sobol/Morris sensitivity wrappers (with mocked SALib).
F8: render_sweep_parallel, render_sweep_tornado, render_sweep_scatter.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from tidal.measurement._sweep_results import SweepResults

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


def _make_sweep_results(
    n_params: int = 3,
    n_runs: int = 12,
    *,
    add_status: bool = False,
) -> SweepResults:
    """Build a synthetic SweepResults with *n_params* swept parameters."""
    rng = np.random.default_rng(42)
    param_names = [f"p{i}" for i in range(n_params)]
    swept = {name: sorted(rng.uniform(0, 1, n_runs).tolist()) for name in param_names}

    rows: list[dict[str, Any]] = []
    run_dirs: list[Path] = []
    for i in range(n_runs):
        row: dict[str, Any] = {}
        for name in param_names:
            row[name] = swept[name][i]
        row["P_max"] = rng.uniform(0, 1)
        row["max_energy_error"] = rng.uniform(1e-6, 1e-3)
        if add_status:
            row["run_status"] = "success"
        rows.append(row)
        run_dirs.append(Path(f"/tmp/run_{i}"))

    return SweepResults(
        swept_params=swept,
        fixed_params={"fixed_a": 1.0},
        sim_settings={"t_end": 10.0},
        rows=rows,
        run_dirs=run_dirs,
        spec_path="spec.json",
        measurements=["conversion"],
    )


# ---------------------------------------------------------------------------
# F4: Sensitivity analysis
# ---------------------------------------------------------------------------


class TestExtractData:
    """Tests for _extract_data helper."""

    def test_basic_extraction(self) -> None:
        from tidal.measurement._sensitivity import _extract_data

        results = _make_sweep_results(n_params=2, n_runs=5)
        param_names, x, y = _extract_data(results, "P_max")
        assert param_names == ["p0", "p1"]
        assert x.shape == (5, 2)
        assert y.shape == (5,)

    def test_no_params_raises(self) -> None:
        from tidal.measurement._sensitivity import _extract_data

        results = SweepResults(
            swept_params={},
            fixed_params={},
            sim_settings={},
            rows=[{"P_max": 0.5}],
            run_dirs=[Path("/tmp/r0")],
            spec_path="spec.json",
            measurements=[],
        )
        with pytest.raises(ValueError, match="No swept parameters"):
            _extract_data(results, "P_max")

    def test_nan_rows_filtered(self) -> None:
        from tidal.measurement._sensitivity import _extract_data

        results = SweepResults(
            swept_params={"g0": [0.1, 0.5]},
            fixed_params={},
            sim_settings={},
            rows=[
                {"g0": 0.1, "P_max": 0.3},
                {"g0": 0.5},  # Missing P_max -> NaN -> filtered
            ],
            run_dirs=[Path("/tmp/r0"), Path("/tmp/r1")],
            spec_path="spec.json",
            measurements=[],
        )
        _, x, y = _extract_data(results, "P_max")
        assert len(y) == 1
        assert x.shape == (1, 1)


class TestRequireSalib:
    """Tests for the SALib import guard."""

    def test_missing_salib_raises(self) -> None:
        from tidal.measurement._sensitivity import _require_salib

        with (
            patch.dict("sys.modules", {"SALib": None}),
            pytest.raises(ImportError, match="SALib is required"),
        ):
            _require_salib()


class TestSobolIndices:
    """Tests for compute_sobol_indices with mocked SALib."""

    def test_insufficient_samples_raises(self) -> None:
        from tidal.measurement._sensitivity import compute_sobol_indices

        results = _make_sweep_results(n_params=3, n_runs=4)
        # Need at least 2*(3+1)=8 samples but only providing 4
        mock_sobol = MagicMock()
        with (
            patch("tidal.measurement._sensitivity._require_salib"),
            patch.dict(
                "sys.modules",
                {
                    "SALib.analyze.sobol": mock_sobol,
                    "SALib.analyze": MagicMock(),
                    "SALib": MagicMock(),
                },
            ),
            pytest.raises(ValueError, match="Sobol analysis needs at least"),
        ):
            compute_sobol_indices(results, "P_max")

    def test_sobol_returns_result(self) -> None:
        from tidal.measurement._sensitivity import compute_sobol_indices

        results = _make_sweep_results(n_params=2, n_runs=12)
        mock_si = {
            "S1": np.array([0.5, 0.3]),
            "ST": np.array([0.6, 0.4]),
            "S1_conf": np.array([0.05, 0.04]),
            "ST_conf": np.array([0.06, 0.05]),
        }

        mock_analyze_fn = MagicMock(return_value=mock_si)
        mock_sobol_mod = MagicMock()
        mock_sobol_mod.analyze = mock_analyze_fn
        mock_salib_analyze = MagicMock()
        mock_salib_analyze.sobol = mock_sobol_mod

        with (
            patch("tidal.measurement._sensitivity._require_salib"),
            patch.dict(
                "sys.modules",
                {
                    "SALib": MagicMock(),
                    "SALib.analyze": mock_salib_analyze,
                    "SALib.analyze.sobol": mock_sobol_mod,
                },
            ),
        ):
            result = compute_sobol_indices(results, "P_max", n_bootstrap=50)

        assert result.method == "sobol"
        assert result.metric == "P_max"
        assert result.param_names == ["p0", "p1"]
        assert result.s1 is not None
        np.testing.assert_allclose(result.s1, [0.5, 0.3])
        assert result.st is not None
        np.testing.assert_allclose(result.st, [0.6, 0.4])
        assert result.mu_star is None


class TestMorrisScreening:
    """Tests for compute_morris_screening with mocked SALib."""

    def test_insufficient_samples_raises(self) -> None:
        from tidal.measurement._sensitivity import compute_morris_screening

        results = _make_sweep_results(n_params=5, n_runs=3)
        # Need at least 5+1=6 but only 3
        mock_morris = MagicMock()
        with (
            patch("tidal.measurement._sensitivity._require_salib"),
            patch.dict(
                "sys.modules",
                {
                    "SALib.analyze.morris": mock_morris,
                    "SALib.analyze": MagicMock(),
                    "SALib": MagicMock(),
                },
            ),
            pytest.raises(ValueError, match="Morris screening needs at least"),
        ):
            compute_morris_screening(results, "P_max")

    def test_morris_returns_result(self) -> None:
        from tidal.measurement._sensitivity import compute_morris_screening

        results = _make_sweep_results(n_params=2, n_runs=10)
        mock_si = {
            "mu_star": np.array([1.2, 0.4]),
            "sigma": np.array([0.3, 0.1]),
        }

        mock_analyze_fn = MagicMock(return_value=mock_si)
        mock_morris_mod = MagicMock()
        mock_morris_mod.analyze = mock_analyze_fn
        mock_salib_analyze = MagicMock()
        mock_salib_analyze.morris = mock_morris_mod

        with (
            patch("tidal.measurement._sensitivity._require_salib"),
            patch.dict(
                "sys.modules",
                {
                    "SALib": MagicMock(),
                    "SALib.analyze": mock_salib_analyze,
                    "SALib.analyze.morris": mock_morris_mod,
                },
            ),
        ):
            result = compute_morris_screening(results, "P_max")

        assert result.method == "morris"
        assert result.mu_star is not None
        np.testing.assert_allclose(result.mu_star, [1.2, 0.4])
        assert result.s1 is None


class TestFormatSensitivityTable:
    """Tests for format_sensitivity_table."""

    def test_sobol_table(self) -> None:
        from tidal.measurement._sensitivity import (
            SensitivityResult,
            format_sensitivity_table,
        )

        result = SensitivityResult(
            method="sobol",
            metric="P_max",
            param_names=["g0", "m2"],
            s1=np.array([0.5, 0.3]),
            st=np.array([0.6, 0.4]),
            s1_conf=np.array([0.05, 0.04]),
            st_conf=np.array([0.06, 0.05]),
            mu_star=None,
            sigma=None,
        )
        table = format_sensitivity_table(result)
        assert "Sobol" in table
        assert "P_max" in table
        assert "g0" in table
        assert "m2" in table

    def test_morris_table(self) -> None:
        from tidal.measurement._sensitivity import (
            SensitivityResult,
            format_sensitivity_table,
        )

        result = SensitivityResult(
            method="morris",
            metric="L_mix",
            param_names=["alpha", "beta"],
            s1=None,
            st=None,
            s1_conf=None,
            st_conf=None,
            mu_star=np.array([1.2, 0.4]),
            sigma=np.array([0.3, 0.1]),
        )
        table = format_sensitivity_table(result)
        assert "Morris" in table
        assert "L_mix" in table
        assert "alpha" in table


# ---------------------------------------------------------------------------
# F4: analyze CLI
# ---------------------------------------------------------------------------


class TestAnalyzeCommand:
    """Tests for tidal analyze CLI dispatch."""

    def test_missing_path(self, tmp_path: Path) -> None:
        from tidal.cli._analyze import analyze_command

        args = MagicMock()
        args.data_path = str(tmp_path / "nonexistent")
        assert analyze_command(args) == 1

    def test_auto_detect_metric(self) -> None:
        from tidal.cli._analyze import _auto_detect_metric

        results = _make_sweep_results(n_params=2, n_runs=5)
        metric = _auto_detect_metric(results)
        assert metric == "P_max"

    def test_auto_detect_metric_none(self) -> None:
        from tidal.cli._analyze import _auto_detect_metric

        results = SweepResults(
            swept_params={"g0": [0.1]},
            fixed_params={},
            sim_settings={},
            rows=[{"g0": 0.1, "custom_metric": 1.0}],
            run_dirs=[Path("/tmp/r0")],
            spec_path="spec.json",
            measurements=[],
        )
        assert _auto_detect_metric(results) is None


# ---------------------------------------------------------------------------
# F8: Advanced visualization
# ---------------------------------------------------------------------------


class TestRenderSweepParallel:
    """Tests for render_sweep_parallel."""

    def _make_ax(self) -> Axes:
        import matplotlib.pyplot as plt

        _, ax = plt.subplots()
        return ax

    def test_too_few_params_raises(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_parallel

        results = _make_sweep_results(n_params=1, n_runs=3)
        ax = self._make_ax()
        with pytest.raises(ValueError, match="2\\+ swept parameters"):
            render_sweep_parallel(ax, results, "P_max")

    def test_smoke_renders(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_parallel

        results = _make_sweep_results(n_params=3, n_runs=8)
        ax = self._make_ax()
        render_sweep_parallel(ax, results, "P_max")
        # Verify polylines were drawn (n_runs lines)
        assert len(ax.lines) == 8

    def test_title_contains_metric(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_parallel

        results = _make_sweep_results(n_params=2, n_runs=4)
        ax = self._make_ax()
        render_sweep_parallel(ax, results, "P_max")
        assert "P_max" in ax.get_title()


class TestRenderSweepTornado:
    """Tests for render_sweep_tornado."""

    def _make_ax(self) -> Axes:
        import matplotlib.pyplot as plt

        _, ax = plt.subplots()
        return ax

    def test_no_params_raises(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_tornado

        results = SweepResults(
            swept_params={},
            fixed_params={},
            sim_settings={},
            rows=[{"P_max": 0.5}],
            run_dirs=[Path("/tmp/r0")],
            spec_path="spec.json",
            measurements=[],
        )
        ax = self._make_ax()
        with pytest.raises(ValueError, match="at least 1"):
            render_sweep_tornado(ax, results, "P_max")

    def test_smoke_renders(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_tornado

        results = _make_sweep_results(n_params=3, n_runs=12)
        ax = self._make_ax()
        render_sweep_tornado(ax, results, "P_max")
        assert "P_max" in ax.get_title()

    def test_single_param(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_tornado

        results = _make_sweep_results(n_params=1, n_runs=5)
        ax = self._make_ax()
        render_sweep_tornado(ax, results, "P_max")
        # Should have exactly 1 horizontal bar
        assert len(ax.patches) == 1


class TestRenderSweepScatter:
    """Tests for render_sweep_scatter."""

    def _make_fig(self) -> Figure:
        import matplotlib.pyplot as plt

        return plt.figure()

    def test_too_few_params_raises(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_scatter

        results = _make_sweep_results(n_params=1, n_runs=3)
        fig = self._make_fig()
        with pytest.raises(ValueError, match="2\\+ swept parameters"):
            render_sweep_scatter(fig, results, "P_max")

    def test_smoke_renders_2_params(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_scatter

        results = _make_sweep_results(n_params=2, n_runs=8)
        fig = self._make_fig()
        render_sweep_scatter(fig, results, "P_max")
        # 2x2 grid of axes
        assert len(fig.axes) >= 4

    def test_smoke_renders_3_params(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_scatter

        results = _make_sweep_results(n_params=3, n_runs=10)
        fig = self._make_fig()
        render_sweep_scatter(fig, results, "P_max")
        # 3x3 grid of axes + colorbar axes
        assert len(fig.axes) >= 9

    def test_title_contains_metric(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_scatter

        results = _make_sweep_results(n_params=2, n_runs=6)
        fig = self._make_fig()
        render_sweep_scatter(fig, results, "max_energy_error")
        # suptitle is on fig, not on individual axes
        assert "max_energy_error" in fig._suptitle.get_text()
