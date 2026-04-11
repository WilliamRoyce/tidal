"""Tests for sweep panel plotting functions and derived metrics.

Verifies that the new plot types (sweep-histogram, sweep-divergence,
campaign), derived metrics (A, C0, log10_A), and helper functions
(_greek_label, _apply_arctan_axis) work correctly.
"""

from __future__ import annotations

import math

import matplotlib as mpl

mpl.use("Agg")
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from tidal.measurement._sweep_results import SweepResults

if TYPE_CHECKING:
    from pathlib import Path


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_sweep_results(  # noqa: C901
    *,
    n_points: int = 10,
    n_params: int = 1,
    include_diverged: bool = False,
) -> SweepResults:
    """Create a mock SweepResults for testing."""
    [f"param{i}" for i in range(n_params)]
    swept_params: dict[str, list[float]] = {}
    rows: list[dict] = []

    if n_params == 1:
        vals = np.linspace(-1.0, 1.0, n_points).tolist()
        swept_params["param0"] = vals
        for i, v in enumerate(vals):
            status = "success"
            if include_diverged and i == n_points - 1:
                status = "diverged"
            rows.append(
                {
                    "param0": v,
                    "P_max": 0.001 * (1 + v**2) if status == "success" else None,
                    "run_status": status,
                    "B0": 0.001,
                    "kappa": 1.0,
                    "t_end": 50.0,
                }
            )
    elif n_params == 2:
        vals0 = np.linspace(-1.0, 1.0, n_points).tolist()
        vals1 = np.linspace(0.0, 2.0, n_points).tolist()
        swept_params["param0"] = vals0
        swept_params["param1"] = vals1
        for v0 in vals0:
            for v1 in vals1:
                status = "success"
                if include_diverged and abs(v0) > 0.8 and v1 > 1.5:
                    status = "diverged"
                rows.append(
                    {
                        "param0": v0,
                        "param1": v1,
                        "P_max": 0.001 * (1 + v0**2 + v1)
                        if status == "success"
                        else None,
                        "run_status": status,
                        "B0": 0.001,
                        "kappa": 1.0,
                        "t_end": 50.0,
                    }
                )
    else:
        # Multi-param LHS-like
        rng = np.random.default_rng(42)
        for i in range(n_params):
            swept_params[f"param{i}"] = np.linspace(-1.0, 1.0, n_points).tolist()
        for _ in range(n_points):
            row: dict = {}
            for i in range(n_params):
                row[f"param{i}"] = float(rng.uniform(-1, 1))
            row["P_max"] = float(rng.uniform(0.0001, 0.01))
            row["run_status"] = "success"
            row["B0"] = 0.001
            row["kappa"] = 1.0
            row["t_end"] = 50.0
            rows.append(row)

    return SweepResults(
        swept_params=swept_params,
        fixed_params={"B0": 0.001, "kappa": 1.0},
        sim_settings={"t_end": 50.0, "grid_shape": 128},
        rows=rows,
        run_dirs=[],
        spec_path="test.json",
        measurements=["conversion"],
        metadata={"sampling_strategy": "grid" if n_params <= 2 else "latin_hypercube"},
    )


# ---------------------------------------------------------------------------
# Greek labels
# ---------------------------------------------------------------------------


class TestGreekLabel:
    def test_known_param(self) -> None:
        from tidal.cli._sweep_panels import _greek_label

        assert _greek_label("delta1") == r"$\delta_1$"
        assert _greek_label("beta2") == r"$\beta_2$"
        assert _greek_label("xi") == r"$\xi$"
        assert _greek_label("B0") == r"$B_0$"

    def test_unknown_param_passthrough(self) -> None:
        from tidal.cli._sweep_panels import _greek_label

        assert _greek_label("my_custom_param") == "my_custom_param"


# ---------------------------------------------------------------------------
# Derived metrics
# ---------------------------------------------------------------------------


class TestDerivedColumns:
    def test_adds_expected_columns(self) -> None:
        results = _make_sweep_results(n_points=5)
        results.add_derived_columns()

        for row in results.rows:
            if row["run_status"] == "success":
                assert "A" in row
                assert "log10_A" in row
                assert "C0" in row
                assert "P_EM" in row
                assert "P_valid" in row

    def test_amplification_baseline(self) -> None:
        results = _make_sweep_results(n_points=5)
        results.add_derived_columns()

        # P_EM should be sin^2(kappa * B0 * t_end / 2)
        expected_p_em = math.sin(1.0 * 0.001 * 50.0 / 2) ** 2
        for row in results.rows:
            if row.get("P_EM") is not None:
                assert abs(row["P_EM"] - expected_p_em) < 1e-10

    def test_validity_check(self) -> None:
        results = _make_sweep_results(n_points=5)
        results.add_derived_columns(validity_threshold=0.1)

        for row in results.rows:
            if row["run_status"] == "success":
                p_max = row["P_max"]
                assert row["P_valid"] == (0 < p_max < 0.1)


# ---------------------------------------------------------------------------
# Render functions (smoke tests — verify they don't crash)
# ---------------------------------------------------------------------------


class TestRenderSweepHistogram:
    def test_renders_without_error(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_histogram

        results = _make_sweep_results(n_points=20)
        fig, ax = plt.subplots()
        render_sweep_histogram(ax, results, "P_max")
        plt.close(fig)

    def test_with_derived_metric(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_histogram

        results = _make_sweep_results(n_points=20)
        results.add_derived_columns()
        fig, ax = plt.subplots()
        render_sweep_histogram(ax, results, "log10_A")
        plt.close(fig)


class TestRenderSweepDivergence:
    def test_renders_without_error(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_divergence

        results = _make_sweep_results(n_points=10, include_diverged=True)
        fig, ax = plt.subplots()
        render_sweep_divergence(ax, results)
        plt.close(fig)


class TestRenderSweep1dUpgrades:
    def test_scatter_mode(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_1d

        results = _make_sweep_results(n_points=10)
        fig, ax = plt.subplots()
        render_sweep_1d(ax, results, "P_max", scatter=True)
        plt.close(fig)

    def test_annotate_extremes(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_1d

        results = _make_sweep_results(n_points=10)
        fig, ax = plt.subplots()
        render_sweep_1d(ax, results, "P_max", annotate_extremes=True)
        plt.close(fig)

    def test_arctan_axes(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_1d

        results = _make_sweep_results(n_points=10)
        fig, ax = plt.subplots()
        render_sweep_1d(ax, results, "P_max", arctan_axes=True)
        plt.close(fig)


class TestRenderSweepScatterUpgrades:
    def test_log_diagonal(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_scatter

        results = _make_sweep_results(n_points=5, n_params=3)
        fig = plt.figure()
        render_sweep_scatter(fig, results, "P_max", log_diagonal=True)
        plt.close(fig)

    def test_params_filter(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_scatter

        results = _make_sweep_results(n_points=5, n_params=3)
        fig = plt.figure()
        render_sweep_scatter(fig, results, "P_max", params=["param0", "param1"])
        plt.close(fig)


class TestRender2dGridUpgrades:
    def test_diverged_hatching(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_2d

        results = _make_sweep_results(n_points=5, n_params=2, include_diverged=True)
        fig, ax = plt.subplots()
        render_sweep_2d(ax, results, "P_max", divergent_center=0.001)
        plt.close(fig)

    def test_clamp_color(self) -> None:
        from tidal.cli._sweep_panels import render_sweep_2d

        results = _make_sweep_results(n_points=5, n_params=2)
        fig, ax = plt.subplots()
        render_sweep_2d(ax, results, "P_max", clamp_color=(0.0, 0.01))
        plt.close(fig)


class TestRenderSweepTornado:
    def test_pvalue_significance(self) -> None:
        """Verify tornado renders with p-value coloring (needs scattered data)."""
        from tidal.cli._sweep_panels import render_sweep_tornado

        results = _make_sweep_results(n_points=30, n_params=3)
        fig, ax = plt.subplots()
        render_sweep_tornado(ax, results, "P_max")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Campaign
# ---------------------------------------------------------------------------


class TestCampaign:
    def test_renders_1d_subdirs(self, tmp_path: Path) -> None:
        """Create mock sweep subdirs and verify campaign renders."""
        from tidal.cli._campaign import render_campaign

        # Create 3 mock sweep subdirectories
        for name in ["param_a", "param_b", "param_c"]:
            subdir = tmp_path / name
            subdir.mkdir()
            results = _make_sweep_results(n_points=5)
            results.to_csv(subdir / "results.csv")  # needed for from_directory
            # Write results.json
            results.to_json(subdir / "results.json")
            results.save_sweep_json(subdir / "sweep.json")

        # Create mock args namespace
        from argparse import Namespace

        args = Namespace(
            data_path=str(tmp_path),
            metric="P_max",
            layout="1x3",
            figsize=None,
            dpi=72,
            scatter=False,
            annotate_extremes=False,
            arctan_axes=False,
            log_y=False,
            hline=None,
            classify=None,
            title=None,
            baseline_formula=None,
            output=str(tmp_path / "campaign.png"),
            divergent_center=None,
            cmap=None,
            clamp_color=None,
        )

        rc = render_campaign(args)
        assert rc == 0
        assert (tmp_path / "campaign.png").exists()

    def test_renders_2d_subdirs(self, tmp_path: Path) -> None:
        """Verify campaign handles 2D sweep subdirectories."""
        from tidal.cli._campaign import render_campaign

        subdir = tmp_path / "pair_ab"
        subdir.mkdir()
        results = _make_sweep_results(n_points=5, n_params=2)
        results.to_json(subdir / "results.json")
        results.save_sweep_json(subdir / "sweep.json")

        from argparse import Namespace

        args = Namespace(
            data_path=str(tmp_path),
            metric="P_max",
            layout=None,
            figsize=None,
            dpi=72,
            scatter=False,
            annotate_extremes=False,
            arctan_axes=False,
            log_y=False,
            hline=None,
            classify=None,
            title=None,
            baseline_formula=None,
            output=str(tmp_path / "campaign_2d.png"),
            divergent_center=None,
            cmap=None,
            clamp_color=None,
        )

        rc = render_campaign(args)
        assert rc == 0
        assert (tmp_path / "campaign_2d.png").exists()
