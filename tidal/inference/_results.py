"""Inference results storage and analysis.

:class:`InferenceResult` stores samples from Monte Carlo or nested
sampling, with conversion to :class:`~tidal.measurement._sweep_results.SweepResults`
for CSV/JSON export.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from tidal.measurement._sweep_results import SweepResults


@dataclass
class InferenceResult:
    """Container for inference samples and diagnostics.

    Parameters
    ----------
    samples : NDArray
        Parameter samples, shape ``(n_samples, n_params)``.
    log_likelihood : NDArray
        Log-likelihood for each sample, shape ``(n_samples,)``.
    log_prior : NDArray
        Log-prior for each sample, shape ``(n_samples,)``.
    param_names : list[str]
        Parameter names, length ``n_params``.
    method : str
        Sampling method: ``"mc"`` or ``"nested"``.
    metrics : dict[str, NDArray] | None
        Additional simulation metrics per sample (P_max, etc.).
    log_evidence : float | None
        Log-evidence (nested sampling only).
    log_evidence_err : float | None
        Uncertainty on log-evidence (nested sampling only).
    weights : NDArray | None
        Importance weights for nested sampling, shape ``(n_samples,)``.
    metadata : dict
        Additional metadata (sampler settings, wall time, etc.).
    """

    samples: NDArray[np.float64]
    log_likelihood: NDArray[np.float64]
    log_prior: NDArray[np.float64]
    param_names: list[str]
    method: str
    metrics: dict[str, NDArray[np.float64]] | None = None
    log_evidence: float | None = None
    log_evidence_err: float | None = None
    weights: NDArray[np.float64] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return self.samples.shape[0]

    @property
    def n_params(self) -> int:
        return self.samples.shape[1]

    @property
    def log_posterior(self) -> NDArray[np.float64]:
        """Log-posterior = log-prior + log-likelihood."""
        return self.log_prior + self.log_likelihood

    def effective_sample_size(self) -> float:
        """Kish effective sample size from weights.

        For unweighted MC samples, returns ``n_samples``.
        For nested sampling, uses importance weights.

        Reference: Kish, L. (1965) *Survey Sampling*, Wiley.
        """
        if self.weights is None:
            return float(self.n_samples)
        w = self.weights
        w_sum = np.sum(w)
        if w_sum == 0:
            return 0.0
        return float(w_sum**2 / np.sum(w**2))

    def best(self) -> dict[str, float]:
        """Return the MAP (maximum a posteriori) parameter values."""
        lp = self.log_posterior
        idx = int(np.argmax(lp))
        return dict(zip(self.param_names, self.samples[idx], strict=False))

    def posterior_mean(self) -> dict[str, float]:
        """Return the posterior mean parameter values.

        For nested sampling, uses importance-weighted mean.
        """
        if self.weights is not None:
            w = self.weights / np.sum(self.weights)
            means = np.average(self.samples, weights=w, axis=0)
        else:
            means = np.mean(self.samples, axis=0)
        return dict(zip(self.param_names, means, strict=False))

    def credible_interval(self, level: float = 0.95) -> dict[str, tuple[float, float]]:
        """Return equal-tailed credible intervals for each parameter.

        Parameters
        ----------
        level : float
            Credible level (default: 0.95 for 95% CI).
        """
        alpha = (1 - level) / 2
        result = {}
        for i, name in enumerate(self.param_names):
            values = self.samples[:, i]
            if self.weights is not None:
                # Weighted quantiles
                w = self.weights / np.sum(self.weights)
                sorted_idx = np.argsort(values)
                sorted_vals = values[sorted_idx]
                cumw = np.cumsum(w[sorted_idx])
                lo = float(sorted_vals[np.searchsorted(cumw, alpha)])
                hi = float(sorted_vals[np.searchsorted(cumw, 1 - alpha)])
            else:
                lo = float(np.quantile(values, alpha))
                hi = float(np.quantile(values, 1 - alpha))
            result[name] = (lo, hi)
        return result

    def to_anesthetic(self) -> Any:
        """Convert to anesthetic ``NestedSamples`` for custom analysis.

        Returns the native anesthetic object with ``.D_KL()``,
        ``.d_G()``, ``.logZ()``, ``.plot_2d()`` etc.
        """
        from tidal.inference._importance import to_anesthetic_samples

        return to_anesthetic_samples(self)

    def parameter_importance(self, n_bootstrap: int = 100) -> Any:
        """Compute parameter importance via KL divergence.

        Uses anesthetic to compute total and per-parameter information
        gain from prior to posterior, plus Bayesian model dimensionality.

        Parameters
        ----------
        n_bootstrap : int
            Number of bootstrap samples for uncertainty estimation.

        Returns
        -------
        ParameterImportanceResult
        """
        from tidal.inference._importance import compute_parameter_importance

        return compute_parameter_importance(self, n_bootstrap)

    @classmethod
    def from_directory(cls, path: Path) -> InferenceResult:
        """Load inference results from a saved directory.

        Expects ``inference.json`` and ``results.csv`` as written by
        :meth:`save`.

        Parameters
        ----------
        path : Path
            Directory containing saved inference results.
        """
        import csv
        import json
        from pathlib import Path as _Path

        path = _Path(path)

        # Load inference metadata
        with (path / "inference.json").open() as f:
            meta = json.load(f)

        # Load samples from CSV
        rows: list[dict[str, str]] = []
        with (path / "results.csv").open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        if not rows:
            msg = f"No rows in {path / 'results.csv'}"
            raise ValueError(msg)

        param_names: list[str] = meta["param_names"]
        n_samples = len(rows)
        n_params = len(param_names)

        samples = np.zeros((n_samples, n_params))
        log_likelihood = np.zeros(n_samples)
        log_prior = np.zeros(n_samples)
        weights_list: list[float] = []

        for i, row in enumerate(rows):
            for j, name in enumerate(param_names):
                samples[i, j] = float(row[name])
            log_likelihood[i] = float(row.get("log_likelihood", 0))
            log_prior[i] = float(row.get("log_prior", 0))
            w = row.get("weight")
            if w is not None:
                weights_list.append(float(w))

        weights = np.array(weights_list) if weights_list else None

        # Reconstruct per-sample metrics from CSV columns that aren't part
        # of the standard schema (params, log_likelihood, log_prior,
        # log_posterior, weight, error_message).  This is what carries
        # run_status, tachyonic_excess, and any per-sample physics metrics
        # back into the InferenceResult so the corner-plot tool sees them.
        standard_cols = {
            *param_names,
            "log_likelihood",
            "log_prior",
            "log_posterior",
            "weight",
            "error_message",
        }
        loaded_metrics: dict[str, np.ndarray] = {}
        if rows:
            for col in rows[0]:
                if col in standard_cols:
                    continue
                # Try to parse as float column; fall back to object/string
                # for run_status and similar string-typed columns.
                try:
                    arr_f = np.asarray(
                        [
                            float(r[col])
                            if r.get(col) not in {None, "", "None"}
                            else np.nan
                            for r in rows
                        ],
                        dtype=float,
                    )
                    loaded_metrics[col] = arr_f
                except (TypeError, ValueError):
                    loaded_metrics[col] = np.asarray(
                        [r.get(col) for r in rows], dtype=object
                    )

        # Reconstruct metadata from saved inference.json fields
        loaded_metadata: dict[str, Any] = {
            "loaded_from": str(path),
        }
        for key in (
            "nlive",
            "sampler",
            "n_iterations",
            "n_calls",
            "priors",
            "rejected_prior_path",
            "parameter_importance",
        ):
            val = meta.get(key)
            if val is not None:
                loaded_metadata[key] = val

        return cls(
            samples=samples,
            log_likelihood=log_likelihood,
            log_prior=log_prior,
            param_names=param_names,
            method=meta.get("method", "nested"),
            log_evidence=meta.get("log_evidence"),
            log_evidence_err=meta.get("log_evidence_err"),
            weights=weights,
            metrics=loaded_metrics or None,
            metadata=loaded_metadata,
        )

    def to_sweep_results(self) -> SweepResults:
        """Convert to a :class:`~tidal.measurement._sweep_results.SweepResults`.

        This allows reusing the existing CSV/JSON serialization infrastructure.
        """
        from tidal.measurement._sweep_results import SweepResults

        rows: list[dict[str, Any]] = []
        for i in range(self.n_samples):
            row: dict[str, Any] = {}
            for j, name in enumerate(self.param_names):
                row[name] = float(self.samples[i, j])
            row["log_likelihood"] = float(self.log_likelihood[i])
            row["log_prior"] = float(self.log_prior[i])
            row["log_posterior"] = float(self.log_prior[i] + self.log_likelihood[i])
            if self.weights is not None:
                row["weight"] = float(self.weights[i])
            row_status = "success"
            row_error: str | None = None
            if self.metrics is not None:
                for key, arr in self.metrics.items():
                    val = arr[i]
                    if key == "run_status":
                        # Stored as object dtype (string); not numeric.
                        row_status = str(val) if val is not None else "success"
                        row[key] = row_status
                        continue
                    # Numeric metrics: NaN -> None for JSON cleanliness.
                    try:
                        row[key] = float(val) if not np.isnan(float(val)) else None
                    except (TypeError, ValueError):
                        row[key] = val if val is not None else None
            row["run_status"] = row_status
            row["error_message"] = row_error
            rows.append(row)

        swept_params = {
            name: list(np.unique(self.samples[:, i]))
            for i, name in enumerate(self.param_names)
        }

        # Filter out non-serializable metadata (e.g. cached anesthetic objects)
        serializable_meta = {
            k: v
            for k, v in self.metadata.items()
            if not k.startswith("_")
            and isinstance(v, (str, int, float, bool, list, dict, type(None)))
        }
        meta: dict[str, Any] = {
            "inference_method": self.method,
            **serializable_meta,
        }
        if self.log_evidence is not None:
            meta["log_evidence"] = self.log_evidence
            meta["log_evidence_err"] = self.log_evidence_err

        return SweepResults(
            swept_params=swept_params,
            fixed_params={},
            sim_settings={},
            rows=rows,
            run_dirs=[],
            spec_path="",
            measurements=["inference"],
            source_fields=None,
            target_fields=None,
            metadata=meta,
            converge_sizes=None,
        )

    def save(self, output_dir: Path) -> None:
        """Save inference results to disk.

        Writes:
        - ``results.csv`` and ``results.json`` via SweepResults
        - ``inference.json`` with evidence, ESS, and summary statistics
        """
        import json
        from datetime import UTC, datetime

        output_dir.mkdir(parents=True, exist_ok=True)

        # Save via SweepResults for CSV/JSON
        sr = self.to_sweep_results()
        sr.to_csv(output_dir / "results.csv")
        sr.to_json(output_dir / "results.json")

        # Save inference-specific metadata
        summary: dict[str, Any] = {
            "method": self.method,
            "n_samples": self.n_samples,
            "n_params": self.n_params,
            "param_names": self.param_names,
            "effective_sample_size": round(self.effective_sample_size(), 1),
            "map_estimate": self.best(),
            "posterior_mean": self.posterior_mean(),
            "credible_interval_95": {
                k: list(v) for k, v in self.credible_interval(0.95).items()
            },
            "timestamp": datetime.now(tz=UTC).isoformat(),
        }
        if self.log_evidence is not None:
            summary["log_evidence"] = self.log_evidence
            summary["log_evidence_err"] = self.log_evidence_err

        # Preserve sampler metadata for round-trip loading
        nlive = self.metadata.get("nlive")
        if nlive is not None:
            summary["nlive"] = nlive
        sampler = self.metadata.get("sampler")
        if sampler is not None:
            summary["sampler"] = sampler
        priors = self.metadata.get("priors")
        if priors is not None:
            summary["priors"] = priors
        rejected_prior_path = self.metadata.get("rejected_prior_path")
        if rejected_prior_path is not None:
            summary["rejected_prior_path"] = str(rejected_prior_path)

        # Compute parameter importance for nested sampling results
        if self.method == "nested":
            try:
                importance = self.parameter_importance()
                summary["parameter_importance"] = {
                    "d_kl": importance.d_kl,
                    "d_kl_err": importance.d_kl_err,
                    "d_g": importance.d_g,
                    "d_g_err": importance.d_g_err,
                    "marginal_d_kl": importance.marginal_d_kl,
                }
                # Also save as standalone file
                with (output_dir / "importance.json").open("w") as f:
                    json.dump(
                        {
                            "param_names": importance.param_names,
                            "d_kl": importance.d_kl,
                            "d_kl_err": importance.d_kl_err,
                            "d_g": importance.d_g,
                            "d_g_err": importance.d_g_err,
                            "marginal_d_kl": importance.marginal_d_kl,
                            "log_evidence": importance.log_evidence,
                            "log_evidence_err": importance.log_evidence_err,
                        },
                        f,
                        indent=2,
                    )
            except Exception:  # noqa: BLE001
                import logging

                logging.getLogger("tidal.inference").debug(
                    "Parameter importance computation skipped during save",
                    exc_info=True,
                )

        with (output_dir / "inference.json").open("w") as f:
            json.dump(summary, f, indent=2)
