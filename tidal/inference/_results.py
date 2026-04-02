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
            if self.metrics is not None:
                for key, arr in self.metrics.items():
                    row[key] = float(arr[i]) if not np.isnan(arr[i]) else None
            row["run_status"] = "success"
            row["error_message"] = None
            rows.append(row)

        swept_params = {
            name: list(np.unique(self.samples[:, i]))
            for i, name in enumerate(self.param_names)
        }

        meta: dict[str, Any] = {
            "inference_method": self.method,
            **self.metadata,
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

        with (output_dir / "inference.json").open("w") as f:
            json.dump(summary, f, indent=2)
