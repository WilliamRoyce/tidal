"""Posterior visualization for inference results.

Uses **anesthetic** (Handley 2019) when available for publication-quality
corner plots that integrate natively with nested sampling output.
Falls back to matplotlib for basic visualization.

References
----------
Handley, W. (2019) "anesthetic: nested sampling visualization",
JOSS 4(37), 1414. https://doi.org/10.21105/joss.01414
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from tidal.inference._results import InferenceResult


def plot_corner(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> None:
    """Generate a corner plot (2D marginals + 1D histograms).

    Uses anesthetic if available, otherwise falls back to matplotlib.

    Parameters
    ----------
    result : InferenceResult
        Inference results with samples and weights.
    output_path : Path | None
        If provided, save figure to this path.
    show : bool
        If True, display the figure interactively.
    """
    try:
        _plot_corner_anesthetic(result, output_path, show=show)
    except ImportError:
        _plot_corner_matplotlib(result, output_path, show=show)


def plot_trace(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> None:
    """Generate trace plots (parameter values vs sample index).

    Parameters
    ----------
    result : InferenceResult
        Inference results.
    output_path : Path | None
        If provided, save figure to this path.
    show : bool
        If True, display the figure interactively.
    """
    import matplotlib.pyplot as plt

    n_params = result.n_params
    fig, axes = plt.subplots(n_params + 1, 1, figsize=(10, 2.5 * (n_params + 1)))
    if n_params == 0:
        return

    axes_list = [axes] if n_params == 0 else list(axes)

    for i, name in enumerate(result.param_names):
        ax = axes_list[i]
        ax.plot(result.samples[:, i], ".", markersize=1, alpha=0.5)
        ax.set_ylabel(name)
        ax.set_xlabel("")

    # Log-likelihood trace
    ax = axes_list[-1]
    ax.plot(result.log_likelihood, ".", markersize=1, alpha=0.5, color="red")
    ax.set_ylabel("log L")
    ax.set_xlabel("Sample index")

    fig.suptitle(f"Trace plot ({result.method}, n={result.n_samples})")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def _plot_corner_anesthetic(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> None:
    """Corner plot using anesthetic (Handley 2019)."""
    import matplotlib.pyplot as plt
    from anesthetic import MCMCSamples, NestedSamples

    if result.method == "nested" and result.weights is not None:
        samples = NestedSamples(
            data=result.samples,
            logL=result.log_likelihood,
            columns=result.param_names,
        )
    else:
        samples = MCMCSamples(
            data=result.samples,
            columns=result.param_names,
        )

    fig, _axes = samples.plot_2d(result.param_names)
    fig.suptitle(
        f"Posterior ({result.method}, n={result.n_samples}, "
        f"ESS={result.effective_sample_size():.0f})",
        y=1.02,
    )

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def _plot_corner_matplotlib(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> None:
    """Basic corner plot using matplotlib (fallback)."""
    import matplotlib.pyplot as plt
    import numpy as np

    n = result.n_params
    if n == 0:
        return

    fig, axes = plt.subplots(n, n, figsize=(3 * n, 3 * n))
    if n == 1:
        axes = np.array([[axes]])

    # Use finite samples only
    mask = np.isfinite(result.log_posterior)
    samples = result.samples[mask]
    weights = result.weights[mask] if result.weights is not None else None

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False)
                continue
            if i == j:
                ax.hist(
                    samples[:, i], bins=30, weights=weights, density=True, alpha=0.7
                )
                ax.set_xlabel(result.param_names[i])
            else:
                ax.scatter(
                    samples[:, j],
                    samples[:, i],
                    c=result.log_likelihood[mask],
                    s=2,
                    alpha=0.5,
                    cmap="viridis",
                )
                if j == 0:
                    ax.set_ylabel(result.param_names[i])
                if i == n - 1:
                    ax.set_xlabel(result.param_names[j])

    fig.suptitle(
        f"Posterior ({result.method}, n={result.n_samples}, "
        f"ESS={result.effective_sample_size():.0f})"
    )
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
