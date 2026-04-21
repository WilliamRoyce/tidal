"""Prior distributions for Bayesian inference.

Each :class:`Prior` defines a 1-D marginal distribution with three
operations required by the sampling pipeline:

- ``sample(rng, n)`` — draw *n* random variates (for Monte Carlo)
- ``log_prob(x)`` — evaluate log-probability density (for posterior weighting)
- ``transform(u)`` — map *u* in [0, 1] to the physical parameter (for
  nested sampling's ``prior_transform`` protocol)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray


# Arctan mapping boundary guard: avoid tan(±pi/2) = ±inf
_ARCTAN_EPS = 0.05


@dataclass(frozen=True)
class Prior:
    """A 1-D marginal prior distribution.

    Parameters
    ----------
    name : str
        Parameter name (must match the JSON spec's parameter key).
    distribution : str
        One of ``"uniform"``, ``"log_uniform"``, ``"normal"``,
        ``"arctan_uniform"``.
    low : float
        Lower bound (for uniform/log_uniform) or mean (for normal).
    high : float
        Upper bound (for uniform/log_uniform) or std (for normal).
    """

    name: str
    distribution: str
    low: float
    high: float

    def __post_init__(self) -> None:
        valid = {"uniform", "log_uniform", "normal", "arctan_uniform"}
        if self.distribution not in valid:
            msg = f"Unknown distribution '{self.distribution}'. Must be one of {sorted(valid)}."
            raise ValueError(msg)
        if self.distribution == "log_uniform" and (self.low <= 0 or self.high <= 0):
            msg = f"log_uniform requires positive bounds, got [{self.low}, {self.high}]"
            raise ValueError(msg)
        if self.distribution in {"uniform", "log_uniform"} and self.low >= self.high:
            msg = f"Prior bounds must satisfy low < high, got [{self.low}, {self.high}]"
            raise ValueError(msg)
        if self.distribution == "normal" and self.high <= 0:
            msg = f"Normal prior std must be positive, got {self.high}"
            raise ValueError(msg)

    # ------------------------------------------------------------------
    # Sampling (Monte Carlo)
    # ------------------------------------------------------------------

    def sample(self, rng: np.random.Generator, n: int) -> NDArray[np.float64]:
        """Draw *n* samples from this prior."""
        if self.distribution == "uniform":
            return rng.uniform(self.low, self.high, size=n)
        if self.distribution == "log_uniform":
            log_lo, log_hi = math.log(self.low), math.log(self.high)
            return np.exp(rng.uniform(log_lo, log_hi, size=n))
        if self.distribution == "normal":
            return rng.normal(self.low, self.high, size=n)
        if self.distribution == "arctan_uniform":
            theta_lo = -math.pi / 2 + _ARCTAN_EPS
            theta_hi = math.pi / 2 - _ARCTAN_EPS
            theta = rng.uniform(theta_lo, theta_hi, size=n)
            return np.tan(theta)
        msg = f"Unsupported distribution: {self.distribution}"  # pragma: no cover
        raise ValueError(msg)

    # ------------------------------------------------------------------
    # Log-probability density (posterior weighting)
    # ------------------------------------------------------------------

    def log_prob(self, x: float) -> float:
        """Evaluate log p(x) under this prior."""
        if self.distribution == "uniform":
            if self.low <= x <= self.high:
                return -math.log(self.high - self.low)
            return -math.inf
        if self.distribution == "log_uniform":
            if self.low <= x <= self.high:
                return -math.log(x) - math.log(math.log(self.high / self.low))
            return -math.inf
        if self.distribution == "normal":
            mean, std = self.low, self.high
            return -0.5 * ((x - mean) / std) ** 2 - math.log(
                std * math.sqrt(2 * math.pi),
            )
        if self.distribution == "arctan_uniform":
            # p(x) = (1/pi) * 1/(1+x^2)  (Cauchy/arctan distribution)
            # Normalized over the truncated range
            theta_lo = -math.pi / 2 + _ARCTAN_EPS
            theta_hi = math.pi / 2 - _ARCTAN_EPS
            if x < math.tan(theta_lo) or x > math.tan(theta_hi):
                return -math.inf
            norm = theta_hi - theta_lo
            return -math.log(1 + x * x) - math.log(norm)
        msg = f"Unsupported distribution: {self.distribution}"  # pragma: no cover
        raise ValueError(msg)

    # ------------------------------------------------------------------
    # Unit-cube transform (nested sampling)
    # ------------------------------------------------------------------

    def transform(self, u: float) -> float:
        """Map *u* in [0, 1] to the physical parameter space.

        This implements the ``prior_transform`` protocol used by PolyChord.
        """
        if self.distribution == "uniform":
            return self.low + u * (self.high - self.low)
        if self.distribution == "log_uniform":
            log_lo, log_hi = math.log(self.low), math.log(self.high)
            return math.exp(log_lo + u * (log_hi - log_lo))
        if self.distribution == "normal":
            from scipy.special import ndtri

            mean, std = self.low, self.high
            return mean + std * float(ndtri(u))
        if self.distribution == "arctan_uniform":
            theta_lo = -math.pi / 2 + _ARCTAN_EPS
            theta_hi = math.pi / 2 - _ARCTAN_EPS
            theta = theta_lo + u * (theta_hi - theta_lo)
            return math.tan(theta)
        msg = f"Unsupported distribution: {self.distribution}"  # pragma: no cover
        raise ValueError(msg)


def parse_prior(spec: str) -> Prior:
    """Parse a CLI prior specification string.

    Format: ``NAME=DISTRIBUTION:ARG1:ARG2``

    Examples::

        "alpha=uniform:0.01:10"
        "xi=log_uniform:0.01:10"
        "delta=normal:0:1"
        "alpha=arctan_uniform:-30:30"

    Raises
    ------
    ValueError
        If the specification string is malformed.
    """
    if "=" not in spec:
        msg = f"Prior spec must contain '=': got '{spec}'"
        raise ValueError(msg)
    name, rest = spec.split("=", 1)
    parts = rest.split(":")
    if len(parts) < 3:
        msg = (
            f"Prior spec needs DIST:ARG1:ARG2, got '{rest}'. "
            f"Example: '{name}=uniform:0.01:10'"
        )
        raise ValueError(msg)
    distribution = parts[0]
    try:
        low = float(parts[1])
        high = float(parts[2])
    except ValueError:
        msg = f"Prior bounds must be numeric, got '{parts[1]}' and '{parts[2]}'"
        raise ValueError(msg) from None
    return Prior(name=name, distribution=distribution, low=low, high=high)


def build_prior_transform(
    priors: list[Prior],
) -> Callable[[NDArray[np.float64]], NDArray[np.float64]]:
    """Build a joint prior_transform for nested sampling.

    Returns a callable that maps a unit-cube vector ``u`` of shape
    ``(ndim,)`` to physical parameters.
    """

    def prior_transform(u: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.array([p.transform(float(u[i])) for i, p in enumerate(priors)])

    return prior_transform
