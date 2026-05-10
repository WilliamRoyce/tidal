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


@dataclass(frozen=True)
class RadialAngularPrior:
    """Joint prior over (r, theta_hat) for a coupling vector ``c = r * theta_hat``.

    The N coupling values form a vector in ``R^N``.  The prior factorises
    into a magnitude ``r`` (sampled ``log_uniform(r_lo, r_hi)``) and a
    direction ``theta_hat`` on the unit sphere ``S^(N-1)`` (sampled
    uniformly within one cubed-sphere sub-tile).  See
    :mod:`tidal.inference._sphere` for the angular chart.

    The ``transform`` method consumes ``N`` cube dimensions:

    - ``u[0]`` -> ``r`` via ``log_uniform(r_lo, r_hi)``
    - ``u[1:N]`` -> ``theta_hat`` via the cubed-sphere chart on the
      ``(face_idx, sub_tile, M, Q)`` cell

    and returns the physical coupling vector ``c = r * theta_hat`` of
    length N.

    Per-coupling signs are encoded in ``theta_hat`` (S^(N-1) is
    sign-symmetric); the prior covers the full real line for each
    coupling component.

    Designed for the nested-sampling code path; ``sample`` and
    ``log_prob`` are not implemented in this session (Monte Carlo
    support deferred — the v3 campaign use case is PolyChord).
    """

    names: tuple[str, ...]
    r_lo: float
    r_hi: float
    face_idx: int
    sub_tile: tuple[int, ...]
    M: int
    Q: NDArray[np.float64]

    def __post_init__(self) -> None:
        if len(self.names) < 2:
            msg = (
                f"RadialAngularPrior requires at least 2 couplings; "
                f"got {len(self.names)}"
            )
            raise ValueError(msg)
        if self.r_lo <= 0 or self.r_hi <= 0:
            msg = f"r bounds must be positive; got [{self.r_lo}, {self.r_hi}]"
            raise ValueError(msg)
        if self.r_lo >= self.r_hi:
            msg = f"r_lo must be < r_hi; got [{self.r_lo}, {self.r_hi}]"
            raise ValueError(msg)
        n = len(self.names)
        if len(self.sub_tile) != n - 1:
            msg = f"sub_tile must have length N - 1 = {n - 1}; got {len(self.sub_tile)}"
            raise ValueError(msg)
        Q = np.asarray(self.Q)  # noqa: N806
        if Q.shape != (n, n):
            msg = f"Q must have shape ({n}, {n}); got {Q.shape}"
            raise ValueError(msg)

    @property
    def n_dims(self) -> int:
        """Number of coupling dimensions covered by this joint prior (= N)."""
        return len(self.names)

    def transform(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        """Map ``u in [0, 1]^N`` to the physical coupling vector ``c``."""
        from tidal.inference._sphere import face_to_direction, tile_bounds

        u = np.asarray(u, dtype=np.float64)
        n = self.n_dims
        if u.shape != (n,):
            msg = f"u must have shape ({n},); got {u.shape}"
            raise ValueError(msg)

        log_lo, log_hi = math.log(self.r_lo), math.log(self.r_hi)
        r = math.exp(log_lo + float(u[0]) * (log_hi - log_lo))

        u_lo, u_hi = tile_bounds(self.sub_tile, self.M)
        u_face = u_lo + u[1:] * (u_hi - u_lo)

        theta_hat = face_to_direction(self.face_idx, u_face, self.Q)
        return r * theta_hat


def parse_joint_prior(spec: str) -> RadialAngularPrior:
    """Parse a CLI joint-prior specification string.

    Format: ``names=N1,N2,...;type=cubed_sphere;M=K;face=F;sub=S1_S2_...;``
    ``r_lo=L;r_hi=H[;Q=identity|random:SEED]``

    Top-level fields are separated by ``;``; each is ``key=value``.
    Optional ``Q=identity`` (default) or ``Q=random:SEED``.

    Examples
    --------
    >>> spec = (
    ...     "names=alpha1,alpha2,alpha3,delta1;"
    ...     "type=cubed_sphere;M=2;face=1;sub=1_1;"
    ...     "r_lo=1e-3;r_hi=1e3"
    ... )
    >>> p = parse_joint_prior(spec)
    >>> p.n_dims
    4

    Raises
    ------
    ValueError
        If the specification string is malformed.
    """
    from tidal.inference._sphere import random_rotation

    fields: dict[str, str] = {}
    for raw in spec.split(";"):
        part = raw.strip()
        if not part:
            continue
        if "=" not in part:
            msg = f"joint-prior field must be 'key=value', got '{part}'"
            raise ValueError(msg)
        k, v = part.split("=", 1)
        fields[k.strip()] = v.strip()

    required = ("names", "type", "M", "face", "sub", "r_lo", "r_hi")
    missing = [r for r in required if r not in fields]
    if missing:
        msg = (
            f"joint-prior missing required field(s): {', '.join(missing)}. "
            f"Format: 'names=N1,N2,...;type=cubed_sphere;M=K;face=F;"
            f"sub=S1_S2_...;r_lo=L;r_hi=H'"
        )
        raise ValueError(msg)

    if fields["type"] != "cubed_sphere":
        msg = f"only type=cubed_sphere is supported, got '{fields['type']}'"
        raise ValueError(msg)

    names = tuple(s.strip() for s in fields["names"].split(",") if s.strip())
    if len(names) < 2:
        msg = f"joint-prior 'names' must list >= 2 couplings, got {names}"
        raise ValueError(msg)

    try:
        m = int(fields["M"])
        face_idx = int(fields["face"])
        sub_tile = tuple(int(s) for s in fields["sub"].split("_") if s)
        r_lo = float(fields["r_lo"])
        r_hi = float(fields["r_hi"])
    except ValueError as e:
        msg = f"joint-prior field parse error: {e}"
        raise ValueError(msg) from None

    n = len(names)
    Q = np.eye(n)  # noqa: N806
    q_spec = fields.get("Q", "identity")
    if q_spec == "identity":
        pass
    elif q_spec.startswith("random:"):
        try:
            seed = int(q_spec.split(":", 1)[1])
        except ValueError:
            msg = f"joint-prior Q=random:SEED expects integer seed, got '{q_spec}'"
            raise ValueError(msg) from None
        Q = random_rotation(n, seed=seed)  # noqa: N806
    else:
        msg = f"joint-prior Q must be 'identity' or 'random:SEED', got '{q_spec}'"
        raise ValueError(msg)

    return RadialAngularPrior(
        names=names,
        r_lo=r_lo,
        r_hi=r_hi,
        face_idx=face_idx,
        sub_tile=sub_tile,
        M=m,
        Q=Q,
    )


def total_prior_ndim(priors: list[Prior | RadialAngularPrior]) -> int:
    """Total cube dimension covered by a (possibly mixed) prior list."""
    return sum(p.n_dims if isinstance(p, RadialAngularPrior) else 1 for p in priors)


def prior_param_names(priors: list[Prior | RadialAngularPrior]) -> list[str]:
    """Flat list of physical-parameter names from a (possibly mixed) prior list.

    Per-coupling Priors contribute their ``name``; RadialAngularPriors
    contribute each entry of ``names`` in order.  Result is the column
    order downstream code (chain matrix, plotting) should use.
    """
    out: list[str] = []
    for p in priors:
        if isinstance(p, RadialAngularPrior):
            out.extend(p.names)
        else:
            out.append(p.name)
    return out


def validate_prior_names(priors: list[Prior | RadialAngularPrior]) -> None:
    """Reject duplicate parameter names in a (possibly mixed) prior list."""
    names = prior_param_names(priors)
    duplicates = {n for n in names if names.count(n) > 1}
    if duplicates:
        msg = (
            f"duplicate parameter name(s) in prior list: "
            f"{sorted(duplicates)}; --prior and --joint-prior name spaces "
            f"must be disjoint"
        )
        raise ValueError(msg)


def build_prior_transform(
    priors: list[Prior] | list[Prior | RadialAngularPrior],
) -> Callable[[NDArray[np.float64]], NDArray[np.float64]]:
    """Build a joint prior_transform for nested sampling.

    Accepts a list mixing :class:`Prior` (per-parameter) and
    :class:`RadialAngularPrior` (joint over N couplings).  Each Prior
    consumes 1 cube dimension; each RadialAngularPrior consumes
    ``n_dims`` (= len(names)).

    Returns a callable that maps a unit-cube vector ``u`` of shape
    ``(total_ndim,)`` to physical parameters in the column order given
    by :func:`prior_param_names`.
    """

    def prior_transform(u: NDArray[np.float64]) -> NDArray[np.float64]:
        out: list[float] = []
        idx = 0
        for p in priors:
            if isinstance(p, RadialAngularPrior):
                n = p.n_dims
                segment = p.transform(u[idx : idx + n])
                out.extend(float(x) for x in segment)
                idx += n
            else:
                out.append(p.transform(float(u[idx])))
                idx += 1
        return np.asarray(out, dtype=np.float64)

    return prior_transform
