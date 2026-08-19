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

#: Every distribution kind ``Prior`` accepts.  Single source of truth:
#: ``Prior.__post_init__`` validates against this set, and the marginal
#: D_KL estimator's kind-coverage test
#: (tests/test_inference.py::TestMarginalDKLPriorTransforms) parametrizes
#: over it — so adding a kind here without teaching
#: ``tidal.inference._importance._uniformizing_transform`` about it FAILS
#: the suite instead of silently falling back to the pre-#420
#: sample-range reference.
VALID_DISTRIBUTIONS = frozenset({"uniform", "log_uniform", "normal", "arctan_uniform"})


def effective_support(
    distribution: str, low: float, high: float
) -> tuple[float, float]:
    """The range a prior actually samples from — not the range that was typed.

    The single source of truth for "what range was sampled", so no
    consumer has to re-derive it from the recorded ``low``/``high``.
    Every place that did so independently got it wrong at least once:
    the marginal D_KL estimator histogrammed a +-20 posterior over
    ``(-89, 89)`` (GH #420) and the ``--full-prior-bounds`` corner plot
    drew ``arctan_uniform`` panels over ``tan(radians(+-89))`` = +-57.3
    (GH #451), because both read bounds that ``arctan_uniform`` ignores
    (GH #425).

    Exposed as a free function as well as :attr:`Prior.effective_support`
    so consumers reading archived prior *records* can ask without
    constructing a :class:`Prior` (which would re-emit the GH #425
    bounds-ignored warning on every replot of an archived chain).

    Returns
    -------
    tuple[float, float]
        ``(low, high)`` for ``uniform``/``log_uniform`` (the recorded
        bounds ARE the support); ``(-S, +S)`` with
        ``S = tan(pi/2 - _ARCTAN_EPS) ~= 19.98`` for ``arctan_uniform``,
        independent of the recorded bounds; ``(-inf, +inf)`` for
        ``normal``, whose support is genuinely unbounded — callers
        needing a finite plotting range must derive one from the mean
        and std themselves rather than be handed a silently truncated
        interval.

    Raises
    ------
    ValueError
        If ``distribution`` is not in :data:`VALID_DISTRIBUTIONS`.
    """
    if distribution not in VALID_DISTRIBUTIONS:
        msg = (
            f"Unknown distribution '{distribution}'. Must be one of "
            f"{sorted(VALID_DISTRIBUTIONS)}."
        )
        raise ValueError(msg)
    if distribution == "arctan_uniform":
        support = math.tan(math.pi / 2 - _ARCTAN_EPS)
        return (-support, support)
    if distribution == "normal":
        return (-math.inf, math.inf)
    return (float(low), float(high))


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
        **Ignored for arctan_uniform** — see below.
    high : float
        Upper bound (for uniform/log_uniform) or std (for normal).
        **Ignored for arctan_uniform** — see below.

    Notes
    -----
    ``arctan_uniform`` does NOT use ``low``/``high`` (GH #425): the angle
    is uniform on the fixed eps-truncated range
    ``(-pi/2 + _ARCTAN_EPS, +pi/2 - _ARCTAN_EPS)``, so the support is
    always ``|x| <= tan(pi/2 - _ARCTAN_EPS) ~= 19.98`` regardless of the
    recorded bounds.  A ``UserWarning`` fires at construction when the
    given bounds differ from that implied support; pass ``0:0`` (the
    sanctioned sentinel used in the docs) to declare the bounds
    deliberately unused without warning.  Honoring the bounds
    would silently redefine the prior for every archived chain whose
    metadata records these unused numbers, so any change must be
    versioned — deliberately deferred, see the options in GH #425.
    """

    name: str
    distribution: str
    low: float
    high: float

    def __post_init__(self) -> None:
        if self.distribution not in VALID_DISTRIBUTIONS:
            msg = (
                f"Unknown distribution '{self.distribution}'. Must be one "
                f"of {sorted(VALID_DISTRIBUTIONS)}."
            )
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
        if self.distribution == "arctan_uniform":
            # GH #425: sample/transform/log_prob use the fixed eps-truncated
            # theta range and ignore low/high entirely.  Tell the user at
            # launch rather than let the recorded bounds masquerade as the
            # support (which is how the #420 marginal D_KL inflation got a
            # (-89, 89) histogram range for a +-20 distribution).
            # ``0:0`` is the sanctioned "bounds unused" sentinel (used by
            # docs/tex/inference.tex examples) and does not warn; bounds
            # matching the implied support do not warn either.
            # stacklevel=3 skips the generated dataclass __init__ so the
            # warning points at the parse_prior / user call site.
            import warnings

            support = math.tan(math.pi / 2 - _ARCTAN_EPS)
            is_sentinel = self.low == 0.0 and self.high == 0.0
            matches_support = math.isclose(
                self.low, -support, rel_tol=1e-9
            ) and math.isclose(self.high, support, rel_tol=1e-9)
            if not (is_sentinel or matches_support):
                warnings.warn(
                    f"arctan_uniform ignores its bounds: requested "
                    f"[{self.low}, {self.high}] but the sampled support is "
                    f"fixed at [-{support:.2f}, {support:.2f}] "
                    f"(theta uniform on +-(pi/2 - {_ARCTAN_EPS})). "
                    f"Use 0:0 to declare the bounds deliberately unused. "
                    f"See GH #425.",
                    UserWarning,
                    stacklevel=3,
                )

    # ------------------------------------------------------------------
    # Support
    # ------------------------------------------------------------------

    @property
    def effective_support(self) -> tuple[float, float]:
        """The range :meth:`sample` actually draws from — not what was typed.

        Thin accessor for :func:`effective_support`; see it for the
        per-distribution contract and the defects that motivated it.
        """
        return effective_support(self.distribution, self.low, self.high)

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
        "alpha=arctan_uniform:-30:30"   # bounds recorded but IGNORED:
                                        # support is fixed at ~+-20 (GH #425)

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

    The N coupling values form a vector in ``R^N``.  The prior factorizes
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
        if self.r_lo > self.r_hi:
            msg = (
                f"r_lo must be <= r_hi; got [{self.r_lo}, {self.r_hi}] "
                f"(equal bounds are allowed as fixed-radius mode)"
            )
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

    @property
    def is_fixed_radius(self) -> bool:
        """True when ``r_lo == r_hi`` — only the angular direction is sampled.

        In fixed-radius mode ``u[0]`` is still consumed (cube dimensionality
        stays at ``N`` for sampler-API stability) but its value is discarded;
        the magnitude is always ``r_lo``.
        """
        return self.r_lo == self.r_hi

    def transform(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        """Map ``u in [0, 1]^N`` to the physical coupling vector ``c``."""
        from tidal.inference._sphere import face_to_direction, tile_bounds

        u = np.asarray(u, dtype=np.float64)
        n = self.n_dims
        if u.shape != (n,):
            msg = f"u must have shape ({n},); got {u.shape}"
            raise ValueError(msg)

        if self.is_fixed_radius:
            r = self.r_lo
        else:
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
    q_arr = np.eye(n)
    q_spec = fields.get("Q", "identity")
    if q_spec == "identity":
        pass
    elif q_spec.startswith("random:"):
        try:
            seed = int(q_spec.split(":", 1)[1])
        except ValueError:
            msg = f"joint-prior Q=random:SEED expects integer seed, got '{q_spec}'"
            raise ValueError(msg) from None
        q_arr = random_rotation(n, seed=seed)
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
        Q=q_arr,
    )


def total_prior_ndim(priors: list[Prior | RadialAngularPrior]) -> int:
    """Total cube dimension covered by a (possibly mixed) prior list."""
    return sum(p.n_dims if isinstance(p, RadialAngularPrior) else 1 for p in priors)


def to_record(prior: Prior | RadialAngularPrior) -> dict[str, object]:
    """Serialize a prior to its on-disk record (``inference.json`` ``priors``).

    One definition of the record shape, shared by the writer in
    :mod:`tidal.cli._sample` and by :meth:`InferenceResult.save`, so a
    result whose ``metadata["priors"]`` holds live prior objects — which
    the importance estimator accepts — reaches disk in the same form as
    one built by the CLI, instead of raising "not JSON serializable"
    there.  Every recorded chain having a faithful priors block is what
    GH #434 was about.

    Scalar records carry ``effective_low``/``effective_high`` alongside
    the requested ``low``/``high`` so post-hoc consumers read what was
    sampled rather than re-deriving it (GH #425/#451).  The keys are
    OMITTED for an unbounded support, since ``Infinity`` is not valid
    strict JSON.
    """
    if isinstance(prior, RadialAngularPrior):
        return {
            "kind": "radial_angular",
            "names": list(prior.names),
            "r_lo": prior.r_lo,
            "r_hi": prior.r_hi,
            "face_idx": prior.face_idx,
            "sub_tile": list(prior.sub_tile),
            "M": prior.M,
            "Q": np.asarray(prior.Q).tolist(),
        }
    record: dict[str, object] = {
        "kind": "scalar",
        "name": prior.name,
        "distribution": prior.distribution,
        "low": prior.low,
        "high": prior.high,
    }
    eff_lo, eff_hi = prior.effective_support
    if math.isfinite(eff_lo) and math.isfinite(eff_hi):
        record["effective_low"] = eff_lo
        record["effective_high"] = eff_hi
    return record


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
