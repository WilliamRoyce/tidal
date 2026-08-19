"""Custom exceptions for TIDAL solvers."""


class SimulationDivergedError(RuntimeError):
    """Raised when simulation fields become non-finite or exceed a norm threshold.

    This indicates a physical or numerical instability — e.g. a coupled mass
    matrix with a negative eigenvalue causing exponentially growing modes.
    """


class KineticEvaluationError(RuntimeError):
    """Raised when a ``kinetic_coefficient_symbolic`` cannot be resolved to a number.

    The modal builders need a concrete ``M`` to populate the mass matrix.
    When evaluation fails — an unbound symbol (a parameter missing from
    ``--param``), an unsupported construct, or a NaN/Inf result — there is
    no correct fallback: proceeding with ``M = 1`` silently changes the
    physics (GH #447). Raised instead of a bare ``ValueError`` on purpose:
    the conversion-stability probe catches ``(LinAlgError, ValueError)``
    and converts it into a "tachyonic" verdict, which would mislabel a
    configuration error as a physics result.
    """
