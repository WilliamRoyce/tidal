"""Custom exceptions for TIDAL solvers."""


class SimulationDivergedError(RuntimeError):
    """Raised when simulation fields become non-finite or exceed a norm threshold.

    This indicates a physical or numerical instability — e.g. a coupled mass
    matrix with a negative eigenvalue causing exponentially growing modes.
    """
