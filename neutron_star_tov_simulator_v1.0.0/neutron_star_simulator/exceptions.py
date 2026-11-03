"""Explicit failures: unsuccessful integrations must never become valid stars."""


class SimulatorError(RuntimeError):
    """Base class for physical or numerical simulation failures."""


class EOSDomainError(SimulatorError, ValueError):
    """EOS input outside a supported physical domain."""


class CausalityError(SimulatorError):
    """Sound speed exceeds c or the EOS has a nonpositive compressibility."""


class IntegrationError(SimulatorError):
    """The integrator failed or encountered an inadmissible geometry."""


class SurfaceNotFoundError(IntegrationError):
    """No zero-pressure surface was reached within max_radius_m."""
