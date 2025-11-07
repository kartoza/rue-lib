"""Custom exceptions for RUE."""


class RUEError(Exception):
    """Base exception for RUE."""

    pass


class ProjectionError(RUEError):
    """Raised when projection operations fail."""

    pass


class GeometryError(RUEError):
    """Raised when geometry operations fail."""

    pass


class ValidationError(RUEError):
    """Raised when input validation fails."""

    pass
