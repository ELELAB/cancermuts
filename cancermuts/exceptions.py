class UnexpectedIsoformError(Exception):
    """Raised when a non-canonical isoform is used in a context that requires the canonical sequence."""
    pass
