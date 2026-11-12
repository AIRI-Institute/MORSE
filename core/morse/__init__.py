"""Full-length mRNA vaccine sequence optimization."""

__version__ = "0.4.0"

__all__ = [
    "load_config", "validate", "score", "optimize", "load_run", "resume",
    "OptimizationResult", "Candidate",
]


def __getattr__(name):
    if name in {"OptimizationResult", "Candidate"}:
        from morse import results
        value = getattr(results, name)
    elif name in __all__:
        from morse import api
        value = getattr(api, name)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
