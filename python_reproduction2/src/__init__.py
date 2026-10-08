"""lpw -- PWRCTV (TGRS 2024) numba-accelerated Python backend.

A self-contained subset synced from ``python_reproduction``: only the
PWRCTV method, only the numba (JIT) implementation.
"""

from .solvers import pwrctv, warm_up

__all__ = ["pwrctv", "warm_up"]

