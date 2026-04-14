from __future__ import annotations

import numpy as np


def compute_transport_k(d31: np.ndarray, d32: np.ndarray) -> np.ndarray:
    d31 = np.asarray(d31, dtype=float)
    d32 = np.asarray(d32, dtype=float)
    if np.any(d31 == 0.0):
        raise ValueError("D31 must be nonzero to compute transport k")
    return 2.5 - d32 / d31
