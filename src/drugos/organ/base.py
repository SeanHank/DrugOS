"""Shared helpers for the Stage-4 organ panels (doc/05 section 4).

Concentration driver for every organ panel is the *free* (unbound) tissue
concentration delivered by Stage-1 PBPK (``PBPKResult.unbound_tissues``),
converted to nanomolar with the compound molecular weight — the same binding
driver used by Stage 2 occupancy — so Stage 2/4 comparers stay on one unit.
"""

from __future__ import annotations

import numpy as np

from drugos.target.occupancy import free_binding_conc_mg_l_to_nm as _to_nm

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]
AUC = float


def free_mg_l_to_nm(c_free_mg_l: NDArray | float, mw: float) -> NDArray:
    """Free tissue concentration (mg/L) -> nanomolar, reusing Stage-2 units."""
    return _to_nm(c_free_mg_l, mw)


def auc_nm_h(t_h: NDArray, c_free_nm: NDArray) -> float:
    """Trapezoidal exposure integral (nM*h) over a free-drug trace."""
    t = np.asarray(t_h, dtype=float)
    c = np.asarray(c_free_nm, dtype=float)
    if t.ndim != 1 or c.ndim != 1 or t.shape[0] != c.shape[0]:
        raise ValueError("t_h and c_free_nm must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")
    return float(np.trapezoid(c, t))


__all__ = ["AUC", "NDArray", "auc_nm_h", "free_mg_l_to_nm"]
