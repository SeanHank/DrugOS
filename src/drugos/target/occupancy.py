"""Target occupancy simulation with turnover (doc/05 section 2.4).

Implements the classic turnover model (Daryaee & Tonge 2019) with reversible
binding plus complex internalization:

    dR/dt   = ksyn - rho*R - kon*D*R + koff*DR          (ksyn = rho*R0)
    dDR/dt  = kon*D*R - koff*DR - kint*DR

driven by the *free* (unbound) tissue concentration time course coming from
Stage 1 PBPK (converted from mg/L to nM using the compound molecular weight).
Fractional occupancy = DR / Rtot(0); the occupancy-integral ("time-at-target")
is the exposure signal consumed by Stage 3.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from drugos.target.targets import Target

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]

_NG_PER_MG = 1e6


def free_binding_conc_mg_l_to_nm(c_free_mg_l: NDArray | float, mw: float) -> NDArray:
    """Convert free tissue concentration (mg/L) to nanomolar binding driver."""
    if mw <= 0:
        raise ValueError("mw must be positive")
    return np.asarray(c_free_mg_l, dtype=float) * _NG_PER_MG / mw


@dataclass(slots=True)
class TargetOccupancyResult:
    """Occupancy trajectory for one binding site.

    ``occupancy`` is the fractional receptor occupancy DR/Rtot(0) on ``t_h``.
    ``time_at_target_h`` is the occupancy integral (area under the occupancy
    curve) — the exposure signal forwarded to Stage 3.
    """

    target: Target
    t_h: NDArray
    free_nm: NDArray
    receptor_nm: NDArray
    complex_nm: NDArray
    occupancy: NDArray
    time_at_target_h: float
    peak_occupancy: float
    solver: str = "LSODA"
    n_eval: int = 0

    def to_series(self) -> dict[str, NDArray]:
        """Stage-2 slice of the pipeline data contract for downstream stages."""
        return {
            "t_h": self.t_h,
            "free_nM": self.free_nm,
            "receptor_nM": self.receptor_nm,
            "complex_nM": self.complex_nm,
            "occupancy": self.occupancy,
        }


def simulate_occupancy(
    t_h: NDArray,
    c_free_mg_l: NDArray,
    target: Target,
    mw: float,
    r0_nm: float | None = None,
    n_eval: int = 601,
    rtol: float = 1e-8,
    atol: float = 1e-9,
) -> TargetOccupancyResult:
    """Integrate target occupancy driven by a free-concentration time course.

    Args:
        t_h: simulation time grid (h), must match ``c_free_mg_l``.
        c_free_mg_l: free (unbound) drug concentration, mg/L.
        target: the binding site to simulate.
        mw: compound molecular weight (g/mol) for nM conversion.
        r0_nm: optional basal abundance override (else ``target.r0_nm``).
        n_eval: number of output samples along the grid for re-evaluation.
    """
    t = np.asarray(t_h, dtype=float)
    free = free_binding_conc_mg_l_to_nm(c_free_mg_l, mw)
    if t.ndim != 1 or free.ndim != 1 or t.shape[0] != free.shape[0]:
        raise ValueError("t_h and c_free_mg_l must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")

    r0 = target.r0_nm if r0_nm is None else r0_nm
    if r0 <= 0:
        raise ValueError("r0_nm must be positive")
    ksyn = target.rho_h * r0
    kon = target.kon_nm_h
    koff = target.koff_1h

    def diff_drive(sol_t: float) -> float:
        return float(np.interp(sol_t, t, free))

    def rhs(sol_t: float, y: NDArray) -> NDArray:
        d = diff_drive(sol_t)
        r, dr = float(y[0]), float(y[1])
        return np.array(
            [
                ksyn - target.rho_h * r - kon * d * r + koff * dr,
                kon * d * r - koff * dr - target.kint_h * dr,
            ],
            dtype=float,
        )

    y0 = np.array([r0, 0.0], dtype=float)
    t_eval = np.linspace(t[0], t[-1], n_eval)
    sol = solve_ivp(
        rhs,
        (t[0], t[-1]),
        y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"occupancy solve failed: {sol.message}")
    receptor = sol.y[0]
    complexed = sol.y[1]
    base = np.full_like(receptor, r0, dtype=float)
    occupancy = np.clip(np.divide(complexed, base, out=np.zeros_like(complexed)), 0.0, 1.0)
    time_at_target = float(np.trapezoid(occupancy, t_eval))
    return TargetOccupancyResult(
        target=target,
        t_h=t_eval,
        free_nm=np.interp(t_eval, t, free),
        receptor_nm=receptor,
        complex_nm=complexed,
        occupancy=occupancy,
        time_at_target_h=time_at_target,
        peak_occupancy=float(np.max(occupancy)),
        n_eval=n_eval,
    )


@dataclass(slots=True)
class PanelEngagement:
    """Aggregated occupancy across the off-target panel."""

    results: dict[str, TargetOccupancyResult] = field(default_factory=dict)

    def ranked_by_time_at_target(self) -> list[tuple[str, TargetOccupancyResult]]:
        rows = sorted(self.results.items(), key=lambda kv: kv[1].time_at_target_h, reverse=True)
        return rows


def simulate_panel(
    t_h: NDArray,
    c_free_mg_l: NDArray,
    panel: tuple[Target, ...],
    mw: float,
    n_eval: int = 601,
) -> PanelEngagement:
    """Simulate every site of the safety panel against one free-drug profile."""
    eng = PanelEngagement()
    for site in panel:
        eng.results[site.name] = simulate_occupancy(t_h, c_free_mg_l, site, mw=mw, n_eval=n_eval)
    return eng


__all__ = [
    "TargetOccupancyResult",
    "PanelEngagement",
    "simulate_occupancy",
    "simulate_panel",
    "free_binding_conc_mg_l_to_nm",
]
