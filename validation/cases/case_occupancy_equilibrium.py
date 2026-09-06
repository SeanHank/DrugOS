"""Occupancy equilibrium (L2 — analytic, doc/08 §1.2, Tier 2).

The occupancy law DR/Rtot = D/(D+Kd) follows from the Stage-2 turnover model
(doc/05 2.4) when target turnover and internalization vanish; the numeric ODE
must reproduce the analytic point value (point-wise match).
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.target import Target, simulate_occupancy


def case_occupancy_equilibrium() -> CaseResult:
    site = Target(name="occupancy-case", kd_nm=100.0, r0_nm=5.0, rho_h=0.0, kint_h=0.0)
    mw = 243.0
    t = np.linspace(0.0, 24.0, 241)
    c = np.full_like(t, 100.0 * mw / 1e6)
    res = simulate_occupancy(t, c, site, mw=mw, n_eval=241)
    pred = float(res.occupancy[-1])
    ok = 0.48 <= pred <= 0.52
    met = MetricResult(
        "occupancy_at_free_equals_Kd",
        pred,
        0.48,
        0.52,
        "fraction",
        "pass" if ok else "FAIL",
    )
    return CaseResult(
        "occupancy equilibrium (drug=Kd)",
        ok,
        [met],
        [f"steady-state occupancy {pred:.3f}; analytic D/(D+Kd) = 0.5"],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_occupancy_equilibrium"]
