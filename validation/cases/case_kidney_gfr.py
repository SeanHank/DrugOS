"""Kidney GFR/AKI panel (L2 — analytic / mechanistic limit, doc/08 §1.2).

Two assertions: (1) the closed-form creatinine balance Scr = P/GFR is
reproduced exactly by the nephron panel, and (2) escalating free-kidney
exposure escalates the KDIGO AKI grade monotonically through the floor-bounded
GFR drop.
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.organ.kidney import KidneyParams, simulate_kidney


def _run(exposure_mg_l: float) -> object:
    t = np.linspace(0.0, 24.0, 25)
    c = np.full_like(t, exposure_mg_l)
    params = KidneyParams(injury_ic50_nm=1.0e3, injury_hill=2.0)
    return simulate_kidney(t, c, 300.0, gfr_base_ml_min=120.0, scr_base_umol_l=80.0, params=params)


def case_kidney_gfr() -> CaseResult:
    zero = _run(0.0)
    low = _run(0.01)
    high = _run(1.0)

    met_closed = MetricResult(
        "scr_at_zero_fn",
        zero.peak_scr_ratio,
        1.0,
        1.0,
        "ratio",
        "pass" if abs(zero.peak_scr_ratio - 1.0) < 1e-9 else "FAIL",
    )
    met_low = MetricResult(
        "low_exposure_scr_ratio",
        low.peak_scr_ratio,
        1.0,
        1.05,
        "ratio",
        "pass" if 1.0 <= low.peak_scr_ratio <= 1.05 else "FAIL",
    )
    met_high = MetricResult(
        "high_exposure_scr_ratio",
        high.peak_scr_ratio,
        2.0,
        8.0,
        "ratio",
        "pass" if 2.0 <= high.peak_scr_ratio <= 8.0 else "FAIL",
    )
    met_escalation = MetricResult(
        "aki_grade_escalation",
        high.aki_grade - low.aki_grade,
        1.0,
        3.0,
        "grade",
        "pass" if high.aki_grade > low.aki_grade else "FAIL",
    )
    ok = all(m.criterion == "pass" for m in (met_closed, met_low, met_high, met_escalation))
    return CaseResult(
        "kidney GFR/AKI escalation (KDIGO)",
        ok,
        [met_closed, met_low, met_high, met_escalation],
        [
            f"closed-form Scr=P/GFR exact at zero exposure; "
            f"1.0 mg/L free kidney exposure -> Scr ratio "
            f"{high.peak_scr_ratio:.2f} (KDIGO stage {high.aki_grade}) vs "
            f"{low.peak_scr_ratio:.3f} (stage {low.aki_grade}); GFR floor "
            f"{high.min_gfr_ml_min:.0f} mL/min.  KDIGO criteria: Scr x2 -> "
            "stage 2, x3 -> stage 3 (or GFR drop)."
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_kidney_gfr"]
