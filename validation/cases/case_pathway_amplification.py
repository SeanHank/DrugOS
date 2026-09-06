"""Pathway amplification (L2 — mechanistic limit, doc/08 §1.2, Tier 3).

With ultrasensitive Hill coupling (doc/05 3.2-3.5) the steady-state readout
EC50 (fit of R(s)=r0+(emax-r0)s^n/(ec50^n+s^n) to constant signal levels)
falls strictly below s=0.5, whereas a one-step mass-action cascade would need
to be driven beyond half-saturation.  Since the occupancy signal is drug/Kd,
EC50<0.5 is the Tier-2 "EC50 below the receptor Kd" amplification criterion.
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pathway import dose_response, mapk_cascade


def case_pathway_amplification() -> CaseResult:
    model = mapk_cascade()
    signals = np.linspace(0.0, 1.0, 21)
    fit = dose_response(model, signals, readout="erk_active", n_eval=200)
    fold = fit.emax / max(fit.r0, 1e-9)
    met_fold = MetricResult(
        "emax_over_baseline_fold",
        fold,
        2.0,
        1e4,
        "ratio",
        "pass" if 2.0 <= fold <= 1e4 else "FAIL",
    )
    met_ec50 = MetricResult(
        "ec50_signal",
        fit.ec50,
        0.0,
        0.5,
        "fraction",
        "pass" if fit.ec50 < 0.5 else "FAIL",
    )
    drug_ec50 = fit.drug_ec50_nm(1.0)
    met_drug = MetricResult(
        "drug_ec50_nm_at_kd_1",
        drug_ec50,
        0.0,
        0.5,
        "nM",
        "pass" if drug_ec50 < 0.5 else "FAIL",
    )
    ok = (
        met_fold.criterion == "pass"
        and met_ec50.criterion == "pass"
        and met_drug.criterion == "pass"
    )
    return CaseResult(
        "pathway amplification (MAPK EC50<Kd equivalent)",
        ok,
        [met_fold, met_ec50, met_drug],
        [
            f"r0={fit.r0:.4g}, emax={fit.emax:.4g}, ec50={fit.ec50:.4g}, "
            f"hill={fit.hill:.3g}; EC50 signal 1e-3 << half-saturation 0.5"
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_pathway_amplification"]
