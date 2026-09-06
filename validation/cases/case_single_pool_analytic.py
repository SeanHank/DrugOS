"""One-compartment analytic limit (L2 — analytic, doc/08 §1.2).

A bolus into a single well-mixed pool: AUC = A0/CL, t1/2 = ln2*V/CL.  The
numeric whole-body ODE must recover the closed form when every tissue has
kp = 1 (no partitioning).
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import (
    CaseResult,
    EvidenceLevel,
    MetricResult,
    _single_pool,
)

from drugos.pk.simulate import simulate_pbpk


def case_single_pool_analytic() -> CaseResult:
    cl = 1.0
    model = _single_pool(cl=cl)
    h = model.physiology
    v = sum(h.organ_volume.values()) + h.arterial_blood_l + h.venous_blood_l
    exp_t12 = np.log(2.0) * v / cl
    res = simulate_pbpk(model, tmax_h=6.0 * exp_t12, n_eval=600)
    m = res.pk_metrics()
    exp_auc = 10.0 / cl
    met_auc = MetricResult(
        "auc_inf_mgh_l",
        m.auc_inf_mgh_l,
        exp_auc * 0.95,
        exp_auc * 1.05,
        "mg.h/L",
        "pass" if 0.95 * exp_auc <= m.auc_inf_mgh_l <= 1.05 * exp_auc else "FAIL",
    )
    met_t12 = MetricResult(
        "t_half_h",
        m.term_half_life_h,
        exp_t12 * 0.9,
        exp_t12 * 1.1,
        "h",
        "pass" if 0.9 * exp_t12 <= m.term_half_life_h <= 1.1 * exp_t12 else "FAIL",
    )
    ok = met_auc.criterion == "pass" and met_t12.criterion == "pass"
    return CaseResult(
        "one-compartment analytic limit",
        ok,
        [met_auc, met_t12],
        [f"V={v:.1f} L; analytic AUC={exp_auc:.2f} mg.h/L, t1/2={exp_t12:.1f} h"],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_single_pool_analytic"]
