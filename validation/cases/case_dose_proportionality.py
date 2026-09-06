"""Dose-proportionality: IV AUC scales ~linearly across 10 and 40 mg.

Cross-stage consistency check (L1, doc/08 §1.3): certifies linear-PBPK
self-consistency, no human data involved.
"""

from __future__ import annotations

from validation.benchmarks import BENCHMARKS
from validation.cases.base import (
    CaseResult,
    EvidenceLevel,
    MetricResult,
    _build_model,
)

from drugos.inputs.models import DosePlan
from drugos.pk.simulate import simulate_pbpk


def case_dose_proportionality() -> CaseResult:
    b = next(x for x in BENCHMARKS if x.name == "midazolam")
    aucs: list[float] = []
    for dose in (10.0, 40.0):
        model = _build_model(b)
        model.dose_plan = DosePlan.iv_bolus(dose_mg=dose)
        res = simulate_pbpk(model, tmax_h=b.tmax_h, n_eval=b.n_eval)
        aucs.append(res.pk_metrics().auc_inf_mgh_l)
    ratio = aucs[1] / aucs[0]
    ok = 3.6 <= ratio <= 4.4
    met = MetricResult("auc_40mg_over_10mg", ratio, 3.6, 4.4, "ratio", "pass" if ok else "FAIL")
    return CaseResult(
        "dose-proportionality (midazolam IV)",
        ok,
        [met],
        [f"AUC10={aucs[0]:.3f}, AUC40={aucs[1]:.3f} mg.h/L (linearity ~{ratio:.2f})"],
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_dose_proportionality"]
