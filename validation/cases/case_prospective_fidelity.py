"""Prospective rerun fidelity (D24 — L1, doc/08 §1.4 runbook).

Phase-6 D24 stands up a prospective re-validation workflow: any pipeline
production run must be reproducible from its fixed-sim seed, and the verdict
must be stable when repeated on an independently-sampled but admissible
profile (the 'held-out subject' surrogate).  This case certifies:

- test-retest determinism: two identical runs produce bit-identical verdicts
  and risks (the runbook's precondition for any prospective claim),
- held-out-subject stability: the same dose on a distinct population profile
  keeps the composite verdict and its driver,
- banded dofetilide risk remains in the QT-high regime with a female elderly
  profile (age/sex tolerance check).

Weakest epistemic weight L1: it certifies measurability/reproducibility, not
human predictivity (that comes from the L3 Tier-1 cases).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from validation.benchmarks import BENCHMARKS
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.inputs.models import HumanProfile
from drugos.pipeline import run_pipeline, spec_from_benchmark_data


def _bench(name: str):
    return next(b for b in BENCHMARKS if b.name == name)


def case_prospective_fidelity() -> CaseResult:
    spec = spec_from_benchmark_data(_bench("dofetilide"))

    r1 = run_pipeline(spec)
    r2 = run_pipeline(spec)
    risks1 = np.asarray([risk.risk for risk in r1.toxicity.risks])
    risks2 = np.asarray([risk.risk for risk in r2.toxicity.risks])
    max_risk_diff = float(np.max(np.abs(risks1 - risks2)))
    deterministic = r1.verdict == r2.verdict and max_risk_diff == 0.0

    held_out = run_pipeline(
        replace(
            spec,
            profile=HumanProfile(sex="female", age_y=70.0, weight_kg=60.0, height_cm=158.0),
        )
    )
    held_qt = next(r.risk for r in held_out.toxicity.risks if r.endpoint.value == "qt")
    verdict_stable = (
        "High composite risk" in r1.verdict
        and r1.toxicity.overall_driver() == "qt"
        and "High composite risk" in held_out.verdict
        and held_out.toxicity.overall_driver() == "qt"
    )
    regime_ok = 0.70 <= held_qt <= 0.99

    metrics = [
        MetricResult(
            "test_retest_max_risk_diff",
            max_risk_diff,
            0.0,
            0.0,
            "P(risk)",
            "pass" if deterministic else "FAIL",
        ),
        MetricResult(
            "held_out_subject_qt_risk",
            held_qt,
            0.70,
            0.99,
            "P(risk)",
            "pass" if regime_ok else "FAIL",
        ),
        MetricResult(
            "held_out_verdict_disagreements",
            0.0 if verdict_stable else 1.0,
            0.0,
            0.0,
            "count",
            "pass" if verdict_stable else "FAIL",
        ),
    ]
    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "D24 prospective rerun fidelity (dofetilide QTc)",
        ok,
        metrics,
        [
            f"Repeated identical runs agree to {max_risk_diff:.1e} in risk and "
            f"keep verdict '{r1.verdict}'; an independent female-70 profile "
            f"also sustains the high-QT regime (dofetilide QT {held_qt:.3f}, "
            "driver qt).  Basis: reproducibility is the precondition of the "
            "runbook; the QTc band itself is anchored by the L3 dofetilide "
            "Tier-1 case (see case_cardiac_qtc)."
        ],
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_prospective_fidelity"]
