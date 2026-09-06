"""Stage-5 clinical fusion — analytic point-match (L2, doc/08 §1.2).

Derives the clinical layer's own equations (doc/05 5.2) as closed-form
limits and checks the shipped code reproduces them exactly:

- the exposure-ratio line `P = sigmoid(slope * (log10(Cmax/IC50) - center))`,
- the log-odds posterior at the pure prior (no evidence, DILI prior 0.25),
- the mechanistic grade->risk table (grade k maps to RULES[].grade_probs[k]),
- the CTCAE-style higher/lower threshold ladders of ``grade_value``,
- time-to-onset / duration from threshold crossing.

No external data: this is the analytic/mechanistic-limit tier for the
clinical scoring layer introduced in Stage 5 (D18-D20).
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.clinical.biomarkers import BiomarkerSpec, crossing_interval, grade_value
from drugos.clinical.toxicity import (
    RULES,
    Endpoint,
    exposure_evidence,
    fuse,
    mechanistic_evidence,
    sigmoid,
)


def _grade_ladder_higher(value: float) -> int:
    return grade_value(value, BiomarkerSpec("x", "u", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0)))


def _grade_ladder_lower(value: float) -> int:
    return grade_value(value, BiomarkerSpec("x", "u", 0.0, 1.0, "lower", (90.0, 60.0, 45.0, 15.0)))


def case_clinical_grading() -> CaseResult:
    qt_rule = RULES[Endpoint.QT]
    expected_line = sigmoid(qt_rule.ratio_slope * (0.0 - qt_rule.ratio_center))
    ev = exposure_evidence(Endpoint.QT, cmax_unbound_nm=100.0, ic50_nm=100.0)
    line_ok = abs(ev.probability - expected_line) < 1e-12

    dili_prior = RULES[Endpoint.DILI].prior
    fused_prior = fuse(dili_prior, [])
    prior_ok = abs(fused_prior[0] - dili_prior) < 1e-12

    mech = mechanistic_evidence(Endpoint.DILI, 3)
    mech_ok = abs(
        mech.probability - RULES[Endpoint.DILI].grade_probs[3]
    ) < 1e-12 and mech.note.startswith("mechanistic")

    hi_grades = [_grade_ladder_higher(v) for v in (0.5, 1.0, 2.0, 3.0, 10.0)]
    hi_ok = hi_grades == [0, 1, 2, 3, 4]
    lo_grades = [_grade_ladder_lower(v) for v in (100.0, 90.0, 60.0, 45.0, 15.0)]
    lo_ok = lo_grades == [0, 1, 2, 3, 4]

    spec = BiomarkerSpec("t", "h", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0))
    cross_interval = crossing_interval(
        np.asarray([0.0, 1.0, 2.0, 3.0, 4.0]),
        np.asarray([0.5, 1.5, 2.5, 0.5, 0.5]),
        spec,
        grade=1,
    )
    onset, duration = cross_interval
    crossing_ok = onset == 1.0 and abs(duration - 1.0) < 1e-12

    metrics = [
        MetricResult(
            "exposure_line_probability",
            ev.probability,
            expected_line,
            expected_line,
            "P(risk)",
            "pass" if line_ok else "FAIL",
        ),
        MetricResult(
            "prior_only_posterior",
            fused_prior[0],
            dili_prior,
            dili_prior,
            "P(risk)",
            "pass" if prior_ok else "FAIL",
        ),
        MetricResult(
            "mechanistic_grade3_probability",
            mech.probability,
            RULES[Endpoint.DILI].grade_probs[3],
            RULES[Endpoint.DILI].grade_probs[3],
            "P(risk)",
            "pass" if mech_ok else "FAIL",
        ),
        MetricResult(
            "grade_ladder_higher",
            float(hi_grades[2]),
            2.0,
            2.0,
            "grade",
            "pass" if hi_ok else "FAIL",
        ),
        MetricResult(
            "grade_ladder_lower",
            float(lo_grades[2]),
            2.0,
            2.0,
            "grade",
            "pass" if lo_ok else "FAIL",
        ),
        MetricResult(
            "crossing_onset_duration",
            onset if onset is not None else float("nan"),
            1.0,
            1.0,
            "h",
            "pass" if crossing_ok else "FAIL",
        ),
    ]
    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "Stage-5 clinical grading: analytic point-matches",
        ok,
        metrics,
        [
            "Exposure ROC line reproduces the closed form "
            f"sigmoid({qt_rule.ratio_slope}*({0.0} - ({qt_rule.ratio_center:.1f}))) = "
            f"{expected_line:.4f}; empty-evidence fusion returns the DILI prior "
            f"{dili_prior:.2f} exactly; grade ladder and crossing windows match "
            "the CTCAE conventions of doc/05 5.1-5.2."
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_clinical_grading"]
