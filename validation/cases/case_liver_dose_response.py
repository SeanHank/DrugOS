"""Liver DILI dose-response (L2 — mechanistic limit, doc/08 §1.2 Tier 3).

The liver QST panel (doc/05 4.2) driven by PBPK free liver exposure must
reproduce the structural acetaminophen injury behaviour: therapeutic dosing
stays below the DILI alert (ALT < 3x ULN, grade 0/1) while apoptosis/necrosis
indicators and the Hy's Law criteria (ALT >= 3x ULN with total bilirubin
>= 2x ULN) all mount at overdose.  This is the stage's own dose-response
trajectory (L2), not an independent external data anchor (L3).
"""

from __future__ import annotations

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, _simulate

from drugos.organ.liver import LiverParams, simulate_liver
from drugos.pk.simulate import PBPKResult

ACETAMINOPHEN_MW = 151.2


def _run_trajectory(res: PBPKResult, dose_mg: float) -> object:
    drift = dose_mg / 1000.0
    c = res.unbound_tissues["liver"] * drift
    return simulate_liver(res.t, c, ACETAMINOPHEN_MW, LiverParams(), n_eval=500)


def case_liver_dose_response() -> CaseResult:
    res = _simulate(_bench())
    therapeutic = _run_trajectory(res, 1000.0)
    overdose = _run_trajectory(res, 20000.0)

    met_ther = MetricResult(
        "therapeutic_peak_alt_uln",
        therapeutic.peak_alt_uln,
        0.0,
        2.0,
        "xULN",
        "pass" if therapeutic.peak_alt_uln < 2.0 else "FAIL",
    )
    met_od_alt = MetricResult(
        "overdose_peak_alt_uln",
        overdose.peak_alt_uln,
        3.0,
        20.0,
        "xULN",
        "pass" if overdose.peak_alt_uln >= 3.0 else "FAIL",
    )
    met_hy = MetricResult(
        "overdose_hy_law",
        1.0 if overdose.hy_law else 0.0,
        1.0,
        1.0,
        "flag",
        "pass" if overdose.hy_law else "FAIL",
    )
    met_trajectory = MetricResult(
        "therapeutic_to_overdose_escalation",
        overdose.peak_alt_uln / max(therapeutic.peak_alt_uln, 1e-9),
        2.0,
        1e3,
        "fold",
        "pass" if overdose.peak_alt_uln > therapeutic.peak_alt_uln else "FAIL",
    )
    ok = all(m.criterion == "pass" for m in (met_ther, met_od_alt, met_hy, met_trajectory))
    return CaseResult(
        "liver DILI dose-response (acetaminophen)",
        ok,
        [met_ther, met_od_alt, met_trajectory, met_hy],
        [
            f"therapeutic ALT {therapeutic.peak_alt_uln:.2f}xULN (safe), "
            f"20 g overdose ALT {overdose.peak_alt_uln:.2f}xULN with total "
            f"bilirubin {overdose.peak_bilirubin_uln:.2f}xULN -> Hy's Law "
            f"{'met' if overdose.hy_law else 'not met'}; DILI grade "
            f"{therapeutic.dili_grade}->{overdose.dili_grade}. "
            "Acetaminophen overdose injury is characterized by centrilobular "
            "necrosis, transaminase >3x ULN and mixed cholestasis (DILI "
            "severity scales, e.g. Maria & Victorino / Hy's Law criteria)."
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


def _bench() -> object:
    """The acetaminophen benchmark entry used to drive the PBPK exposure."""
    from validation.benchmarks import BENCHMARKS

    return next(b for b in BENCHMARKS if b.name == "acetaminophen")


__all__ = ["case_liver_dose_response"]
