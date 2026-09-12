"""Multi-ionic ORd-CiPA-v1-2017 QTw lane (doc/12, L19 realized, L3 evidence).

Runs the vendored O'Hara-Rudy CiPA-v1-2017 human ventricular action-potential
model (Dutta et al., Front Physiol 2017; BSD Myokit encoding at
``data/models/ohara-cipa-v1-2017.mmt``) under **concurrent fractional block of
up to six current systems** (IKr, IKs, IK1, Ito, INaL, ICaL) and reports the
net APD90 effect — the realized multi-ionic upgrade of the R-3 cross-check,
which only ever scaled the single IKr conductance of the 2011 model.

Asserted physics (net effect computed by the action-potential model, never
assumed from the algebraic encoder):

- the un-blocked APD90 at 1 Hz is physiological (published endocardial band),
- a strong hERG block prolongs APD90 strongly (as R-2/R-3 establish for
  dofetilide), so the CiPA lane reproduces the same direction as ORd-2011 and
  the calibrated encoder,
- late-INa co-block *relieves* part of the hERG-driven prolongation (the CiPA
  multi-current net-signal logic: blocking the plateau late-sodium influx
  shortens APD), so the co-blocked delta is strictly below the hERG-only
  delta,
- the combined IKr+INaL+ICaL profile still prolongs (net positive risk
  signal), strictly below the hERG-only value but unambiguously proarrhythmic.

Runtime requirements: myokit + SUNDIALS headers (doc/06, doc/12).
"""

from __future__ import annotations

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.organ.cardiac_ap import cipa_apd90

_BASE_LO_MS = 200.0
_BASE_HI_MS = 350.0


def case_ord_multi_ionic_qt() -> CaseResult:
    base = cipa_apd90({})
    ikr = cipa_apd90({"IKr": 0.5})
    ikr_inal = cipa_apd90({"IKr": 0.5, "INaL": 0.3})
    combined = cipa_apd90({"IKr": 0.5, "INaL": 0.15, "ICaL": 0.15})

    base_ms = base.apd90_base_ms
    ikr_delta = ikr.delta_apd90_ms
    inal_delta = ikr_inal.delta_apd90_ms
    comb_delta = combined.delta_apd90_ms

    base_ok = _BASE_LO_MS <= base_ms <= _BASE_HI_MS
    ikr_strong = ikr_delta >= 60.0
    late_na_relieves = inal_delta < ikr_delta - 5.0
    combined_net = 30.0 <= comb_delta < ikr_delta - 5.0

    metrics = [
        MetricResult(
            "cipa_apd90_base_ms",
            base_ms,
            _BASE_LO_MS,
            _BASE_HI_MS,
            "ms",
            "pass" if base_ok else "FAIL",
        ),
        MetricResult(
            "IKr_0.5_delta_apd90_ms",
            ikr_delta,
            60.0,
            1000.0,
            "ms",
            "pass" if ikr_strong else "FAIL",
        ),
        MetricResult(
            "IKr0.5_INaL0.3_delta_apd90_ms",
            inal_delta,
            0.0,
            ikr_delta - 5.0,
            "ms",
            "pass" if late_na_relieves else "FAIL",
        ),
        MetricResult(
            "IKr0.5_INaL0.15_ICaL0.15_delta_apd90_ms",
            comb_delta,
            30.0,
            ikr_delta - 5.0,
            "ms",
            "pass" if combined_net else "FAIL",
        ),
    ]
    notes = [
        f"ORd-CiPA-v1-2017 (myokit, endo, 30 s pre-pace @1 Hz): baseline "
        f"APD90={base_ms:.1f} ms; "
        f"delta-APD90 hERG-only (50% IKr block)={ikr_delta:.1f} ms; "
        f"with 30% late-INa co-block={inal_delta:.1f} ms; "
        f"combined IKr+INaL+ICaL block={comb_delta:.1f} ms; "
        "late-INa/ICaL co-block relieves part of the hERG-driven prolongation "
        "(net multi-current effect computed by the CiPA action-potential "
        "model, mirroring the CiPA score logic); direction matches the "
        "algebraic encoder and the R-3 ORd-2011 anchor"
    ]
    return CaseResult(
        benchmark="multi-ionic QTw cross-check (ORd-CiPA-v1-2017)",
        passed=all(m.criterion == "pass" for m in metrics),
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L3_EMPIRICAL,
    )
