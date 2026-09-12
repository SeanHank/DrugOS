"""Sympathetic-suppression (beta-like) cardiac branch (L2, doc/05 4.3).

An earlier roadmap item in doc/05 4.3 asked for exactly this: the shipped
ERK/MAPK readout encodes only the *stimulatory* pathway tone, so a drug that
suppresses sympathetic drive (bradycardia + negative inotropy, e.g.
beta-adrenergic site blockade) was not representable.  That branch is now
``CardiacParams.sympathetic_tone`` (model-layer default-off at 1.0;
auto-anchored by the pipeline in full fidelity): a residual
tone in (0, 1] computed from a saturable Emax axis
``tone = IC50/(IC50 + C_free)`` (``sympathetic_tone_from_emax``) and
applied multiplicatively to heart rate *and* stroke volume in
``simulate_hemodynamics``, so a beta-like blockade that halves the tone
quarters cardiac output.

The checks pin, on the analytic Windkessel with a fixed systemic resistance:

- at IC50 free exposure (C = IC50) tone = 1/2, HR halves, SV halves, and
  CO drops to exactly Q/4 with MAP = CVP + (Q/4)/(Q) * (MAP - CVP);
- zero exposure (tone 1.0) preserves the baseline load (Q, MAP 93) and the
  branch is model-layer default-off (null-effect anchor auto-assigned
  in full fidelity);
- a saturating exposure drives tone -> 0 and CO/MAP toward the CVP floor;
- the QST output is monotone in the exposure/(IC50) ratio;
- degenerate inputs (IC50 <= 0, tone out of (0, 1]) raise instead of
  corrupting the ODE.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.organ.cardiac import (
    CardiacParams,
    simulate_hemodynamics,
    sympathetic_tone_from_emax,
)


def _q(ton: float, base: SimpleNamespace) -> float:
    """Cardiac output under residual sympathetic tone ``ton`` (CO ~ tone^2)."""
    return base.co_l_min * ton * ton


def case_cardiac_sympathetic_suppression() -> CaseResult:
    phy = profile()
    base = simulate_hemodynamics(phy, CardiacParams())
    cvp = 5.0
    metrics: list[MetricResult] = []
    notes: list[str] = []

    # 1. IC50 exposure halves the tone and quarters CO, exactly.
    ton_ic50 = sympathetic_tone_from_emax(50.0, 50.0)
    ok1 = math.isclose(ton_ic50, 0.5, abs_tol=1e-12)
    r_ic50 = simulate_hemodynamics(phy, CardiacParams(sympathetic_tone=ton_ic50))
    ok1b = math.isclose(r_ic50.hr_bpm, base.hr_bpm * 0.5, rel_tol=1e-9) and math.isclose(
        r_ic50.sv_ml, base.sv_ml * 0.5, rel_tol=1e-9
    )
    ok1c = math.isclose(r_ic50.co_l_min, _q(ton_ic50, base), rel_tol=1e-9) and math.isclose(
        r_ic50.pa_mmhg[-1] - r_ic50.pv_mmhg[-1], 0.25 * (base.map_mmhg - cvp), rel_tol=1e-9
    )
    metrics.append(
        MetricResult(
            "ic50_exposure_quarters_cardiac_output",
            r_ic50.co_l_min / base.co_l_min,
            0.25,
            0.25,
            "CO/CO_base @ IC50",
            "pass" if (ok1 and ok1b and ok1c) else "FAIL",
        )
    )

    # 2. Zero exposure is baseline and the branch is off by default.
    ok2a = sympathetic_tone_from_emax(0.0, 50.0) == 1.0
    ok2b = math.isclose(base.co_l_min * 1000.0, phy.cardiac_output_ml_min, rel_tol=1e-9)
    ok2c = CardiacParams().sympathetic_tone == 1.0
    metrics.append(
        MetricResult(
            "off_by_default_preserves_baseline",
            base.map_mmhg,
            93.0,
            93.0,
            "MAP mmHg (no anchor)",
            "pass" if (ok2a and ok2b and ok2c) else "FAIL",
        )
    )

    # 3. Saturating exposure drives tone -> 0 and the pressure drop to zero
    #    (pa falls toward CVP on the fixed-resistance manifold).
    ton_sat = sympathetic_tone_from_emax(1.0e6, 1.0)
    r_sat = simulate_hemodynamics(phy, CardiacParams(sympathetic_tone=ton_sat))
    ok3a = ton_sat < 1.0e-4
    ok3b = r_sat.co_l_min < 1.0e-3 * base.co_l_min
    ok3c = (r_sat.pa_mmhg[-1] - r_sat.pv_mmhg[-1]) < 0.1
    metrics.append(
        MetricResult(
            "saturating_exposure_cvp_floor",
            r_sat.pa_mmhg[-1],
            cvp,
            cvp + 1.0,
            "arterial pressure mmHg @ saturating blockade",
            "pass" if (ok3a and ok3b and ok3c) else "FAIL",
        )
    )

    # 4. Monotone dose-response in the exposure/(IC50) ratio.
    ratios = [0.0, 0.1, 0.5, 1.0, 2.0, 10.0]
    co_seq = [
        simulate_hemodynamics(
            phy, CardiacParams(sympathetic_tone=sympathetic_tone_from_emax(r, 1.0))
        ).co_l_min
        for r in ratios
    ]
    ok4 = (
        all(co_seq[i] < co_seq[i - 1] for i in range(1, len(co_seq))) and co_seq[0] == base.co_l_min
    )
    metrics.append(
        MetricResult(
            "monotone_exposure_response",
            co_seq[-1],
            0.0,
            base.co_l_min * (1.0 / 11.0) ** 2,
            "CO @ C/IC50 = 10",
            "pass" if ok4 else "FAIL",
        )
    )

    # 5. Degenerate inputs raise instead of corrupting the ODE.
    raised = [0.0, 0.0]
    try:
        sympathetic_tone_from_emax(1.0, 0.0)
    except ValueError:
        raised[0] = 1.0
    try:
        simulate_hemodynamics(phy, CardiacParams(sympathetic_tone=1.5))
    except ValueError:
        raised[1] = 1.0
    metrics.append(
        MetricResult(
            "degenerate_sympathetic_inputs_rejected",
            sum(raised),
            2.0,
            2.0,
            "count of rejected probes",
            "pass" if sum(raised) == 2.0 else "FAIL",
        )
    )

    ok = all(m.criterion == "pass" for m in metrics)
    notes.append(
        f"C=IC50 -> CO/CO_base={r_ic50.co_l_min / base.co_l_min:.4f} "
        f"(target 0.25); saturating MAP={r_sat.map_mmhg:.2f} mmHg (floor {cvp}); "
        f"monotone over C/IC50 in {ratios}"
    )
    return CaseResult(
        "Sympathetic suppression branch (beta-like Emax on HR and SV)",
        ok,
        metrics,
        notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_cardiac_sympathetic_suppression"]
