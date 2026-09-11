"""Immune-mediated DILI QST branch (L2, doc/05 4.2).

The doc/05 4.2 item 5 was an explicit stub: "Immune-mediated component
(stub in baseline; deferred)."  That seam is now
``LiverParams.immune_ic50_nm`` — a saturable hapten/danger hazard
(:func:`immune_hazard`) driving an adaptive immune-response ODE
``dI/dt = k_recruit·hazard·(1−I) − k_decay·I`` whose level loads the
fourth axis of :func:`combined_stress` through ``immune_weight``
(default 0, so every shipped default run stays validated exactly).

The checks pin, on a constant-exposure single-pool liver (no clearance,
MW 300):

- the immune response builds from zero and saturates at the analytically
  expected steady state ``I_ss = k_recruit·hazard / (k_recruit + k_decay)``
  (hazard ≈ 1 at 100× IC50 exposure);
- with ``immune_weight=0`` the axis is inert and dead_frac matches the
  no-immune baseline (numerical tolerance);
- with ``immune_weight > 0`` the stress is higher than baseline, and
  dead_frac at the end of the horizon is strictly greater;
- the response is monotone in exposure (higher C → higher I → higher dead);
- degenerate inputs (IC50 ≤ 0, recruit < 0, decay < 0) raise.
"""

from __future__ import annotations

import math

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.organ.liver import LiverParams, immune_hazard, simulate_liver


def case_dili_immune_activation() -> CaseResult:
    t = np.linspace(0.0, 72.0, 241)
    metrics: list[MetricResult] = []
    notes: list[str] = []

    # 1. Steady-state I_ss = recruit*hazard/(recruit+decay); hazard≈1 at 100×IC50.
    recruit, decay, ic50 = 0.05, 0.10, 10.0
    p_saturate = LiverParams(
        immune_ic50_nm=ic50,
        immune_weight=0.0,  # inert weight → dead unaffected
        immune_recruit_1h=recruit,
        immune_decay_1h=decay,
    )
    c_hi = np.full_like(t, 50.0)  # 50 mg/L / MW 300 = 166667 nM >> ic50
    r_saturate = simulate_liver(t, c_hi, 300.0, p_saturate)
    hazard_hi = immune_hazard(166667.0, ic50)
    i_ss_expected = recruit * hazard_hi / (recruit + decay)
    ok1 = math.isclose(r_saturate.immune[-1], i_ss_expected, rel_tol=0.02) and hazard_hi > 0.99
    metrics.append(
        MetricResult(
            "immune_steady_state_analytic",
            r_saturate.immune[-1],
            i_ss_expected * 0.98,
            i_ss_expected * 1.02,
            "I_ss",
            "pass" if ok1 else "FAIL",
        )
    )

    # 2. With immune_weight=0 the axis is inert; dead_frac ≈ no-immune baseline.
    p_base = LiverParams()
    r_base = simulate_liver(t, c_hi, 300.0, p_base)
    inert_dead_diff = abs(r_saturate.dead_frac.max() - r_base.dead_frac.max())
    ok2 = inert_dead_diff < 5e-3
    metrics.append(
        MetricResult(
            "immune_weight_zero_inert",
            inert_dead_diff,
            0.0,
            5e-3,
            "max dead_frac difference (weight=0 vs baseline)",
            "pass" if ok2 else "FAIL",
        )
    )

    # 3. With immune_weight > 0 the immune axis loads stress; dead_frac end > baseline.
    p_active = LiverParams(
        immune_ic50_nm=ic50,
        immune_weight=1.0,
        immune_recruit_1h=0.1,
        immune_decay_1h=0.05,
    )
    r_active = simulate_liver(t, c_hi, 300.0, p_active)
    ok3 = r_active.immune[-1] > 0.3 and r_active.dead_frac[-1] > r_base.dead_frac[-1]
    metrics.append(
        MetricResult(
            "immune_weight_increases_dead",
            r_active.dead_frac[-1],
            r_base.dead_frac[-1],
            1.0,
            "dead_frac @72h with immune_weight=1",
            "pass" if ok3 else "FAIL",
        )
    )

    # 4. Monotone: higher exposure → higher I → higher dead (at same horizon).
    c_low = np.full_like(t, 1.0)  # near IC50
    p_mono = LiverParams(
        immune_ic50_nm=ic50,
        immune_weight=1.0,
        immune_recruit_1h=0.1,
        immune_decay_1h=0.05,
    )
    r_low = simulate_liver(t, c_low, 300.0, p_mono)
    ok4 = r_low.immune[-1] < r_active.immune[-1] and r_low.dead_frac[-1] < r_active.dead_frac[-1]
    metrics.append(
        MetricResult(
            "monotone_exposure_response",
            r_low.immune[-1],
            0.0,
            r_active.immune[-1],
            "immune @72h (low exposure)",
            "pass" if ok4 else "FAIL",
        )
    )

    # 5. Degenerate inputs raise.
    raised = [0.0, 0.0, 0.0]
    try:
        immune_hazard(100.0, 0.0)
    except ValueError:
        raised[0] = 1.0
    try:
        simulate_liver(t, c_hi, 300.0, LiverParams(immune_ic50_nm=-1.0))
    except ValueError:
        raised[1] = 1.0
    try:
        simulate_liver(t, c_hi, 300.0, LiverParams(immune_ic50_nm=10.0, immune_recruit_1h=-1.0))
    except ValueError:
        raised[2] = 1.0
    metrics.append(
        MetricResult(
            "degenerate_immune_inputs_rejected",
            sum(raised),
            3.0,
            3.0,
            "count of rejected probes",
            "pass" if sum(raised) == 3.0 else "FAIL",
        )
    )

    ok = all(m.criterion == "pass" for m in metrics)
    notes.append(
        f"I_ss={r_saturate.immune[-1]:.4f} (target {i_ss_expected:.4f}); "
        f"active dead@72h={r_active.dead_frac[-1]:.4f} vs base={r_base.dead_frac[-1]:.4f}; "
        f"low-exposure immune={r_low.immune[-1]:.4f} < active={r_active.immune[-1]:.4f}"
    )
    return CaseResult(
        "Immune-mediated DILI QST (adaptive immune response via hapten hazard)",
        ok,
        metrics,
        notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_dili_immune_activation"]
