"""Mass budget: no-elimination bolus conserves mass (L1, doc/08 §1.4).

Weakest evidence level: verifies only that the ODE solver and mass
bookkeeping neither create nor destroy drug (self-consistency).
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


def case_mass_balance() -> CaseResult:
    b = next(x for x in BENCHMARKS if x.name == "midazolam")
    model = _build_model(b)
    model.cl_hep_l_h = 0.0
    model.cl_renal_l_h = 0.0
    model.dose_plan = DosePlan.iv_bolus(dose_mg=b.dose_mg)
    res = simulate_pbpk(model, tmax_h=b.tmax_h, n_eval=b.n_eval)
    if res.final_state is None:
        rel = 1.0
    else:
        body = model.state_total_mass(res.final_state)
        rel = abs(body - b.dose_mg) / b.dose_mg
    notes = [
        f"no-elimination IV bolus; end-state body mass {body:.6g} mg vs "
        f"dose {b.dose_mg:.6g} mg (rel. err {rel:.2e})"
    ]
    met = MetricResult(
        name="mass_budget_rel_error",
        predicted=rel,
        lo=0.0,
        hi=1e-6,
        unit="fraction",
        criterion="pass" if rel <= 1e-6 else "FAIL",
    )
    return CaseResult(
        "mass-balance (midazolam IV, no elimination)",
        rel <= 1e-6,
        [met],
        notes,
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_mass_balance"]
