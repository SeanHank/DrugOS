"""Bioavailability F reporting (L1 — internal consistency, doc/08 §1.3).

doc/05 §2.1: the pipeline must report a route-dependent systemic
bioavailability.  IV is 100%; a depot route reports its depot availability;
oral dosing loses the unabsorbed fraction (colon-feces sink) and then suffers
well-stirred first-pass hepatic extraction ``F = Fa * (1 - Eh)`` with
``Eh = CL_h/(Q_h + CL_h)``.  In a single-pool (kp = 1) model the systemic
AUC ratio between equal oral and IV doses measures F directly, so the
reported value and the simulated exposure must agree to floating point —
the bookkeeping/ODE self-consistency that L1 certifies.
"""

from __future__ import annotations

from types import SimpleNamespace

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pk.pbpk_build import TISSUE_LIST, AbsorptionParams, PBPKModel, absorption_rate_from_fa
from drugos.pk.simulate import simulate_pbpk

DOSE = 10.0
CL = 0.5  # L/h


def _model(dose_plan: object, **abs_kwargs: object) -> PBPKModel:
    h = profile()
    names: list[str] = [*TISSUE_LIST]
    part = SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})
    return PBPKModel(
        physiology=h,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=CL,
        cl_renal_l_h=0.0,
        dose_plan=dose_plan,
        absorption=AbsorptionParams(**abs_kwargs),
    )


def case_bioavailability_f() -> CaseResult:
    iv = simulate_pbpk(_model(build_dose_plan("iv_bolus", DOSE)), tmax_h=500.0, n_eval=800)
    oral = simulate_pbpk(_model(build_dose_plan("oral", DOSE)), tmax_h=500.0, n_eval=800)
    depot = simulate_pbpk(
        _model(build_dose_plan("subcutaneous", DOSE), depot_bioavailability=0.7),
        tmax_h=500.0,
        n_eval=800,
    )

    f_oral = oral.bioavailability_f
    f_depot = depot.bioavailability_f
    met_iv = MetricResult(
        "iv_bioavailability",
        iv.bioavailability_f or 0.0,
        0.9999,
        1.0001,
        "fraction",
        "pass" if iv.bioavailability_f == 1.0 else "FAIL",
    )
    met_oral_range = MetricResult(
        "oral_bioavailability_range",
        f_oral or 0.0,
        0.0,
        1.0,
        "fraction",
        "pass" if f_oral is not None and 0.0 < f_oral <= 1.0 else "FAIL",
    )
    met_first_pass = MetricResult(
        "first_pass_splits_iv",
        (iv.bioavailability_f or 1.0) - (f_oral or 0.0),
        1e-6,
        1.0,
        "fraction",
        "pass" if (iv.bioavailability_f or 1.0) > (f_oral or 0.0) else "FAIL",
    )
    met_depot = MetricResult(
        "depot_bioavailability",
        f_depot or 0.0,
        0.6999,
        0.7001,
        "fraction",
        "pass" if f_depot is not None and abs(f_depot - 0.7) <= 1e-3 else "FAIL",
    )
    ratio = oral.pk_metrics().auc_inf_mgh_l / max(iv.pk_metrics().auc_inf_mgh_l, 1e-12)
    met_auc_ratio = MetricResult(
        "auc_ratio_matches_reported_F",
        ratio,
        (f_oral or 0.0) * 0.98,
        (f_oral or 0.0) * 1.02,
        "fraction",
        "pass" if f_oral is not None and abs(ratio - f_oral) / f_oral <= 0.02 else "FAIL",
    )
    met_metrics_f = MetricResult(
        "metrics_f_abs_reported",
        oral.pk_metrics().f_abs,
        (f_oral or 0.0) - 1e-9,
        (f_oral or 0.0) + 1e-9,
        "fraction",
        "pass" if f_oral is not None and oral.pk_metrics().f_abs == f_oral else "FAIL",
    )

    # Permeability-gated absorption (doc/05 §1.4): the small-intestine rate is
    # tuned to the fraction-absorbed target by absorption_rate_from_fa, so a
    # low-fa compound peaks lower, later and loses more to the colon-feces
    # sink than an equimolar high-fa one (same F after first-pass is not
    # asserted — only the absorption-axis ordering the pipeline now wires).
    fast = simulate_pbpk(
        _model(
            build_dose_plan("oral", DOSE),
            k_si_absorption=absorption_rate_from_fa(0.95),
        ),
        tmax_h=500.0,
        n_eval=800,
    )
    slow = simulate_pbpk(
        _model(
            build_dose_plan("oral", DOSE),
            k_si_absorption=absorption_rate_from_fa(0.3),
        ),
        tmax_h=500.0,
        n_eval=800,
    )
    met_gated_cmax = MetricResult(
        "permeability_gated_cmax_ordering",
        float(fast.plasma_total.max()) - float(slow.plasma_total.max()),
        1e-6,
        1e6,
        "mg/L",
        "pass" if float(fast.plasma_total.max()) > float(slow.plasma_total.max()) else "FAIL",
    )
    met_gated_feces = MetricResult(
        "permeability_gated_feces_ordering",
        float(slow.feces_cum_mg[-1]) - float(fast.feces_cum_mg[-1]),
        -1e6,
        -1e-6,
        "mg",
        "pass" if float(slow.feces_cum_mg[-1]) > float(fast.feces_cum_mg[-1]) else "FAIL",
    )

    metrics = (
        met_iv,
        met_oral_range,
        met_first_pass,
        met_depot,
        met_auc_ratio,
        met_metrics_f,
        met_gated_cmax,
        met_gated_feces,
    )
    ok = sum(1 for m in metrics if m.criterion == "pass") == len(metrics)
    return CaseResult(
        "bioavailability F reporting (IV/depot/oral first-pass)",
        ok,
        list(metrics),
        [
            f"I F={iv.bioavailability_f:.3f}, depot F={f_depot:.3f}, "
            f"oral F={f_oral:.4f}; oral/IV AUC ratio={ratio:.4f} while CL={CL} L/h. "
            "The unabsorbed colon-transit fraction leaves the oral body and "
            "first-pass hepatic extraction is inside the reported F. "
            "Permeability-gated absorption: fa 0.95 peaks above fa 0.3 with "
            "the low-fa molecules losing more to the colon/feces sink."
        ],
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_bioavailability_f"]
