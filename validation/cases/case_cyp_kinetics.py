"""Per-isoform CYP saturable hepatic kinetics (L2, doc/08 doc/05 1.4/4.1).

The per-CYP Michaelis-Menten/Hill hepatic term (``PBPKModel.cyp_terms``,
off by default) replaces the lumped linear ``cl_hep`` with the sum of
isoform fluxes on the *unbound* liver concentration.  The abundance-scaled
Vmax is built from the physiology CYP-abundance table
(``HumanPhysiology.hepatic_cyp_content_nmol``, Barter et al. 2013) via
``cyp_vmax_mg_h``, so the resource the stage was deferred on is now the
actual input.

Checks pin, in a kp=1 single pool:

- the resource pipe: Vmax comes from the abundance table (helper units);
- two-isoform additivity: below Km the low-dose plasma clearance equals the
  lumped twin ``CL_0 = sum(Vmax_i/Km_i)``;
- abundance scaling: doubling per-isoform liver content doubles the low-dose
  slope;
- Hill (n>1) sigmoid: at sub-Km unbound liver concentrations the cooperative
  form carries *less* flux than Michaelis-Menten (sharper saturation region)
  and above Km it approaches Vmax faster;
- parameter validation: a degenerate term raises instead of silently corrupting
  the ODE.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pk.pbpk_build import (
    TISSUE_LIST,
    AbsorptionParams,
    CypTerm,
    PBPKModel,
    cyp_vmax_mg_h,
)
from drugos.pk.simulate import simulate_pbpk


def _one_pool_model(cyp_terms: tuple[CypTerm, ...], dose_mg: float) -> tuple[PBPKModel, float]:
    h = profile()
    v = sum(h.organ_volume.values()) + h.arterial_blood_l + h.venous_blood_l
    names: list[str] = [*TISSUE_LIST]
    part = SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})
    model = PBPKModel(
        physiology=h,
        partition=part,
        bp=1.0,
        fup=1.0,
        cyp_terms=cyp_terms,
        dose_plan=build_dose_plan("iv_bolus", dose_mg),
        absorption=AbsorptionParams(),
    )
    return model, v


def _app_cl(cyp_terms: tuple[CypTerm, ...], dose_mg: float) -> float:
    m, _ = _one_pool_model(cyp_terms, dose_mg)
    res = simulate_pbpk(m, tmax_h=48.0, n_eval=300)
    return dose_mg / float(res.pk_metrics().auc_inf_mgh_l)


def case_cyp_kinetics() -> CaseResult:
    metrics: list[MetricResult] = []
    h = profile()

    # 1. Resource pipe: Vmax = kcat * E_total * MW from the abundance table.
    content = h.hepatic_cyp_content_nmol["CYP3A4"]
    kcat, mw = 2.0, 500.0
    vmax = cyp_vmax_mg_h(content, kcat, mw)
    metrics.append(
        MetricResult(
            "abundance_scaled_vmax_from_table",
            vmax,
            0.999 * content * 1.0e-6 * kcat * mw,
            1.001 * content * 1.0e-6 * kcat * mw,
            "mg/h",
            "pass" if math.isclose(vmax, content * 1.0e-6 * kcat * mw, rel_tol=1e-3) else "FAIL",
        )
    )

    # 2. Two-isoform additivity: below Km the low-dose clearance is the sum of
    #    Vmax/Km slopes, matching the lumped cl_hep twin.
    terms = (
        CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=12.0),
        CypTerm(isoform="CYP2D6", km_mg_l=1.0, vmax_mg_h=8.0),
    )
    cl0 = sum(t.vmax_mg_h / t.km_mg_l for t in terms)
    dose = 5.0
    cl_sum = _app_cl(terms, dose)
    m_lin, _ = _one_pool_model((), dose)
    m_lin.cl_hep_l_h = cl0
    res = simulate_pbpk(m_lin, tmax_h=48.0, n_eval=300)
    cl_lin = dose / float(res.pk_metrics().auc_inf_mgh_l)
    metrics.append(
        MetricResult(
            "isoform_additivity_low_dose_slope",
            cl_sum,
            cl0 * 0.8,
            cl0 * 1.2,
            "L/h",
            "pass" if 0.8 * cl0 <= cl_sum <= 1.2 * cl0 else "FAIL",
        )
    )
    metrics.append(
        MetricResult(
            "matches_lumped_linear_twin",
            cl_sum,
            cl_lin * 0.8,
            cl_lin * 1.2,
            "L/h",
            "pass" if 0.8 * cl_lin <= cl_sum <= 1.2 * cl_lin else "FAIL",
        )
    )

    # 3. Abundance scaling: doubling per-isoform content doubles the low-dose
    #    (linear) slope, since Vmax is linear in hepatic CYP content.
    content2 = h.hepatic_cyp_content_nmol["CYP2D6"]
    v1 = cyp_vmax_mg_h(content2, kcat, mw)
    v2 = cyp_vmax_mg_h(2.0 * content2, kcat, mw)
    t1 = (CypTerm(isoform="CYP2D6", km_mg_l=1.0, vmax_mg_h=v1),)
    t2 = (CypTerm(isoform="CYP2D6", km_mg_l=1.0, vmax_mg_h=v2),)
    cl_1, cl_2 = _app_cl(t1, dose), _app_cl(t2, dose)
    metrics.append(
        MetricResult(
            "abundance_doubling_doubles_slope",
            cl_2 / cl_1,
            1.8,
            2.2,
            "ratio",
            "pass" if 1.8 <= cl_2 / cl_1 <= 2.2 else "FAIL",
        )
    )

    # 4. Hill sigmoid vs Michaelis-Menten (same Vmax/Km): below Km the
    #    cooperative term saturates harder (less flux), above Km it approaches
    #    Vmax faster (more flux) — the sharper transition region.
    km = 1.0
    mm = CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=1.0)
    hill = CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=1.0, hill=2.0)
    metric_below = hill.rate(0.5 * km) / mm.rate(0.5 * km)
    metric_above = hill.rate(2.0 * km) / mm.rate(2.0 * km)
    metrics.append(
        MetricResult(
            "hill_saturates_more_below_km",
            metric_below,
            0.0,
            0.85,
            "ratio vs MM",
            "pass" if metric_below < 0.85 and metric_below > 0.0 else "FAIL",
        )
    )
    metrics.append(
        MetricResult(
            "hill_approaches_vmax_faster_above_km",
            metric_above,
            1.05,
            1.5,
            "ratio vs MM",
            "pass" if metric_above > 1.05 else "FAIL",
        )
    )

    # 5. Degenerate parameters raise instead of corrupting the ODE.
    raised = 0.0
    try:
        CypTerm(isoform="CYP3A4", km_mg_l=0.0, vmax_mg_h=1.0)
    except ValueError:
        raised = 1.0
    metrics.append(
        MetricResult(
            "invalid_parameters_rejected",
            raised,
            1.0,
            1.0,
            "flag",
            "pass" if raised == 1.0 else "FAIL",
        )
    )

    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "Per-CYP hepatic kinetics (MM/Hill, abundance-scaled Vmax)",
        ok,
        metrics,
        [
            f"source Vmax(CYP3A4)={vmax:.3f} mg/h from {content:.0f} nmol content; "
            f"low-dose CL={cl_sum:.2f} vs twin {cl_lin:.2f} L/h; "
            f"doubling ratio={cl_2 / cl_1:.2f}; "
            f"Hill/MM flux at 0.5Km={metric_below:.2f}, at 2Km={metric_above:.2f}"
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_cyp_kinetics"]
