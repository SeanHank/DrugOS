"""SC/IM depot analytic limit (L2 — analytic, doc/08 §1.2).

A subcutaneous dose enters a first-order-available depot: with kp = 1
everywhere the whole body is a single pool of volume V, so the blood
concentration follows Bateman:

    C(t) = D*F*ka / (V*(ka - Ke)) * (e^{-Ke*t} - e^{-ka*t}),
    Ke = CL/V,  tmax = ln(ka/Ke) / (ka - Ke).

The exercise pins Ka, F and CL: the simulated Cmax, Tmax, AUC and the
unabsorbed fraction (F < 1 must reappear in feces) must match the closed
form.  This is the depot-equivalent of the IV one-compartment limit.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pk.pbpk_build import TISSUE_LIST, AbsorptionParams, PBPKModel
from drugos.pk.simulate import simulate_pbpk


def case_sc_im_depot() -> CaseResult:
    cl = 0.5  # L/h
    ka = 0.3  # 1/h
    f = 0.9
    dose = 10.0
    h = profile()
    v = sum(h.organ_volume.values()) + h.arterial_blood_l + h.venous_blood_l
    ke = cl / v
    names: list[str] = [*TISSUE_LIST]
    part = SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})
    model = PBPKModel(
        physiology=h,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=cl,
        cl_renal_l_h=0.0,
        dose_plan=build_dose_plan("subcutaneous", dose),
        absorption=AbsorptionParams(k_depot_absorption=ka, depot_bioavailability=f),
    )
    res = simulate_pbpk(model, tmax_h=500.0, n_eval=1000)
    m = res.pk_metrics()

    exp_auc = dose * f / cl
    tmax = math.log(ka / ke) / (ka - ke)
    cmax_expr = ka * (math.exp(-ke * tmax) - math.exp(-ka * tmax)) / (ka - ke)
    exp_cmax = dose * f * cmax_expr / v
    exp_feces = dose * (1.0 - f)

    met_cmax = MetricResult(
        "cmax_mg_l",
        m.cmax_mg_l,
        exp_cmax * 0.95,
        exp_cmax * 1.05,
        "mg/L",
        "pass" if 0.95 * exp_cmax <= m.cmax_mg_l <= 1.05 * exp_cmax else "FAIL",
    )
    met_tmax = MetricResult(
        "tmax_h",
        m.tmax_h,
        tmax * 0.9,
        tmax * 1.1,
        "h",
        "pass" if 0.9 * tmax <= m.tmax_h <= 1.1 * tmax else "FAIL",
    )
    met_auc = MetricResult(
        "auc_inf_mgh_l",
        m.auc_inf_mgh_l,
        exp_auc * 0.95,
        exp_auc * 1.05,
        "mg.h/L",
        "pass" if 0.95 * exp_auc <= m.auc_inf_mgh_l <= 1.05 * exp_auc else "FAIL",
    )
    feces = float(res.feces_cum_mg[-1]) if res.feces_cum_mg is not None else 0.0
    met_feces = MetricResult(
        "unabsorbed_feces_mg",
        feces,
        exp_feces * 0.9,
        exp_feces * 1.1,
        "mg",
        "pass" if 0.9 * exp_feces <= feces <= 1.1 * exp_feces else "FAIL",
    )
    ok = (
        met_cmax.criterion == "pass"
        and met_tmax.criterion == "pass"
        and met_auc.criterion == "pass"
        and met_feces.criterion == "pass"
    )
    return CaseResult(
        "SC/IM depot analytic (Bateman single pool)",
        ok,
        [met_cmax, met_tmax, met_auc, met_feces],
        [
            f"V={v:.1f} L, ka={ka}/h, F={f}; analytic Cmax={exp_cmax:.3f} mg/L, "
            f"Tmax={tmax:.1f} h, AUC={exp_auc:.1f} mg.h/L"
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_sc_im_depot"]
