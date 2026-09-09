"""Extended Stage-1 mechanisms: tubular secretion, gut-wall extraction,
saturable hepatic clearance, and enterohepatic recirculation (L2, doc/08).

Each pair check pins one of the off-by-default realism terms in
``pk/pbpk_build.py`` in a kp=1 single-pool setting, which turns every
compartment into one well-mixed pool of volume ``V`` and gives closed-form
references:

- renal tubular secretion adds a first-order term to the filtration flux, so
  the single-pool urine fraction in [0, T] is ``1 - exp(-(CLr+CLs)/V*T)``;
- first-pass gut-wall extraction multiplies the absorbed fraction ``Fa`` by
  ``(1 - Eg)`` before the hepatic pass, so ``F = (1-Fa_loss)*(1-Eg)*(1-Eh)``;
- the Michaelis-Menten hepatic term reproduces the linear slope
  ``CL_0 = Vmax/Km`` below Km and bends apparent clearance down once the
  unbound liver concentration approaches Km (super-proportional AUC at dose);
- biliary secretion ``cl_bil`` moves parent through a bile pool that empties
  into the SI lumen: with high SI reabsorption the faecal spill is small
  (recirculation) and system mass is conserved exactly.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pk.pbpk_build import TISSUE_LIST, AbsorptionParams, PBPKModel
from drugos.pk.simulate import simulate_pbpk


def _one_pool_model(
    *,
    route: str,
    dose_mg: float,
    cl_hep: float = 0.0,
    cl_renal: float = 0.0,
    cl_sec: float = 0.0,
    vmax: float | None = None,
    km: float | None = None,
    cl_bil: float = 0.0,
    k_bile: float = 0.5,
    k_si_abs: float = 0.55,
    eg: float = 0.0,
) -> tuple[PBPKModel, float]:
    h = profile()
    v = sum(h.organ_volume.values()) + h.arterial_blood_l + h.venous_blood_l
    names: list[str] = [*TISSUE_LIST]
    part = SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})
    model = PBPKModel(
        physiology=h,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=cl_hep,
        cl_renal_l_h=cl_renal,
        cl_sec_l_h=cl_sec,
        hepatic_vmax_mg_h=vmax,
        hepatic_km_mg_l=km,
        cl_bil_l_h=cl_bil,
        k_bile_emptying_1h=k_bile,
        dose_plan=build_dose_plan(route, dose_mg),
        absorption=AbsorptionParams(k_si_absorption=k_si_abs, gut_extraction_eg=eg),
    )
    return model, v


def _fraction(a: float, b: float) -> float:
    return a / b if b else 0.0


def case_clearance_mechanisms() -> CaseResult:
    metrics: list[MetricResult] = []

    # 1. Renal tubular secretion: single-pool urine fraction in [0, 6 h].
    #    Analytic: 1 - exp(-(CLr + CLs)/V * 6).
    dose = 100.0
    m_sec, v = _one_pool_model(route="iv_bolus", dose_mg=dose, cl_renal=0.4, cl_sec=0.6)
    r_sec = simulate_pbpk(m_sec, tmax_h=6.0, n_eval=300)
    urine_frac = _fraction(float(r_sec.urine_cum_mg[-1]), dose)
    exp_frac = 1.0 - math.exp(-(0.4 + 0.6) / v * 6.0)
    metrics.append(
        MetricResult(
            "secretion_urine_fraction",
            urine_frac,
            exp_frac * 0.9,
            exp_frac * 1.1,
            "fraction",
            "pass" if 0.9 * exp_frac <= urine_frac <= 1.1 * exp_frac else "FAIL",
        )
    )
    m_filt, _ = _one_pool_model(route="iv_bolus", dose_mg=dose, cl_renal=0.4)
    r_filt = simulate_pbpk(m_filt, tmax_h=6.0, n_eval=300)
    frac_filt = _fraction(float(r_filt.urine_cum_mg[-1]), dose)
    metrics.append(
        MetricResult(
            "secretion_lifts_urine",
            urine_frac - frac_filt,
            1e-6,
            1.0,
            "fraction",
            "pass" if urine_frac > frac_filt else "FAIL",
        )
    )

    # 2. Gut-wall extraction: oral F drops by (1-Eg) at fixed Eh and landfill.
    odose = 100.0
    m0, _ = _one_pool_model(route="oral", dose_mg=odose, cl_hep=0.5, k_si_abs=0.55, eg=0.0)
    mg, _ = _one_pool_model(route="oral", dose_mg=odose, cl_hep=0.5, k_si_abs=0.55, eg=0.5)
    r0 = simulate_pbpk(m0, tmax_h=72.0, n_eval=600)
    rg = simulate_pbpk(mg, tmax_h=72.0, n_eval=600)
    f0, fg = r0.bioavailability_f, rg.bioavailability_f
    metrics.append(
        MetricResult(
            "gut_wall_reduces_f",
            fg / f0,
            0.45,
            0.55,
            "ratio",
            "pass" if f0 is not None and fg is not None and 0.45 <= fg / f0 <= 0.55 else "FAIL",
        )
    )

    # 3. Saturable (MM) hepatic clearance: below Km the slope is Vmax/Km;
    #    high-dose apparent clearance collapses (super-proportional AUC).
    low, v = _one_pool_model(route="iv_bolus", dose_mg=5.0, vmax=20.0, km=1.0)
    high, _ = _one_pool_model(route="iv_bolus", dose_mg=500.0, vmax=20.0, km=1.0)
    r_low = simulate_pbpk(low, tmax_h=96.0, n_eval=600)
    r_high = simulate_pbpk(high, tmax_h=96.0, n_eval=600)
    cl_low = 5.0 / r_low.pk_metrics().auc_inf_mgh_l
    cl_high = 500.0 / r_high.pk_metrics().auc_inf_mgh_l
    metrics.append(
        MetricResult(
            "mm_low_dose_slope",
            cl_low,
            16.0,
            24.0,
            "L/h",
            "pass" if 16.0 <= cl_low <= 24.0 else "FAIL",
        )
    )
    metrics.append(
        MetricResult(
            "mm_saturation_drops_cl",
            cl_high,
            0.0,
            cl_low * 0.5,
            "L/h",
            "pass" if cl_high < cl_low * 0.5 else "FAIL",
        )
    )

    # 4. Biliary excretion + enterohepatic recirculation: mass conserved, slow
    #    SI reabsorption spills bile to feces, fast reabsorption recirculates.
    bdose = 100.0

    def bile_run(k_si_abs: float) -> tuple[PBPKModel, object]:
        m, _ = _one_pool_model(
            route="iv_bolus",
            dose_mg=bdose,
            cl_bil=2.0,
            k_bile=0.5,
            k_si_abs=k_si_abs,
        )
        return m, simulate_pbpk(m, tmax_h=96.0, n_eval=600)

    m_fast, r_fast = bile_run(0.9)
    m_slow, r_slow = bile_run(0.05)
    conserved = all(
        math.isclose(m.state_total_mass(r.final_state), bdose, rel_tol=1e-3)
        for m, r in ((m_fast, r_fast), (m_slow, r_slow))
        if r.final_state is not None
    )
    feces_fast = float(r_fast.feces_cum_mg[-1])
    feces_slow = float(r_slow.feces_cum_mg[-1])
    metrics.append(
        MetricResult(
            "ehc_mass_conservative",
            1.0 if conserved else 0.0,
            1.0,
            1.0,
            "flag",
            "pass" if conserved else "FAIL",
        )
    )
    metrics.append(
        MetricResult(
            "ehc_reabsorption_vs_feces",
            feces_fast,
            0.0,
            feces_slow,
            "mg",
            "pass" if feces_fast < feces_slow else "FAIL",
        )
    )
    metrics.append(
        MetricResult(
            "ehc_lifts_plasma",
            float(r_fast.plasma_total.sum()) - float(r_slow.plasma_total.sum()),
            1e-6,
            1e6,
            "mg.h/L",
            "pass"
            if float(r_fast.plasma_total.sum()) > float(r_slow.plasma_total.sum())
            else "FAIL",
        )
    )

    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "Clearance & absorption realism (secretion / gut-wall / MM / EHC)",
        ok,
        metrics,
        [
            f"V={v:.1f} L single pool; secretion urine_frac={urine_frac:.3f} "
            f"vs {exp_frac:.3f} analytic; gut-wall F ratio={_fraction(fg or 0, f0 or 0):.3f}; "
            f"MM CL low={cl_low:.2f}/high={cl_high:.2f} L/h; "
            f"EHC feces fast={feces_fast:.2f}/slow={feces_slow:.2f} mg"
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_clearance_mechanisms"]
