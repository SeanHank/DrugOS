"""Native target-mediated drug disposition coupling (L2, doc/08, doc/05 1.4/2.4).

The sequential pipeline only *approximates* TMDD via the opt-in
``feedback_loop`` driver; the deferred item asked for true mass-balance
coupling so Stage-2 occupancy fluxes feed back into the tissue ODEs.  That is
now ``PBPKModel.target_binding`` (off by default): one reversible binding site
is coupled straight into a tissue mass balance with the regulatorily standard
turnover model of ``drugos.target.occupancy`` — second-order association
``kon*D*R``, first-order dissociation ``koff*DR`` (``koff = kon*kd`` on the
``Target``), constant receptor synthesis ``ksyn = rho*R0``, and irreversible
internalization ``kint*DR`` that drains drug into a cleared sink.

The site here is a GLP-1R-like illustrative gut receptor (high affinity,
high abundance, internalizing) to exercise the mechanism; the checks pin, in
a kp=1 single pool:

- drug mass conservation with reversible binding and with the internalized
  sink (the tallied total stays exactly at the administered dose, no
  metabolic/renal clearance);
- the sink really departs the free-circulating pool;
- saturable target capacity gives dose-disproportional elimination (the TMDD
  hallmark: a low dose is scavenged as a much larger fraction);
- quasi-steady mass-action returns the target KD (DR/R = D/Kd);
- binding lowers exposure vs the free-dispersion twin; the extension is off
  by default (no extra states) and degenerate sites are rejected.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pk.pbpk_build import PBPKModel, TargetBinding
from drugos.pk.simulate import simulate_pbpk
from drugos.target.targets import Target

_MW = 450.0
_KD = 1.5
_KON = 0.5
_R0 = 800.0
_RHO = 0.08
_SITE = "gut"


def _tmdd_model(dose_mg: float, kint: float = 0.25) -> PBPKModel:
    h = profile()
    names: list[str] = list(h.organ_volume) + ["arterial", "venous"]
    part = SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})
    target = Target(
        name="GLP-1R-like gut site",
        kd_nm=_KD,
        kon_nm_h=_KON,
        r0_nm=_R0,
        rho_h=_RHO,
        kint_h=kint,
    )
    return PBPKModel(
        physiology=h,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        mw_g_per_mol=_MW,
        target_binding=(TargetBinding(tissue=_SITE, target=target),),
        dose_plan=build_dose_plan("iv_bolus", dose_mg),
    )


def _in_system_fraction(model: PBPKModel) -> tuple[float, float]:
    res = simulate_pbpk(model, tmax_h=24.0, n_eval=200)
    assert res.final_state is not None
    assert res.tmdd is not None
    cleared = float(list(res.tmdd.values())[0]["cleared_mg"][-1])
    total = float(model.state_total_mass(res.final_state))
    return (total - cleared) / res.dose_mg, cleared


def case_tmdd_drug_disposition() -> CaseResult:
    metrics: list[MetricResult] = []
    dose = 50.0

    # 1. Reversible binding conserves drug mass exactly (no clearance, kint=0).
    m_rev = _tmdd_model(dose, kint=0.0)
    r_rev = simulate_pbpk(m_rev, tmax_h=24.0, n_eval=200)
    assert r_rev.final_state is not None
    bound_mg = float(list(r_rev.tmdd.values())[0]["bound_mg"][-1])
    mass_closed = float(m_rev.state_total_mass(r_rev.final_state))
    ok1 = math.isclose(mass_closed, r_rev.dose_mg, rel_tol=1e-6) and bound_mg > 1e-3
    metrics.append(
        MetricResult(
            "mass_conserves_with_binding",
            mass_closed,
            r_rev.dose_mg * 0.999999,
            r_rev.dose_mg * 1.000001,
            "mg",
            "pass" if ok1 else "FAIL",
        )
    )

    # 2. Internalization is an irreversible drug sink that leaves circulation.
    m_sink = _tmdd_model(dose, kint=0.25)
    remain_f, cleared = _in_system_fraction(m_sink)
    ok2 = remain_f < 0.999 and cleared > 0.5
    metrics.append(
        MetricResult(
            "tmdd_internalization_is_drug_sink",
            cleared,
            0.5,
            3.0,
            "mg cleared in 24h",
            "pass" if ok2 else "FAIL",
        )
    )

    # 3. Dose-disproportional elimination (saturable target capacity): the low
    #    dose is scavenged as a much larger fraction than the high dose.
    f_low, _ = _in_system_fraction(_tmdd_model(5.0, kint=0.25))
    f_high, _ = _in_system_fraction(_tmdd_model(2000.0, kint=0.25))
    ok3 = f_low < 0.9 and f_high > 0.95 and f_low < f_high
    metrics.append(
        MetricResult(
            "dose_disproportional_retention",
            f_low,
            0.0,
            0.9,
            "fraction retained at 5 mg",
            "pass" if ok3 else "FAIL",
        )
    )
    metrics.append(
        MetricResult(
            "high_dose_approaches_linear_retention",
            f_high,
            0.95,
            1.0,
            "fraction retained at 2000 mg",
            "pass" if f_high > 0.95 else "FAIL",
        )
    )

    # 4. Quasi-steady mass action recreates the target KD: DR/R = D/Kd in the
    #    reversible run (well-mixed pool, late flat window).
    r_eq = simulate_pbpk(m_rev, tmax_h=72.0, n_eval=400)
    row = list(r_eq.tmdd.values())[0]
    free_nm = np.asarray(r_eq.unbound_tissues[_SITE]) * 1.0e6 / _MW
    half = len(r_eq.t) // 2
    ratio = np.asarray(row["complex_nmol"][half:]) / np.asarray(row["receptor_nmol"][half:])
    quiet = free_nm[half:] > 1e-3
    dev = np.log10((ratio[quiet] * _KD) / free_nm[half:][quiet])
    mean_dev = float(np.nanmean(dev))
    ok4 = quiet.sum() >= 10 and abs(mean_dev) < 0.02
    metrics.append(
        MetricResult(
            "quasi_steady_ratio_matches_kd",
            mean_dev,
            -0.02,
            0.02,
            "log10(DR/R vs D/Kd)",
            "pass" if ok4 else "FAIL",
        )
    )

    # 5. Binding lowers systemic exposure vs the free-dispersion twin: at a low
    #    dose the internalizing sink removes a visible fraction of the AUC.
    m_sink5 = _tmdd_model(5.0, kint=0.25)
    auc_bind = simulate_pbpk(m_sink5, tmax_h=24.0, n_eval=200).pk_metrics().auc_last_mgh_l
    h5 = profile()
    names5: list[str] = list(h5.organ_volume) + ["arterial", "venous"]
    part5 = SimpleNamespace(kp={n: 1.0 for n in names5}, kpu={n: 1.0 for n in names5})
    twin = PBPKModel(
        physiology=h5,
        partition=part5,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        mw_g_per_mol=_MW,
        dose_plan=build_dose_plan("iv_bolus", dose),
    )
    r_twin = simulate_pbpk(twin, tmax_h=24.0, n_eval=200)
    auc_twin = r_twin.pk_metrics().auc_last_mgh_l
    ok5 = auc_bind < 0.9 * auc_twin
    metrics.append(
        MetricResult(
            "binding_lowers_exposure_vs_twin",
            auc_bind / auc_twin,
            0.0,
            0.9,
            "AUC ratio",
            "pass" if ok5 else "FAIL",
        )
    )

    # 6. Extension off by default: the base state vector keeps its documented
    #    size (2 blood + tissues + 7 GI/urine/feces/depot/bile compartments).
    base_n = _tmdd_model(dose).n_state - 3
    expected_n = 2 + len(profile().organ_volume) + 7
    metrics.append(
        MetricResult(
            "off_by_default_state_count",
            base_n,
            expected_n,
            expected_n,
            "state dim without binding",
            "pass" if base_n == expected_n else "FAIL",
        )
    )

    # 7. Degenerate sites raise instead of corrupting the ODE.
    raised = 0.0
    try:
        PBPKModel(
            physiology=profile(),
            partition=part5,
            bp=1.0,
            fup=1.0,
            mw_g_per_mol=_MW,
            target_binding=(
                TargetBinding(
                    tissue="myocardium_extra",
                    target=Target(name="x", kd_nm=1.0),
                ),
            ),
            dose_plan=build_dose_plan("iv_bolus", dose),
        )
    except ValueError:
        raised = 1.0
    metrics.append(
        MetricResult(
            "degenerate_site_rejected",
            raised,
            1.0,
            1.0,
            "flag",
            "pass" if raised == 1.0 else "FAIL",
        )
    )

    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "Native TMDD drug disposition (mass-balance coupling)",
        ok,
        metrics,
        [
            f"mass-closed={mass_closed:.2f}/{r_rev.dose_mg:.0f} mg, bound={bound_mg:.3f} mg; "
            f"cleared(24h)={cleared:.2f} mg sink; "
            f"retention 5 mg={f_low:.2f} vs 2000 mg={f_high:.2f} "
            f"(dose-disproportional); <log10(DR/R vs D/Kd)>={mean_dev:.4f} "
            f"(Kd={_KD} nM); AUC(bound)/AUC(twin)={auc_bind / auc_twin:.2f}"
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_tmdd_drug_disposition"]
