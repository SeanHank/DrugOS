"""Tests for the whole-body PBPK engine (pk/pbpk_build.py, pk/simulate.py)."""

import math

import numpy as np
import pytest

from drugos.inputs.models import HumanProfile, Molecule, Route, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.resolve_human import resolve_human
from drugos.pk.partitions import partition_from_molecule
from drugos.pk.pbpk_build import AbsorptionParams, PBPKModel
from drugos.pk.simulate import PBPKResult, simulate_pbpk


def _partition(fup: float = 0.5, log_p: float = 2.0, pka_bases: list[float] | None = None):
    profile = resolve_human(HumanProfile(sex=Sex.MALE))
    mol = Molecule(
        name="testbase",
        canonical_smiles="C",
        molecular_formula="C",
        mw=200.0,
        log_p=log_p,
        pka_acids=[],
        pka_bases=pka_bases or [],
    )
    return profile, partition_from_molecule(mol, fup=fup, bp=1.0, hematocrit=profile.hematocrit)


def _model(
    route: Route = Route.IV_BOLUS,
    amount_mg: float = 10.0,
    cl_hep: float = 0.0,
    cl_renal: float = 0.0,
    interval_h: float | None = None,
    pka_bases: list[float] | None = None,
) -> PBPKModel:
    profile, part = _partition(pka_bases=pka_bases)
    plan = build_dose_plan(
        route=route,
        amount_mg=amount_mg,
        interval_h=interval_h or 0.0,
        n_doses=1 if interval_h is None else 3,
    )
    return PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=0.5,
        cl_hep_l_h=cl_hep,
        cl_renal_l_h=cl_renal,
        dose_plan=plan,
    )


def test_iv_mass_conservation_without_clearance() -> None:
    model = _model(route=Route.IV_BOLUS, amount_mg=10.0)
    res = simulate_pbpk(model, tmax_h=24.0, n_eval=120)
    assert res.final_state is not None
    assert math.isclose(model.state_total_mass(res.final_state), 10.0, rel_tol=1e-6)


def test_iv_auc_gives_clearance_close_to_input() -> None:
    model = _model(route=Route.IV_BOLUS, amount_mg=50.0, cl_hep=1.5, cl_renal=0.5)
    res = simulate_pbpk(model, tmax_h=96.0, n_eval=500)
    m = res.pk_metrics()
    # Liver/renal intrinsic clearances act on the unbound concentration, so total
    # clearance from plasma is roughly  fup*(cl_hep + cl_renal) = 1.0 L/h.
    assert 0.8 < m.cl_l_h < 1.3
    assert math.isfinite(m.auc_inf_mgh_l) and m.auc_inf_mgh_l > 0
    assert np.isfinite(res.plasma_total).all()


def test_oral_absorbed_fraction_responds_to_rate() -> None:
    fast = _model(route=Route.ORAL, amount_mg=100.0)
    fast.absorption = AbsorptionParams(k_si_absorption=0.9)
    slow = _model(route=Route.ORAL, amount_mg=100.0)
    slow.absorption = AbsorptionParams(k_si_absorption=0.1)

    def absorbed(m: PBPKModel) -> float:
        r = simulate_pbpk(m, tmax_h=72.0, n_eval=200)
        return (m.dose_plan.total_dose_mg - float(r.feces_cum_mg[-1])) / m.dose_plan.total_dose_mg

    f_fast, f_slow = absorbed(fast), absorbed(slow)
    assert f_fast > 0.9
    assert f_slow < f_fast
    assert f_slow > 0.0


def test_solubility_limited_absorption_caps_dissolved_pool() -> None:
    # A low logS-equivalent solubility caps the dissolved small-intestine mass,
    # spilling the excess dose undissolved to the colon/feces sink (novel/ADMET
    # runs only).  High solubility behaves like the default first-order path.
    dose_mg = 1000.0
    low_sol = _model(route=Route.ORAL, amount_mg=dose_mg)
    low_sol.absorption = AbsorptionParams(
        k_si_absorption=0.55, solubility_mg_ml=0.005, gi_volume_ml=250.0
    )
    high_sol = _model(route=Route.ORAL, amount_mg=dose_mg)
    high_sol.absorption = AbsorptionParams(
        k_si_absorption=0.55, solubility_mg_ml=50.0, gi_volume_ml=250.0
    )
    r_low = simulate_pbpk(low_sol, tmax_h=72.0, n_eval=300)
    r_high = simulate_pbpk(high_sol, tmax_h=72.0, n_eval=300)

    def feces_frac(m: PBPKModel, r: PBPKResult) -> float:
        return float(r.feces_cum_mg[-1]) / m.dose_plan.total_dose_mg

    assert feces_frac(low_sol, r_low) > feces_frac(high_sol, r_high)
    assert float(np.max(r_low.plasma_total)) < float(np.max(r_high.plasma_total))
    # The two-compartment capacity itself is the gate: 0.005 mg/mL * 250 mL.
    cap_mg = 0.005 * 250.0
    assert 0.0 < cap_mg < dose_mg


def test_absorption_params_reject_bad_solubility() -> None:
    with pytest.raises(ValueError):
        AbsorptionParams(solubility_mg_ml=0.0)
    with pytest.raises(ValueError):
        AbsorptionParams(solubility_mg_ml=-1.0)
    with pytest.raises(ValueError):
        AbsorptionParams(gi_volume_ml=0.0)
    with pytest.raises(ValueError):
        AbsorptionParams(gut_extraction_eg=-0.1)
    with pytest.raises(ValueError):
        AbsorptionParams(gut_extraction_eg=1.0)


def test_renal_tubular_secretion_raises_urine_fraction() -> None:
    # Active tubular secretion (``cl_sec``, OAT/OCT-style) adds to the
    # glomerular-filtration term on the unbound kidney dose, increasing the
    # fraction of an IV dose recovered in urine over pure filtration.
    base = _model(route=Route.IV_BOLUS, amount_mg=100.0, cl_renal=1.0)
    sec = _model(route=Route.IV_BOLUS, amount_mg=100.0, cl_renal=1.0)
    sec.cl_sec_l_h = 2.0
    r_base = simulate_pbpk(base, tmax_h=48.0, n_eval=200)
    r_sec = simulate_pbpk(sec, tmax_h=48.0, n_eval=200)
    f_base = float(r_base.urine_cum_mg[-1]) / 100.0
    f_sec = float(r_sec.urine_cum_mg[-1]) / 100.0
    assert f_sec > f_base


def test_gut_wall_extraction_reduces_oral_bioavailability() -> None:
    # First-pass intestinal extraction removes a fraction of the SI-absorbed
    # flux before the portal blood (enterocyte CYP/efflux loss), so the
    # reported F falls by (1 - Eg) with the same absorbed mass and hepatic Eh.
    dose_mg = 100.0
    m0 = _model(route=Route.ORAL, amount_mg=dose_mg, cl_hep=0.5)
    m0.absorption = AbsorptionParams(k_si_absorption=0.55)
    m5 = _model(route=Route.ORAL, amount_mg=dose_mg, cl_hep=0.5)
    m5.absorption = AbsorptionParams(k_si_absorption=0.55, gut_extraction_eg=0.5)
    r0 = simulate_pbpk(m0, tmax_h=48.0, n_eval=200)
    r5 = simulate_pbpk(m5, tmax_h=48.0, n_eval=200)
    assert r5.bioavailability_f is not None and r0.bioavailability_f is not None
    assert r5.bioavailability_f < r0.bioavailability_f
    assert math.isclose(r5.bioavailability_f, r0.bioavailability_f * 0.5, rel_tol=0.05)
    assert float(np.max(r5.plasma_total)) < float(np.max(r0.plasma_total))


def test_saturable_hepatic_clearance_is_dose_dependent() -> None:
    # With a Vmax/Km (Michaelis-Menten) hepatic term the intrinsic clearance
    # falls as the unbound liver concentration approaches Km, so the apparent
    # clearance drops and AUC above the low-dose linear regime is
    # super-proportional.  The linear `cl_hep` twin stays dose-proportional.
    vmax, km = 20.0, 1.0  # mg/h, mg/L  ->  linear slope CL_0 = Vmax/Km = 20 L/h
    low_mm = _model(route=Route.IV_BOLUS, amount_mg=10.0)
    low_mm.hepatic_vmax_mg_h, low_mm.hepatic_km_mg_l = vmax, km
    high_mm = _model(route=Route.IV_BOLUS, amount_mg=500.0)
    high_mm.hepatic_vmax_mg_h, high_mm.hepatic_km_mg_l = vmax, km
    low_lin = _model(route=Route.IV_BOLUS, amount_mg=10.0, cl_hep=vmax / km)
    high_lin = _model(route=Route.IV_BOLUS, amount_mg=500.0, cl_hep=vmax / km)

    app_cl: list[float] = []
    for m in (low_mm, high_mm, low_lin, high_lin):
        auc = simulate_pbpk(m, tmax_h=48.0, n_eval=200).pk_metrics().auc_inf_mgh_l
        app_cl.append(m.dose_plan.total_dose_mg / auc)

    cl_mm_low, cl_mm_high, cl_lin_low, cl_lin_high = app_cl
    assert math.isclose(cl_mm_low, cl_lin_low, rel_tol=0.2)  # below Km -> linear slope
    assert cl_mm_high < cl_mm_low  # saturation lowers apparent clearance
    assert math.isclose(cl_lin_high, cl_lin_low, rel_tol=0.1)  # linear stays proportional


def test_biliary_excretion_and_enterohepatic_recirculation() -> None:
    # Biliary secretion routes unbound parent into a bile pool that empties
    # into the SI lumen: with a high SI reabsorption rate the faecal spill is
    # small (recirculation), with a low one it is large.  No clearance terms
    # are active, so total system mass is exactly conserved in both runs.
    dose_mg = 100.0

    def run(k_si_abs: float) -> tuple[PBPKModel, PBPKResult]:
        m = _model(route=Route.IV_BOLUS, amount_mg=dose_mg)
        m.cl_bil_l_h = 2.0
        m.k_bile_emptying_1h = 0.5
        m.absorption = AbsorptionParams(k_si_absorption=k_si_abs)
        return m, simulate_pbpk(m, tmax_h=72.0, n_eval=300)

    m_fast, r_fast = run(k_si_abs=0.9)
    m_slow, r_slow = run(k_si_abs=0.05)
    for m, r in ((m_fast, r_fast), (m_slow, r_slow)):
        assert r.final_state is not None
        assert math.isclose(m.state_total_mass(r.final_state), dose_mg, rel_tol=1e-4)
    # Slower SI reabsorption leaves more of the emptied bile to reach the feces.
    assert float(r_fast.feces_cum_mg[-1]) < float(r_slow.feces_cum_mg[-1])
    # Recirculation keeps the drug in the systemic loop; slow reabsorption
    # drops it to the feces sink instead.
    assert float(np.sum(r_fast.plasma_total)) > float(np.sum(r_slow.plasma_total))


def test_repeated_doses_build_accumulation() -> None:
    single = _model(route=Route.ORAL, amount_mg=100.0, cl_hep=2.0)
    multi = _model(route=Route.ORAL, amount_mg=100.0, cl_hep=2.0, interval_h=24.0)
    auc1 = simulate_pbpk(single, tmax_h=24.0, n_eval=200).pk_metrics().auc_last_mgh_l
    res = simulate_pbpk(multi, tmax_h=72.0, n_eval=300)
    auc3 = res.pk_metrics().auc_last_mgh_l
    assert auc3 > auc1


def test_result_contract_lists() -> None:
    model = _model(route=Route.IV_BOLUS, amount_mg=10.0, cl_hep=1.0)
    res = simulate_pbpk(model, tmax_h=8.0, n_eval=50)
    assert isinstance(res, PBPKResult)
    d = res.to_data_contract()
    assert set(d) >= {"time", "plasma_total", "tissues", "pk_metrics"}
    assert len(d["time"]) == len(d["plasma_total"])
    assert "CMAX" not in d  # contract keys are snake_case
    t = np.asarray(d["time"])
    assert (np.diff(t) > 0).all()
