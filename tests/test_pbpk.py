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
