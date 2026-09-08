"""Boundary/edge coverage for the PK engine: partition helpers, absorption,
event application across all routes, and non-compartmental metrics edge cases."""

import math

import numpy as np
import pytest

from drugos.inputs.models import DosePlan, HumanProfile, Molecule, Route, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.resolve_human import resolve_human
from drugos.pk.partitions import (
    ka_ap_strong_base,
    partition_from_molecule,
    rodgers_rowland_partition,
)
from drugos.pk.pbpk_build import AbsorptionParams, PBPKModel, absorption_rate_from_fa
from drugos.pk.simulate import (
    bioavailable_fraction,
    compute_pk_metrics,
    simulate_pbpk,
)


# ---------------------------------------------------------------------------
# partitions: parameter validation and helper branches
# ---------------------------------------------------------------------------
def test_ka_ap_strong_base_validation() -> None:
    with pytest.raises(ValueError, match="fup"):
        ka_ap_strong_base(2.0, 9.0, 0.0, 1.0, 0.45)
    with pytest.raises(ValueError, match="fup"):
        ka_ap_strong_base(2.0, 9.0, 1.5, 1.0, 0.45)
    with pytest.raises(ValueError, match="blood-to-plasma"):
        ka_ap_strong_base(2.0, 9.0, 0.5, -1.0, 0.45)


def test_rodgers_invalid_fup() -> None:
    with pytest.raises(ValueError, match="fup"):
        rodgers_rowland_partition(fup=2.0, log_p=2.0)


def test_partition_protein_ratio_modes() -> None:
    # albumin mode -> comp.ar; lipoprotein mode -> comp.lr.
    al = rodgers_rowland_partition(fup=0.2, log_p=2.0, pka_acids=[3.0], protein_ratio="albumin")
    lp = rodgers_rowland_partition(fup=0.2, log_p=2.0, pka_acids=[3.0], protein_ratio="lipoprotein")
    assert al.kpu["liver"] != lp.kpu["liver"]
    # auto: ionizable compounds use the albumin ratio.
    auto = rodgers_rowland_partition(fup=0.2, log_p=2.0, pka_acids=[3.0])
    assert auto.kpu["liver"] == pytest.approx(al.kpu["liver"])


def test_partition_from_molecule_wrapper() -> None:
    from drugos.inputs.models import RrClass

    profile = resolve_human(HumanProfile(sex=Sex.MALE))
    mol = Molecule(name="x", log_p=2.0, pka_bases=[9.0], rr_class=RrClass.STRONG_BASE)
    part = partition_from_molecule(mol, fup=0.5, bp=1.0, hematocrit=profile.hematocrit)
    assert part.rr_class is RrClass.STRONG_BASE
    assert part.kpu_bc is not None

    untyped = Molecule(name="x", log_p=None)
    with pytest.raises(ValueError, match="log_p"):
        partition_from_molecule(untyped, fup=0.5)


# ---------------------------------------------------------------------------
# absorption rate helper
# ---------------------------------------------------------------------------
def test_absorption_rate_from_fa_clamps() -> None:
    assert absorption_rate_from_fa(0.0) == pytest.approx(0.55 * (0.05 / 0.85))
    assert absorption_rate_from_fa(0.99) == pytest.approx(0.55 * (0.98 / 0.85))
    assert absorption_rate_from_fa(0.85) == pytest.approx(0.55)
    assert absorption_rate_from_fa(0.5, base=1.0) == pytest.approx(1.0 * (0.5 / 0.85))


# ---------------------------------------------------------------------------
# PBPK: event application across routes, infusion dynamics
# ---------------------------------------------------------------------------
def _model(route: Route, amount_mg: float = 10.0, **kwargs) -> PBPKModel:
    profile = resolve_human(HumanProfile(sex=Sex.MALE))
    mol = Molecule(name="t", log_p=2.0, pka_bases=[9.0])
    part = partition_from_molecule(mol, fup=0.5, bp=1.0, hematocrit=profile.hematocrit)
    plan = build_dose_plan(route=route, amount_mg=amount_mg, **kwargs)
    return PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=0.5,
        dose_plan=plan,
    )


@pytest.mark.parametrize(
    "route",
    [Route.ORAL, Route.SUBCUTANEOUS, Route.INTRAMUSCULAR, Route.TRANSDERMAL],
)
def test_apply_event_injects_into_correct_compartment(route: Route) -> None:
    model = _model(route, amount_mg=7.0)
    state = model.initial_state()
    model.apply_event(state, 0.0)
    if route is Route.ORAL:
        assert state[model.state_index("stomach")] == 7.0
    else:
        assert state[model.state_index("depot")] == 7.0
    assert abs(model.state_total_mass(state) - 7.0) < 1e-12


def test_apply_event_ignores_non_matching_time() -> None:
    model = _model(Route.IV_BOLUS, amount_mg=5.0)
    state = model.initial_state()
    model.apply_event(state, 1.0)
    assert abs(model.state_total_mass(state)) < 1e-12


def test_multiple_routes_in_one_plan() -> None:
    plan = DosePlan(events=[])
    plan.add(DosePlan.iv_bolus(1.0, time_h=0.0).events[0])
    plan.add(DosePlan.oral(2.0, time_h=0.0).events[0])
    plan.add(build_dose_plan("subcutaneous", 3.0, start_h=0.0).events[0])
    model = _model(Route.IV_BOLUS)
    model.dose_plan = plan
    state = model.initial_state()
    model.apply_event(state, 0.0)
    assert state[model.state_index("venous")] == pytest.approx(1.0)
    assert state[model.state_index("stomach")] == pytest.approx(2.0)
    assert state[model.state_index("depot")] == pytest.approx(3.0)


def test_infusion_raises_plasma_more_than_bolus_early() -> None:
    infusion = _model(Route.IV_INFUSION, amount_mg=50.0, duration_h=4.0)
    bolus = _model(Route.IV_BOLUS, amount_mg=50.0)
    r_inf = simulate_pbpk(infusion, tmax_h=2.0, n_eval=60)
    r_bolus = simulate_pbpk(bolus, tmax_h=2.0, n_eval=60)
    # During the infusion the plasma level is still building up (lower AUC than
    # the instantaneous bolus) but it must be rising over the interval.
    assert np.isfinite(r_inf.plasma_total).all()
    assert r_inf.plasma_total[-1] < r_bolus.plasma_total[-1]
    assert float(np.diff(r_inf.plasma_total).max()) > 0


def test_sc_depot_flattens_and_delays_peak_but_preserves_auc() -> None:
    bolus = _model(Route.IV_BOLUS, amount_mg=10.0)
    sc = _model(Route.SUBCUTANEOUS, amount_mg=10.0)
    r_bolus = simulate_pbpk(bolus, tmax_h=200.0, n_eval=400)
    r_sc = simulate_pbpk(sc, tmax_h=200.0, n_eval=400)
    m_bolus = r_bolus.pk_metrics()
    m_sc = r_sc.pk_metrics()
    # Depot absorption blunts and delays the peak relative to an instant bolus.
    assert m_sc.cmax_mg_l < m_bolus.cmax_mg_l
    assert m_sc.tmax_h > m_bolus.tmax_h
    # With F=1 the systemic exposure is conserved (same dose). Compare over a
    # horizon long enough that the depot has fully emptied, so the absorption
    # lag no longer truncates the SC profile.
    assert m_sc.auc_last_mgh_l == pytest.approx(m_bolus.auc_last_mgh_l, rel=0.05)
    # The depot is drained by the end of the horizon.
    assert r_sc.final_state is not None
    assert r_sc.final_state[sc.state_index("depot")] < 1e-6


def test_absorption_params_reject_invalid_depot_values() -> None:
    with pytest.raises(ValueError):
        AbsorptionParams(k_depot_absorption=0.0)
    with pytest.raises(ValueError):
        AbsorptionParams(depot_bioavailability=0.0)
    with pytest.raises(ValueError):
        AbsorptionParams(depot_bioavailability=1.5)


def test_empty_dose_plan_simulates_to_zero() -> None:
    profile = resolve_human(HumanProfile(sex=Sex.MALE))
    mol = Molecule(name="t", log_p=2.0)
    part = partition_from_molecule(mol, fup=0.5)
    model = PBPKModel(physiology=profile, partition=part, bp=1.0, fup=0.5, dose_plan=DosePlan())
    res = simulate_pbpk(model, tmax_h=2.0, n_eval=40)
    assert np.allclose(res.plasma_total, 0.0)
    assert res.route is Route.ORAL  # default route taken


def test_n_state_matches_state_vector() -> None:
    model = _model(Route.IV_BOLUS)
    assert model.n_state == len(model.initial_state())


def test_solver_failure_raises(monkeypatch) -> None:
    import types as _types

    fake = _types.SimpleNamespace(success=False, message="forced failure")
    monkeypatch.setattr("drugos.pk.simulate.solve_ivp", lambda *a, **k: fake)
    with pytest.raises(RuntimeError, match="forced failure"):
        simulate_pbpk(_model(Route.IV_BOLUS), tmax_h=1.0)


def test_simulate_rejects_nonpositive_tmax() -> None:
    with pytest.raises(ValueError, match="tmax_h"):
        simulate_pbpk(_model(Route.IV_BOLUS), tmax_h=0.0)


# ---------------------------------------------------------------------------
# bioavailability F (doc/05 2.1): route-dependent reporting, first-pass
# ---------------------------------------------------------------------------
def _model_with_abs(
    route: Route, amount_mg: float, cl_hep_l_h: float, **abs_kwargs: object
) -> PBPKModel:
    profile = resolve_human(HumanProfile(sex=Sex.MALE))
    mol = Molecule(name="t", log_p=2.0, pka_bases=[9.0])
    part = partition_from_molecule(mol, fup=0.5, bp=1.0, hematocrit=profile.hematocrit)
    plan = build_dose_plan(route=route, amount_mg=amount_mg)
    return PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=0.5,
        cl_hep_l_h=cl_hep_l_h,
        cl_renal_l_h=0.2,
        dose_plan=plan,
        absorption=AbsorptionParams(**abs_kwargs),  # type: ignore[arg-type]
    )


def test_bioavailable_fraction_iv_is_unity() -> None:
    m = _model_with_abs(Route.IV_BOLUS, 10.0, 0.0)
    assert bioavailable_fraction(m, Route.IV_BOLUS, 10.0, None) == 1.0
    assert bioavailable_fraction(m, Route.IV_INFUSION, 10.0, np.array([5.0])) == 1.0


def test_bioavailable_fraction_depot_uses_availability() -> None:
    from types import SimpleNamespace

    sc = _model_with_abs(Route.SUBCUTANEOUS, 10.0, 0.0, depot_bioavailability=0.64)
    assert bioavailable_fraction(sc, Route.SUBCUTANEOUS, 10.0, None) == pytest.approx(0.64)
    td = _model_with_abs(Route.TRANSDERMAL, 10.0, 0.0, depot_bioavailability=0.1)
    assert bioavailable_fraction(td, Route.TRANSDERMAL, 10.0, np.array([0.0])) == pytest.approx(0.1)
    # A (theoretically) out-of-range availability is clamped onto [0, 1].
    stub = SimpleNamespace(
        absorption=SimpleNamespace(depot_bioavailability=1.3),
        physiology=SimpleNamespace(organ_flow={"liver": 1.4}),
        cl_hep_l_h=0.0,
    )
    assert bioavailable_fraction(stub, Route.INTRAMUSCULAR, 10.0, None) == pytest.approx(1.0)
    stub_low = SimpleNamespace(
        absorption=SimpleNamespace(depot_bioavailability=-0.5),
        physiology=SimpleNamespace(organ_flow={"liver": 1.4}),
        cl_hep_l_h=0.0,
    )
    got = bioavailable_fraction(stub_low, Route.INTRAMUSCULAR, 10.0, np.array([0.0]))
    assert got == pytest.approx(0.0)


def test_bioavailable_fraction_oral_first_pass() -> None:
    m = _model_with_abs(Route.ORAL, 100.0, cl_hep_l_h=0.5)
    # No feces sink -> Fa = 1 and F = 1 - Eh < 1 because of hepatic extraction.
    f_no_sink = bioavailable_fraction(m, Route.ORAL, 100.0, None)
    f_clean = bioavailable_fraction(m, Route.ORAL, 100.0, np.array([0.0, 0.0]))
    assert 0.0 < f_clean < 1.0
    assert f_no_sink == pytest.approx(f_clean)
    # All of the dose overflowing into the feces sink -> nothing absorbed.
    assert bioavailable_fraction(m, Route.ORAL, 100.0, np.array([0.0, 100.0])) == 0.0
    # Zero dose with a computed (finite) feces sink: absorbed term is skipped.
    assert bioavailable_fraction(m, Route.ORAL, 0.0, np.array([0.0, 1.0])) == pytest.approx(
        bioavailable_fraction(m, Route.ORAL, 0.0, None)
    )


def test_simulate_reports_bioavailability_f_and_passes_to_metrics() -> None:
    m = _model_with_abs(Route.ORAL, 100.0, cl_hep_l_h=0.5)
    r = simulate_pbpk(m, tmax_h=24.0, n_eval=120)
    assert r.bioavailability_f is not None
    assert 0.0 < r.bioavailability_f <= 1.0
    assert r.pk_metrics().f_abs == pytest.approx(r.bioavailability_f)
    iv = simulate_pbpk(_model_with_abs(Route.IV_BOLUS, 10.0, 0.0), tmax_h=24.0, n_eval=120)
    assert iv.bioavailability_f == 1.0


# ---------------------------------------------------------------------------
# compute_pk_metrics edge cases
# ---------------------------------------------------------------------------
def _metrics(conc, dose=100.0, route=Route.IV_BOLUS):
    t = np.linspace(0.0, 24.0, len(conc))
    return compute_pk_metrics(np.asarray(t), np.asarray(conc, dtype=float), dose, route)


def test_pk_metrics_declining_profile() -> None:
    c = np.exp(-0.2 * np.linspace(0, 24, 50)) * 10.0
    m = _metrics(c)
    assert m.tmax_h == 0.0
    assert m.cl_l_h > 0
    assert math.isfinite(m.term_half_life_h)
    assert m.vss_l > 0
    assert math.isclose(m.term_half_life_h, math.log(2) / 0.2, rel_tol=0.1)
    assert m.auc_inf_mgh_l > m.auc_last_mgh_l


def test_pk_metrics_increasing_profile_no_negative_lambda() -> None:
    # A rising tail has a positive slope -> lambda_z <= 0 -> fallback branch.
    c = np.linspace(0.0, 5.0, 20)
    m = _metrics(c)
    assert m.lambda_z_1h <= 0.0
    assert m.term_half_life_h == math.inf
    assert m.auc_inf_mgh_l == pytest.approx(m.auc_last_mgh_l)


def test_pk_metrics_few_points_rejected() -> None:
    with pytest.raises(ValueError, match="3 time points"):
        compute_pk_metrics(np.array([0.0, 1.0]), np.array([1.0, 1.0]), 10.0, Route.ORAL)


def test_pk_metrics_tail_too_short_uses_last_points() -> None:
    # After Cmax the floor (5% of Cmax) trims all but two points -> last-3 branch.
    c = np.array([0.0, 8.0, 5.0, 2.0, 0.1, 0.05])
    m = _metrics(c, dose=80.0)
    assert m.lambda_z_1h > 0
    assert math.isfinite(m.term_half_life_h)


def test_pk_metrics_zero_profile_gives_zero_metrics() -> None:
    m = _metrics(np.zeros(10), dose=100.0)
    assert m.auc_last_mgh_l == 0.0
    assert m.cl_l_h == 0.0
    assert m.mrt_h == 0.0
    assert m.vss_l == 0.0


def test_pk_metrics_zero_dose_gives_zero_metrics() -> None:
    c = np.exp(-0.2 * np.linspace(0, 24, 50)) * 10.0
    m = _metrics(c, dose=0.0)
    assert m.cl_l_h == 0.0
    assert m.mrt_h == 0.0
    assert m.vss_l == 0.0


def test_pk_metrics_f_abs_passthrough_and_validation() -> None:
    t = np.linspace(0.0, 24.0, 50)
    c = np.exp(-0.2 * t) * 10.0
    m = compute_pk_metrics(t, c, 100.0, Route.ORAL, f_abs=0.7)
    assert m.f_abs == pytest.approx(0.7)
    assert m.as_dict()["bioavailability_f"] == pytest.approx(0.7)
    assert compute_pk_metrics(t, c, 100.0, Route.ORAL).f_abs == 1.0
    with pytest.raises(ValueError, match="f_abs"):
        compute_pk_metrics(t, c, 100.0, Route.ORAL, f_abs=-0.1)
    with pytest.raises(ValueError, match="f_abs"):
        compute_pk_metrics(t, c, 100.0, Route.ORAL, f_abs=1.5)
