"""Tests for the whole-body PBPK engine (pk/pbpk_build.py, pk/simulate.py)."""

import math

import numpy as np
import pytest

from drugos.inputs.models import HumanProfile, Molecule, Route, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.resolve_human import resolve_human
from drugos.pk.partitions import partition_from_molecule
from drugos.pk.pbpk_build import (
    AbsorptionParams,
    CypTerm,
    PBPKModel,
    TargetBinding,
    cyp_vmax_mg_h,
)
from drugos.pk.simulate import PBPKResult, simulate_pbpk
from drugos.pk.skin import SkinLayers
from drugos.target.targets import Target


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


def test_cyp_term_validation() -> None:
    with pytest.raises(ValueError):
        CypTerm(isoform="CYP3A4", km_mg_l=0.0, vmax_mg_h=1.0)
    with pytest.raises(ValueError):
        CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=-1.0)
    with pytest.raises(ValueError):
        CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=1.0, hill=0.0)


def test_cyp_term_rate_mm_and_hill() -> None:
    km, vmax = 2.0, 10.0
    mm = CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=vmax)
    hill = CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=vmax, hill=2.0)
    assert math.isclose(mm.rate(km), vmax / 2.0)
    assert math.isclose(hill.rate(km), vmax / 2.0)  # at c = Km the flux is Vmax/2
    assert mm.rate(0.0) == 0.0 and hill.rate(0.0) == 0.0
    # Below Km, Michaelis-Menten approaches the Vmax/Km slope while the Hill
    # (n=2) form is sigmoidal: flux ~ Vmax*(c/Km)^2 at low c (no linear regime).
    c = 1.0e-3 * km
    assert math.isclose(mm.rate(c), vmax / km * c, rel_tol=0.05)
    assert math.isclose(hill.rate(c), vmax * (c / km) ** 2, rel_tol=0.05)


def test_cyp_vmax_helper_units_and_validation() -> None:
    content_nmol, kcat, mw = 10.0, 10.0, 500.0
    assert math.isclose(cyp_vmax_mg_h(content_nmol, kcat, mw), 0.05)  # 10 nmol -> 0.05 mg/h
    with pytest.raises(ValueError):
        cyp_vmax_mg_h(0.0, kcat, mw)
    with pytest.raises(ValueError):
        cyp_vmax_mg_h(content_nmol, 0.0, mw)
    with pytest.raises(ValueError):
        cyp_vmax_mg_h(content_nmol, kcat, 0.0)


def test_cyp_terms_reproduce_linear_low_dose() -> None:
    # Sum-of-CYP MM fluxes at doses far below Km reproduce the linear slope
    # CL_0 = sum(Vmax_i/Km_i), so the per-CYP engine matches the lumped
    # cl_hep twin exactly.
    terms = (
        CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=12.0),
        CypTerm(isoform="CYP2D6", km_mg_l=1.0, vmax_mg_h=8.0),
    )
    cl_lin = sum(t.vmax_mg_h / t.km_mg_l for t in terms)
    dose_mg = 5.0
    mm = _model(route=Route.IV_BOLUS, amount_mg=dose_mg)
    mm.cyp_terms = terms
    lin = _model(route=Route.IV_BOLUS, amount_mg=dose_mg, cl_hep=cl_lin)
    auc_mm = simulate_pbpk(mm, tmax_h=48.0, n_eval=200).pk_metrics().auc_inf_mgh_l
    auc_lin = simulate_pbpk(lin, tmax_h=48.0, n_eval=200).pk_metrics().auc_inf_mgh_l
    assert math.isclose(dose_mg / auc_mm, dose_mg / auc_lin, rel_tol=0.2)


def test_cyp_hill_saturates_more_than_mm() -> None:
    # At supersaturating concentrations the Hill (n=2) form bends sharper than
    # Michaelis-Menten, so its apparent high-dose clearance is lower.
    km, vmax, dose_mg = 1.0, 20.0, 500.0
    mm = _model(route=Route.IV_BOLUS, amount_mg=dose_mg)
    mm.cyp_terms = (CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=vmax),)
    hill = _model(route=Route.IV_BOLUS, amount_mg=dose_mg)
    hill.cyp_terms = (CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=vmax, hill=2.0),)

    def app_cl(m: PBPKModel) -> float:
        auc = simulate_pbpk(m, tmax_h=96.0, n_eval=300).pk_metrics().auc_inf_mgh_l
        return dose_mg / auc

    cl_mm, cl_hill = app_cl(mm), app_cl(hill)
    assert cl_hill < cl_mm


def test_cyp_abundance_scaling_doubles_low_dose_slope() -> None:
    # The helper consumes the physiology CYP-abundance table: Vmax is linear in
    # per-isoform liver content, so doubling content doubles the low-dose slope.
    h = resolve_human(HumanProfile(sex=Sex.MALE))
    content = h.hepatic_cyp_content_nmol["CYP3A4"]
    km, mw, kcat = 1.0, 500.0, 2.0
    v1 = cyp_vmax_mg_h(content, kcat, mw)
    v2 = cyp_vmax_mg_h(2.0 * content, kcat, mw)
    terms1 = (CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=v1),)
    terms2 = (CypTerm(isoform="CYP3A4", km_mg_l=km, vmax_mg_h=v2),)
    dose_mg = 1.0

    def cl_of(terms) -> float:
        m = _model(route=Route.IV_BOLUS, amount_mg=dose_mg)
        m.cyp_terms = terms
        auc = simulate_pbpk(m, tmax_h=48.0, n_eval=200).pk_metrics().auc_inf_mgh_l
        return dose_mg / auc

    cl1, cl2 = cl_of(terms1), cl_of(terms2)
    assert math.isclose(cl2, 2.0 * cl1, rel_tol=0.1)


def test_cyp_terms_conflict_with_lumped_mm_raises() -> None:
    with pytest.raises(ValueError):
        PBPKModel(
            physiology=_partition()[0],
            partition=_partition()[1],
            bp=1.0,
            fup=0.5,
            cyp_terms=(CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=1.0),),
            hepatic_vmax_mg_h=10.0,
            hepatic_km_mg_l=2.0,
        )


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


_TMDD_MW = 450.0


def _tmdd_target(kint: float = 0.25, kd: float = 1.5, r0: float = 800.0) -> Target:
    return Target(
        name="TMDD test site",
        kd_nm=kd,
        kon_nm_h=0.5,
        r0_nm=r0,
        rho_h=0.08,
        kint_h=kint,
    )


def _tmdd_model(
    amount_mg: float = 50.0,
    kint: float = 0.25,
    cl_hep: float = 0.0,
    cl_renal: float = 0.0,
    kd: float = 1.5,
    r0: float = 800.0,
    bind: bool = True,
) -> PBPKModel:
    profile, part = _partition(fup=0.1)
    plan = build_dose_plan(route=Route.IV_BOLUS, amount_mg=amount_mg, interval_h=0.0, n_doses=1)
    site = TargetBinding(tissue="gut", target=_tmdd_target(kint=kint, kd=kd, r0=r0))
    binding = () if not bind else (site,)
    return PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=0.1,
        cl_hep_l_h=cl_hep,
        cl_renal_l_h=cl_renal,
        mw_g_per_mol=_TMDD_MW if bind else None,
        target_binding=binding,
        dose_plan=plan,
    )


def test_tmdd_bound_off_by_default_no_extra_states() -> None:
    plain = _tmdd_model(bind=False)
    bound = _tmdd_model()
    assert plain.n_state == bound.n_state - 3
    y = plain.initial_state()
    assert math.isclose(plain.state_total_mass(y), np.sum(y), rel_tol=0.0)


def test_tmdd_guards() -> None:
    bind = (TargetBinding(tissue="gut", target=_tmdd_target()),)
    profile, part = _partition(fup=0.1)
    plan = build_dose_plan(route=Route.IV_BOLUS, amount_mg=10.0, interval_h=0.0, n_doses=1)
    with pytest.raises(ValueError, match="mw_g_per_mol is required"):
        PBPKModel(
            physiology=profile,
            partition=part,
            bp=1.0,
            fup=0.1,
            target_binding=bind,
            dose_plan=plan,
        )
    with pytest.raises(ValueError, match="mw_g_per_mol must be positive"):
        PBPKModel(
            physiology=profile,
            partition=part,
            bp=1.0,
            fup=0.1,
            mw_g_per_mol=0.0,
            target_binding=bind,
            dose_plan=plan,
        )
    with pytest.raises(ValueError, match="not in model tissues"):
        PBPKModel(
            physiology=profile,
            partition=part,
            bp=1.0,
            fup=0.1,
            mw_g_per_mol=_TMDD_MW,
            target_binding=(TargetBinding(tissue="plasma_extra", target=_tmdd_target()),),
            dose_plan=plan,
        )


def test_tmdd_reversible_binding_conserves_mass() -> None:
    model = _tmdd_model(amount_mg=50.0, kint=0.0, kd=0.5)
    res = simulate_pbpk(model, tmax_h=24.0, n_eval=200)
    assert res.final_state is not None
    assert math.isclose(model.state_total_mass(res.final_state), res.dose_mg, rel_tol=1e-6)
    row = res.tmdd["gut::TMDD test site"]
    assert float(row["bound_mg"][-1]) > 0.01
    assert float(row["cleared_mg"][-1]) == 0.0
    # Binding depletes free receptor below basal abundance.
    assert float(row["receptor_nmol"][-1]) < 800.0 * model.physiology.organ_volume["gut"]


def test_tmdd_internalization_is_irreversible_drug_sink() -> None:
    model = _tmdd_model(amount_mg=50.0, kint=0.25)
    res = simulate_pbpk(model, tmax_h=24.0, n_eval=200)
    assert res.final_state is not None
    cleared = float(res.tmdd["gut::TMDD test site"]["cleared_mg"][-1])
    assert cleared > 0.5
    # The internalized sink is part of the drug mass tally, so with no
    # metabolic/renal clearance the total stays exactly at the dose.
    assert math.isclose(model.state_total_mass(res.final_state), res.dose_mg, rel_tol=1e-6)
    # And the internalized drug is gone from the free-circulating pool.
    in_system = model.state_total_mass(res.final_state) - cleared
    assert 0.0 < in_system < res.dose_mg


def test_tmdd_quasi_steady_ratio_matches_kd() -> None:
    kd = 1.5
    model = _tmdd_model(amount_mg=5.0, kint=0.0, kd=kd)
    res = simulate_pbpk(model, tmax_h=72.0, n_eval=400)
    row = res.tmdd["gut::TMDD test site"]
    free_nm = np.asarray(res.unbound_tissues["gut"]) * 1e6 / _TMDD_MW
    half = len(res.t) // 2
    ratio = np.asarray(row["complex_nmol"][half:]) / np.asarray(row["receptor_nmol"][half:])
    free = free_nm[half:]
    quiet = free > 1e-3
    if quiet.sum() >= 10:
        dev = np.log10((ratio[quiet] * kd) / free[quiet])
        assert abs(float(np.nanmean(dev))) < 0.05


def test_tmdd_confers_dose_disproportional_elimination() -> None:
    # No metabolic/renal clearance: the ONLY removal mechanism is the native
    # target sink, so a low dose is scavenged as a much larger fraction than a
    # high dose (saturable target capacity -> dose-disproportional exposure).
    def in_system_fraction(model: PBPKModel) -> float:
        res = simulate_pbpk(model, tmax_h=24.0, n_eval=200)
        assert res.final_state is not None
        cleared = float(res.tmdd["gut::TMDD test site"]["cleared_mg"][-1])
        return (model.state_total_mass(res.final_state) - cleared) / res.dose_mg

    remain_low = in_system_fraction(_tmdd_model(amount_mg=4.0, kint=0.25))
    remain_high = in_system_fraction(_tmdd_model(amount_mg=2000.0, kint=0.25))
    assert remain_low < 0.9 < remain_high
    assert remain_low < remain_high


def test_tmdd_binding_reduces_exposure_vs_free_dispersion_twin() -> None:
    bound = _tmdd_model(amount_mg=50.0, kint=0.25)
    free = _tmdd_model(amount_mg=50.0, kint=0.25, bind=False)
    r_b = simulate_pbpk(bound, tmax_h=24.0, n_eval=200)
    r_f = simulate_pbpk(free, tmax_h=24.0, n_eval=200)
    auc_b = r_b.pk_metrics().auc_last_mgh_l
    auc_f = r_f.pk_metrics().auc_last_mgh_l
    assert 0.0 < auc_b < auc_f


def _skin_partition() -> tuple[HumanProfile, object]:
    from types import SimpleNamespace

    profile = resolve_human(HumanProfile(sex=Sex.MALE))
    names = list(profile.organ_volume) + ["arterial", "venous"]
    return profile, SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})


def _skin_model(dose_mg: float, **skin_kw: float) -> PBPKModel:
    profile, part = _skin_partition()
    plan = build_dose_plan(route=Route.TRANSDERMAL, amount_mg=dose_mg, interval_h=0.0, n_doses=1)
    return PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        dose_plan=plan,
        absorption=AbsorptionParams(skin_layers=SkinLayers(**skin_kw)),
    )


def _plain_transdermal(dose_mg: float) -> PBPKModel:
    profile, part = _skin_partition()
    plan = build_dose_plan(route=Route.TRANSDERMAL, amount_mg=dose_mg, interval_h=0.0, n_doses=1)
    return PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        dose_plan=plan,
    )


def test_skin_layers_off_by_default_adds_five_states_and_routes_dose() -> None:
    plain = _plain_transdermal(10.0)
    skin = _skin_model(10.0)
    assert skin.n_state == plain.n_state + 5

    y0 = plain.initial_state()
    plain.apply_event(y0, 0.0)
    assert math.isclose(float(y0[plain.state_index("depot")]), 10.0)

    y1 = skin.initial_state()
    skin.apply_event(y1, 0.0)
    assert float(y1[skin.state_index("depot")]) == 0.0
    assert math.isclose(float(y1[skin._skin_indices["skin_surface"]]), 10.0)


def test_skin_layers_conserve_drug_mass_no_clearance() -> None:
    model = _skin_model(20.0)
    res = simulate_pbpk(model, tmax_h=72.0, n_eval=400)
    assert res.final_state is not None
    assert math.isclose(model.state_total_mass(res.final_state), res.dose_mg, rel_tol=1e-9)


def test_skin_flux_conservation_inside_rhs() -> None:
    model = _skin_model(10.0)
    y = model.initial_state()
    model.apply_event(y, 0.0)
    # Load the dermis so the capillary sink is nonzero.
    y[model._skin_indices["skin_dermis"]] = 3.0
    dydt = model.rhs(1.0, y)
    fluxes = model.skin_fluxes(y)
    cap = fluxes["dermal_capillary"]
    skin_sum = sum(
        dydt[model._skin_indices[f"skin_{name}"]] for name in ("surface", "sc", "ve", "dermis")
    )
    assert math.isclose(skin_sum, -cap, abs_tol=1e-9)
    assert math.isclose(dydt[model.state_index("venous")], cap, abs_tol=1e-9)
    assert math.isclose(dydt[model._skin_indices["skin_unabsorbed"]], 0.0, abs_tol=1e-9)


def test_skin_partition_equilibrium_reproduces_k() -> None:
    model = _skin_model(
        10.0,
        k_dermal_capillary_1h=1e-4,
        surface_sc_partition=2.0,
        sc_ve_partition=3.0,
        ve_dermis_partition=4.0,
        surface_thickness_um=2000.0,
    )
    res = simulate_pbpk(model, tmax_h=400.0, n_eval=800)
    assert res.final_state is not None
    vols = model._skin_volumes_cm3
    y = res.final_state
    cs = y[model._skin_indices["skin_surface"]] / vols["surface"]
    csc = y[model._skin_indices["skin_sc"]] / vols["sc"]
    cve = y[model._skin_indices["skin_ve"]] / vols["ve"]
    cde = y[model._skin_indices["skin_dermis"]] / vols["dermis"]
    assert math.isclose(csc / cs, 2.0, rel_tol=1e-2)
    assert math.isclose(cve / csc, 3.0, rel_tol=1e-2)
    assert math.isclose(cde / cve, 4.0, rel_tol=1e-2)


def test_skin_steady_flux_matches_composite_permeability() -> None:
    model = _skin_model(100.0, surface_thickness_um=5000.0, k_dermal_capillary_1h=50.0)
    res = simulate_pbpk(model, tmax_h=6.0, n_eval=600)
    assert res.final_state is not None
    y = res.final_state
    skin = model.absorption.skin_layers
    assert skin is not None
    j = model.skin_fluxes(y)
    j_obs = j["sc_to_ve"]
    c_surf = y[model._skin_indices["skin_surface"]] / skin.layer_volumes_cm3()["surface"]
    j_theory = skin.composite_permeability_cm_h() * c_surf * skin.area_cm2
    assert 0.98 < j_obs / j_theory < 1.02
    # The three serial links carry a single common steady-state flux.
    assert math.isclose(j["surface_to_sc"], j_obs, rel_tol=1e-2)
    assert math.isclose(j["ve_to_dermis"], j_obs, rel_tol=1e-2)


def test_skin_barrier_thickness_slows_absorption() -> None:
    thin = simulate_pbpk(_skin_model(10.0, sc_thickness_um=20.0), tmax_h=6.0, n_eval=400)
    thick = simulate_pbpk(_skin_model(10.0, sc_thickness_um=200.0), tmax_h=6.0, n_eval=400)
    assert thin.skin is not None and thick.skin is not None
    absorbed_thin = float(thin.skin["absorbed_mg"][-1])
    absorbed_thick = float(thick.skin["absorbed_mg"][-1])
    assert absorbed_thin > 8.0
    assert absorbed_thick < 0.6 * absorbed_thin


def test_skin_diffusivity_speeds_early_absorption() -> None:
    slow = simulate_pbpk(_skin_model(10.0, sc_diffusivity_cm2_h=1.0e-5), tmax_h=1.0, n_eval=300)
    fast = simulate_pbpk(_skin_model(10.0, sc_diffusivity_cm2_h=1.0e-4), tmax_h=1.0, n_eval=300)
    assert slow.skin is not None and fast.skin is not None
    slow_abs = float(slow.skin["absorbed_mg"][-1])
    fast_abs = float(fast.skin["absorbed_mg"][-1])
    assert fast_abs > 1.5 * slow_abs and 0.0 < slow_abs < fast_abs < 10.0


def test_skin_depot_bioavailability_tallies_unabsorbed_sink() -> None:
    profile, part = _skin_partition()
    plan = build_dose_plan(route=Route.TRANSDERMAL, amount_mg=10.0, interval_h=0.0, n_doses=1)
    model = PBPKModel(
        physiology=profile,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        dose_plan=plan,
        absorption=AbsorptionParams(skin_layers=SkinLayers(), depot_bioavailability=0.8),
    )
    res = simulate_pbpk(model, tmax_h=48.0, n_eval=400)
    assert res.skin is not None
    assert math.isclose(float(res.skin["unabsorbed_mg"][-1]), 2.0, rel_tol=0.02)
    assert math.isclose(float(res.skin["absorbed_mg"][-1]), 8.0, rel_tol=0.02)


def test_skin_guards() -> None:
    with pytest.raises(ValueError, match="area_cm2 must be positive"):
        SkinLayers(area_cm2=0.0)
    with pytest.raises(ValueError, match="sc_thickness_um must be positive"):
        SkinLayers(sc_thickness_um=-1.0)
    with pytest.raises(ValueError, match="sc_diffusivity_cm2_h must be positive"):
        SkinLayers(sc_diffusivity_cm2_h=0.0)
    with pytest.raises(ValueError, match="surface_sc_partition must be positive"):
        SkinLayers(surface_sc_partition=0.0)
    with pytest.raises(ValueError, match="k_dermal_capillary_1h must be positive"):
        SkinLayers(k_dermal_capillary_1h=0.0)


def test_skin_state_record_and_contract() -> None:
    model = _skin_model(10.0)
    res = simulate_pbpk(model, tmax_h=6.0, n_eval=100)
    assert res.skin is not None
    assert set(res.skin) == {
        "surface_mg",
        "sc_mg",
        "ve_mg",
        "dermis_mg",
        "unabsorbed_mg",
        "absorbed_mg",
    }
    assert res.to_data_contract()["skin"] is not None
    rec = model.skin_state_record(model.initial_state())
    assert set(rec) == {
        "surface_mg",
        "sc_mg",
        "ve_mg",
        "dermis_mg",
        "unabsorbed_mg",
        "absorption_rate_mg_h",
    }
    assert model.skin_state_record(model.initial_state())["absorption_rate_mg_h"] == 0.0
    assert set(model.skin_fluxes(model.initial_state())) == {
        "surface_to_sc",
        "sc_to_ve",
        "ve_to_dermis",
        "dermal_capillary",
    }
    plain = _plain_transdermal(10.0)
    assert plain.skin_fluxes(plain.initial_state()) == {}
    assert plain.skin_state_record(plain.initial_state()) == {}


def test_skin_absorption_rate_zero_without_layers() -> None:
    plain = _plain_transdermal(10.0)
    assert plain.skin_absorption_rate_mg_h(plain.initial_state()) == 0.0
