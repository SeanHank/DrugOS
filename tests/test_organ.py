"""Stage-4 organ panel tests (liver QST, cardiac, kidney, feedback).

Covers the `drugos.organ` package to 100% branch coverage and locks the
dose-response behaviour that the Stage-4 validation cases assert.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from drugos.inputs.models import DosePlan, HumanProfile, Sex
from drugos.inputs.resolve_human import resolve_human
from drugos.organ import feedback as fb
from drugos.organ.base import auc_nm_h, free_mg_l_to_nm
from drugos.organ.cardiac import (
    CardiacParams,
    ikr_fraction,
    predict_qtc,
    simulate_cardiac,
    simulate_hemodynamics,
    tdpr_band,
)
from drugos.organ.kidney import (
    KidneyParams,
    aki_grade,
    ckdepi_2021_egfr,
    gfr_trajectory,
    nephron_injury,
    scr_from_gfr,
    scr_mg_dl_to_umol_l,
    simulate_kidney,
)
from drugos.organ.liver import (
    BileAcidParams,
    LiverParams,
    aten_floor_factor,
    bile_acid_stress,
    bsep_ki_from_ic50_nm,
    combined_stress,
    dili_grade,
    inhibition,
    liver_params_from_panel,
    mitochondrial_block,
    redox_state,
    simulate_gcdca_pbk,
    simulate_liver,
)
from drugos.pk.partitions import partition_from_molecule
from drugos.pk.pbpk_build import PBPKModel
from drugos.pk.physiology import build_human
from drugos.target.targets import safety_panel


# ---------------------------------------------------------------------------
# organ.base
# ---------------------------------------------------------------------------
def test_auc_nm_h_matches_trapezoid() -> None:
    t = np.linspace(0.0, 1.0, 11)
    c = np.linspace(0.0, 10.0, 11)
    assert auc_nm_h(t, c) == pytest.approx(np.trapezoid(c, t))


def test_auc_nm_h_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        auc_nm_h(np.array([0.0, 1.0]), np.array([1.0, 2.0, 3.0]))
    with pytest.raises(ValueError):
        auc_nm_h(np.array([1.0]), np.array([2.0]))


def test_free_mg_l_to_nm_converts() -> None:
    out = free_mg_l_to_nm(np.array([1.0, 2.0]), 100.0)
    assert out == pytest.approx(np.array([1.0e4, 2.0e4]))


# ---------------------------------------------------------------------------
# liver: unit levels
# ---------------------------------------------------------------------------
def test_inhibition_saturates() -> None:
    assert inhibition(0.0, 1.0) == 0.0
    assert inhibition(1.0, 1.0) == pytest.approx(0.5)
    assert inhibition(1.0e6, 1.0) == pytest.approx(1.0, abs=1e-6)


def test_inhibition_rejects_nonpositive_ic50() -> None:
    with pytest.raises(ValueError):
        inhibition(1.0, 0.0)


def test_mitochondrial_block_empty_tuple() -> None:
    assert mitochondrial_block(1.0e4, ()) == 0.0


def test_mitochondrial_block_worst_case() -> None:
    b = mitochondrial_block(100.0, (1.0e3, 100.0, 1.0e6))
    assert b == pytest.approx(100.0 / 200.0)


def test_redox_state_zero_exposure() -> None:
    ros, gsh, depletion = redox_state(0.0, 1.0e3)
    assert ros == 0.0
    assert gsh == 1.0
    assert depletion == 0.0


def test_redox_state_saturated_gsh_floor() -> None:
    ros, gsh, depletion = redox_state(1.0e6, 1.0e3)
    assert ros == pytest.approx(1.0, abs=1e-3)
    assert gsh == pytest.approx(0.15)
    assert depletion == pytest.approx(ros * (1.0 - ros))


def test_redox_state_rejects_nonpositive_ic50() -> None:
    with pytest.raises(ValueError):
        redox_state(1.0, 0.0)


def test_combined_stress_boundaries() -> None:
    assert combined_stress(1.0, 1.0, 1.0) == pytest.approx(0.5)
    assert combined_stress(1.0, 0.0, 0.0) == pytest.approx(1.0)
    assert combined_stress(0.0, 1.0, 1.0) == pytest.approx(0.0)


def test_kill_rate_zero_and_positive() -> None:
    from drugos.organ.liver import _kill_rate

    assert _kill_rate(0.0, 0.03, 0.25, 3.0) == 0.0
    k = _kill_rate(0.5, 0.03, 0.25, 3.0)
    assert 0.0 < k < 0.03
    with pytest.raises(ValueError):
        _kill_rate(0.5, 0.0, 0.25, 3.0)


def test_aten_floor_factor() -> None:
    assert aten_floor_factor(0.0, 0.2, 0.05) == pytest.approx(1.0)
    assert aten_floor_factor(0.5, 0.2, 0.05) == pytest.approx(0.2 + 0.8 * 0.5 * 1.05)
    assert aten_floor_factor(1.0, 0.2, 0.05) == pytest.approx(0.2)
    assert aten_floor_factor(0.0, 1.0, 0.05) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        aten_floor_factor(1.0, 1.5, 0.0)


def test_dili_grade_bands() -> None:
    assert dili_grade(1.9, 1.0) == (0, False)
    assert dili_grade(2.5, 1.0) == (1, False)
    assert dili_grade(3.2, 1.5) == (2, False)
    assert dili_grade(3.2, 2.2) == (3, True)


def test_liver_params_from_panel() -> None:
    panel = safety_panel()
    p = liver_params_from_panel(panel)
    bsep = next(t for t in panel if t.name.startswith("BSEP")).kd_nm
    assert p.bsep_ic50_nm == bsep
    assert 500.0 in p.mito_ic50_nm  # complex IV 0.5 uM
    empty = liver_params_from_panel(())
    assert empty.bsep_ic50_nm == 1.0e6
    assert empty.mito_ic50_nm == (1.0e6,)


# ---------------------------------------------------------------------------
# liver: integration
# ---------------------------------------------------------------------------
def _exposure(dose_mg: float, mw: float = 151.2) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian liver free-exposure bump; 1 g gives a ~5.85 mg/L peak."""
    t = np.linspace(0.0, 24.0, 241)
    peak = 5.85 * dose_mg / 1000.0
    c = peak * np.exp(-((t - 1.0) ** 2) / (2.0 * 3.0**2))
    return t, c


def test_liver_healthy_zero_exposure() -> None:
    t = np.linspace(0.0, 24.0, 61)
    r = simulate_liver(t, np.zeros(61), 151.2)
    assert r.dili_grade == 0
    assert r.hy_law is False
    assert r.peak_alt_uln < 2.0
    assert r.max_dead_frac < 0.01
    assert "alt_U_L" in r.to_series()


def test_liver_dose_response_therapeutic_safe() -> None:
    t, c = _exposure(1000.0)
    r = simulate_liver(t, c, 151.2)
    assert r.dili_grade == 0
    assert r.peak_alt_uln < 2.0
    assert r.peak_bilirubin_mg_dl < 1.5


def test_liver_overdose_meets_hy_law() -> None:
    t, c = _exposure(20000.0)
    r = simulate_liver(t, c, 151.2)
    assert r.dili_grade == 3
    assert r.hy_law is True
    assert r.peak_alt_uln >= 3.0
    assert r.peak_bilirubin_uln >= 2.0


def test_liver_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        simulate_liver(np.array([0.0, 1.0]), np.array([1.0, 2.0, 3.0]), 10.0)
    with pytest.raises(ValueError):
        simulate_liver(np.array([1.0]), np.array([2.0]), 10.0)


def test_liver_accepts_explicit_params() -> None:
    t, c = _exposure(1000.0)
    r = simulate_liver(t, c, 151.2, params=LiverParams(kill_max_1h=0.001))
    assert r.dili_grade == 0


def test_liver_solve_failure_raises() -> None:
    with patch("drugos.organ.liver.solve_ivp") as mock:
        mock.return_value = SimpleNamespace(success=False, message="boom")
        with pytest.raises(RuntimeError, match="boom"):
            simulate_liver(np.array([0.0, 1.0]), np.array([1.0, 1.0]), 10.0)


# ---------------------------------------------------------------------------
# liver: bile-acid PBK cholestasis anchor (R-7)
# ---------------------------------------------------------------------------
def test_bsep_ki_conversion() -> None:
    assert bsep_ki_from_ic50_nm(100.0) == pytest.approx(0.05)
    with pytest.raises(ValueError):
        bsep_ki_from_ic50_nm(0.0)


def test_bile_acid_stress_mapping() -> None:
    assert bile_acid_stress(1.0) == 0.0
    assert bile_acid_stress(0.5) == 0.0
    # At the validated 1.5x risk threshold the stress reaches 0.5.
    assert bile_acid_stress(1.5) == pytest.approx(0.5)
    s = bile_acid_stress(3.0)
    assert 0.0 < s < 1.0
    assert s > 0.9
    with pytest.raises(ValueError):
        bile_acid_stress(1.0, risk_fold=1.0)


def test_gcdca_pbk_no_drug_baseline() -> None:
    t = np.linspace(0.0, 24.0, 121)
    fold = simulate_gcdca_pbk(t, np.zeros_like(t), 45.0)
    # No drug: the intrahepatic pool stays at its daily-cycling baseline, so
    # the fold never exceeds unity (measured against the baseline peak).
    assert float(fold.max()) == pytest.approx(1.0, abs=1e-3)
    assert float(fold.min()) > 0.0


def test_gcdca_pbk_strong_inhibitor_accumulates() -> None:
    t = np.linspace(0.0, 24.0, 121)
    c = 1.0 * np.ones_like(t)  # 1 uM free hepatic, Ki = 0.1 uM (ritonavir-class)
    fold = simulate_gcdca_pbk(t, c, 0.1)
    assert float(fold.max()) > 3.0
    assert bile_acid_stress(float(fold.max())) > 0.9
    # A weak inhibitor at low exposure stays near baseline.
    weak = simulate_gcdca_pbk(t, 0.05 * np.ones_like(t), 45.0)
    assert float(weak.max()) < 1.15


def test_gcdca_pbk_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        simulate_gcdca_pbk(np.array([0.0, 1.0, 2.0]), np.array([1.0, 2.0]), 1.0)
    with pytest.raises(ValueError):
        simulate_gcdca_pbk(np.array([1.0]), np.array([2.0]), 1.0)
    with pytest.raises(ValueError):
        simulate_gcdca_pbk(np.array([0.0, 1.0]), np.array([0.0, 0.0]), 0.0)
    with pytest.raises(ValueError):
        simulate_gcdca_pbk(np.array([1.0, 0.0]), np.array([0.0, 0.0]), 1.0)


def test_gcdca_pbk_solve_failure_raises() -> None:
    with patch("drugos.organ.liver.solve_ivp") as mock:
        mock.return_value = SimpleNamespace(success=False, message="bang")
        with pytest.raises(RuntimeError, match="bang"):
            simulate_gcdca_pbk(
                np.array([0.0, 24.0]),
                np.array([0.1, 0.1]),
                45.0,
                params=BileAcidParams(),
            )


def test_gcdca_pbk_baseline_collapse_raises() -> None:
    with patch("drugos.organ.liver.solve_ivp") as mock:
        # A converged baseline whose intrahepatic pool is identically zero
        # must be rejected, not silently divided by.
        mock.return_value = SimpleNamespace(success=True, y=np.zeros((12, 250)))
        with pytest.raises(RuntimeError, match="non-positive baseline pool"):
            simulate_gcdca_pbk(
                np.array([0.0, 24.0]),
                np.array([1.0, 1.0]),
                45.0,
                params=BileAcidParams(),
            )


def test_liver_death_solve_failure_raises() -> None:
    import scipy.integrate

    real = scipy.integrate.solve_ivp
    calls: dict[str, int] = {"n": 0}

    def fake(*args: object, **kwargs: object) -> object:
        calls["n"] += 1
        if calls["n"] <= 5:
            return real(*args, **kwargs)
        return SimpleNamespace(success=False, message="death-boom")

    with patch("drugos.organ.liver.solve_ivp", side_effect=fake):
        with pytest.raises(RuntimeError, match="death-boom"):
            simulate_liver(
                np.array([0.0, 24.0]),
                np.array([0.1, 0.1]),
                151.2,
            )


# ---------------------------------------------------------------------------
# cardiac
# ---------------------------------------------------------------------------
def _physiology() -> object:
    return build_human(HumanProfile(sex=Sex.MALE, age=35, body_weight_kg=70, height_cm=175))


def test_ikr_fraction_emax() -> None:
    assert float(ikr_fraction(0.0)) == pytest.approx(1.0)
    assert float(ikr_fraction(0.5, b50=0.5)) == pytest.approx(0.5)
    assert float(ikr_fraction(np.array([1.0]))[0]) == pytest.approx(1.0 - 1.0 / 1.5)


def test_ikr_fraction_rejects_nonpositive_b50() -> None:
    with pytest.raises(ValueError):
        ikr_fraction(0.5, b50=0.0)


def test_predict_qtc_clips_blockade() -> None:
    qtc, prolong = predict_qtc(np.array([-1.0, 0.0, 2.0]), 415.0, 40.0)
    assert qtc[0] == 415.0
    assert prolong[0] == 0.0
    assert prolong[2] == pytest.approx(40.0 * (1.0 - (1.0 - 1.0 / 1.5)))
    with pytest.raises(ValueError):
        predict_qtc(np.array([0.0]), 415.0, -1.0)


def test_tdpr_band_thresholds() -> None:
    assert tdpr_band(np.array([410.0, 420.0])) == ("none", 0)
    assert tdpr_band(np.array([455.0])) == ("low", 1)
    assert tdpr_band(np.array([485.0])) == ("moderate", 2)
    assert tdpr_band(np.array([505.0])) == ("high", 3)


def test_simulate_hemodynamics_default_and_explicit() -> None:
    phy = _physiology()
    h = simulate_hemodynamics(phy, CardiacParams())
    assert h.map_mmhg == pytest.approx(93.0, rel=1e-3)
    assert h.sv_ml == pytest.approx(phy.cardiac_output_ml_min / 70.0)
    h2 = simulate_hemodynamics(
        phy,
        CardiacParams(sv_ml=90.0, inotropy=1.1, chronotropy=0.9, r_sys_mmhg_min_l=16.0),
    )
    assert h2.sv_ml == pytest.approx(99.0)
    assert h2.hr_bpm == pytest.approx(63.0)
    assert h2.co_l_min == pytest.approx(99.0 * 63.0 / 1000.0)


def test_hemodynamics_rejects_bad_constants() -> None:
    phy = _physiology()
    import drugos.organ.cardiac as cc

    with pytest.raises(ValueError):
        cc._systemic_resistance(0.0, 93.0, 5.0)
    with pytest.raises(ValueError):
        simulate_hemodynamics(phy, CardiacParams(bpm=0.0))
    with pytest.raises(ValueError):
        simulate_hemodynamics(phy, CardiacParams(sv_ml=-1.0))
    with pytest.raises(ValueError):
        simulate_hemodynamics(phy, CardiacParams(r_sys_mmhg_min_l=-1.0))
    with pytest.raises(ValueError):
        simulate_hemodynamics(phy, CardiacParams(c_art_l_mmhg=0.0))


def test_hemodynamics_solve_failure_raises() -> None:
    with patch("drugos.organ.cardiac.solve_ivp") as mock:
        mock.return_value = SimpleNamespace(success=False, message="pump")
        with pytest.raises(RuntimeError, match="pump"):
            simulate_hemodynamics(_physiology(), CardiacParams())


def test_simulate_cardiac_end_to_end() -> None:
    phy = _physiology()
    t = np.linspace(0.0, 6.0, 61)
    b = 0.5 * np.exp(-((t - 0.5) ** 2) / (2.0 * 0.3**2))
    r = simulate_cardiac(t, b, phy, n_eval=120)
    assert r.qtc_ms.max() > 415.0
    assert r.qtc_ms.min() == pytest.approx(415.0, rel=1e-3)
    assert r.delta_qtc_ms.max() > 0.0
    assert r.delta_qtc_ms.min() == pytest.approx(0.0, abs=1e-3)
    assert r.tdpr_grade >= 0
    series = r.to_series()
    assert "qtc_ms" in series and "pa_mmHg" in series


def test_simulate_cardiac_default_params() -> None:
    phy = _physiology()
    t = np.linspace(0.0, 2.0, 21)
    r = simulate_cardiac(t, np.zeros_like(t), phy)
    assert r.delta_qtc_ms.max() == pytest.approx(0.0)
    assert r.tdpr_grade == 0


def test_simulate_cardiac_rejects_bad_input() -> None:
    phy = _physiology()
    with pytest.raises(ValueError):
        simulate_cardiac(np.array([0.0, 1.0]), np.array([0.1]), phy)
    with pytest.raises(ValueError):
        simulate_cardiac(np.array([1.0]), np.array([0.1]), phy)


# ---------------------------------------------------------------------------
# kidney
# ---------------------------------------------------------------------------
def test_nephron_injury_sigmoid() -> None:
    assert float(nephron_injury(0.0, 1.0e3, 2.0)) == 0.0
    assert float(nephron_injury(-5.0, 1.0e3, 2.0)) == 0.0
    assert float(nephron_injury(1.0e3, 1.0e3, 2.0)) == pytest.approx(0.5)
    assert float(nephron_injury(1.0e6, 1.0e3, 2.0)) > 0.99


def test_nephron_injury_rejects_bad_constants() -> None:
    with pytest.raises(ValueError):
        nephron_injury(1.0, 0.0, 2.0)
    with pytest.raises(ValueError):
        nephron_injury(1.0, 1.0, 0.0)


def test_ckdepi_2021_reference_points() -> None:
    # A 60-year-old with Scr = 1.0 mg/dL is grade CKD-2 territory (~ eGFR 85-90
    # per 1.73 m2 under the 2021 race-free equation).
    egfr_m = ckdepi_2021_egfr(1.0, 60, female=False)
    assert egfr_m > 85.0
    egfr_f = ckdepi_2021_egfr(1.0, 60, female=True)
    assert egfr_f < egfr_m
    # Same sex, higher Scr (1.5) => lower eGFR, but lower Scr (0.8) => higher.
    egfr_high = ckdepi_2021_egfr(1.5, 60, female=False)
    assert egfr_high < egfr_m
    assert ckdepi_2021_egfr(0.8, 60, female=False) > egfr_m
    # BSA scaling: > 1.73 m2 raises the absolute GFR, < 1.73 lowers it.
    assert ckdepi_2021_egfr(1.0, 60, False, bsa_m2=2.0) > egfr_m
    assert ckdepi_2021_egfr(1.0, 60, False, bsa_m2=1.5) < egfr_m


def test_ckdepi_2021_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        ckdepi_2021_egfr(0.0, 60, False)
    with pytest.raises(ValueError):
        ckdepi_2021_egfr(1.0, 200, False)
    with pytest.raises(ValueError):
        ckdepi_2021_egfr(1.0, 60, False, bsa_m2=0.0)


def test_scr_units_conversion() -> None:
    assert scr_mg_dl_to_umol_l(1.0) == pytest.approx(88.4)
    with pytest.raises(ValueError):
        scr_mg_dl_to_umol_l(0.0)


def test_gfr_trajectory_floor() -> None:
    inj = np.array([0.0, 1.0])
    assert gfr_trajectory(120.0, inj) == pytest.approx(np.array([120.0, 18.0]))
    assert gfr_trajectory(120.0, inj, floor=0.0) == pytest.approx(np.array([120.0, 0.0]))
    with pytest.raises(ValueError):
        gfr_trajectory(0.0, inj)
    with pytest.raises(ValueError):
        gfr_trajectory(120.0, inj, floor=1.0)


def test_scr_from_gfr_inverse() -> None:
    g = np.array([120.0, 60.0, 30.0])
    scr = scr_from_gfr(g, 9.6)
    assert scr[0] == pytest.approx(80.0)
    assert scr[1] == pytest.approx(160.0)
    assert scr[2] == pytest.approx(320.0)
    with pytest.raises(ValueError):
        scr_from_gfr(g, 0.0)
    with pytest.raises(ValueError):
        scr_from_gfr(np.array([0.0]), 9.6)


def test_aki_grade_branches() -> None:
    assert aki_grade(3.2, 0.9) == 3
    assert aki_grade(1.2, 0.3) == 3
    assert aki_grade(2.2, 0.9) == 2
    assert aki_grade(1.6, 0.9) == 1
    assert aki_grade(1.2, 0.45) == 1
    assert aki_grade(1.2, 0.6) == 0


def test_simulate_kidney_end_to_end() -> None:
    t = np.linspace(0.0, 24.0, 49)
    c = np.full_like(t, 1.0)  # mg/L -> (nM) exposure
    r = simulate_kidney(t, c, 300.0, gfr_base_ml_min=120.0, scr_base_umol_l=80.0)
    assert r.scr_umol_l[0] == pytest.approx(80.0, rel=1e-6)
    assert r.peak_scr_ratio >= 1.0
    assert "scr_umol_L" in r.to_series()
    with pytest.raises(ValueError):
        simulate_kidney(np.array([0.0, 1.0]), np.array([1.0, 2.0, 3.0]), 10.0, 120.0, 80.0)
    with pytest.raises(ValueError):
        simulate_kidney(np.array([1.0]), np.array([2.0]), 10.0, 120.0, 80.0)
    with pytest.raises(ValueError):
        simulate_kidney(t, c, 10.0, gfr_base_ml_min=0.0, scr_base_umol_l=80.0)
    with pytest.raises(ValueError):
        simulate_kidney(t, c, 10.0, gfr_base_ml_min=120.0, scr_base_umol_l=0.0)


def test_simulate_kidney_explicit_params() -> None:
    t = np.linspace(0.0, 24.0, 49)
    c = np.full_like(t, 1.0)
    r = simulate_kidney(t, c, 300.0, 120.0, 80.0, params=KidneyParams(injury_ic50_nm=1.0e3))
    assert r.aki_grade >= 1


# ---------------------------------------------------------------------------
# feedback
# ---------------------------------------------------------------------------
def _pbpk_model() -> PBPKModel:
    phy = resolve_human(HumanProfile(sex=Sex.MALE))
    part = partition_from_molecule(
        _fake_molecule(),
        fup=0.5,
        bp=1.0,
        hematocrit=phy.hematocrit,
    )
    return PBPKModel(
        physiology=phy,
        partition=part,
        bp=1.0,
        fup=0.5,
        cl_hep_l_h=10.0,
        cl_renal_l_h=2.0,
        dose_plan=DosePlan(),
    )


def _fake_molecule() -> object:
    from drugos.inputs.models import Molecule

    return Molecule(name="t", log_p=2.0, pka_bases=[9.0])


def test_clo01_clamps() -> None:
    assert fb.clo_01(0.5) == 0.5
    assert fb.clo_01(1.5) == 1.0
    assert fb.clo_01(-1.0) == 0.0


def test_organ_feedback_factory() -> None:
    f = fb.organ_feedback(0.8, 0.6, 0.4)
    assert f.co_scale == pytest.approx(0.8)
    assert f.hepatic_cl_scale == pytest.approx(0.6)
    assert f.gfr_scale == pytest.approx(0.4)


def test_feedback_from_results_scales() -> None:
    f = fb.feedback_from_results(0.2, 0.5, 1.0)
    assert f.hepatic_cl_scale == pytest.approx(0.9)
    assert f.gfr_scale == pytest.approx(0.5)
    assert f.co_scale == pytest.approx(1.0)


def test_apply_pk_scaling_keeps_original() -> None:
    model = _pbpk_model()
    f = fb.OrganFeedback(co_scale=0.5, hepatic_cl_scale=0.8, gfr_scale=0.6)
    new = fb.apply_pk_scaling(model, f)
    assert new.cl_hep_l_h == pytest.approx(8.0)
    assert new.cl_renal_l_h == pytest.approx(1.2)
    assert new.physiology.cardiac_output_l_min == pytest.approx(
        model.physiology.cardiac_output_l_min * 0.5
    )
    assert model.cl_hep_l_h == 10.0
    assert model.physiology.cardiac_output_l_min == pytest.approx(
        new.physiology.cardiac_output_l_min * 2.0
    )
    assert new.physiology.organ_flow["liver"] == pytest.approx(
        model.physiology.organ_flow["liver"] * 0.5
    )


def test_organ_package_re_exports() -> None:
    from drugos import organ

    for name in organ.__all__:
        assert hasattr(organ, name)
