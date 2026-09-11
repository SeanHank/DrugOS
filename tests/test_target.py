"""Tests for Stage 2 (target binding): targets, panel, occupancy model."""

from __future__ import annotations

import numpy as np
import pytest

from drugos.target import (
    PanelEngagement,
    Target,
    free_binding_conc_mg_l_to_nm,
    safety_panel,
    simulate_occupancy,
    simulate_panel,
)

MW = 243.3


def _sustained_free(level_nm: float, mw: float = MW, tmax_h: float = 24.0, n: int = 601):
    t = np.linspace(0.0, tmax_h, n)
    c = np.full_like(t, level_nm * mw / 1e6)
    return t, c


def test_target_defaults_and_koff() -> None:
    site = Target(name="x", kd_nm=10.0)
    assert site.koff_1h == pytest.approx(site.kon_nm_h * site.kd_nm)
    assert site.low_confidence is True
    with pytest.raises(ValueError):
        Target(name="bad", kd_nm=0.0)
    with pytest.raises(ValueError):
        Target(name="bad", kd_nm=1.0, kon_nm_h=-1.0)
    with pytest.raises(ValueError):
        Target(name="bad", kd_nm=1.0, r0_nm=0.0)


def test_kd_from_ic50_cheng_prusoff() -> None:
    base = Target(name="BSEP (cholestasis)", kd_nm=90000.0)
    # Assay convention [S]/Km = 1 -> Ki = IC50/2 (conservative).
    meas = base.kd_from_ic50_um(90.0)
    assert not meas.low_confidence
    assert meas.kd_nm == pytest.approx(45000.0)
    assert meas.name == base.name
    # Ratio 0 -> Ki = IC50 exactly; larger ratios sharpen the affinity.
    assert base.kd_from_ic50_um(90.0, substrate_ratio=0.0).kd_nm == pytest.approx(90000.0)
    assert base.kd_from_ic50_um(90.0, substrate_ratio=3.0).kd_nm == pytest.approx(22500.0)
    with pytest.raises(ValueError):
        base.kd_from_ic50_um(0.0)
    with pytest.raises(ValueError):
        base.kd_from_ic50_um(90.0, substrate_ratio=-1.0)


def test_cheng_prusoff_helper() -> None:
    from drugos.target.targets import cheng_prusoff_ki_nm

    assert cheng_prusoff_ki_nm(1.0, 0.0) == pytest.approx(1000.0)
    assert cheng_prusoff_ki_nm(1.0, 1.0) == pytest.approx(500.0)
    assert cheng_prusoff_ki_nm(1.0, 9.0) == pytest.approx(100.0)
    with pytest.raises(ValueError):
        cheng_prusoff_ki_nm(-1.0, 1.0)
    with pytest.raises(ValueError):
        cheng_prusoff_ki_nm(1.0, -0.5)


def test_safety_panel_composition() -> None:
    panel = safety_panel()
    names = [t.name for t in panel]
    for required in (
        "hERG (Kv11.1)",
        "CYP3A4 inhibition",
        "BSEP (cholestasis)",
        "Mitochondrial complex I",
        "Glucocorticoid receptor",
    ):
        assert required in names
    assert all(t.low_confidence for t in panel)
    assert all(hasattr(t, "koff_1h") for t in panel)


def test_free_conversion_units() -> None:
    c = free_binding_conc_mg_l_to_nm(np.array([0.2433, 0.4866]), mw=243.3)
    assert c[0] == pytest.approx(1000.0)
    assert c[1] == pytest.approx(2000.0)
    with pytest.raises(ValueError):
        free_binding_conc_mg_l_to_nm(np.array([1.0]), mw=0.0)


def test_occupancy_equilibrium_half_at_kd() -> None:
    site = Target(name="eq", kd_nm=100.0, r0_nm=10.0, rho_h=0.0, kint_h=0.0)
    t, c = _sustained_free(level_nm=100.0)
    res = simulate_occupancy(t, c, site, mw=MW)
    assert res.peak_occupancy == pytest.approx(0.5, rel=0.02)
    assert res.occupancy[-1] == pytest.approx(0.5, rel=0.02)


def test_occupancy_bound_and_monotonic() -> None:
    site = Target(name="x", kd_nm=50.0, rho_h=0.0, kint_h=0.0)
    t, c = _sustained_free(level_nm=25.0)
    res = simulate_occupancy(t, c, site, mw=MW)
    assert float(np.max(res.occupancy)) <= 1.0
    assert float(np.min(res.occupancy)) >= 0.0
    tail = res.occupancy[-6:]
    assert tail[-1] >= tail[0]


def test_occupancy_turnover_and_internalization_reduce() -> None:
    turnover = Target(name="t", kd_nm=100.0, r0_nm=10.0, rho_h=0.5, kint_h=1.0)
    t, c = _sustained_free(level_nm=100.0)
    res = simulate_occupancy(t, c, turnover, mw=MW)
    assert res.peak_occupancy < 0.5


def test_time_at_target_integral() -> None:
    site = Target(name="ta", kd_nm=100.0, r0_nm=10.0, rho_h=0.0, kint_h=0.0)
    t, c = _sustained_free(level_nm=100.0, tmax_h=24.0)
    res = simulate_occupancy(t, c, site, mw=MW)
    # Sustained D=Kd => occupancy ~0.5 over the full 24 h horizon.
    assert res.time_at_target_h == pytest.approx(0.5 * 24.0, rel=0.03)
    assert res.peak_occupancy == pytest.approx(0.5, rel=0.02)


def test_washout_recovers_baseline() -> None:
    site = Target(name="w", kd_nm=100.0, r0_nm=10.0, rho_h=0.1, kint_h=0.2)
    t, c = _sustained_free(level_nm=50.0, tmax_h=48.0)
    t = np.concatenate([t, t[-1] + 1.0 + np.arange(200)])
    c = np.concatenate([c, np.zeros(200)])
    res = simulate_occupancy(t, c, site, mw=MW)
    assert res.occupancy[-1] < 0.05
    assert float(np.mean(res.occupancy[:100])) > 0.2


def test_occupancy_extra_high_affinity_full() -> None:
    site = Target(name="h", kd_nm=0.1, r0_nm=0.05, rho_h=0.0, kint_h=0.0)
    t, c = _sustained_free(level_nm=20.0)
    res = simulate_occupancy(t, c, site, mw=MW)
    assert res.peak_occupancy > 0.98


def test_occupancy_validation_of_inputs() -> None:
    site = Target(name="v", kd_nm=10.0)
    with pytest.raises(ValueError):
        simulate_occupancy(np.array([0.0]), np.array([1.0]), site, mw=MW)
    with pytest.raises(ValueError):
        simulate_occupancy(np.array([0.0, 1.0]), np.full(3, 1.0), site, mw=MW)
    with pytest.raises(ValueError):
        simulate_occupancy(np.array([0.0, 1.0]), np.array([1.0, 1.0]), site, mw=MW, r0_nm=-1.0)


def test_r0_override() -> None:
    site = Target(name="o", kd_nm=100.0, r0_nm=10.0, rho_h=0.0, kint_h=0.0)
    t, c = _sustained_free(level_nm=100.0)
    res = simulate_occupancy(t, c, site, mw=MW, r0_nm=100.0)
    assert res.receptor_nm[0] == pytest.approx(100.0)
    # Equilibrium occupancy D/(D+Kd) is R0-independent, so still ~0.5.
    assert res.peak_occupancy == pytest.approx(0.5, rel=0.02)
    # But the bound mass scales with the larger receptor pool.
    assert res.complex_nm[-1] == pytest.approx(50.0, rel=0.1)


def test_solver_failure_raises(monkeypatch) -> None:
    import types as _types

    import drugos.target.occupancy as occ

    fake = _types.SimpleNamespace(success=False, message="forced occupancy failure")
    monkeypatch.setattr(occ, "solve_ivp", lambda *a, **k: fake)
    site = Target(name="f", kd_nm=10.0)
    t, c = _sustained_free(level_nm=10.0)
    with pytest.raises(RuntimeError, match="forced occupancy failure"):
        simulate_occupancy(t, c, site, mw=MW)


def test_panel_simulation_and_ranking() -> None:
    panel = safety_panel()
    t, c = _sustained_free(level_nm=2000.0)
    eng = simulate_panel(t, c, panel, mw=118.0)
    assert isinstance(eng, PanelEngagement)
    assert len(eng.results) == len(panel)
    ranked = eng.ranked_by_time_at_target()
    assert ranked[0][1].time_at_target_h >= ranked[-1][1].time_at_target_h
    # Low-kd sites dominate the exposure ranking at equal concentration.
    hERG = eng.results["hERG (Kv11.1)"]
    cyp = eng.results["CYP3A4 inhibition"]
    assert hERG.time_at_target_h > cyp.time_at_target_h
    for row in ranked:
        assert row[1].to_series()["occupancy"].shape[0] > 1


def test_panel_to_series_contract() -> None:
    site = Target(name="s", kd_nm=10.0, rho_h=0.0, kint_h=0.0)
    t, c = _sustained_free(level_nm=10.0)
    res = simulate_occupancy(t, c, site, mw=MW)
    series = res.to_series()
    assert set(series) == {"t_h", "free_nM", "receptor_nM", "complex_nM", "occupancy"}
    delta = np.array(series["free_nM"])
    assert delta[0] == pytest.approx(10.0, rel=0.2)


def test_occupancy_long_free_drug_low_affinity() -> None:
    site = Target(name="z", kd_nm=1e6, r0_nm=1.0, rho_h=0.1, kint_h=0.1)
    t, c = _sustained_free(level_nm=1e4, tmax_h=6.0)
    res = simulate_occupancy(t, c, site, mw=MW)
    assert res.peak_occupancy < 0.05
    assert res.time_at_target_h >= 0.0
