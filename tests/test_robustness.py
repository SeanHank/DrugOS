"""Phase-6 robustness tests (drugos.robustness: D21-D23)."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pytest

from drugos.pipeline import RunSpec, run_pipeline
from drugos.robustness.population import (
    PopulationConfig,
    run_population,
    sample_profiles,
)
from drugos.robustness.sensitivity import (
    SensitivityConfig,
    first_total_indices,
    local_sensitivity,
    run_sobol_sensitivity,
    saltelli_design,
)
from drugos.robustness.uncertainty import (
    CURVE_NAMES,
    DEFAULT_KEYS,
    EnsembleConfig,
    EnsembleResult,
    apply_multipliers,
    curves_from_result,
    points_from_result,
    run_uncertainty,
    sample_multipliers,
)


# ---------------------------------------------------------------- D21
def test_ensemble_config_validation() -> None:
    with pytest.raises(ValueError):
        EnsembleConfig(n_runs=0)
    with pytest.raises(ValueError):
        EnsembleConfig(cv=-1.0)


def test_sample_multipliers() -> None:
    config = EnsembleConfig(n_runs=3, seed=1, cv=0.1)
    samples = sample_multipliers(config)
    assert len(samples) == 3
    assert set(samples[0]) == set(DEFAULT_KEYS)
    assert all(v > 0 for v in samples[0].values())
    again = sample_multipliers(config)
    assert samples == again


def test_apply_multipliers(fast_warfarin: RunSpec) -> None:
    scaled = apply_multipliers(fast_warfarin, {"fup": 0.001, "cl_hep_l_h": 2.0})
    assert scaled.fup == 1.0e-4  # clamped up
    scaled = apply_multipliers(fast_warfarin, {"fup": 1.0e6, "qt_ic50_nm": 2.0, "cns_ic50_nm": 2.0})
    assert scaled.fup == 1.0  # clamped down
    assert scaled.qt_ic50_nm == pytest.approx(2.0 * fast_warfarin.qt_ic50_nm)
    # None potency keys are skipped untouched (warfarin cns is None)
    assert scaled.cns_ic50_nm is None
    assert scaled.dose_plan.events[0].dose_mg == fast_warfarin.dose_plan.events[0].dose_mg
    with pytest.raises(ValueError):
        apply_multipliers(fast_warfarin, {"cl_hep_l_h": 0.0})
    with pytest.raises(ValueError):
        apply_multipliers(fast_warfarin, {"cl_hep_l_h": float("nan")})


def test_curves_and_points_from_result(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    curves = curves_from_result(result)
    assert set(curves) == set(CURVE_NAMES)
    points = points_from_result(result)
    for name in (
        "peak_qtc_ms",
        "peak_alt_uln",
        "overall_risk",
        "dili",
        "qt",
        "aki",
        "cns",
        "dili_grade",
    ):
        assert name in points


def test_run_uncertainty_default_config(fast_warfarin: RunSpec) -> None:
    config = EnsembleConfig(n_runs=2, seed=3, cv=0.1)
    ens = run_uncertainty(fast_warfarin, config)
    assert len(ens.verdicts) == 2
    assert ens.t_h.shape == ens.curves["plasma_total_mg_l"].shape[1:]
    band = ens.band("qtc_ms")
    assert set(band) == {"t_h", "lo", "med", "hi"}
    assert np.all(band["lo"] <= band["med"]) and np.all(band["med"] <= band["hi"])
    q = ens.point_quantiles("overall_risk")
    assert q[0] <= q[1] <= q[2]
    assert ens.verdict_counts()
    assert ens.config.n_runs == 2
    # default config path
    assert run_uncertainty(fast_warfarin, EnsembleConfig(n_runs=1)).config.n_runs == 1


def test_ensemble_result_to_dict_decimation() -> None:
    n, m = 3, 130
    rng = np.random.default_rng(5)
    curves = {
        "plasma_total_mg_l": rng.random((n, m)),
        "qtc_ms": rng.random((n, m)),
    }
    points = {"overall_risk": rng.random(n)}
    t_h = np.arange(m, dtype=float)
    ens = EnsembleResult(
        config=EnsembleConfig(n_runs=n),
        t_h=t_h,
        curves=curves,
        points=points,
        verdicts=["A", "B", "A"],
    )
    d = ens.to_dict()
    step = max(1, m // 60)
    assert len(d["band_90"]["plasma_total_mg_l"]["t_h"]) == 1 + (m - 1) // step
    assert d["verdict_counts"]["A"] == 2
    assert d["points_q5_q50_q95"]["overall_risk"]


# ---------------------------------------------------------------- D22
def test_population_config_validation() -> None:
    with pytest.raises(ValueError):
        PopulationConfig(n_individuals=0)
    with pytest.raises(ValueError):
        PopulationConfig(age_range=(90.0, 80.0))
    with pytest.raises(ValueError):
        PopulationConfig(age_range=(18.0, 101.0))
    with pytest.raises(ValueError):
        PopulationConfig(bmi_mean=10.0)
    with pytest.raises(ValueError):
        PopulationConfig(bmi_sd=0.0)


def test_sample_profiles_diversity() -> None:
    config = PopulationConfig(n_individuals=32, seed=11)
    profiles = sample_profiles(config)
    assert len(profiles) == 32
    sexes = {p.sex.value for p in profiles}
    assert sexes == {"male", "female"}
    for p in profiles:
        assert config.age_range[0] <= p.age_y <= config.age_range[1]
        assert 145.0 <= p.height_cm <= 205.0
        assert 40.0 <= p.weight_kg <= 160.0


def test_run_population(fast_warfarin: RunSpec) -> None:
    config = PopulationConfig(n_individuals=4, seed=2)
    result = run_population(fast_warfarin, config)
    assert result.overall_risk.shape == (4,)
    assert len(result.verdicts) == 4
    assert set(result.risks) == {"dili", "qt", "aki", "cns"}
    q = result.risk_quantiles("qt")
    assert q[0] <= q[1] <= q[2]
    counts = result.incidence(1)
    assert set(counts) == set(result.grades)
    assert result.incidence(2).keys() == counts.keys()
    assert result.verdict_counts()
    d = result.to_dict()
    assert d["n_individuals"] == 4
    assert len(d["body"]["weight_kg"]) == 4
    assert d["incidence_grade_ge_1"]
    assert len(result.sexes) == 4


# ---------------------------------------------------------------- D23
def test_sensitivity_config_validation() -> None:
    with pytest.raises(ValueError):
        SensitivityConfig(delta=0.0)
    with pytest.raises(ValueError):
        SensitivityConfig(delta=1.0)
    with pytest.raises(ValueError):
        SensitivityConfig(sobol_n=0)
    with pytest.raises(ValueError):
        SensitivityConfig(sobol_outputs=())


def test_local_sensitivity_default_evaluator(fast_warfarin: RunSpec) -> None:
    config = SensitivityConfig(delta=0.2, keys=("cl_hep_l_h", "fup"))
    res = local_sensitivity(fast_warfarin, output="overall_risk", config=config)
    assert res.baseline > 0.0
    assert set(res.values) == {"cl_hep_l_h", "fup"}
    ranked = res.drivers
    assert abs(ranked[0][1]) >= abs(ranked[-1][1])


def test_local_sensitivity_missing_output_is_explicit_error(fast_warfarin: RunSpec) -> None:
    def evaluator(_: RunSpec) -> Mapping[str, float]:
        return {"other": 1.0}

    with pytest.raises(KeyError):
        local_sensitivity(
            fast_warfarin,
            output="nonexistent",
            evaluator=evaluator,  # type: ignore[arg-type]
        )


def test_local_sensitivity_custom_evaluator(fast_warfarin: RunSpec) -> None:
    def evaluator(spec: RunSpec) -> Mapping[str, float]:
        return {"overall_risk": spec.cl_hep_l_h}

    config = SensitivityConfig(delta=0.1, keys=("cl_hep_l_h",))
    res = local_sensitivity(
        fast_warfarin, output="overall_risk", config=config, evaluator=evaluator
    )
    assert res.values["cl_hep_l_h"] == pytest.approx(1.0)


def test_saltelli_design() -> None:
    a, b = saltelli_design(8, 3, seed=9)
    assert a.shape == (8, 3) and b.shape == (8, 3)
    assert np.all((a >= 0) & (a < 1)) and np.all((b >= 0) & (b < 1))
    with pytest.raises(ValueError):
        saltelli_design(0, 3, 1)
    with pytest.raises(ValueError):
        saltelli_design(8, 0, 1)


def test_first_total_indices_normal() -> None:
    rng = np.random.default_rng(4)
    n, d = 16, 3
    ya = rng.random(n)
    yb = rng.random(n)
    yc = rng.random((d, n))
    s1, st = first_total_indices(ya, yb, yc)
    assert s1.shape == (d,) and st.shape == (d,)
    assert np.all(s1 >= -1.0) and np.all(s1 <= 1.0)
    assert np.all(st >= 0.0) and np.all(st <= 1.0)


def test_first_total_indices_edges() -> None:
    n, d = 8, 2
    flat = np.zeros(n)
    s1, st = first_total_indices(flat, flat, np.zeros((d, n)))
    assert np.all(s1 == 0.0) and np.all(st == 0.0)
    s1, st = first_total_indices(np.ones(n), np.ones(n), np.zeros((2, n, 1)))
    assert np.all(s1 == 0.0)


def test_first_total_indices_clip() -> None:
    n, d = 3, 2
    ya = np.zeros(n)
    yb = np.ones(n)
    yc = -1.0 * np.ones((d, n))
    s1, st = first_total_indices(ya, yb, yc)
    assert np.all(s1 == 1.0)
    assert np.all(st == 1.0)


def test_run_sobol_with_injected_evaluator(fast_warfarin: RunSpec) -> None:
    config = SensitivityConfig(
        keys=("cl_hep_l_h", "fup"), sobol_n=4, sobol_outputs=("qt", "overall_risk")
    )

    def evaluator(multipliers: dict[str, float]) -> Mapping[str, float]:
        total = sum(multipliers.values())
        return {"qt": total, "overall_risk": total}

    result = run_sobol_sensitivity(fast_warfarin, config, evaluator)  # type: ignore[arg-type]
    assert result.design_n == 4
    assert set(result.indices) == {"qt", "overall_risk"}
    assert result.to_dict()["keys"] == ["cl_hep_l_h", "fup"]
    drivers = result.first_drivers("qt")
    assert len(drivers) == 2
    assert abs(drivers[0][1]) >= abs(drivers[1][1])
    assert len(result.total_drivers("overall_risk")) == 2


def test_run_sobol_with_default_evaluator(fast_warfarin: RunSpec) -> None:
    config = SensitivityConfig(
        keys=("cl_hep_l_h", "fup"), sobol_n=2, sobol_outputs=("overall_risk",)
    )
    result = run_sobol_sensitivity(fast_warfarin, config)
    s1, st = result.indices["overall_risk"]
    assert np.all(np.isfinite(s1))
    assert np.all(np.isfinite(st))
