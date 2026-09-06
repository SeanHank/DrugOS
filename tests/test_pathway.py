"""Tests for Stage 3 (pathway): DSL, ODE engine, dose-response, amplification."""

from __future__ import annotations

import numpy as np
import pytest

from drugos.pathway import (
    Activation,
    DoseResponseFit,
    FirstOrderDegradation,
    PathwayModel,
    PathwayResult,
    ReversibleBinding,
    SourceProduction,
    compile_model,
    dose_response,
    mapk_cascade,
    pathway_steady_state,
    simulate_pathway,
)


def _linear_two_step() -> PathwayModel:
    """Basal + signal arms, degradation, and source-maintained pools."""
    model = PathwayModel(
        name="linear",
        species={"x": 1.0, "y": 0.0, "z": 0.0},
        input_node="signal",
        readout="z",
    )
    model.add_activation(
        Activation(
            substrate="x",
            product="y",
            driver="x",
            vmax=0.05,
            km=0.5,
            driver_half=5.0,
        )
    )
    model.add_activation(Activation(substrate="x", product="y", vmax=1.0, km=0.5, input_drive=True))
    model.add_activation(
        Activation(substrate="y", product="z", vmax=1.5, km=0.5, driver="y", driver_half=1.0)
    )
    model.add_reaction(SourceProduction(species="x", rate=0.2))
    model.add_reaction(FirstOrderDegradation(species="x", rate=0.1))
    model.add_reaction(FirstOrderDegradation(species="y", rate=0.2))
    model.add_reaction(FirstOrderDegradation(species="z", rate=0.2))
    return model


def test_reaction_validation() -> None:
    with pytest.raises(ValueError):
        Activation(substrate="x", product="y", vmax=0.0, km=1.0, input_drive=True)
    with pytest.raises(ValueError):
        Activation(substrate="x", product="y", vmax=1.0, km=1.0)
    with pytest.raises(ValueError):
        Activation(substrate="x", product="y", vmax=1.0, km=1.0, driver="e", input_drive=True)
    with pytest.raises(ValueError):
        ReversibleBinding(a="a", b="b", complex="ab", kon=0.0)
    with pytest.raises(ValueError):
        FirstOrderDegradation(species="x", rate=0.0)
    with pytest.raises(ValueError):
        SourceProduction(species="x", rate=-1.0)


def test_compile_ordering() -> None:
    model = _linear_two_step()
    compiled = compile_model(model)
    assert compiled.state_names() == sorted(model.species)
    assert compiled.n_species == 3


def test_steady_state_zero_signal_baseline() -> None:
    model = _linear_two_step()
    ss = pathway_steady_state(model, horizon_h=200.0, n_eval=80)
    assert ss["y"] > 0.0
    assert ss["z"] > 0.0
    ss2 = pathway_steady_state(model, horizon_h=200.0, n_eval=80)
    assert ss["z"] == pytest.approx(ss2["z"], rel=1e-3)


def test_signal_drives_readout_and_fold_change() -> None:
    model = _linear_two_step()
    t = np.linspace(0.0, 24.0, 100)
    signal = np.concatenate([np.full(50, 0.0), np.full(50, 1.0)])
    res = simulate_pathway(model, t, signal, n_eval=150)
    assert res.t_h.shape[0] == 150
    fold = res.readout_fold_change()
    assert float(fold[-1]) > float(fold[0])
    assert float(fold[-1]) > 1.0


def test_no_readout_raises() -> None:
    model = PathwayModel(name="plain", species={"x": 1.0})
    res = PathwayResult(
        model=model,
        t_h=np.array([0.0]),
        concentrations={"x": np.array([0.5])},
        signal=np.array([0.0]),
        baseline={"x": 1.0},
    )
    with pytest.raises(ValueError):
        res.readout_fold_change()
    with pytest.raises(ValueError):
        dose_response(model, signals=np.array([0.0, 1.0]))


def test_fold_change_rejects_zero_baseline() -> None:
    model = PathwayModel(name="plain", species={"x": 1.0}, readout="x")
    res = PathwayResult(
        model=model,
        t_h=np.array([0.0]),
        concentrations={"x": np.array([0.5])},
        signal=np.array([0.0]),
        baseline={"x": 0.0},
    )
    with pytest.raises(ValueError, match="baseline"):
        res.readout_fold_change()


def test_reversible_binding_conserves_mass() -> None:
    model = PathwayModel(
        name="binding",
        species={"a": 5.0, "b": 5.0, "ab": 0.0},
        readout="ab",
    )
    model.add_reaction(ReversibleBinding(a="a", b="b", complex="ab", kon=1.0, koff=0.2))
    model.add_reaction(SourceProduction(species="a", rate=1e-9))
    res = simulate_pathway(model, np.linspace(0.0, 5.0, 50), np.zeros(50), n_eval=100)
    a = res.concentrations["a"]
    b = res.concentrations["b"]
    ab = res.concentrations["ab"]
    # Each binding pool is conserved: A: a+ab = 5, B: b+ab = 5.
    assert float((a + ab)[-1]) == pytest.approx(5.0, rel=1e-6)
    assert float((b + ab)[-1]) == pytest.approx(5.0, rel=1e-6)
    eq = float(ab[-1]) / (float(a[-1]) * float(b[-1]))
    assert eq == pytest.approx(1.0 / 0.2, rel=1e-2)


def test_binding_loose_complex_mostly_free() -> None:
    model = PathwayModel(
        name="loose",
        species={"a": 5.0, "b": 5.0, "ab": 0.0},
        readout="ab",
    )
    model.add_reaction(ReversibleBinding(a="a", b="b", complex="ab", kon=0.1, koff=100.0))
    t = np.linspace(0.0, 3.0, 40)
    res = simulate_pathway(model, t, np.zeros(40), n_eval=80)
    assert float(res.concentrations["ab"][-1]) < 0.1


def test_source_and_degradation_balance() -> None:
    model = PathwayModel(name="bal", species={"p": 0.0, "q": 0.0})
    model.add_reaction(SourceProduction(species="p", rate=0.5))
    model.add_reaction(FirstOrderDegradation(species="p", rate=0.1))
    model.add_reaction(SourceProduction(species="q", rate=0.25))
    ss = pathway_steady_state(model, horizon_h=100.0, n_eval=60)
    assert ss["p"] == pytest.approx(5.0, rel=0.02)
    assert ss["q"] > 10.0


def test_inhibition_blocks_readout() -> None:
    def run(with_inhibitor: bool) -> float:
        model = PathwayModel(
            name="inh",
            species={"s": 400.0, "i": 200.0 if with_inhibitor else 0.0, "z": 0.0},
            readout="z",
        )
        model.add_activation(
            Activation(
                substrate="s",
                product="z",
                vmax=2.0,
                km=0.3,
                input_drive=True,
                inhibitor="i",
                ki=10.0,
            )
        )
        model.add_reaction(FirstOrderDegradation(species="z", rate=0.05))
        res = simulate_pathway(model, np.linspace(0.0, 24.0, 100), np.ones(100), n_eval=120)
        return float(res.concentrations["z"][-1])

    assert run(with_inhibitor=False) > 3.0 * run(with_inhibitor=True)


def test_input_signal_clamped_to_unit_support() -> None:
    model = _linear_two_step()
    t = np.linspace(0.0, 4.0, 40)
    r1 = simulate_pathway(model, t, np.full(40, 1.0), n_eval=80)
    r5 = simulate_pathway(model, t, np.full(40, 5.0), n_eval=80)
    rn = simulate_pathway(model, t, np.full(40, -3.0), n_eval=80)
    assert float(r5.signal[-1]) == pytest.approx(1.0)
    assert float(rn.signal[-1]) == pytest.approx(0.0)
    assert np.allclose(r1.concentrations["z"], r5.concentrations["z"], rtol=1e-9)


def test_integration_input_mismatch_raises() -> None:
    model = _linear_two_step()
    with pytest.raises(ValueError):
        simulate_pathway(model, np.linspace(0.0, 1.0, 10), np.zeros(12))
    with pytest.raises(ValueError):
        simulate_pathway(model, np.linspace(0.0, 1.0, 10), np.zeros((2, 5)))


def test_solver_failure_raises(monkeypatch) -> None:
    import types as _types

    import drugos.pathway.simulator as sim

    fake = _types.SimpleNamespace(success=False, message="forced pathway failure")
    monkeypatch.setattr(sim, "solve_ivp", lambda *a, **k: fake)
    model = _linear_two_step()
    with pytest.raises(RuntimeError, match="forced pathway failure"):
        simulate_pathway(model, np.linspace(0.0, 1.0, 10), np.zeros(10))


def test_dose_response_fit_prediction() -> None:
    x = np.linspace(0.02, 0.98, 25)
    truth = DoseResponseFit(r0=1.0, emax=9.0, ec50=0.3, hill=1.5)
    y = truth.predict(x)
    refit = DoseResponseFit(r0=1.0, emax=9.0, ec50=0.3, hill=1.5)
    assert np.allclose(refit.predict(x), y)


def test_dose_response_on_model_amplification() -> None:
    model = mapk_cascade()
    x = np.linspace(0.01, 0.99, 21)
    fit = dose_response(model, signals=x, horizon_h=300.0, n_eval=120)
    assert 0.0 < fit.ec50 < 0.5
    assert fit.hill > 0.0
    pred = fit.predict(x)
    assert pred[-1] > pred[0]
    # Amplification: half-max readout at occupancy below 0.5 == drug below Kd.
    d50 = fit.drug_ec50_nm(kd_nm=100.0)
    assert d50 < 100.0


def test_mapk_cascade_dose_response_points_increase() -> None:
    model = mapk_cascade()
    x = np.linspace(0.05, 0.95, 9)
    fit = dose_response(model, signals=x, horizon_h=200.0, n_eval=100)
    assert fit.emax > fit.r0
    assert fit.predict(np.array(x)[-1]) > fit.predict(np.array(x)[0])


def test_drug_ec50_requires_bounded_ec50() -> None:
    fit = DoseResponseFit(r0=1.0, emax=9.0, ec50=0.5, hill=1.0)
    assert fit.drug_ec50_nm(100.0) == pytest.approx(100.0)
    bad = DoseResponseFit(r0=1.0, emax=9.0, ec50=1.5, hill=1.0)
    with pytest.raises(ValueError):
        bad.drug_ec50_nm(100.0)


def test_mapk_cascade_perturbation_reaches_plateau() -> None:
    model = mapk_cascade()
    t = np.linspace(0.0, 300.0, 200)
    signal = np.clip(np.linspace(0.0, 1.5, 200), 0.0, 1.0)
    res = simulate_pathway(model, t, signal, n_eval=120)
    erk = res.concentrations["erk_active"]
    assert float(erk[-1]) > float(erk[0])
    fold = res.readout_fold_change()
    assert float(fold[-1]) > 2.0


def test_empty_pathway_initial_condition() -> None:
    model = PathwayModel(name="empty")
    compiled = compile_model(model)
    assert compiled.initial_condition(model).shape == (0,)


def test_source_only_model_and_signal_arg_checks() -> None:
    # SourceProduction as the final reaction exercises the loop-exit path.
    model = PathwayModel(name="src", species={"p": 0.0}, readout="p")
    model.add_reaction(SourceProduction(species="p", rate=0.7))
    ss = pathway_steady_state(model, horizon_h=10.0, n_eval=20)
    assert ss["p"] > 5.0
    with pytest.raises(ValueError, match="1-D"):
        dose_response(model, signals=np.zeros((2, 3)))
