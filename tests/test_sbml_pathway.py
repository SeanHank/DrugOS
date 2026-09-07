"""Tests for Stage 3 (pathway): production SBML MAPK cascade integration.

Validates parsing of the vendored Huang/Levchenko ultrasensitive MAPK
cascade (BioModels BIOMD0000000009, CC0) and its signal-driven inhibition
response used by the pipeline default.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from drugos.pathway import parse_sbml, simulate_sbml_pathway
from drugos.pathway.sbml_pathway import _safe_formula, model_file


def _time_axis(hours: float = 72.0, n: int = 200) -> np.ndarray:
    return np.linspace(0.0, hours, n)


def test_parse_vendored_model() -> None:
    spec = parse_sbml(model_file())
    assert "PP_K" in spec.species
    assert "E1" in spec.species
    # Huang/Ferrell: 20 reactions, 22 fitted species (26 incl. outputs).
    assert len(spec.reactions) == 20
    assert spec.compartment > 0


def test_drug_free_baseline_matches_steady_state() -> None:
    t = _time_axis()
    res = simulate_sbml_pathway(t, np.zeros_like(t), n_eval=120)
    assert res.model.readout == "PP_K"
    assert res.model.species["E1"] == pytest.approx(3e-5, rel=1e-3)
    # Drug-free system should relax to its (fully activated) steady state.
    assert res.baseline["PP_K"] > 0.5
    assert res.concentrations["PP_K"][-1] == pytest.approx(res.baseline["PP_K"], rel=1e-2)


def test_monotone_inhibition_with_signal() -> None:
    t = _time_axis()
    res0 = simulate_sbml_pathway(t, np.zeros_like(t), n_eval=160)
    res1 = simulate_sbml_pathway(t, np.full_like(t, 1.0), n_eval=160)
    # Full occupancy (E1 -> 0) shuts off the cascade readout.
    assert res1.baseline is not None
    assert float(np.max(res1.concentrations["PP_K"])) < 0.05
    fc_full = res1.readout_fold_change("PP_K")
    assert float(np.min(fc_full)) < 0.05

    # Partial occupancy gives an intermediate, monotone response.
    res_mid = simulate_sbml_pathway(t, np.full_like(t, 0.9), n_eval=160)
    mid_peak = float(np.max(res_mid.concentrations["PP_K"]))
    full_peak = float(np.max(res1.concentrations["PP_K"]))
    zero_peak = float(np.max(res0.concentrations["PP_K"]))
    assert full_peak <= mid_peak <= zero_peak


def test_missing_model_file_raises() -> None:
    t = _time_axis()
    with pytest.raises(FileNotFoundError):
        simulate_sbml_pathway(t, np.zeros_like(t), path=Path("/no/such/file.xml"))


def test_mismatched_time_signal_lengths_raise() -> None:
    with pytest.raises(ValueError, match="equal-length"):
        simulate_sbml_pathway(np.linspace(0.0, 72.0, 100), np.zeros(80), n_eval=20)


def test_safe_formula_rejects_unknown_identifier() -> None:
    with pytest.raises(ValueError, match="unknown identifiers"):
        _safe_formula("k_mystery * X", allowed={"X"}, local_params=set())


def test_safe_formula_rewrites_math_functions() -> None:
    out = _safe_formula("sqrt(X) + log(Y)", allowed={"X", "Y"}, local_params=set())
    assert out == "np.sqrt(X) + np.log(Y)"


def test_parse_invalid_sbml_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.xml"
    bad.write_text("<not-sbml>", encoding="utf-8")
    with pytest.raises((ValueError, RuntimeError)):
        parse_sbml(bad)


def test_parse_skips_reaction_without_kinetic_law(tmp_path: Path) -> None:
    src = tmp_path / "min.xml"
    src.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model>
    <listOfCompartments>
      <compartment id="c" size="1"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="X" compartment="c" initialConcentration="1"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="r0">
        <listOfReactants><speciesReference species="X" stoichiometry="1"/></listOfReactants>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
""",
        encoding="utf-8",
    )
    spec = parse_sbml(src)
    assert spec.species == {"X": 1.0}
    assert len(spec.reactions) == 0
    assert spec.compartment == 1.0


def test_parse_compartment_without_size_defaults_to_one(tmp_path: Path) -> None:
    src = tmp_path / "nosize.xml"
    src.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model>
    <listOfCompartments>
      <compartment id="c"/>
    </listOfCompartments>
  </model>
</sbml>
""",
        encoding="utf-8",
    )
    spec = parse_sbml(src)
    assert spec.compartment == 1.0


def test_steady_state_solver_failure_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    class _FakeSol:
        success = False
        message = "boom"

    import drugos.pathway.sbml_pathway as mod

    t = _time_axis()
    monkeypatch.setattr(mod, "solve_ivp", lambda *a, **k: _FakeSol())
    with pytest.raises(RuntimeError, match="steady-state solve failed"):
        simulate_sbml_pathway(t, np.zeros_like(t), n_eval=20)


def test_dynamic_solver_failure_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """After the drug-free baseline converges, force the driven solve to fail."""

    class _FakeSol:
        success = False
        message = "boom"

    import drugos.pathway.sbml_pathway as mod

    real = mod.solve_ivp
    t = _time_axis()

    def flaky(*args, **kwargs):
        if kwargs.get("t_eval") is None:
            return real(*args, **kwargs)  # steady-state solve keeps working
        return _FakeSol()

    monkeypatch.setattr(mod, "solve_ivp", flaky)
    with pytest.raises(RuntimeError, match="pathway solve failed"):
        simulate_sbml_pathway(t, np.zeros_like(t), n_eval=20)
