"""ADMET-AI wrapper coverage tests with a structurally-typed fake model.

The heavy torch model is never loaded: ``admet_ai.ADMETModel`` is replaced by a
fake whose ``predict`` returns protocol-compatible row access so that the
mapping/validation logic in :mod:`drugos.pk.admet` is exercised end to end.
"""

import sys
import types

import pytest

from drugos.pk.admet import AdmetOutput, ADMETPredictor, predict_admet


class _FakeCol:
    def __init__(self, values) -> None:
        self._values = values

    @property
    def iloc(self):
        return self._values


class _FakeFrame:
    def __init__(self, table: dict[str, list]) -> None:
        self._table = table

    def __getitem__(self, key: str) -> _FakeCol:
        if key not in self._table:
            raise KeyError(key)
        return _FakeCol(self._table[key])


class _FakeModel:
    table: dict[str, list] = {}

    def predict(self, smiles):
        t = {k: v * len(smiles) for k, v in self.table.items()}
        return _FakeFrame(t)


_GOOD_TABLE = {
    "PPBR_AZ": [90.0],  # -> fup = 0.10
    "Solubility_AqSolDB": [-3.5],
    "Clearance_Hepatocyte_AZ": [1.2],
    "Clearance_Microsome_AZ": [-0.5],  # non-positive -> None
    "Half_Life_Obach": [4.0],
    "VDss_Lombardo": [0.5],
    "hERG": [1.5],  # clamped to 1.0
    "DILI": [-0.2],  # clamped to 0.0
    "AMES": [0.5],
    "BBB_Martins": [0.8],
    "Pgp_Broccatelli": [0.3],
    "HIA_Hou": [0.9],
    "Bioavailability_Ma": [0.7],
    "Caco2_Wang": [1.1],
    "Lipophilicity_AstraZeneca": [2.15],
    "CYP2D6_Veith": [0.1],
    "CYP3A4_Veith": [0.2],
    "CYP2C9_Veith": [0.3],
    "CYP1A2_Veith": [0.4],
    "CYP2C19_Veith": [0.5],
}


def _install(monkeypatch, table: dict[str, list]) -> None:
    monkeypatch.setattr(_FakeModel, "table", table)
    mod = types.ModuleType("admet_ai")
    mod.ADMETModel = _FakeModel
    monkeypatch.setitem(sys.modules, "admet_ai", mod)


def test_admet_predict_single_and_batch(monkeypatch) -> None:
    _install(monkeypatch, _GOOD_TABLE)

    pred = ADMETPredictor()
    out_single = pred.predict("c1ccccc1")
    assert isinstance(out_single, AdmetOutput)
    assert out_single.fup_plasma == pytest.approx(0.10)
    assert out_single.cl_int_hep_ml_min_kg == pytest.approx(1.2)
    assert out_single.cl_int_mic_ml_min_kg is None  # non-positive
    assert out_single.hERG == pytest.approx(1.0)  # clamped
    assert out_single.DILI == pytest.approx(0.0)  # clamped
    assert out_single.log_p_pred == pytest.approx(2.15)

    outs = pred.predict(["c1ccccc1", "CCO"])
    assert isinstance(outs, list) and len(outs) == 2
    assert outs[0].smiles == "c1ccccc1"
    assert outs[1].smiles == "CCO"


def test_predict_admet_one_shot(monkeypatch) -> None:
    _install(monkeypatch, _GOOD_TABLE)
    out = predict_admet("c1ccccc1")
    assert isinstance(out, AdmetOutput)


def test_predict_admet_reuses_predictor(monkeypatch) -> None:
    _install(monkeypatch, _GOOD_TABLE)
    pred = ADMETPredictor()
    assert pred._model is None
    out = predict_admet("C", predictor=pred)
    assert isinstance(out, AdmetOutput)
    # Internally created predictor is torn down; a handed-in one is kept.
    assert pred._model is not None


def _alt_table(ppbr) -> dict[str, list]:
    """Good table with a replaced PPBR cell (fup = 1 - ppbr/100 is mapped)."""
    base = {k: v for k, v in _GOOD_TABLE.items() if k != "PPBR_AZ"}
    return {**base, "PPBR_AZ": ppbr}


def test_missing_column_raises(monkeypatch) -> None:
    _install(monkeypatch, {"PPBR_AZ": [50.0]})  # every other schema column is absent
    pred = ADMETPredictor()
    with pytest.raises(KeyError):
        pred.predict("C")


def test_malformed_numeric_value_raises(monkeypatch) -> None:
    _install(monkeypatch, _alt_table(["not-a-number"]))
    pred = ADMETPredictor()
    with pytest.raises(ValueError):
        pred.predict("C")


def test_nan_cell_yields_none(monkeypatch) -> None:
    _install(monkeypatch, _alt_table([float("nan")]))
    pred = ADMETPredictor()
    out = pred.predict("C")
    assert out.fup_plasma is None
    assert out.hERG == pytest.approx(1.0)


def test_nan_prob_cell_yields_none(monkeypatch) -> None:
    base = {k: v for k, v in _GOOD_TABLE.items() if k != "hERG"}
    _install(monkeypatch, {**base, "hERG": [float("nan")]})
    out = ADMETPredictor().predict("C")
    assert out.hERG is None


def test_to_dict_excludes_raw() -> None:
    out = AdmetOutput(smiles="C", fup_plasma=0.2)
    d = out.to_dict()
    assert d["fup_plasma"] == 0.2
    assert "raw" not in d
