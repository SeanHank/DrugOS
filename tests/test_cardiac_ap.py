"""Coverage for the ORd cardiac APD90 cross-check lane (G3, 100 % branch).

The myokit-coupled runner is exercised in unit tests with a mocked myokit
(short-circuited via ``sys.modules`` and ``_import_myokit``); the end-to-end
(real compiled) simulation runs in the R-3 validation case
(validation/cases/case_cardiac_ap_ord.py), not in pytest.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

import drugos.organ.cardiac_ap as cap
from drugos.organ.cardiac_ap import (
    Apd90Result,
    MultiIonicResult,
    apd90_from_trace,
    cipa_apd90,
    ord_apd90,
)


def _beat_trace(apd90_ms: float) -> tuple[np.ndarray, np.ndarray]:
    t = np.arange(501, dtype=float)
    v = np.full_like(t, -85.0)
    rise_end = 100
    v[rise_end:] = 40.0
    v[rise_end + int(apd90_ms) :] = -85.0
    return t, v


def _runner(result: Apd90Result):
    t, v = _beat_trace(result.apd90_base_ms)

    def run(block_frac: float) -> tuple[np.ndarray, np.ndarray]:
        if block_frac == 0.0:
            return t, v
        return _beat_trace(result.apd90_block_ms)

    return run


def test_apd90_from_trace_nominal() -> None:
    t, v = _beat_trace(250.0)
    assert apd90_from_trace(t, v) == pytest.approx(250.0, rel=0.01)


def test_apd90_from_trace_no_crossing_returns_nan() -> None:
    t = np.linspace(0.0, 1000.0, 100)
    v = np.full_like(t, 40.0)
    assert np.isnan(apd90_from_trace(t, v))


def test_apd90_from_trace_shape_mismatch_raises() -> None:
    with pytest.raises(ValueError, match="equal-length"):
        apd90_from_trace(np.array([0.0, 1.0]), np.array([0.0]))


def test_apd90_from_trace_too_short_raises() -> None:
    with pytest.raises(ValueError, match="at least two"):
        apd90_from_trace(np.array([0.0]), np.array([0.0]))


def test_apd90_from_trace_two_samples_recovers_half_width() -> None:
    t = np.array([0.0, 1.0])
    v = np.array([40.0, -85.0])
    assert apd90_from_trace(t, v) == pytest.approx(1.0)


def test_import_myokit_short_circuits_sys_modules(monkeypatch) -> None:
    class FakeMod:
        pass

    monkeypatch.setitem(sys.modules, "myokit", FakeMod)
    assert cap._import_myokit() is FakeMod


def test_import_myokit_missing_propagates(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "myokit", None)
    with pytest.raises(ImportError):
        cap._import_myokit()


def test_load_and_run_scales_ikr_conductance(monkeypatch) -> None:
    pre_calls: list[float] = []

    class FakeSim:
        def __init__(self) -> None:
            self.constants: dict[str, float] = {}

        def set_constant(self, name: str, value: float) -> None:
            self.constants[name] = value

        def pre(self, duration: float) -> None:
            pre_calls.append(duration)

        def run(self, duration: float, log: list[str], log_interval: float) -> dict[str, object]:
            t = np.linspace(0.0, 500.0, 500)
            v = np.full_like(t, -85.0)
            return {"engine.time": t, "membrane.V": v}

    class FakeVariable:
        def value(self) -> float:
            return 0.046

    class FakeModel:
        def get(self, name: str) -> FakeVariable:
            assert name == "ikr.gKr"
            return FakeVariable()

    class FakeMyokit:
        def __init__(self) -> None:
            self.sims: list[FakeSim] = []

        def load(self, path: str) -> tuple[FakeModel, str, str]:
            assert path.endswith("ohara-2011.mmt")
            return FakeModel(), "protocol", "script"

        def Simulation(self, model: FakeModel, protocol: str) -> FakeSim:
            sim = FakeSim()
            self.sims.append(sim)
            return sim

    fake = FakeMyokit()
    monkeypatch.setattr(cap, "_import_myokit", lambda: fake)

    t, v = cap._load_and_run(Path("x/ohara-2011.mmt"), 0.3, 0, 50, 1000.0, 0.1)
    assert t.shape == v.shape == (500,)
    assert fake.sims[0].constants == {"cell.mode": 0, "ikr.gKr": 0.046 * 0.7}
    assert pre_calls == [50_000.0]

    cap._load_and_run(Path("x/ohara-2011.mmt"), 0.0, 0, 50, 1000.0, 0.1)
    assert fake.sims[1].constants == {"cell.mode": 0}


def test_ord_apd90_block_frac_out_of_range_raises() -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        ord_apd90(1.5)


def test_ord_apd90_uses_provided_runner() -> None:
    res = ord_apd90(0.5, runner=_runner(Apd90Result(0.5, 260.0, 330.0, 70.0)))
    assert res.apd90_base_ms == pytest.approx(260.0, rel=0.01)
    assert res.apd90_block_ms == pytest.approx(330.0, rel=0.01)
    assert res.delta_apd90_ms == pytest.approx(70.0, rel=0.02)
    assert res.to_dict()["block_frac"] == 0.5


def test_ord_apd90_default_runner_builds_model_lane(monkeypatch) -> None:
    calls: list[float] = []

    def fake_load_and_run(
        model_path: Path,
        block_frac: float,
        cell_mode: int,
        pre_paces: int,
        run_ms: float,
        log_interval_ms: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        assert model_path.name == "ohara-2011.mmt"
        assert cell_mode == 0
        calls.append(block_frac)
        return _beat_trace(260.0 if block_frac == 0.0 else 330.0)

    monkeypatch.setattr(cap, "_load_and_run", fake_load_and_run)
    res = ord_apd90(0.4)
    assert calls == [0.0, 0.4]
    assert res.block_frac == 0.4
    assert res.delta_apd90_ms == pytest.approx(70.0, rel=0.02)


def test_validate_block_map_rejects_unknown_current() -> None:
    with pytest.raises(ValueError, match="unknown current system"):
        cap._validate_block_map({"IKr": 0.5, "IRyR": 0.1})
    assert "IKr" in str(cap.CIPA_BLOCKABLE_CONSTANTS)


def test_validate_block_map_rejects_out_of_range() -> None:
    with pytest.raises(ValueError, match="within \\[0, 1\\]"):
        cap._validate_block_map({"INaL": 1.5})
    with pytest.raises(ValueError, match="within \\[0, 1\\]"):
        cap._validate_block_map({"ICaL": -0.1})


def test_validate_block_map_accepts_subset() -> None:
    cap._validate_block_map({"IKr": 0.0, "IKs": 1.0})
    assert set(cap.CIPA_BLOCKABLE_CONSTANTS) == {"IKr", "IKs", "IK1", "Ito", "INaL", "ICaL"}


def test_load_and_run_cipa_scales_literals(monkeypatch) -> None:
    class FakeSim:
        def __init__(self) -> None:
            self.constants: dict[str, float] = {}
            self.pre_calls: list[float] = []

        def set_constant(self, name: str, value: float) -> None:
            self.constants[name] = value

        def pre(self, duration: float) -> None:
            self.pre_calls.append(duration)

        def run(self, duration: float, log: list[str], log_interval: float) -> dict[str, object]:
            t = np.linspace(0.0, 500.0, 500)
            return {"engine.time": t, "membrane.V": np.full_like(t, -85.0)}

    values = {
        "ikr.gKr": 0.046,
        "iks.gKs": 0.0034,
        "ik1.gK1": 0.1908,
        "ito.gto": 0.02,
        "inal.gNaL": 0.0075,
        "ical.PCa_base": 0.0001,
    }

    class FakeVariable:
        def __init__(self, v: float) -> None:
            self.v = v

        def value(self) -> float:
            return self.v

    class FakeModel:
        def __init__(self) -> None:
            self.gotten: list[str] = []

        def get(self, name: str) -> FakeVariable:
            self.gotten.append(name)
            return FakeVariable(values[name])

    class FakeMyokit:
        def __init__(self) -> None:
            self.sims: list[FakeSim] = []

        def load(self, path: str) -> tuple[FakeModel, str, str]:
            assert path.endswith("ohara-cipa-v1-2017.mmt")
            return FakeModel(), "protocol", "script"

        def Simulation(self, model: FakeModel, protocol: str) -> FakeSim:
            sim = FakeSim()
            self.sims.append(sim)
            return sim

    fake = FakeMyokit()
    monkeypatch.setattr(cap, "_import_myokit", lambda: fake)
    t, v = cap._load_and_run_cipa(
        Path("x/ohara-cipa-v1-2017.mmt"),
        {"IKr": 0.5, "INaL": 0.25, "ICaL": 0.1},
        0,
        30_000.0,
        1500.0,
        0.05,
    )
    assert t.shape == v.shape == (500,)
    assert fake.sims[0].constants == {
        "cell.mode": 0,
        "ikr.gKr": 0.046 * 0.5,
        "inal.gNaL": 0.0075 * 0.75,
        "ical.PCa_base": 0.0001 * 0.9,
    }
    assert fake.sims[0].pre_calls == [30_000.0]

    cap._load_and_run_cipa(Path("x/ohara-cipa-v1-2017.mmt"), {"IKs": 0.0}, 0, 1.0, 2.0, 0.1)
    assert fake.sims[1].constants == {"cell.mode": 0}


def test_cipa_apd90_uses_provided_runner() -> None:
    t, v = _beat_trace(260.0)

    def run(block_map: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
        if not block_map:
            return t, v
        return _beat_trace(345.0)

    res = cipa_apd90({"IKr": 0.5, "INaL": 0.25}, runner=run)
    assert isinstance(res, MultiIonicResult)
    assert res.block_map == {"IKr": 0.5, "INaL": 0.25}
    assert res.apd90_base_ms == pytest.approx(260.0, rel=0.01)
    assert res.apd90_block_ms == pytest.approx(345.0, rel=0.01)
    assert res.delta_apd90_ms == pytest.approx(85.0, rel=0.02)
    assert res.to_dict()["delta_apd90_ms"] == pytest.approx(85.0, rel=0.02)


def test_cipa_apd90_empty_map_control() -> None:
    t, v = _beat_trace(300.0)

    def run(block_map: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
        assert block_map == {}
        return t, v

    res = cipa_apd90({}, runner=run)
    assert abs(res.delta_apd90_ms) < 1e-6


def test_cipa_apd90_default_runner_builds_model_lane(monkeypatch) -> None:
    calls: list[dict[str, float]] = []

    def fake_load_and_run_cipa(
        model_path: Path,
        block_map: dict[str, float],
        cell_mode: int,
        pre_ms: float,
        run_ms: float,
        log_interval_ms: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        assert model_path.name == "ohara-cipa-v1-2017.mmt"
        assert cell_mode == 0
        assert pre_ms == 30_000.0
        calls.append(block_map)
        return _beat_trace(260.0) if not block_map else _beat_trace(350.0)

    monkeypatch.setattr(cap, "_load_and_run_cipa", fake_load_and_run_cipa)
    res = cipa_apd90({"IKr": 0.5})
    assert calls == [{}, {"IKr": 0.5}]
    assert res.delta_apd90_ms == pytest.approx(90.0, rel=0.02)


def test_cipa_apd90_rejects_bad_map() -> None:
    with pytest.raises(ValueError, match="unknown current system"):
        cipa_apd90({"QX": 0.5}, runner=_beat_trace)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="within \\[0, 1\\]"):
        cipa_apd90({"IKr": 2.0}, runner=_beat_trace)  # type: ignore[arg-type]
