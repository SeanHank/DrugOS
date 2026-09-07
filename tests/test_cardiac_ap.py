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
from drugos.organ.cardiac_ap import Apd90Result, apd90_from_trace, ord_apd90


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
