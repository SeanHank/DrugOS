"""R bridge tests: literature-PK cross-check (required-R path, G5-clean)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from drugos.rbridge import RVerify, _estimator_py, verify_pk

_T = np.linspace(0.0, 96.0, 400, dtype=float)
_C = 2.0 * np.exp(-0.3 * _T) + 0.9 * np.exp(-0.015 * _T)


def test_estimator_py_matches_closed_form() -> None:
    auc, cl, lam = _estimator_py(_T, _C, 2.0)
    assert math.isfinite(auc) and auc > 0
    assert math.isfinite(cl) and cl > 0
    assert abs(cl - 2.0 / auc) < 1e-12
    assert lam > 0


def test_estimator_py_no_terminal_tail() -> None:
    t = np.linspace(0.0, 10.0, 60, dtype=float)
    c = np.concatenate([np.exp(-t[:-12]), np.zeros(12)])
    auc, cl, lam = _estimator_py(t, c, 5.0)
    assert math.isfinite(auc) and math.isfinite(cl)
    assert math.isnan(lam)  # tail falls below the 5% Cmax gate -> no lambda_z


def test_estimator_py_rejects_bad_input() -> None:
    with pytest.raises(ValueError, match="dose must be positive"):
        _estimator_py(_T, _C, 0.0)
    with pytest.raises(ValueError, match="at least 4"):
        _estimator_py(np.array([0.0, 0.1], dtype=float), np.array([1.0, 0.5], dtype=float), 1.0)


def test_verify_pk_agreement_loglinear() -> None:
    v = verify_pk(_T, _C, 2.0)
    assert isinstance(v, RVerify)
    assert v.verdict == "r:agree"
    assert v.agreement_frac <= 0.02
    assert v.fit in ("loglinear", "two_comp")
    d = v.to_dict()
    assert d["verdict"] == "r:agree"
    assert set(d) == {
        "r_cl_l_h",
        "py_cl_l_h",
        "r_auc_inf_mg_h_l",
        "r_t_half_beta_h",
        "agreement_frac",
        "fit",
        "verdict",
        "note",
    }


def test_verify_pk_two_comp_fit_label() -> None:
    t = np.linspace(0.0, 96.0, 500, dtype=float)
    c = 2.0 * np.exp(-0.9 * t) + 0.9 * np.exp(-0.015 * t)
    v = verify_pk(t, c, 5.0)
    assert v.verdict == "r:agree"
    assert v.fit == "two_comp"
    assert math.isfinite(v.r_t_half_beta_h)


def test_verify_pk_mismatch_when_r_estimator_diverges(monkeypatch) -> None:
    from drugos import rbridge

    real = rbridge._run_r_fit

    def shifted(_t, _c, dose):
        out = real(_t, _c, dose)
        return {**out, "cl_l_h": out["cl_l_h"] * 1.5}

    monkeypatch.setattr(rbridge, "_run_r_fit", shifted)
    v = verify_pk(_T, _C, 2.0)
    assert v.verdict == "r:mismatch"
    assert v.agreement_frac > 0.02


def test_verify_pk_propagates_runtime_error_when_r_missing(monkeypatch) -> None:
    from drugos import rbridge

    def no_r() -> object:
        raise RuntimeError("R bridge unavailable")

    monkeypatch.setattr(rbridge, "_import_r", no_r)
    with pytest.raises(RuntimeError, match="R bridge unavailable"):
        verify_pk(_T, _C, 2.0)


def test_verify_pk_nonfinite_r_result(monkeypatch) -> None:
    from drugos import rbridge

    def nan_fit(_t, _c, _dose) -> dict[str, float]:
        return {
            "n_obs": 4.0,
            "auc_inf_mg_hl": math.nan,
            "cl_l_h": math.nan,
            "t_half_beta_h": math.nan,
            "two_comp": 0.0,
        }

    monkeypatch.setattr(rbridge, "_run_r_fit", nan_fit)
    with pytest.raises(RuntimeError, match="non-finite"):
        verify_pk(_T, _C, 2.0)


def test_verify_pk_input_validation() -> None:
    with pytest.raises(ValueError, match="dose must be positive"):
        verify_pk(_T, _C, -1.0)
    with pytest.raises(ValueError, match="at least 4"):
        verify_pk(_T[:3], _C[:3], 1.0)


def test_verify_pk_rising_short_window_tracks_numpy() -> None:
    """A window shorter than t1/2 (still-rising curve) must agree, not crash.

    Both estimators produce the same non-positive CL (a window-coverage
    artifact, not an R/numpy disagreement); the bridge reports agreement
    instead of hard-failing the pipeline.
    """
    t = np.linspace(0.0, 6.0, 62, dtype=float)
    c = 0.1 * np.exp(0.26 * t)  # monotonic rise, no peak inside the window
    v = verify_pk(t, c, 5.0)
    assert v.verdict == "r:agree"
    assert v.agreement_frac <= 0.02
    assert v.r_cl_l_h == pytest.approx(v.py_cl_l_h, rel=1e-6)


def test_verify_pk_nan_tail_not_serialized() -> None:
    t = np.linspace(0.0, 10.0, 60, dtype=float)
    c = np.concatenate([np.exp(-t[:-12]), np.zeros(12)])
    v = verify_pk(t, c, 5.0)
    assert v.verdict == "r:agree"
    assert math.isnan(v.r_t_half_beta_h)
    d = v.to_dict()
    assert d["r_t_half_beta_h"] is None  # NaN must never reach JSON
    assert isinstance(d["r_cl_l_h"], float)
    assert isinstance(d["py_cl_l_h"], float)


def test_cached_r_session_branches(monkeypatch) -> None:
    from drugos import rbridge

    # First import in this process populated the cache; force re-probe to
    # cover the None-cache branch of _import_r.
    original_imported, original_loaded = rbridge._r_imported, rbridge._r_loaded
    rbridge._r_imported = None
    rbridge._r_loaded = False
    try:
        v = verify_pk(_T, _C, 2.0)
        assert v.verdict == "r:agree"
        assert rbridge._r_imported is not None
        assert rbridge._r_loaded is True
    finally:
        rbridge._r_imported, rbridge._r_loaded = original_imported, original_loaded


def test_import_r_raises_guidance(monkeypatch) -> None:
    from drugos import rbridge

    def no_rpy2(name, *args, **kwargs):
        if name.split(".")[0] == "rpy2":
            raise ImportError("no R")
        return original_import(name, *args, **kwargs)

    original_import = __import__
    original_imported = rbridge._r_imported
    rbridge._r_imported = None
    try:
        monkeypatch.setattr("builtins.__import__", no_rpy2)
        with pytest.raises(RuntimeError, match="install R"):
            rbridge._import_r()
    finally:
        monkeypatch.undo()
        rbridge._r_imported = original_imported


def test_publishable_roundtrip() -> None:
    v = verify_pk(_T, _C, 2.0)
    d = v.to_dict()
    assert isinstance(d["py_cl_l_h"], float)
    assert 0.0 <= d["agreement_frac"] <= 1.0


def test_is_finite_flags_only_central_estimates() -> None:
    good = RVerify(1.0, 1.0, 2.0, math.nan, 0.0, "loglinear", "r:agree", "")
    assert good.is_finite() is True  # NaN t½ is acceptable; CL/AUC are what matter
    bad_cl = RVerify(math.nan, 1.0, 2.0, 1.0, 0.0, "loglinear", "r:agree", "")
    assert bad_cl.is_finite() is False
    bad_auc = RVerify(1.0, 1.0, math.inf, 1.0, 0.0, "loglinear", "r:agree", "")
    assert bad_auc.is_finite() is False
