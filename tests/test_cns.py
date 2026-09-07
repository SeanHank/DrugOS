"""Phase-7 CNS organ-panel tests (drugos.organ.cns)."""

from __future__ import annotations

import numpy as np
import pytest

from drugos.organ.cns import CnsParams, cns_grade, kpu_brain_from_bbb, simulate_cns


def test_cns_params_validation() -> None:
    with pytest.raises(ValueError):
        CnsParams(ic50_nm=0)
    with pytest.raises(ValueError):
        CnsParams(kpu_brain=-1)
    with pytest.raises(ValueError):
        CnsParams(ratio_thresholds=(0.1, 0.2, 0.3))
    with pytest.raises(ValueError):
        CnsParams(ratio_thresholds=(0.1, 0.3, 0.2, 0.4))


def test_cns_grade_ladder() -> None:
    thresholds = (0.1, 0.32, 1.0, 3.2)
    assert cns_grade(0.0, thresholds) == 0
    assert cns_grade(0.099, thresholds) == 0
    assert cns_grade(0.1, thresholds) == 1
    assert cns_grade(0.5, thresholds) == 2
    assert cns_grade(1.0, thresholds) == 3
    assert cns_grade(3.2, thresholds) == 4
    assert cns_grade(1e6, thresholds) == 4


def test_simulate_cns_free_drug_hypothesis() -> None:
    t = np.array([0.0, 1.0, 2.0])
    free_mg_l = np.array([0.0, 0.151, 0.0])  # for mw=151 -> 1 uM at t=1
    result = simulate_cns(t, free_mg_l, mw=151.0)
    assert result.t_h.shape == t.shape
    assert result.brain_free_nm.shape == t.shape
    assert result.peak_brain_free_nm == pytest.approx(1e3 * 0.8, rel=1e-9)
    ratio = result.exposure_ratio
    assert ratio == pytest.approx(800.0 / 1.0e5, rel=1e-9)
    assert result.grade == cns_grade(ratio, CnsParams().ratio_thresholds)


def test_simulate_cns_errors() -> None:
    with pytest.raises(ValueError):
        simulate_cns(np.array([0.0, 1.0]), np.array([0.0, 1.0, 2.0]), mw=100.0)
    with pytest.raises(ValueError):
        simulate_cns(np.array([]), np.array([]), mw=100.0)
    with pytest.raises(ValueError):
        simulate_cns(np.array([0.0]), np.array([1.0]), mw=0.0)


def test_simulate_cns_custom_params() -> None:
    t = np.array([0.0, 1.0])
    free = np.array([0.0, 0.151])
    params = CnsParams(ic50_nm=5.0e3, kpu_brain=1.0)
    result = simulate_cns(t, free, mw=151.0, params=params)
    assert result.peak_brain_free_nm == pytest.approx(1000.0, rel=1e-9)
    assert result.exposure_ratio == pytest.approx(0.2, rel=1e-9)
    assert result.grade == 1


def test_kpu_brain_from_bbb_mapping() -> None:
    assert kpu_brain_from_bbb(None) == 1.0
    assert kpu_brain_from_bbb(0.0) == 0.2
    assert kpu_brain_from_bbb(0.49) == 0.2
    assert kpu_brain_from_bbb(0.5) == 1.0
    assert kpu_brain_from_bbb(0.9) == 1.0
    assert kpu_brain_from_bbb(1.0) == 1.0
    assert kpu_brain_from_bbb(-0.5) == 0.2
    assert kpu_brain_from_bbb(2.0) == 1.0
