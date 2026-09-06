"""Tests for Rodgers-Rowland partition coefficients (pk/partitions.py).

Golden references for strong bases were produced by an independent
transcription of the Metrum Research PBPK_PC ``CalcKp_R&R.R`` implementation,
which for strong bases agrees exactly with the published Rodgers & Rowland
equations (R&R 2005; 2007 erratum; Korzekwa et al. DMD 2019).
"""

import math

import pytest

from drugos.pk.partitions import ka_ap_strong_base, kpu_blood_cell, rodgers_rowland_partition

TOL = 1.2e-2

METOPROLOL_REF_KP = {
    "adipose": 1.9942,
    "bone_rest": 3.8694,
    "brain": 3.3281,
    "gut": 12.1429,
    "heart": 11.4008,
    "kidney": 23.9053,
    "liver": 21.8954,
    "lung": 18.8533,
    "muscle": 8.3694,
    "skin": 6.9426,
    "spleen": 15.7450,
}

CAFFEINE_REF_KP = {
    "adipose": 0.5648,
    "bone_rest": 1.4029,
    "brain": 1.6140,
    "gut": 3.6777,
    "heart": 3.4936,
    "kidney": 6.5918,
    "liver": 6.1480,
    "lung": 5.3289,
    "muscle": 2.8547,
    "skin": 2.2220,
    "spleen": 4.6586,
}


def test_metoprolol_strong_base() -> None:
    p = rodgers_rowland_partition(fup=0.879, log_p=2.15, pka_bases=[9.7], bp=1.52, hematocrit=0.45)
    assert math.isclose(p.kpu_bc, 2.4523, abs_tol=TOL)
    assert math.isclose(p.ka_ap, 2.0442, abs_tol=TOL)
    for tissue, ref in METOPROLOL_REF_KP.items():
        assert math.isclose(p.kp[tissue], ref, abs_tol=TOL), (tissue, p.kp[tissue], ref)


def test_caffeine_strong_base() -> None:
    p = rodgers_rowland_partition(
        fup=0.681, log_p=-0.07, pka_bases=[10.4], bp=0.98, hematocrit=0.45
    )
    assert math.isclose(p.kpu_bc, 1.4032, abs_tol=TOL)
    assert math.isclose(p.ka_ap, 0.6492, abs_tol=TOL)
    for tissue, ref in CAFFEINE_REF_KP.items():
        assert math.isclose(p.kp[tissue], ref, abs_tol=TOL), (tissue, p.kp[tissue], ref)


def test_ka_ap_monotonic() -> None:
    # Qualitative sanity: Ka_AP must rise with lipophilicity and fall as fup
    # rises (the exact reference numbers of Assmus 2017 could not be reused:
    # that report is not open-access and its Ka_AP unit convention differs).
    lipo = ka_ap_strong_base(4.80, 9.5, 0.05, 1.06, 0.45)
    polar = ka_ap_strong_base(3.65, 9.5, 0.05, 1.00, 0.45)
    assert lipo > polar > 0
    bound = ka_ap_strong_base(4.80, 9.5, 0.05, 1.06, 0.45)
    less_bound = ka_ap_strong_base(4.80, 9.5, 0.20, 1.06, 0.45)
    assert bound > less_bound


def test_strong_base_requires_bp() -> None:
    with pytest.raises(ValueError):
        rodgers_rowland_partition(fup=0.5, log_p=2.0, pka_bases=[9.0])


def test_kpu_blood_cell_formula() -> None:
    # (BP + H - 1) / (H * fup), H = 0.45
    assert math.isclose(kpu_blood_cell(0.6, 0.05, 0.45), 2.2222, abs_tol=1e-4)


def test_neutral_kpu_increases_with_binding() -> None:
    lo = rodgers_rowland_partition(fup=1.0, log_p=1.5).kpu["brain"]
    hi = rodgers_rowland_partition(fup=0.2, log_p=1.5).kpu["brain"]
    assert hi > lo
    # Neutral case: kp == kpu (fup=1).
    neutral = rodgers_rowland_partition(fup=1.0, log_p=1.5)
    assert math.isclose(lo, neutral.kp["brain"], rel_tol=1e-12)


def test_acid_kpu_corner_cases() -> None:
    # Stronger base branch disabled; a weak acid with high pKa behaves like neutral.
    weak = rodgers_rowland_partition(fup=0.5, log_p=2.0, pka_acids=[10.0])
    strong = rodgers_rowland_partition(fup=0.5, log_p=2.0, pka_acids=[3.0])
    assert strong.kpu["liver"] < weak.kpu["liver"]  # ionized acids partition less
