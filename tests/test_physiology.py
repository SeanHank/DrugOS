"""Physiology-layer coverage tests: HumanPhysiology accessors, organ scaling,
impairment modifiers and the circulation-closure invariant."""

import math

import pytest

from drugos.inputs.models import HumanProfile, Sex
from drugos.pk.physiology import (
    PORTAL_DRAINING_TISSUES,
    build_human,
    glom_filtration_clearance,
    mosteller_bsa,
)

MALE = Sex.MALE
FEMALE = Sex.FEMALE


def _phys(**kwargs) -> object:
    prof = HumanProfile(sex=kwargs.pop("sex", MALE), **kwargs)
    return build_human(prof)


def test_accessors() -> None:
    phys = _phys()
    assert math.isclose(phys.cardiac_output_ml_min, phys.cardiac_output_l_min * 1000.0)
    assert math.isclose(phys.gfr_ml_min, phys.gfr_l_min * 1000.0)
    assert sorted(phys.organ_names) == sorted(phys.organ_volume)
    comp = phys.tissue_composition()
    assert comp["brain"].water_fraction == pytest.approx(0.162 + 0.620)
    art = phys.liver_flow_art_l_min
    portal = phys.hepatic_portal_flow_l_min
    assert portal > 0
    assert art > 0
    assert phys.organ_flow["liver"] == pytest.approx(portal + art)
    d = phys.to_dict()
    assert d["weight_kg"] == 70.0


def test_circulation_closed_exactly() -> None:
    # Cardiac output equals vena-cava return; the lung receives it.  Sum of
    # venous returns must equal the arterial demand exactly.
    phys = _phys()
    assert phys.cardiac_output_l_min == pytest.approx(
        sum(
            f
            for t, f in phys.organ_flow.items()
            if t != "lung" and t not in PORTAL_DRAINING_TISSUES
        )
    )
    assert phys.organ_flow["lung"] == pytest.approx(phys.cardiac_output_l_min)


def test_mosteller_bsa() -> None:
    # (170*70/3600)^0.5 = 1.818
    assert mosteller_bsa(170.0, 70.0) == pytest.approx(1.8185, rel=1e-3)


def test_overrides_take_effect() -> None:
    phys = _phys(
        hematocrit=0.40,
        albumin_g_l=38.0,
        agp_g_l=1.1,
        gfr_ml_min=90.0,
        cardiac_index_l_min_m2=3.5,
    )
    assert phys.hematocrit == 0.40
    assert phys.albumin_g_l == 38.0
    assert phys.agp_g_l == 1.1
    assert math.isclose(phys.gfr_ml_min, 90.0)
    # Cardiac output is the circulation-closed vena-cava return: within the
    # tissue-fraction closure factor (~0.998) of CI*BSA.
    assert math.isclose(phys.cardiac_output_l_min, 3.5 * mosteller_bsa(170.0, 70.0), rel_tol=1e-2)


def test_female_defaults() -> None:
    phys = build_human(HumanProfile(sex=FEMALE))
    assert phys.hematocrit == 0.40
    assert phys.albumin_g_l == 43.0
    assert phys.agp_g_l == 0.7


def test_age_adjustments() -> None:
    phys60 = _phys(age_y=60.0)
    assert phys60.albumin_g_l == 41.0
    assert phys60.agp_g_l == 0.9
    phys75 = _phys(age_y=75.0)
    assert phys75.albumin_g_l == 41.0
    assert phys75.gfr_ml_min == pytest.approx(125.0 * 0.90)


def test_renal_impairment() -> None:
    mild = _phys(mild_renal_impairment=True)
    mod = _phys(moderate_renal_impairment=True)
    base = _phys()
    assert mild.gfr_l_min < base.gfr_l_min
    assert mod.gfr_l_min < mild.gfr_l_min


def test_hepatic_impairment_and_heart_failure() -> None:
    base = _phys()
    mild_liver = _phys(mild_hepatic_impairment=True)
    mod_liver = _phys(moderate_hepatic_impairment=True)
    hf = _phys(heart_failure=True)
    assert mild_liver.organ_volume["liver"] == pytest.approx(base.organ_volume["liver"] * 0.80)
    assert mod_liver.organ_volume["liver"] == pytest.approx(base.organ_volume["liver"] * 0.65)
    assert mod_liver.organ_flow["liver"] == pytest.approx(base.organ_flow["liver"] * 0.85)
    assert hf.cardiac_output_l_min < base.cardiac_output_l_min


def test_scaling_tracks_weight() -> None:
    heavy = _phys(weight_kg=90.0)
    assert heavy.organ_volume["liver"] > _phys().organ_volume["liver"]


def test_glom_filtration_clearance() -> None:
    # 125 mL/min, fup 0.4, pure filtration -> 0.4*0.125 L/min*60 = 3.0 L/h.
    assert glom_filtration_clearance(125.0, 0.4) == pytest.approx(3.0)
    # fe folds in net secretion/reabsorption.
    assert glom_filtration_clearance(125.0, 0.4, fe_unchanged=0.5) == pytest.approx(1.5)
    assert glom_filtration_clearance(125.0, 0.4, fe_unchanged=0.0) == pytest.approx(0.0)
    with pytest.raises(ValueError):
        glom_filtration_clearance(0.0, 0.4)
    with pytest.raises(ValueError):
        glom_filtration_clearance(125.0, 0.0)
    with pytest.raises(ValueError):
        glom_filtration_clearance(125.0, 1.1)
    with pytest.raises(ValueError):
        glom_filtration_clearance(125.0, 0.4, fe_unchanged=-0.1)
    with pytest.raises(ValueError):
        glom_filtration_clearance(125.0, 0.4, fe_unchanged=1.1)
