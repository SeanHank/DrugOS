"""Input-layer coverage tests: models, dosing schedule, human resolution,
chemical structure parsing."""

import os

import pytest
from pydantic import ValidationError

from drugos.inputs.models import (
    DoseEvent,
    DosePlan,
    HumanProfile,
    IonClass,
    Molecule,
    Route,
    RrClass,
    Sex,
)
from drugos.inputs.parse_dosing import build_dose_plan, dose_plan_to_events
from drugos.inputs.parse_structure import StructureParseError, fraction_neutral, parse_structure
from drugos.inputs.resolve_human import human_profile, resolve_human
from drugos.pk.physiology import HumanPhysiology


# --------------------------------------------------------------------------
# Sex coercion (_missing_)
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw,expected",
    [("m", Sex.MALE), ("M", Sex.MALE), ("male ", Sex.MALE), ("f", Sex.FEMALE), ("FEM", Sex.FEMALE)],
)
def test_sex_missing_fuzzy_matches(raw: str, expected: Sex) -> None:
    assert Sex(raw) is expected


@pytest.mark.parametrize("raw", ["unknown", "x", ""])
def test_sex_missing_raises_for_unknown(raw: str) -> None:
    with pytest.raises(ValueError):
        Sex(raw)


# --------------------------------------------------------------------------
# Molecule descriptor accessors
# --------------------------------------------------------------------------
def test_molecule_pka_accessors() -> None:
    mol = Molecule(pka_acids=[3.0, 5.0], pka_bases=[9.0, 7.0])
    assert mol.strongest_acid_pka == 3.0
    assert mol.strongest_base_pka == 9.0
    assert mol.is_ionizable

    empty = Molecule()
    assert empty.strongest_acid_pka is None
    assert empty.strongest_base_pka is None
    assert not empty.is_ionizable


# --------------------------------------------------------------------------
# HumanProfile validation
# --------------------------------------------------------------------------
def test_human_profile_coerces_sex_string() -> None:
    p = HumanProfile(sex="female")
    assert p.sex is Sex.FEMALE


@pytest.mark.parametrize(
    "field,kwargs",
    [
        ("age", {"age_y": -1.0}),
        ("age", {"age_y": 121.0}),
        ("height", {"height_cm": 0.0}),
        ("weight", {"weight_kg": -5.0}),
        ("hematocrit", {"hematocrit": 0.1}),
        ("albumin", {"albumin_g_l": 200.0}),
        ("agp", {"agp_g_l": 20.0}),
        ("gfr", {"gfr_ml_min": 0.0}),
        ("cardiac_index", {"cardiac_index_l_min_m2": 11.0}),
        ("creatinine", {"serum_creatinine_mg_dl": -1.0}),
    ],
)
def test_human_profile_validation_rejects(field: str, kwargs: dict) -> None:
    with pytest.raises(ValidationError):
        HumanProfile(sex=Sex.MALE, **kwargs)


# --------------------------------------------------------------------------
# DoseEvent / DosePlan
# --------------------------------------------------------------------------
def test_dose_event_rejects_negative_dose() -> None:
    with pytest.raises(ValueError, match="dose_mg"):
        DoseEvent(dose_mg=-1.0)


def test_dose_event_requires_infusion_duration() -> None:
    with pytest.raises(ValueError, match="infusion"):
        DoseEvent(dose_mg=10.0, route=Route.IV_INFUSION)


@pytest.mark.parametrize("duration", [None, -1.0, 0.0])
def test_dose_event_infusion_duration_must_be_positive(duration) -> None:
    with pytest.raises(ValueError, match="infusion"):
        DoseEvent(dose_mg=10.0, route=Route.IV_INFUSION, infusion_duration_h=duration)


def test_dose_event_infusion_ok() -> None:
    DoseEvent(dose_mg=10.0, route=Route.IV_INFUSION, infusion_duration_h=2.0)


def test_dose_plan_factories() -> None:
    assert DosePlan.iv_bolus(5.0).total_dose_mg == 5.0
    assert DosePlan.iv_bolus(5.0).routes == [Route.IV_BOLUS]
    assert DosePlan.iv_infusion(10.0, 2.0).events[0].infusion_duration_h == 2.0
    oral = DosePlan.oral(20.0, food_state="fed")
    assert oral.events[0].route is Route.ORAL
    assert oral.events[0].food_state == "fed"


def test_dose_plan_add_sorts_events() -> None:
    plan = DosePlan.iv_bolus(1.0, time_h=4.0)
    plan.add(DoseEvent(time_h=1.0, dose_mg=2.0, route=Route.ORAL))
    assert [e.time_h for e in plan.events] == [1.0, 4.0]


def test_dose_plan_repeat_validates() -> None:
    with pytest.raises(ValueError, match="repeat"):
        DosePlan.iv_bolus(1.0).repeat(24.0, -1)
    with pytest.raises(ValueError, match="repeat"):
        DosePlan.iv_bolus(1.0).repeat(0.0, 2)


def test_dose_plan_repeat_builds_schedule() -> None:
    plan = DosePlan.oral(100.0).repeat(24.0, 3)
    assert len(plan.events) == 3
    assert [e.time_h for e in plan.events] == [0.0, 24.0, 48.0]


# --------------------------------------------------------------------------
# build_dose_plan / dose_plan_to_events
# --------------------------------------------------------------------------
def test_build_dose_plan_validations() -> None:
    with pytest.raises(ValueError, match="amount_mg"):
        build_dose_plan("iv_bolus", 0.0)
    with pytest.raises(ValueError, match="duration_h"):
        build_dose_plan(Route.IV_INFUSION, 10.0, duration_h=None)
    with pytest.raises(ValueError, match="n_doses"):
        build_dose_plan("oral", 10.0, n_doses=0)


def test_build_dose_plan_repeated_regimen() -> None:
    plan = build_dose_plan("iv_bolus", 15.0, interval_h=12.0, n_doses=4)
    assert len(plan.events) == 4
    events = dose_plan_to_events(plan)
    assert [e.time_h for e in events] == [0.0, 12.0, 24.0, 36.0]


# --------------------------------------------------------------------------
# resolve_human
# --------------------------------------------------------------------------
def test_resolve_human_default_male() -> None:
    phys: HumanPhysiology = resolve_human(HumanProfile(sex=Sex.MALE))
    assert phys.weight_kg == 70.0
    assert phys.hematocrit == 0.45
    assert phys.cardiac_output_l_min > 0


def test_human_profile_factory_with_overrides() -> None:
    prof = human_profile("male", age_y=40.0, weight_kg=80.0, heart_failure=True)
    assert prof.weight_kg == 80.0
    assert prof.heart_failure


# --------------------------------------------------------------------------
# parse_structure
# --------------------------------------------------------------------------
def test_parse_smiles_neutral() -> None:
    mol = parse_structure("c1ccccc1", name="benzene")
    assert mol.ion_class is IonClass.NEUTRAL
    assert mol.rr_class is RrClass.NEUTRAL
    assert mol.canonical_smiles  # canonicalized
    assert mol.inchikey and len(mol.inchikey) == 27
    assert mol.molecular_formula
    assert mol.mw and mol.mw > 0
    assert mol.log_p is not None
    assert mol.hbd == 0
    assert mol.aromatic_rings == 1


def test_parse_smiles_acid() -> None:
    mol = parse_structure("CC(=O)O")
    assert mol.rr_class is RrClass.ACID
    assert mol.strongest_acid_pka == pytest.approx(4.2)
    assert mol.ion_class is IonClass.ANION


def test_parse_smiles_base() -> None:
    mol = parse_structure("CCN")
    assert mol.rr_class is RrClass.STRONG_BASE
    assert mol.strongest_base_pka == pytest.approx(10.3)
    assert mol.ion_class is IonClass.CATION


def test_parse_smiles_zwitterion() -> None:
    mol = parse_structure("NCC(=O)O")
    assert mol.ion_class is IonClass.ZWITTERION
    assert mol.rr_class is RrClass.STRONG_BASE


def test_parse_smiles_weak_base() -> None:
    mol = parse_structure("Nc1ccccc1")
    assert mol.rr_class is RrClass.WEAK_BASE
    assert mol.ion_class is IonClass.NEUTRAL


def test_parse_inchi() -> None:
    mol = parse_structure("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)")
    assert mol.rr_class is RrClass.ACID


def test_parse_sdf_file(tmp_path) -> None:
    from rdkit import Chem

    path = tmp_path / "mol.sdf"
    writer = Chem.SDWriter(os.fspath(path))
    writer.write(Chem.MolFromSmiles("CCO"))
    writer.close()
    mol = parse_structure(os.fspath(path))
    assert mol.molecular_formula == "C2H6O"


def test_parse_sdf_missing_file() -> None:
    with pytest.raises(StructureParseError, match="not found"):
        parse_structure("does_not_exist.sdf")


def test_parse_invalid_smiles() -> None:
    with pytest.raises(StructureParseError, match="Could not parse"):
        parse_structure("not a molecule((((")


def test_fraction_neutral_function() -> None:
    mol = Molecule(pka_acids=[4.2], pka_bases=[])
    assert fraction_neutral(mol, 4.2) == pytest.approx(0.5, rel=0.01)


# --------------------------------------------------------------------------
# heuristic helpers (typed, private API)
# --------------------------------------------------------------------------
def test_macro_pka_ignores_invalid_smarts() -> None:
    from rdkit import Chem

    from drugos.inputs.parse_structure import _macro_pka

    mol = Chem.MolFromSmiles("CC(=O)O")
    values = _macro_pka(mol, table=[("[not a pattern@", 1.0)])
    assert values == []


def test_rr_class_zwitterion() -> None:
    from drugos.inputs.parse_structure import _rr_class

    assert _rr_class([4.2], [6.0]) is RrClass.ZWITTERION
