"""End-to-end full-chain trust tests (SMILES -> ADMET-AI -> PBPK -> report).

Fast, deterministic legs of the trust mechanism (doc/08 full-chain E2E): the
heavy torch model is never loaded — a structurally-typed fake ADMET predictor
feeds ``predict_admet`` — and the exact production entry used by the CLI and
web playground (``parse_structure`` -> ``predict_admet`` -> ``spec_from_admet``
-> ``run_pipeline`` -> ``to_contract`` -> ``render_all``/``write_report``) is
exercised and asserted for hand-off integrity, determinism, physical
plausibility and fidelity-provenance disclosure.  The heavy real-model leg of
the same chain lives in ``case_full_chain_admet_to_report``.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from drugos.inputs.models import HumanProfile, Route
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.parse_structure import parse_structure
from drugos.pipeline import (
    EmpiricalObservations,
    ExposureProfile,
    Measurements,
    RunSpec,
    fidelity_provenance,
    run_pipeline,
    spec_from_admet,
)
from drugos.pk.admet import AdmetOutput
from drugos.pk.pbpk_build import CypTerm, TargetBinding
from drugos.pk.skin import SkinLayers
from drugos.report import render_all, write_report
from drugos.target.targets import Target

_SUBSTRATES = (
    "CC(=O)Nc1ccc(O)cc1",  # acetaminophen (oral analgesic)
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine (CNS-penetrant stimulant)
)


class _FakePredictor:
    """Protocol-compatible predict_admet stand-in returning AdmetOutput rows."""

    def __init__(self, outputs: list[AdmetOutput]) -> None:
        self._outputs = outputs

    def predict(self, smiles: str | list[str]) -> AdmetOutput | list[AdmetOutput]:
        if isinstance(smiles, str):
            return self._outputs[0]
        return list(self._outputs)


def _admet(smiles: str) -> AdmetOutput:
    return AdmetOutput(
        smiles=smiles,
        fup_plasma=0.10,
        cl_int_hep_ml_min_kg=1.2,
        log_s=-3.5,
        hERG=0.9,
        DILI=0.1,
        AMES=0.5,
        BBB=0.8,
        Pgp=0.3,
        HIA=0.9,
        bioavailable_Ma=0.7,
        caco2_log=1.1,
        log_p_pred=2.15,
        cyp3a4_inhibitor=0.2,
        cyp2d6_inhibitor=0.1,
        cyp2c9_inhibitor=0.3,
        cyp1a2_inhibitor=0.4,
        cyp2c19_inhibitor=0.5,
    )


def _spec(smiles: str, route: Route = Route.ORAL) -> RunSpec:
    out = _admet(smiles)
    molecule = parse_structure(smiles, name=smiles)
    profile = HumanProfile(sex="male", age_y=35.0, height_cm=175.0, weight_kg=70.0)
    dose_plan = build_dose_plan(route=route, amount_mg=10.0, interval_h=0.0, n_doses=1)
    return spec_from_admet(molecule, out, profile, dose_plan)


def test_full_chain_handoff_integrity_and_trust() -> None:
    """SMILES -> ADMET -> PBPK -> organ/clinical -> contract -> report files."""
    spec = _spec(_SUBSTRATES[0])
    result = run_pipeline(spec)
    contract = result.to_contract()

    assert set(contract) == {
        "manifest",
        "pk",
        "occupancy",
        "pathway",
        "organ",
        "clinical",
        "trust",
        "r_verify",
    }
    assert contract["pk"]["cmax_mg_l"] > 0.0
    assert 0.0 < contract["pk"]["bioavailability_f"] <= 1.0
    assert contract["clinical"]["verdict"]

    trust = contract["trust"]
    assert trust["fidelity"] == "full"
    assert "G5" in trust["policy"]
    assert any("cholestasis" in a for a in trust["anchors_wired"])
    assert any("CKD-EPI" in a for a in trust["anchors_wired"])
    assert any("BBB_Martins" in a for a in trust["anchors_wired"])
    assert len(trust["no_public_data_sites"]) >= 1
    assert "raw" not in (trust["admet_ml_estimates"] or [])
    engaged = trust["mechanism_terms_engaged"]
    # Full-fidelity auto-engages the route-applicable ADME/physiology terms.
    assert len(engaged) >= 5
    assert trust["mechanism_terms_degraded"] == []
    estimates = trust["estimates"]
    assert estimates and isinstance(estimates[0], dict)
    assert all({"term", "basis"} <= set(e) for e in estimates)
    # No runtime inventory of unwired seams survives (G5/G6): every term that
    # contributes is engaged in this record, so the record carries no seam key.
    assert not [k for k in trust if "seam" in k]

    # DISCLAIMER §2: every run carries an explicit predictive-regime disclosure.
    reliability = trust["reliability"]
    for key in ("regime", "label", "reliability", "basis", "band_cv", "disclaimer"):
        assert key in reliability
    assert reliability["band_cv"] > 0.0
    assert reliability["regime"] in {
        "novel_molecule",
        "partial_evidence",
        "validated_offlabel_route",
        "validated_extrapolated_dose",
        "validated_in_range_on_label",
        "measured_in_range_on_label",
    }

    payloads = render_all(result)
    assert all(payloads.values())
    re_json = render_all(result)
    assert re_json == payloads


def test_full_chain_measured_override_upgrades_regime() -> None:
    """All-true-parameters enter the same chain and are honoured end to end."""
    spec = _spec(_SUBSTRATES[0])
    novel = run_pipeline(spec).to_contract()["trust"]["reliability"]["regime"]
    assert novel == "novel_molecule"

    measured = spec_from_admet(
        parse_structure(_SUBSTRATES[0], name=_SUBSTRATES[0]),
        _admet(_SUBSTRATES[0]),
        HumanProfile(sex="male", age_y=35.0, height_cm=175.0, weight_kg=70.0),
        build_dose_plan(route=Route.ORAL, amount_mg=10.0, interval_h=0.0, n_doses=1),
        measurements=Measurements(
            fup=0.10,
            cl_hep_l_h=0.31,
            cl_renal_l_h=0.02,
            qt_ic50_nm=1.0e5,
        ),
        empirical=EmpiricalObservations(plasma_cmax_mg_l=0.5),
    )
    result = run_pipeline(measured)
    contract = result.to_contract()
    assert contract["trust"]["reliability"]["regime"] == "measured_in_range_on_label"
    assert contract["trust"]["empirical_agreement"] is not None
    # the run used the measured PK rather than the ADMET-AI prediction
    assert result.spec.cl_hep_l_h == pytest.approx(0.31)
    assert result.spec.fup == pytest.approx(0.10)


def test_full_chain_two_compounds_generic() -> None:
    """The chain is generic over chemically distinct molecules."""
    c1 = run_pipeline(_spec(_SUBSTRATES[0])).to_contract()
    c2 = run_pipeline(_spec(_SUBSTRATES[1])).to_contract()
    assert c1["manifest"]["name"] != c2["manifest"]["name"]
    for c in (c1, c2):
        assert c["pk"]["cmax_mg_l"] > 0.0
        assert 0.0 < c["pk"]["bioavailability_f"] <= 1.0
        assert c["trust"]["no_public_data_sites"]
        assert len(c["trust"]["mechanism_terms_engaged"]) >= 5
        assert c["trust"]["mechanism_terms_degraded"] == []


def test_full_chain_deterministic_repeat() -> None:
    spec = _spec(_SUBSTRATES[0])
    a = run_pipeline(spec).to_contract()
    b = run_pipeline(spec).to_contract()
    assert a == b


def test_full_chain_report_writes_to_disk(tmp_path) -> None:
    result = run_pipeline(_spec(_SUBSTRATES[0]))
    paths = write_report(result, directory=str(tmp_path))
    assert all(Path(p).exists() for p in paths.values())
    md = Path(paths["markdown"]).read_text(encoding="utf-8")
    assert result.name in md
    assert "composite risk" in md
    js = Path(paths["json"]).read_text(encoding="utf-8")
    assert '"trust"' in js
    assert "G5" in js


def test_fidelity_provenance_default_disclosures() -> None:
    spec = _spec(_SUBSTRATES[0])
    exposure = ExposureProfile(
        plasma_cmax_unbound_nm=1.0,
        liver_cmax_free_nm=1.0,
        kidney_cmax_free_nm=1.0,
        qt_ic50_nm=None,
        dili_ic50_nm=None,
        cns_ic50_nm=None,
    )
    prov = fidelity_provenance(spec, exposure)
    # BBB kpu wiring and the hERG sieve are disclosed when the ADMET heads
    # are present, even though CNS *endpoint* grading is not anchored.
    assert prov["fidelity"] == "full"
    assert len(prov["anchors_wired"]) == 4  # CKD-EPI + cholestasis + BBB + hERG
    assert prov["cns_grading_anchored"] is False
    assert prov["admet_ml_estimates"]
    assert prov["mechanism_terms_degraded"] == []
    assert len(prov["mechanism_terms_engaged"]) >= 5
    estimates = prov["estimates"]
    assert isinstance(estimates, list) and len(estimates) >= 1
    assert all({"term", "basis"} <= set(e) for e in estimates)
    assert "corpus-calibrated hERG P->KD sieve (R-8)" in prov["anchors_wired"]
    assert "BBB_Martins" in " ".join(prov["anchors_wired"])


def test_fidelity_provenance_anchors_and_realism_terms() -> None:
    spec = _spec(_SUBSTRATES[0])
    # Base spec already engages MM, secretion, biliary, etc. from auto-estimates.
    # Replace with per-CYP kinetics (clear the lumped MM pair to avoid mutual
    # exclusion in __post_init__) and add all optional realism overrides.
    full = replace(
        spec,
        hepatic_vmax_mg_h=None,
        hepatic_km_mg_l=None,
        cyp_terms=(CypTerm(isoform="CYP3A4", km_mg_l=5.0, vmax_mg_h=1.0),),
        cl_sec_l_h=0.1,
        cl_bil_l_h=0.05,
        gut_extraction_eg=0.3,
        target_binding=(TargetBinding(tissue="liver", target=Target(name="x", kd_nm=2.0)),),
        skin_layers=SkinLayers(),
        dili_immune_ic50_nm=50.0,
        dili_immune_weight=1.0,
        beta_block_ic50_nm=30.0,
        feedback_loop=1,
        sc_im_ka_per_h=0.5,
        estimate_basis={},
    )
    exposure = ExposureProfile(
        plasma_cmax_unbound_nm=1.0,
        liver_cmax_free_nm=2.0,
        kidney_cmax_free_nm=1.0,
        qt_ic50_nm=100.0,
        dili_ic50_nm=200.0,
        cns_ic50_nm=1.0e5,
    )
    anchored = replace(exposure, cns_anchored=True)
    prov = fidelity_provenance(full, anchored)
    assert prov["fidelity"] == "full"
    assert prov["cns_grading_anchored"] is True
    assert any("BBB_Martins" in a for a in prov["anchors_wired"])
    engaged = prov["mechanism_terms_engaged"]
    for term in (
        "tubular secretion",
        "biliary excretion + enterohepatic recirculation",
        "first-pass gut-wall extraction",
        "per-CYP abundance-scaled kinetics",
        "TMDD target binding (native mass balance)",
        "multi-layer transdermal skin permeation",
        "immune-mediated DILI axis",
        "sympathetic-suppression cardiac branch",
        "organ-feedback loop (coupled clearance)",
        "SC/IM depot absorption",
    ):
        assert term in engaged
    assert "saturable (Michaelis-Menten) hepatic clearance" not in engaged
    assert prov["mechanism_terms_degraded"] == []


def test_fidelity_provenance_mixed_confidence_panel() -> None:
    spec = _spec(_SUBSTRATES[0])
    low = Target(name="low", kd_nm=1000.0, low_confidence=True)
    high = Target(name="high", kd_nm=5.0, low_confidence=False)
    spec = replace(spec, panel=(low, high))
    exposure = ExposureProfile(
        plasma_cmax_unbound_nm=1.0,
        liver_cmax_free_nm=1.0,
        kidney_cmax_free_nm=1.0,
        qt_ic50_nm=None,
        dili_ic50_nm=None,
        cns_ic50_nm=None,
    )
    prov = fidelity_provenance(spec, exposure)
    assert prov["fidelity"] == "full"
    assert prov["no_public_data_sites"] == ["low"]


def test_fidelity_provenance_admet_none() -> None:
    spec = replace(_spec(_SUBSTRATES[0]), admet=None)
    prov = fidelity_provenance(
        spec,
        ExposureProfile(
            plasma_cmax_unbound_nm=1.0,
            liver_cmax_free_nm=1.0,
            kidney_cmax_free_nm=1.0,
            qt_ic50_nm=None,
            dili_ic50_nm=None,
            cns_ic50_nm=None,
        ),
    )
    assert prov["fidelity"] == "full"
    assert prov["admet_ml_estimates"] is None
    assert not prov["cns_grading_anchored"]
    assert "BBB_Martins" not in " ".join(prov["anchors_wired"])
    assert len(prov["anchors_wired"]) == 2
    # Auto-anchors still engage all route-applicable terms.
    assert len(prov["mechanism_terms_engaged"]) >= 5
    assert prov["mechanism_terms_degraded"] == []


def test_realism_error_for_incomplete_full_spec() -> None:
    """A manually built full-fidelity spec with missing realism data raises."""
    from drugos.inputs.parse_dosing import build_dose_plan

    mol = parse_structure("CC(=O)NC1=CC=C(O)C=C1", name="tiny")
    profile = HumanProfile(sex="male", age_y=35.0, height_cm=175.0, weight_kg=70.0)
    plan = build_dose_plan(route="oral", amount_mg=250.0, interval_h=0.0, n_doses=1)
    incomplete = RunSpec(
        molecule=mol,
        profile=profile,
        dose_plan=plan,
        cl_hep_l_h=2.0,
        cl_renal_l_h=0.4,
        mw=151.0,
        fup=0.8,
    )
    assert incomplete.fidelity == "full"
    with pytest.raises(Exception, match="full-fidelity run refused"):
        run_pipeline(incomplete)


def test_baseline_optout_bypasses_gate_and_discloses() -> None:
    """An explicitly opted-in baseline run bypasses the gate and discloses."""
    mol = parse_structure("CC(=O)NC1=CC=C(O)C=C1", name="tiny")
    profile = HumanProfile(sex="male", age_y=35.0, height_cm=175.0, weight_kg=70.0)
    plan = build_dose_plan(route="oral", amount_mg=250.0, interval_h=0.0, n_doses=1)
    baseline = RunSpec(
        molecule=mol,
        profile=profile,
        dose_plan=plan,
        cl_hep_l_h=2.0,
        cl_renal_l_h=0.4,
        mw=151.0,
        fup=0.8,
        fidelity="baseline",
    )
    # Full-fidelity terms off → should raise if full; baseline passes.
    result = run_pipeline(baseline)
    contract = result.to_contract()
    assert contract["trust"]["fidelity"] == "baseline"
    assert len(contract["trust"]["mechanism_terms_degraded"]) >= 3
    assert contract["trust"]["estimates"] is None
