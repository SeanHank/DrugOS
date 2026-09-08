"""Stage-5 clinical grading and composite toxicity tests (drugos.clinical)."""

from __future__ import annotations

import numpy as np
import pytest

from drugos.clinical.biomarkers import (
    BIOMARKERS,
    BiomarkerSpec,
    crossing_interval,
    grade_absolute,
    grade_timeseries,
    grade_value,
    peak_fraction,
    summarize_cardiac,
    summarize_kidney,
    summarize_liver,
)
from drugos.clinical.toxicity import (
    RULES,
    Endpoint,
    EndpointRisk,
    Evidence,
    EvidenceKind,
    ToxicityReport,
    clamp_prob,
    exposure_evidence,
    fuse,
    logit,
    mechanistic_evidence,
    risk_grade,
    score_toxicity,
    sigmoid,
    structural_evidence,
)
from drugos.organ.cns import CnsParams
from drugos.pk.admet import AdmetOutput

T = np.array([0.0, 1.0, 2.0, 3.0, 4.0])


# ----------------------------------------------------------------- biomarker
def test_spec_validation() -> None:
    with pytest.raises(ValueError):
        BiomarkerSpec("x", "u", 0.0, 1.0, worse="sideways")
    with pytest.raises(ValueError):
        BiomarkerSpec("x", "u", 0.0, 1.0, thresholds=(1.0, 2.0, 3.0))


def test_grade_value_directions() -> None:
    hi = BiomarkerSpec("x", "u", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0))
    lo = BiomarkerSpec("x", "u", 0.0, 1.0, "lower", (70.0, 65.0, 60.0, 50.0))
    assert grade_value(0.5, hi) == 0
    assert grade_value(1.0, hi) == 1
    assert grade_value(50.0, hi) == 4
    assert grade_value(71.0, lo) == 0
    assert grade_value(70.0, lo) == 1
    assert grade_value(20.0, lo) == 4


def test_peak_fraction_directions() -> None:
    hi = BiomarkerSpec("x", "u", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0))
    lo = BiomarkerSpec("x", "u", 0.0, 1.0, "lower", (90.0, 60.0, 45.0, 15.0))
    v, t = peak_fraction(T, np.array([1.0, 3.0, 2.0, 1.0, 0.5]), hi)
    assert (v, t) == (3.0, 1.0)
    v, t = peak_fraction(T, np.array([100.0, 80.0, 90.0, 120.0, 110.0]), lo)
    assert (v, t) == (80.0, 1.0)


def test_crossing_interval() -> None:
    hi = BiomarkerSpec("x", "u", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0))
    lo = BiomarkerSpec("x", "u", 0.0, 1.0, "lower", (90.0, 60.0, 45.0, 15.0))
    onset, dur = crossing_interval(T, np.array([0.5, 1.5, 2.5, 0.9, 0.4]), hi, grade=1)
    assert onset == 1.0
    assert dur == pytest.approx(1.0)
    onset, dur = crossing_interval(T, np.array([0.2, 0.2, 0.2, 1.2, 1.5]), hi)
    assert onset == 3.0
    assert dur == pytest.approx(1.0)
    onset, dur = crossing_interval(T, np.array([95.0, 60.0, 55.0, 92.0, 91.0]), lo, grade=2)
    assert onset == 1.0
    assert dur == pytest.approx(1.0)
    assert crossing_interval(T, np.array([0.0, 0.0, 0.0, 0.0, 0.0]), hi) == (
        None,
        None,
    )


def test_grade_timeseries_full() -> None:
    g = grade_timeseries(T, np.array([0.5, 15.0, 20.0, 20.0, 0.5]), BIOMARKERS["ALT"])
    assert g.grade == 4
    assert g.severity == "life-threatening"
    assert g.value == pytest.approx(20.0)
    assert g.onset_h == 1.0
    assert g.duration_h == pytest.approx(2.0)
    assert g.to_dict()["name"] == "ALT (alanine transaminase)"
    assert g.ref_hi == 1.0


def test_grade_absolute() -> None:
    g = grade_absolute(200.0, BIOMARKERS["heart_rate"])
    assert g.grade == 4
    assert g.onset_h is None
    assert g.duration_h is None
    assert grade_absolute(50.0, BIOMARKERS["map"]).grade == 4
    assert grade_absolute(90.0, BIOMARKERS["map"]).grade == 0


def test_summary_helpers() -> None:
    alt = np.array([1.0, 2.0, 0.5, 0.5, 0.5])
    ast = np.array([0.3, 0.3, 0.3, 0.3, 0.3])
    tbi = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
    liver = summarize_liver(T, alt, ast, tbi)
    assert [b.name for b in liver] == [
        "ALT (alanine transaminase)",
        "AST (aspartate transaminase)",
        "total bilirubin",
    ]
    assert liver[0].grade == 2  # ALT peak 2.0 xULN crosses two thresholds
    qtc = np.array([400.0, 470.0, 470.0, 380.0, 380.0])
    cardiac = summarize_cardiac(T, qtc, heart_rate_bpm=60.0, map_mmhg=90.0)
    assert [b.grade for b in cardiac] == [1, 0, 0]
    gfr = np.array([100.0, 55.0, 50.0, 50.0, 90.0])
    scr = np.array([0.8, 1.0, 1.0, 0.9, 0.9])
    kidney = summarize_kidney(T, gfr, scr)
    assert kidney[0].grade == 2  # gfr <= 60
    assert kidney[1].grade == 0
    kim1 = np.array([1.0, 2.5, 2.5, 1.0, 1.0])
    kidney_k = summarize_kidney(T, gfr, scr, kim1_xunl=kim1)
    assert [b.name for b in kidney_k] == ["GFR", "serum creatinine", "KIM-1 (urinary)"]
    assert kidney_k[2].grade == 1  # peak 2.5 xUNL crosses 1.5


def test_kim1_catalogue_thresholds() -> None:
    assert BIOMARKERS["KIM_1"].thresholds == (1.5, 3.0, 5.0, 10.0)
    assert BIOMARKERS["KIM_1"].unit == "xUNL"


def test_biomarkers_catalogue_has_cns() -> None:
    assert BIOMARKERS["CNS_exposure"].thresholds == CnsParams().ratio_thresholds


# ------------------------------------------------------------------ toxicity
def test_enums() -> None:
    assert Endpoint.DILI.value == "dili"
    assert list(Endpoint) == [
        Endpoint.DILI,
        Endpoint.QT,
        Endpoint.AKI,
        Endpoint.CNS,
    ]
    assert EvidenceKind.EXPOSURE.value == "exposure_ratio"
    assert RULES[Endpoint.AKI].structural_field is None


def test_clamp_logit_sigmoid() -> None:
    assert clamp_prob(0.0) == 1e-6
    assert clamp_prob(1.0) == 1.0 - 1e-6
    assert clamp_prob(0.5) == 0.5
    assert logit(0.5) == 0.0
    assert sigmoid(0.0) == 0.5
    assert sigmoid(logit(0.3)) == pytest.approx(0.3)


def test_mechanistic_evidence_clamps() -> None:
    e = mechanistic_evidence(Endpoint.DILI, grade=-3)
    assert e.probability == RULES[Endpoint.DILI].grade_probs[0]
    e = mechanistic_evidence(Endpoint.QT, grade=99)
    assert e.probability == RULES[Endpoint.QT].grade_probs[4]
    assert EvidenceKind.MECHANISTIC in e.to_dict().values()


def test_exposure_evidence() -> None:
    e = exposure_evidence(Endpoint.QT, cmax_unbound_nm=100.0, ic50_nm=1.0)
    assert e.probability > 0.5
    with pytest.raises(ValueError):
        exposure_evidence(Endpoint.QT, 100.0, 0.0)
    with pytest.raises(ValueError):
        exposure_evidence(Endpoint.QT, -1.0, 10.0)


def test_structural_evidence() -> None:
    admet = AdmetOutput(smiles="c1ccccc1", DILI=0.9, hERG=0.1, BBB=0.5)
    e = structural_evidence(Endpoint.DILI, admet)
    assert e.kind is EvidenceKind.STRUCTURAL
    assert e.probability == 0.9
    with pytest.raises(ValueError):
        structural_evidence(Endpoint.AKI, admet)
    with pytest.raises(ValueError):
        structural_evidence(Endpoint.QT, AdmetOutput(smiles="C", hERG=None))


def test_fuse_and_risk_grade() -> None:
    risk, lo, hi = fuse(0.25, [])
    assert risk == 0.25
    assert 0.0 <= lo <= risk <= hi <= 1.0
    ev = Evidence(
        kind=EvidenceKind.EXPOSURE,
        endpoint=Endpoint.QT,
        probability=0.99,
        weight=1.0,
        note="n",
    )
    risk, _, _ = fuse(0.2, [ev])
    assert risk > 0.5
    risk2, _, _ = fuse(0.2, [ev, mechanistic_evidence(Endpoint.QT, 4)])
    assert risk2 > risk
    assert risk_grade(0.0) == 0
    assert risk_grade(0.36) == 1
    assert risk_grade(0.51) == 2
    assert risk_grade(0.71) == 3
    assert risk_grade(0.91) == 4


def test_score_toxicity_full_path() -> None:
    admet = AdmetOutput(smiles="CC(=O)Nc1ccc(O)cc1", DILI=0.8, hERG=0.02, BBB=0.4)
    report = score_toxicity(
        {Endpoint.DILI: 3, Endpoint.QT: 0, Endpoint.AKI: 1, Endpoint.CNS: 0},
        cmax_unbound_nm=1000.0,
        ic50_nm={
            Endpoint.DILI: 5.0e5,
            Endpoint.QT: 1.0e5,
            Endpoint.AKI: None,
            Endpoint.CNS: 1.0e5,
        },
        admet=admet,
    )
    assert report.by_endpoint(Endpoint.DILI) is not None
    assert report.by_endpoint(Endpoint.CNS) is not None
    dili = report.by_endpoint(Endpoint.DILI)
    assert dili is not None
    assert dili.risk > 0.5
    assert dili.risk >= dili.ci_lo <= dili.ci_hi
    assert dili.driver in {EvidenceKind.MECHANISTIC, EvidenceKind.EXPOSURE}
    assert any(e.kind is EvidenceKind.STRUCTURAL for e in dili.evidence)
    assert report.overall_driver() == "dili"


def test_score_toxicity_empty_evidence() -> None:
    report = score_toxicity({}, cmax_unbound_nm=None, ic50_nm={}, admet=None)
    assert report.overall_risk() == pytest.approx(0.25)  # DILI prior
    assert report.overall_driver() == "dili"
    qt = report.by_endpoint(Endpoint.QT)
    assert qt is not None
    assert qt.driver is EvidenceKind.MECHANISTIC
    assert qt.grade == 0
    assert qt.risk == pytest.approx(0.20)


def test_score_toxicity_unanchored_cns() -> None:
    report = score_toxicity(
        {}, cmax_unbound_nm=None, ic50_nm={}, admet=AdmetOutput(smiles="C", BBB=0.9)
    )
    cns = report.by_endpoint(Endpoint.CNS)
    assert cns is not None
    assert cns.risk > 0.20  # structural BBB line (=0.9) engaged
    assert any(e.kind is EvidenceKind.STRUCTURAL for e in cns.evidence)


def test_score_toxicity_structural_mask_suppresses_cns() -> None:
    # doc/05 5.2: the CNS line is only fused when the molecule is CNS-anchored;
    # the pipeline passes a structural mask for unanchored molecules.
    admet = AdmetOutput(smiles="C", BBB=0.9)
    raw = score_toxicity({}, cmax_unbound_nm=None, ic50_nm={}, admet=admet)
    assert raw.by_endpoint(Endpoint.CNS).risk > 0.20
    masked = score_toxicity(
        {}, cmax_unbound_nm=None, ic50_nm={}, admet=admet, structural_mask={Endpoint.CNS}
    )
    cns = masked.by_endpoint(Endpoint.CNS)
    assert cns is not None
    assert cns.risk == pytest.approx(0.20)  # class prior only
    assert not any(e.kind is EvidenceKind.STRUCTURAL for e in cns.evidence)
    # other endpoints are unaffected by the CNS mask
    assert masked.by_endpoint(Endpoint.DILI).risk == pytest.approx(
        raw.by_endpoint(Endpoint.DILI).risk
    )


def test_toxicity_report_methods() -> None:
    Empty = tuple[EndpointRisk, ...]
    report = ToxicityReport(risks=Empty())
    assert report.by_endpoint(Endpoint.QT) is None
    assert report.overall_risk() == 0.0
    assert report.overall_driver() == "no_evidence"
    assert report.to_dict()["overall_driver"] == "no_evidence"

    ev = Evidence(
        kind=EvidenceKind.EXPOSURE,
        endpoint=Endpoint.QT,
        probability=0.8,
        weight=1.0,
    )
    r = EndpointRisk(
        endpoint=Endpoint.QT,
        label="QT",
        risk=0.9,
        ci_lo=0.1,
        ci_hi=0.99,
        grade=4,
        severity="life-threatening",
        driver=EvidenceKind.EXPOSURE,
        evidence=(ev,),
    )
    single = ToxicityReport(risks=(r,))
    assert single.by_endpoint(Endpoint.QT) is r
    assert single.overall_risk() == 0.9
    assert single.overall_driver() == "qt"
    assert single.to_dict()["endpoints"][0]["driver"] == "exposure_ratio"


def test_endpointrisk_to_dict_clamps() -> None:
    ev = Evidence(
        kind=EvidenceKind.STRUCTURAL,
        endpoint=Endpoint.DILI,
        probability=0.5,
        weight=0.7,
        baseline_probability=0.5,
        note="x",
    )
    r = EndpointRisk(
        endpoint=Endpoint.DILI,
        label="x",
        risk=0.0,
        ci_lo=0.0,
        ci_hi=1.0,
        grade=0,
        severity="normal",
        driver=EvidenceKind.STRUCTURAL,
        evidence=(ev,),
    )
    d = r.to_dict()
    assert d["ci"] == [0.0, 1.0]
    assert d["evidence"][0]["note"] == "x"
