"""Stage-0..5 pipeline driver tests (drugos.pipeline)."""

from __future__ import annotations

from dataclasses import replace

import pytest

import drugos.pipeline as pl
from drugos.clinical.toxicity import Endpoint, EndpointRisk, Evidence, EvidenceKind, ToxicityReport
from drugos.pipeline import (
    RunSpec,
    _BenchmarkLike,
    _other_benchmark_affinities,
    _score_clinical,
    _verdict,
    benchmark_names,
    run_pipeline,
    spec_from_benchmark_data,
)


def _bench(**kw: object) -> _BenchmarkLike:
    base: dict[str, object] = dict(
        name="warfarin",
        smiles="CC(=O)c1ccc(cc1)C(C(=O)O)",
        fup=0.01,
        bp=1.15,
        route="oral",
        dose_mg=5.0,
        tmax_h=6.0,
        n_eval=61,
        log_p=3.1,
        pka_acids=[5.0],
        pka_bases=[],
        published={"cl_plasma_l_h": (0.15, 0.2)},
    )
    base.update(kw)
    return _BenchmarkLike(**base)


def test_benchmark_registry_self_contained() -> None:
    """The built-in registry mirrors validation.benchmarks; the wheel must not
    depend on the repo-local validation tree for CLI/playground runs."""
    names = pl.benchmark_names()
    assert names == ["acetaminophen", "warfarin", "midazolam", "ciprofloxacin", "dofetilide"]
    for name in names:
        data = pl.benchmark_data(name)
        assert data is not None
        assert data.name == name
        spec = pl.spec_from_benchmark_data(data)
        assert spec.name == name
        assert spec.dose_plan is not None
    assert pl.benchmark_data("nope") is None


def test_run_spec_validation() -> None:
    spec = _bench()
    good = spec_from_benchmark_data(spec)
    with pytest.raises(ValueError):
        replace(good, mw=-1.0)
    with pytest.raises(ValueError):
        replace(good, fup=0.0)
    with pytest.raises(ValueError):
        replace(good, fup=1.5)
    with pytest.raises(ValueError):
        replace(good, qt_ic50_nm=0.0)
    with pytest.raises(ValueError):
        replace(good, dili_ic50_nm=-5.0)
    with pytest.raises(ValueError):
        replace(good, cns_ic50_nm=0.0)
    with pytest.raises(ValueError):
        replace(good, feedback_loop=-1)
    with pytest.raises(ValueError):
        replace(good, sc_im_ka_per_h=0.0)


def test_sc_im_ka_override_reaches_model_and_blunts_peak(fast_warfarin: RunSpec) -> None:
    from drugos.inputs.parse_dosing import build_dose_plan

    dose = fast_warfarin.dose_plan.total_dose_mg
    subcut = replace(fast_warfarin, dose_plan=build_dose_plan("subcutaneous", dose))
    slow = replace(subcut, sc_im_ka_per_h=0.02)
    fast = replace(subcut, sc_im_ka_per_h=5.0)
    r_slow = run_pipeline(slow)
    r_fast = run_pipeline(fast)
    # A slow depot release reaches a lower peak and later Tmax than a fast one.
    assert float(r_slow.pk.plasma_total.max()) < float(r_fast.pk.plasma_total.max())
    assert r_fast.pk.pk_metrics().tmax_h < r_slow.pk.pk_metrics().tmax_h


def test_run_pipeline_warfarin(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    assert result.verdict == "No elevated composite risk detected"
    contract = result.to_contract()
    assert contract["pathway"]["available"] is True
    assert contract["pathway"]["readout_peak"] > 0.0
    assert contract["occupancy"]["primary"] > 0.0
    assert contract["clinical"]["toxicity"]["endpoints"][0]["risk"] < 0.3
    assert "cns" in contract["organ"]
    assert "brain_free_nm" in contract["organ"]["trajectories"]
    assert contract["clinical"]["exposure"]["cns_ic50_nm"] == 1.0e5
    assert result.exposure.cns_anchored is False


def test_run_pipeline_apap_od_high_dili(fast_apap_od: RunSpec) -> None:
    result = run_pipeline(fast_apap_od)
    assert "High composite risk" in result.verdict
    assert "dili" in result.verdict
    dili = result.toxicity.by_endpoint(Endpoint.DILI)
    assert dili is not None
    assert dili.risk >= 0.5
    assert dili.ci_lo <= dili.risk <= dili.ci_hi


def test_run_pipeline_cns_anchored(cns_anchored_apap: RunSpec) -> None:
    result = run_pipeline(cns_anchored_apap)
    assert result.exposure.cns_anchored is True
    cns = result.toxicity.by_endpoint(Endpoint.CNS)
    assert cns is not None
    assert cns.risk != 0.20  # mechanistic/exposure CNS lines engaged
    contract = result.to_contract()
    assert contract["clinical"]["exposure"]["cns_ic50_nm"] == 2.0e4


def test_no_pathway(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(replace(fast_warfarin, include_pathway=False))
    assert result.pathway is None
    assert result.primary_signal is not None
    contract = result.to_contract()
    assert contract["pathway"]["available"] is False
    assert contract["pathway"]["readout_peak"] == 0.0


def test_empty_safety_panel(no_safety_panel: RunSpec) -> None:
    result = run_pipeline(no_safety_panel)
    assert result.primary_signal is None
    assert result.pathway is None
    assert result.exposure.qt_ic50_nm is None
    assert result.exposure.dili_ic50_nm is None
    contract = result.to_contract()
    assert contract["occupancy"]["primary"] == 0.0
    assert contract["occupancy"]["ranked"] == []
    assert contract["pathway"]["available"] is False
    assert contract["clinical"]["exposure"]["qt_ic50_nm"] is None


def test_non_herg_panel(non_herg_panel: RunSpec) -> None:
    result = run_pipeline(non_herg_panel)
    assert result.primary_signal is not None
    assert result.pathway is not None
    assert result.exposure.qt_ic50_nm is None
    contract = result.to_contract()
    assert contract["organ"]["cardiac"]["qtc_peak_ms"] < 500.0


def test_feedback_loop_default_and_neutral(fast_warfarin: RunSpec) -> None:
    # Zero loops (default) and a single loop with no organ injury reproduce the
    # same PK: warfarin at this dose does not injure the liver or kidney.
    base = run_pipeline(fast_warfarin)
    looped = run_pipeline(replace(fast_warfarin, feedback_loop=1))
    assert looped.metrics.auc_inf_mgh_l == pytest.approx(base.metrics.auc_inf_mgh_l, rel=1e-6)
    assert looped.metrics.cmax_mg_l == pytest.approx(base.metrics.cmax_mg_l, rel=1e-6)
    assert looped.organ.liver.dili_grade == base.organ.liver.dili_grade
    assert looped.to_contract()["clinical"]["verdict"] == base.to_contract()["clinical"]["verdict"]


def test_feedback_loop_wiring_applies_injury_scaling(
    fast_apap_od: RunSpec, monkeypatch: pytest.MonkeyPatch
) -> None:
    # doc/05 4.5 outer loop: organ injury must produce a hepatic/Kidney PK
    # scaling and that scaling must be passed through apply_pk_scaling.
    captured: dict[str, object] = {}

    def spy(model: object, feedback: object) -> object:
        captured["feedback"] = feedback
        return model

    monkeypatch.setattr(pl, "apply_pk_scaling", spy)
    toxic = replace(
        fast_apap_od,
        feedback_loop=1,
        liver_params=pl.LiverParams(kill_max_1h=0.6, kill_ec50=0.05, kill_hill=2.0),
    )
    result = run_pipeline(toxic)
    assert result.organ.liver.max_dead_frac > 0.0
    feedback = captured["feedback"]
    assert feedback is not None
    assert 0.0 < feedback.hepatic_cl_scale < 1.0


def test_pipeline_cns_structural_line_gated_by_anchoring(fast_warfarin: RunSpec) -> None:
    from drugos.pk.admet import AdmetOutput

    spec = replace(fast_warfarin, admet=AdmetOutput(smiles="CC(=O)c1ccc(cc1)C(C(=O)O)", BBB=0.9))
    cns = run_pipeline(spec).toxicity.by_endpoint(Endpoint.CNS)
    assert cns is not None
    assert cns.risk == pytest.approx(0.20)  # unanchored -> class prior only
    assert not any(e.kind is EvidenceKind.STRUCTURAL for e in cns.evidence)


def test_default_panel_herg_branch(default_panel_no_override: RunSpec) -> None:
    result = run_pipeline(default_panel_no_override)
    assert result.primary_signal is not None
    assert result.exposure.qt_ic50_nm is not None  # read from the hERG panel prior
    assert result.toxicity.by_endpoint(Endpoint.QT) is not None


def test_run_pipeline_tmax_error(fast_warfarin: RunSpec) -> None:
    with pytest.raises(ValueError):
        run_pipeline(replace(fast_warfarin, tmax_h=0.0))


def test_summarize_clinical_helpers(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    liver = result.organ.liver
    cardiac = result.organ.cardiac
    kidney = result.organ.kidney
    cns = result.organ.cns
    assert [b.name for b in pl.summarize_clinical_liver(liver)] == [
        "ALT (alanine transaminase)",
        "AST (aspartate transaminase)",
        "total bilirubin",
    ]
    assert [b.name for b in pl.summarize_clinical_cardiac(cardiac)] == [
        "QTc (Fridericia)",
        "Delta QTc",
        "heart rate",
        "mean arterial pressure",
    ]
    assert [b.name for b in pl.summarize_clinical_kidney(kidney)] == [
        "GFR",
        "serum creatinine",
    ]
    assert [b.name for b in pl.summarize_clinical_cns(cns)] == ["CNS exposure ratio"]
    assert result.organ.grade_by_name("CNS exposure ratio") is not None
    assert result.organ.grade_by_name("does-not-exist") is None


def test_score_clinical_defaults_when_biomarkers_missing() -> None:
    result = run_pipeline(_bench_spec())
    org = replace(result.organ, biomarkers={})
    rep = _score_clinical(org, result.exposure, result.spec.admet)
    qt = rep.by_endpoint(Endpoint.QT)
    assert qt is not None
    assert qt.grade == 0


def _bench_spec() -> RunSpec:
    return spec_from_benchmark_data(_bench())


def test_verdict_thresholds() -> None:
    low = _report(0.1)
    assert _verdict(low) == "No elevated composite risk detected"
    mid = _report(0.4)
    assert _verdict(mid) == "Elevated composite risk: monitor on the flagged endpoint(s)"
    high = _report(0.9)
    assert _verdict(high) == "High composite risk (90%, driver dili)"


def test_benchmark_names() -> None:
    assert benchmark_names() == [
        "acetaminophen",
        "warfarin",
        "midazolam",
        "ciprofloxacin",
        "dofetilide",
    ]


def test_affinities() -> None:
    assert _other_benchmark_affinities("dofetilide") == (2.0, None)
    assert _other_benchmark_affinities("warfarin") == (1e6, None)
    assert _other_benchmark_affinities("acetaminophen") == (1e6, 5.0e5)
    assert _other_benchmark_affinities("unknown-compound") == (None, None)


def test_spec_from_benchmark_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    bad = _bench(name="boom", smiles="CC")
    monkeypatch.setattr(
        "drugos.pipeline.parse_structure",
        lambda smiles, name: type("Mol", (), {"mw": 0.0})(),
    )
    with pytest.raises(ValueError):
        spec_from_benchmark_data(bad)


def test_benchmark_plan_override(fast_warfarin: RunSpec) -> None:
    spec = spec_from_benchmark_data(_bench(), dose_override_mg=25.0)
    assert spec.dose_plan.events[0].dose_mg == 25.0
    spec2 = spec_from_benchmark_data(_bench())
    assert spec2.dose_plan.events[0].dose_mg == 5.0


def test_organ_stage_grade_by_name(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    g = result.organ.grade_by_name("QTc (Fridericia)")
    assert g is not None and g.name == "QTc (Fridericia)"
    assert result.organ.grade_by_name("nonexistent") is None


def _report(risk: float) -> ToxicityReport:
    ev = Evidence(
        kind=EvidenceKind.EXPOSURE,
        endpoint=Endpoint.DILI,
        probability=0.99,
        weight=1.0,
    )
    r = EndpointRisk(
        endpoint=Endpoint.DILI,
        label="Drug-induced liver injury",
        risk=risk,
        ci_lo=max(0.0, risk - 0.1),
        ci_hi=min(1.0, risk + 0.1),
        grade=4,
        severity="life-threatening",
        driver=EvidenceKind.EXPOSURE,
        evidence=(ev,),
    )
    return ToxicityReport(risks=(r,))
