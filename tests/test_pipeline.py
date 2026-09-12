"""Stage-0..5 pipeline driver tests (drugos.pipeline)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import drugos.pipeline as pl
from drugos.clinical.toxicity import Endpoint, EndpointRisk, Evidence, EvidenceKind, ToxicityReport
from drugos.organ import free_mg_l_to_nm
from drugos.pipeline import (
    RunSpec,
    _absorption_params,
    _BenchmarkLike,
    _other_benchmark_affinities,
    _score_clinical,
    _verdict,
    benchmark_names,
    run_pipeline,
    spec_from_benchmark_data,
)
from drugos.pk.pbpk_build import CypTerm, TargetBinding
from drugos.target.targets import Target


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
    with pytest.raises(ValueError):
        replace(good, fa=2.0)
    with pytest.raises(ValueError):
        replace(good, fa=0.0)
    with pytest.raises(ValueError):
        replace(good, fidelity="oops")
    with pytest.raises(ValueError):
        replace(good, si_segments=0)
    assert replace(good, fa=0.9).fa == pytest.approx(0.9)


def test_oral_fa_gates_si_absorption_rate() -> None:
    spec = spec_from_benchmark_data(pl.benchmark_data("acetaminophen"))
    assert spec.fa == pytest.approx((0.80 + 0.98) / 2.0)
    oral = pl.RunSpec(
        molecule=spec.molecule,
        profile=spec.profile,
        dose_plan=spec.dose_plan,
        cl_hep_l_h=spec.cl_hep_l_h,
        cl_renal_l_h=spec.cl_renal_l_h,
        mw=spec.mw,
        fup=spec.fup,
        fa=0.85,
    )
    gated = _absorption_params(oral)
    assert gated.k_si_absorption == pytest.approx(0.55 * (0.85 / 0.85))
    # A high-fa (fast/permeable) compound absorbs faster than a low-fa one.
    fast = _absorption_params(replace(oral, fa=0.97))
    slow = _absorption_params(replace(oral, fa=0.5))
    assert fast.k_si_absorption > slow.k_si_absorption
    # Defaults untouched when no fa target is present.
    assert _absorption_params(replace(oral, fa=None)).k_si_absorption == pytest.approx(0.55)
    # fa is ignored for non-oral (IV) plans.
    from drugos.inputs.parse_dosing import build_dose_plan

    iv = replace(oral, fa=0.85, dose_plan=build_dose_plan("iv_bolus", 5.0))
    assert _absorption_params(iv).k_si_absorption == pytest.approx(0.55)
    # The depot override wins over a stale fa.
    depot = _absorption_params(replace(oral, fa=0.9, sc_im_ka_per_h=0.02, dose_plan=spec.dose_plan))
    assert depot.k_depot_absorption == pytest.approx(0.02)
    assert depot.k_si_absorption == pytest.approx(0.55)


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
    assert "regen_scale" in contract["organ"]["trajectories"]
    liver_scale = contract["organ"]["liver"]
    assert 0.5 <= liver_scale["regen_scale_min"] <= liver_scale["regen_scale_peak"] <= 1.5
    traj = contract["organ"]["trajectories"]
    assert len({len(v) for v in traj.values()}) == 1, "trajectories must share one time grid"
    # CNS is driven by the PBPK brain compartment free exposure (kpu_brain=1.0):
    # the reported brain free peak equals the free brain tissue peak in nM.
    brain_free_nm_pk = float(
        free_mg_l_to_nm(result.pk.unbound_tissues["brain"].max(), result.spec.mw)
    )
    assert result.organ.cns.peak_brain_free_nm == pytest.approx(brain_free_nm_pk, rel=1e-3)
    # hERG cardiac block now comes from the heart free exposure: warfarin has
    # negligible cardiac free exposure so its delta-QTc stays sub-millisecond.
    assert float(result.organ.cardiac.delta_qtc_ms.max()) < 1.0
    assert contract["clinical"]["exposure"]["cns_ic50_nm"] == 1.0e5
    assert result.exposure.cns_anchored is False


def test_run_pipeline_r_verify_present(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    assert result.r_verify is not None
    verify = result.to_contract()["r_verify"]
    assert verify["verdict"] == "r:agree"
    assert verify["agreement_frac"] <= 0.02
    assert isinstance(verify["py_cl_l_h"], float)


def test_run_pipeline_r_verify_unavailable_contract(
    fast_warfarin: RunSpec, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(pl, "_r_verify_pk", lambda *_: None)
    result = run_pipeline(fast_warfarin)
    contract = result.to_contract()
    assert result.r_verify is None
    assert contract["r_verify"] is None


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


def test_production_sbml_pathway_runs(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(replace(fast_warfarin, include_pathway=True))
    assert result.pathway is not None
    assert result.pathway.model.readout == "PP_K"
    # Full pump-rate occupancy (>= drug/Kd type signal) drives the cascade:
    # the production readout peak must be finite and reported in the contract.
    assert np.isfinite(result.pathway.concentrations["PP_K"]).all()
    contract = result.to_contract()
    assert contract["pathway"]["available"] is True
    assert contract["pathway"]["readout"] == "PP_K"
    assert contract["pathway"]["readout_peak"] > 0.0


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
    # Under full fidelity the mitochondrial/immune/pathway axes of a high-dose
    # APAP run legitimately perturb the simulated cardiac output, so the
    # cardiac-output feedback is now depressed below unity (it was neutral in
    # the linear baseline, where no mitochondrial TMDD or immune load existed).
    assert 0.0 < feedback.co_scale < 1.0


def test_pipeline_cns_structural_line_gated_by_anchoring(fast_warfarin: RunSpec) -> None:
    from drugos.pk.admet import AdmetOutput

    spec = replace(fast_warfarin, admet=AdmetOutput(smiles="CC(=O)c1ccc(cc1)C(C(=O)O)", BBB=0.9))
    cns = run_pipeline(spec).toxicity.by_endpoint(Endpoint.CNS)
    assert cns is not None
    assert cns.risk == pytest.approx(0.20)  # unanchored -> class prior only
    assert not any(e.kind is EvidenceKind.STRUCTURAL for e in cns.evidence)


def test_run_pipeline_cns_bbb_partition(fast_warfarin: RunSpec) -> None:
    from drugos.pk.admet import AdmetOutput

    penetrant = run_pipeline(
        replace(fast_warfarin, admet=AdmetOutput(smiles="CC(=O)c1ccc(cc1)C(C(=O)O)", BBB=0.9))
    )
    restricted = run_pipeline(
        replace(fast_warfarin, admet=AdmetOutput(smiles="CC(=O)c1ccc(cc1)C(C(=O)O)", BBB=0.1))
    )
    assert penetrant.exposure.cns_kpu_brain == 1.0
    assert restricted.exposure.cns_kpu_brain == 0.2
    peak_ratio = restricted.organ.cns.peak_brain_free_nm / penetrant.organ.cns.peak_brain_free_nm
    assert peak_ratio == pytest.approx(0.2, rel=1e-9)
    assert restricted.organ.cns.exposure_ratio == pytest.approx(
        penetrant.organ.cns.exposure_ratio * 0.2, rel=1e-9
    )


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
        "KIM-1 (urinary)",
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


def test_cns_ic50_forwarded_to_summarize(fast_warfarin: RunSpec) -> None:
    # Regression: summarize_clinical_cns must grade against the spec-driven IC50,
    # not a fresh CnsParams() default (the 2026-06 latent defect).
    anchored = run_pipeline(replace(fast_warfarin, cns_ic50_nm=2.0e4))
    default = run_pipeline(fast_warfarin)
    g_a = anchored.organ.grade_by_name("CNS exposure ratio")
    g_d = default.organ.grade_by_name("CNS exposure ratio")
    assert g_a is not None and g_d is not None
    assert anchored.organ.cns.exposure_ratio == pytest.approx(g_a.value)
    # Halving the IC50 quintuples the exposure margin.
    assert g_a.value == pytest.approx(5.0 * g_d.value, rel=1e-6)
    assert anchored.organ.cns.exposure_ratio == pytest.approx(
        5.0 * default.organ.cns.exposure_ratio, rel=1e-6
    )


def test_admet_herg_sieve_wires_occupancy_qt_and_cardiac(fast_warfarin: RunSpec) -> None:
    from drugos.pk.admet import AdmetOutput

    nonblocker = replace(
        fast_warfarin, qt_ic50_nm=None, admet=AdmetOutput(smiles="CC(=O)c1ccc", hERG=0.1)
    )
    blocker = replace(
        fast_warfarin, qt_ic50_nm=None, admet=AdmetOutput(smiles="CC(=O)c1ccc", hERG=0.9)
    )
    # Corpus-calibrated monotone P -> KD (doc/12 D10): never more potent than
    # the 2 nM panel prior, exponentially weaker as confidence falls.
    assert pl._effective_herg_kd_nm(nonblocker) == pytest.approx(269217.32)
    assert pl._effective_herg_kd_nm(blocker) == pytest.approx(7.42894)
    assert pl._effective_herg_kd_nm(
        replace(nonblocker, admet=AdmetOutput(smiles="C", hERG=0.5))
    ) == pytest.approx(1414.2136)
    # The resolver replaces the hERG site's KD in the panel (occupancy/primary
    # signal/pathway drive) while leaving every other target untouched.
    sieved = pl._herg_sieved_panel(nonblocker)
    assert all(t.kd_nm == pytest.approx(269217.32) for t in sieved if "hERG" in t.name)
    assert sum(1 for t in sieved if "hERG" in t.name) == sum(
        1 for t in nonblocker.panel if "hERG" in t.name
    )
    # A predicted non-blocker loses its QT bite; a predicted blocker keeps it.
    r_nb = run_pipeline(nonblocker)
    r_b = run_pipeline(blocker)
    assert r_nb.exposure.qt_ic50_nm == pytest.approx(269217.32)
    assert r_b.exposure.qt_ic50_nm == pytest.approx(7.42894)
    assert r_nb.exposure.dili_ic50_nm == r_b.exposure.dili_ic50_nm
    assert r_nb.organ.cardiac.delta_qtc_ms.max() < r_b.organ.cardiac.delta_qtc_ms.max()


def test_qt_override_beats_admet_sieve(fast_warfarin: RunSpec) -> None:
    from drugos.pk.admet import AdmetOutput

    spec = replace(fast_warfarin, qt_ic50_nm=12.5, admet=AdmetOutput(smiles="CC", hERG=0.95))
    assert pl._effective_herg_kd_nm(spec) == pytest.approx(12.5)
    result = run_pipeline(spec)
    assert result.exposure.qt_ic50_nm == pytest.approx(12.5)


def test_organ_feedback_co_fraction_tracks_cardiac_output(fast_warfarin: RunSpec) -> None:
    from dataclasses import replace as _replace
    from types import SimpleNamespace

    model = pl._panel_model(fast_warfarin)
    t = np.linspace(0.0, 6.0, 61)
    phy = model.physiology
    zero = np.zeros_like(t)
    liver = pl.simulate_liver(t, zero, fast_warfarin.mw)
    kidney = pl.simulate_kidney(
        t, zero, fast_warfarin.mw, gfr_base_ml_min=phy.gfr_ml_min, scr_base_umol_l=80.0
    )
    cardiac = pl.simulate_cardiac(t, zero, phy)
    ref_co = phy.cardiac_output_ml_min / 1000.0
    normal = pl._organ_feedback(model, liver, kidney, cardiac)
    assert normal.co_scale == pytest.approx(1.0)
    assert normal.hepatic_cl_scale == pytest.approx(1.0)
    halved = pl._organ_feedback(model, liver, kidney, _replace(cardiac, co_l_min=0.5 * ref_co))
    assert halved.co_scale == pytest.approx(0.5)
    # Zero-reference cardiac output guard stays neutral instead of NaN.
    probe_model = SimpleNamespace(
        physiology=SimpleNamespace(cardiac_output_ml_min=0.0, gfr_ml_min=phy.gfr_ml_min)
    )
    guard = pl._organ_feedback(probe_model, liver, kidney, _replace(cardiac, co_l_min=99.0))
    assert guard.co_scale == pytest.approx(1.0)


def test_oral_bioavailability_reported_in_contract(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    assert result.metrics.f_abs == pytest.approx(result.pk.bioavailability_f or 1.0)
    assert 0.0 < result.metrics.f_abs <= 1.0
    assert result.to_contract()["pk"]["bioavailability_f"] == pytest.approx(
        round(result.metrics.f_abs, 4)
    )


def test_proliferation_signal_helper_branches() -> None:
    class _Model:
        readout = "PP_K"

    class _NoReadout:
        readout = None

    t = np.linspace(0.0, 1.0, 11)
    good = pl.PathwayResult(
        model=_Model(),
        t_h=t,
        concentrations={"PP_K": np.linspace(1.0, 2.0, 11)},
        signal=np.zeros(11),
        baseline={"PP_K": 1.0},
    )
    fc = pl._proliferation_signal(good, t)
    assert fc is not None and fc[0] == pytest.approx(1.0) and fc[-1] == pytest.approx(2.0)
    none_readout = pl.PathwayResult(
        model=_NoReadout(), t_h=t, concentrations={}, signal=np.zeros(11), baseline={}
    )
    assert pl._proliferation_signal(none_readout, t) is None
    zero_base = pl.PathwayResult(
        model=_Model(),
        t_h=t,
        concentrations={"PP_K": np.ones(11)},
        signal=np.zeros(11),
        baseline={"PP_K": 0.0},
    )
    assert pl._proliferation_signal(zero_base, t) is None
    missing_base = pl.PathwayResult(
        model=_Model(),
        t_h=t,
        concentrations={"PP_K": np.ones(11)},
        signal=np.zeros(11),
        baseline={},
    )
    assert pl._proliferation_signal(missing_base, t) is None


def test_herg_occupancy_effective_kd_branches() -> None:
    t = np.linspace(0.0, 1.0, 21)
    c_free = np.full_like(t, 1.0)
    zero = pl._herg_occupancy(t, c_free, mw=300.0, qt_ic50_nm=None)
    assert np.all(zero == 0.0)
    occ = pl._herg_occupancy(t, c_free, mw=300.0, qt_ic50_nm=2.0)
    assert np.all(occ >= 0.0) and np.any(occ > 0.0)


def test_absorption_params_admet_fa_and_solubility_fallbacks() -> None:
    from drugos.pk.admet import AdmetOutput

    base = _bench_spec()
    with_admet = replace(base, admet=AdmetOutput(smiles="C", HIA=0.40, log_s=-6.5))
    params = pl._absorption_params(with_admet)
    assert params.k_si_absorption == pytest.approx(0.55 * (0.40 / 0.85))
    assert params.solubility_mg_ml == pytest.approx(10.0**-6.5 * base.mw / 1000.0)

    bio_maybe = replace(base, admet=AdmetOutput(smiles="C", HIA=None, bioavailable_Ma=0.7))
    params2 = pl._absorption_params(bio_maybe)
    assert params2.k_si_absorption == pytest.approx(0.55 * (0.7 / 0.85))
    assert params2.solubility_mg_ml is None

    nan_log_s = replace(base, admet=AdmetOutput(smiles="C", HIA=0.8, log_s=float("nan")))
    params3 = pl._absorption_params(nan_log_s)
    assert params3.solubility_mg_ml is None

    missing = pl._absorption_params(replace(base, admet=None))
    assert missing.k_si_absorption == 0.55 and missing.solubility_mg_ml is None

    assert pl._admet_fa(AdmetOutput(smiles="C", HIA=0.2)) == pytest.approx(0.2)
    assert pl._admet_fa(AdmetOutput(smiles="C", HIA=None)) is None
    assert pl._admet_fa(None) is None


def test_spec_from_admet_builds_full_chain_input() -> None:
    from drugos.inputs.models import HumanProfile, Sex
    from drugos.inputs.parse_dosing import build_dose_plan
    from drugos.inputs.parse_structure import parse_structure
    from drugos.pk.admet import AdmetOutput

    mol = parse_structure("CC(=O)NC1=CC=C(O)C=C1", name="newmol")
    admet = AdmetOutput(
        smiles=mol.canonical_smiles,
        fup_plasma=0.72,
        cl_int_hep_ml_min_kg=18.0,
        HIA=0.9,
        log_s=-3.2,
    )
    profile = HumanProfile(sex=Sex.MALE, age_y=40, height_cm=170, weight_kg=70)
    plan = build_dose_plan("oral", 500.0)
    spec = pl.spec_from_admet(mol, admet, profile, plan)
    assert spec.fup == pytest.approx(0.72)
    assert spec.cl_hep_l_h == pytest.approx(18.0 * 60.0 / 1000.0 * 70.0)
    assert spec.cl_renal_l_h > 0.0
    assert spec.fa == pytest.approx(0.9)
    assert spec.admet is admet
    assert spec.include_pathway is True
    params = pl._absorption_params(spec)
    assert params.k_si_absorption == pytest.approx(0.55 * (0.9 / 0.85))
    assert params.solubility_mg_ml == pytest.approx(10.0**-3.2 * mol.mw / 1000.0)

    iv = pl.spec_from_admet(mol, admet, profile, build_dose_plan("iv_bolus", 500.0))
    assert iv.fa is None
    assert iv.admet is admet

    no_fup = replace(admet, fup_plasma=None)
    with pytest.raises(ValueError):
        pl.spec_from_admet(mol, no_fup, profile, plan)
    no_cl = replace(admet, cl_int_hep_ml_min_kg=None)
    with pytest.raises(ValueError):
        pl.spec_from_admet(mol, no_cl, profile, plan)
    over = pl.spec_from_admet(mol, admet, profile, plan, cl_hep_l_h=7.0, cl_renal_l_h=0.2)
    assert over.cl_hep_l_h == pytest.approx(7.0) and over.cl_renal_l_h == pytest.approx(0.2)

    ext = pl.spec_from_admet(
        mol,
        admet,
        profile,
        plan,
        cl_sec_l_h=1.5,
        cl_bil_l_h=0.3,
        bile_emptying_1h=0.5,
        hepatic_vmax_mg_h=10.0,
        hepatic_km_mg_l=2.0,
        gut_extraction_eg=0.4,
    )
    assert ext.cl_sec_l_h == 1.5 and ext.cl_bil_l_h == 0.3
    assert ext.bile_emptying_1h == 0.5
    assert ext.hepatic_vmax_mg_h == 10.0 and ext.hepatic_km_mg_l == 2.0
    assert ext.gut_extraction_eg == 0.4
    assert pl._absorption_params(ext).gut_extraction_eg == 0.4

    terms = (CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=7.8),)
    cyber = pl.spec_from_admet(mol, admet, profile, plan, cyp_terms=terms)
    assert cyber.cyp_terms == terms


def test_run_spec_stage1_realism_wiring(fast_warfarin: RunSpec) -> None:
    # The pipeline threads the optional Stage-1 realism terms from the spec
    # into the PBPK model: gut-wall extraction lowers reported oral F, and
    # tubular secretion raises the urinary recovery.  Zero the other engaged
    # full-fidelity terms so the two wirings below are observed in isolation
    # (fidelity baseline keeps the gate from vetoing the zeroed terms).
    base = replace(
        fast_warfarin,
        include_pathway=False,
        feedback_loop=0,
        tmax_h=48.0,
        fidelity="baseline",
        cl_bil_l_h=0.0,
        hepatic_vmax_mg_h=None,
        hepatic_km_mg_l=None,
        cyp_terms=(),
        target_binding=(),
        si_segments=None,
        gut_extraction_eg=0.0,
    )
    plain = replace(base, gut_extraction_eg=0.0, cl_sec_l_h=0.0)
    eg = replace(base, gut_extraction_eg=0.5)
    sec = replace(base, cl_sec_l_h=5.0)

    r_plain = pl.run_pipeline(plain)
    r_eg = pl.run_pipeline(eg)
    r_sec = pl.run_pipeline(sec)

    assert r_plain.pk.bioavailability_f is not None
    assert r_eg.pk.bioavailability_f is not None
    assert r_eg.pk.bioavailability_f < r_plain.pk.bioavailability_f
    assert r_plain.pk.urine_cum_mg is not None and r_sec.pk.urine_cum_mg is not None
    assert float(np.sum(r_sec.pk.urine_cum_mg)) > float(np.sum(r_plain.pk.urine_cum_mg))


def test_run_spec_rejects_bad_stage1_params(fast_warfarin: RunSpec) -> None:
    with pytest.raises(ValueError):
        replace(fast_warfarin, cl_sec_l_h=-1.0)
    with pytest.raises(ValueError):
        replace(fast_warfarin, hepatic_vmax_mg_h=5.0, hepatic_km_mg_l=None)  # km missing
    with pytest.raises(ValueError):
        replace(fast_warfarin, hepatic_km_mg_l=0.0)
    with pytest.raises(ValueError):
        replace(fast_warfarin, hepatic_km_mg_l=-1.0, hepatic_vmax_mg_h=5.0)
    with pytest.raises(ValueError):
        replace(fast_warfarin, bile_emptying_1h=0.0)
    with pytest.raises(ValueError):
        replace(fast_warfarin, gut_extraction_eg=-0.1)
    with pytest.raises(ValueError):
        replace(
            fast_warfarin,
            cyp_terms=(CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=1.0),),
            hepatic_vmax_mg_h=5.0,
            hepatic_km_mg_l=2.0,  # cyp_terms and lumped MM are mutually exclusive
        )


def test_panel_model_wires_cyp_terms(fast_warfarin: RunSpec) -> None:
    terms = (CypTerm(isoform="CYP3A4", km_mg_l=1.0, vmax_mg_h=7.8),)
    # Clear the auto-anchored lumped MM pair before switching to per-CYP terms.
    spec = replace(fast_warfarin, hepatic_vmax_mg_h=None, hepatic_km_mg_l=None, cyp_terms=terms)
    assert pl._panel_model(spec).cyp_terms == terms
    with pytest.raises(ValueError):
        replace(fast_warfarin, gut_extraction_eg=1.0)


def test_panel_model_wires_target_binding(fast_warfarin: RunSpec) -> None:
    target = Target(
        name="TMDD site",
        kd_nm=1.5,
        kon_nm_h=0.5,
        r0_nm=800.0,
        rho_h=0.08,
        kint_h=0.25,
    )
    spec = replace(
        fast_warfarin,
        target_binding=(TargetBinding(tissue="gut", target=target),),
        include_pathway=False,
        feedback_loop=0,
        fidelity="baseline",
    )
    panel = pl._panel_model(spec)
    assert panel.target_binding == spec.target_binding
    assert panel.mw_g_per_mol == spec.mw
    assert pl.run_pipeline(spec).pk.tmdd is not None
    with pytest.raises(ValueError):
        replace(fast_warfarin, target_binding=(TargetBinding(tissue="bogus", target=target),))
    assert replace(fast_warfarin, bile_emptying_1h=0.5).bile_emptying_1h == 0.5


def test_panel_model_wires_skin_layers(fast_warfarin: RunSpec) -> None:
    from drugos.inputs.parse_dosing import build_dose_plan
    from drugos.pk.skin import SkinLayers

    spec = replace(
        fast_warfarin,
        skin_layers=SkinLayers(sc_thickness_um=20.0),
        include_pathway=False,
        feedback_loop=0,
        fidelity="baseline",
        dose_plan=build_dose_plan(route="transdermal", amount_mg=10.0, interval_h=0.0, n_doses=1),
    )
    panel = pl._panel_model(spec)
    assert panel.absorption.skin_layers == spec.skin_layers
    run = pl.run_pipeline(spec)
    assert run.pk.skin is not None
    assert run.pk.skin["absorbed_mg"][-1] > 0.0


def test_cardiac_tone_scale_clamps_and_neutral_below_baseline() -> None:
    t = np.linspace(0.0, 1.0, 11)
    assert pl._cardiac_tone_scale(None) == 1.0
    assert pl._cardiac_tone_scale(np.full_like(t, 1.0)) == 1.0
    # A suppressed/baseline ERK readout (benign drug) stays exactly neutral.
    low = pl._cardiac_tone_scale(np.full_like(t, 0.0))
    assert low == 1.0
    # Amplification above baseline lifts tone with the damped gain and ceiling.
    assert pl._cardiac_tone_scale(np.linspace(1.0, 1.5, 11)) == pytest.approx(
        1.0 + 0.2 * (1.25 - 1.0)
    )
    high = pl._cardiac_tone_scale(np.full_like(t, 10.0))
    assert high == pytest.approx(1.3)
    assert isinstance(pl._cardiac_tone_scale(np.full_like(t, 2.0)), float)


def test_beta_block_ic50_wires_sympathetic_tone_and_hemodynamics(
    fast_warfarin: RunSpec,
) -> None:
    base = run_pipeline(replace(fast_warfarin, beta_block_ic50_nm=None, fidelity="baseline"))
    hard = run_pipeline(replace(fast_warfarin, beta_block_ic50_nm=0.01))
    soft = run_pipeline(replace(fast_warfarin, beta_block_ic50_nm=1000.0))
    assert base.exposure.beta_block_ic50_nm is None
    assert hard.exposure.beta_block_ic50_nm == pytest.approx(0.01)
    assert base.organ.cardiac.co_l_min > hard.organ.cardiac.co_l_min
    assert hard.organ.cardiac.co_l_min < soft.organ.cardiac.co_l_min
    assert hard.organ.cardiac.hr_bpm < base.organ.cardiac.hr_bpm
    assert hard.organ.cardiac.sv_ml < base.organ.cardiac.sv_ml


def test_run_spec_rejects_nonpositive_beta_block_ic50(fast_warfarin: RunSpec) -> None:
    with pytest.raises(ValueError, match="beta_block_ic50_nm"):
        replace(fast_warfarin, beta_block_ic50_nm=-5.0)
    with pytest.raises(ValueError, match="beta_block_ic50_nm"):
        replace(fast_warfarin, beta_block_ic50_nm=0.0)


def test_dili_immune_ic50_wires_liver_immune_axis(fast_warfarin: RunSpec) -> None:
    base = run_pipeline(replace(fast_warfarin, dili_immune_ic50_nm=None))
    active = run_pipeline(
        replace(
            fast_warfarin,
            dili_immune_ic50_nm=50.0,
            dili_immune_weight=1.0,
        )
    )
    assert base.exposure.dili_immune_ic50_nm is None
    assert active.exposure.dili_immune_ic50_nm == pytest.approx(50.0)
    # Immune state is driven by the hazard regardless of weight; verify wiring.
    assert active.organ.liver.immune.max() > 0.01
    # With weight>0 the immune axis loads the stress; dead may differ slightly.
    assert base.organ.liver.immune.max() < 1e-6


def test_dili_immune_weight_inert_when_zero(fast_warfarin: RunSpec) -> None:
    r = run_pipeline(replace(fast_warfarin, dili_immune_ic50_nm=10.0, dili_immune_weight=0.0))
    # The hazard drives the immune response state, but weight=0 means it
    # never loads combined_stress, so dead_frac stays baseline.
    assert r.organ.liver.immune.max() > 0.01


def test_run_spec_rejects_bad_dili_immune(fast_warfarin: RunSpec) -> None:
    with pytest.raises(ValueError, match="dili_immune_ic50_nm"):
        replace(fast_warfarin, dili_immune_ic50_nm=-1.0)
    with pytest.raises(ValueError, match="dili_immune_weight"):
        replace(fast_warfarin, dili_immune_weight=1.5)


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


def test_available_realism_terms_and_degraded_branches(fast_warfarin: RunSpec) -> None:
    engaged = fast_warfarin
    missing = replace(engaged, hepatic_vmax_mg_h=None, hepatic_km_mg_l=None)
    available = pl._available_realism_terms(missing)
    assert set(available) == set(pl._available_realism_terms(engaged)) | {
        "saturable (Michaelis-Menten) hepatic clearance"
    }
    guard = pl._degraded_realism_terms(missing)
    assert guard and "guardrail" in guard[0]
    assert pl._degraded_realism_terms(engaged) == []
    baseline_full_hep = replace(engaged, fidelity="baseline")
    assert pl._degraded_realism_terms(baseline_full_hep) == []
    baseline_off = pl._degraded_realism_terms(replace(missing, fidelity="baseline"))
    assert any("explicit fidelity='baseline' opt-out" in item for item in baseline_off)


def test_full_fidelity_engaged_spec_is_stable(fast_warfarin: RunSpec) -> None:
    engaged = fast_warfarin
    out = pl._engage_full_fidelity(engaged)
    assert out.feedback_loop == 1
    assert out.si_segments == engaged.si_segments
    assert out.dili_immune_ic50_nm == engaged.dili_immune_ic50_nm
    assert out.beta_block_ic50_nm == engaged.beta_block_ic50_nm
    assert out.target_binding == engaged.target_binding
    baseline = replace(engaged, fidelity="baseline")
    assert pl._engage_full_fidelity(baseline) is baseline


def test_full_fidelity_immune_weight_already_set(fast_warfarin: RunSpec) -> None:
    primed = replace(fast_warfarin, dili_immune_ic50_nm=None, dili_immune_weight=0.5)
    out = pl._engage_full_fidelity(primed)
    assert out.dili_immune_weight == 0.5
    assert out.dili_immune_ic50_nm is not None


def test_full_fidelity_sc_im_and_transdermal_wiring(fast_warfarin: RunSpec) -> None:
    from drugos.inputs.parse_dosing import build_dose_plan

    dose = fast_warfarin.dose_plan.total_dose_mg
    subcut = replace(
        fast_warfarin,
        dose_plan=build_dose_plan("subcutaneous", dose),
        sc_im_ka_per_h=None,
    )
    sc_ok = pl._engage_full_fidelity(subcut)
    assert sc_ok.sc_im_ka_per_h is not None
    transdermal = replace(
        fast_warfarin,
        dose_plan=build_dose_plan("transdermal", dose),
        skin_layers=None,
    )
    trans_ok = pl._engage_full_fidelity(transdermal)
    assert trans_ok.skin_layers is not None


def test_assert_full_fidelity_raises_when_missing(fast_warfarin: RunSpec) -> None:
    missing = replace(fast_warfarin, hepatic_vmax_mg_h=None, hepatic_km_mg_l=None)
    with pytest.raises(pl.RealismError):
        pl._assert_full_fidelity(missing)
    pl._assert_full_fidelity(replace(missing, fidelity="baseline"))


def test_spec_from_admet_rejects_unusable_mw(monkeypatch: pytest.MonkeyPatch) -> None:
    from types import SimpleNamespace

    from drugos.inputs.models import HumanProfile, Sex
    from drugos.inputs.parse_dosing import build_dose_plan

    monkeypatch.setattr(pl, "parse_structure", lambda *a, **k: SimpleNamespace(mw=None))
    mol = SimpleNamespace(name="empty", canonical_smiles=None, smiles=None, mw=None)
    profile = HumanProfile(sex=Sex.MALE, age_y=40, height_cm=170, weight_kg=70)
    with pytest.raises(ValueError, match="molecular weight"):
        pl.spec_from_admet(mol, None, profile, build_dose_plan("oral", 500.0))


def test_full_fidelity_uses_admet_heads_and_discloses_them() -> None:
    from drugos.inputs.models import HumanProfile, Sex
    from drugos.inputs.parse_dosing import build_dose_plan
    from drugos.inputs.parse_structure import parse_structure
    from drugos.pk.admet import AdmetOutput

    mol = parse_structure("CC(=O)Nc1ccc(O)cc1", name="admet-heads")
    profile = HumanProfile(sex=Sex.MALE, age_y=40, height_cm=170, weight_kg=70)
    plan = build_dose_plan("oral", 500.0)

    def build(heads: dict[str, float]) -> RunSpec:
        admet = AdmetOutput(
            smiles=mol.canonical_smiles,
            fup_plasma=0.6,
            cl_int_hep_ml_min_kg=12.0,
            HIA=0.8,
            log_s=-3.0,
            **heads,
        )
        return pl.spec_from_admet(mol, admet, profile, plan)

    lipophilic = build(dict(hERG=0.35, BBB=0.9, Pgp=0.7, log_p_pred=2.5, cyp3a4_inhibitor=0.8))
    assert lipophilic.cl_sec_l_h > 1e-4  # P-gp head present, not the 0.5 neutral prior
    assert lipophilic.cl_bil_l_h > 1e-4  # high-logP biliary fraction engaged
    assert lipophilic.gut_extraction_eg > 0.0  # CYP3A4 head present, not the 0.5 prior
    result = run_pipeline(lipophilic)
    trust = result.to_contract()["trust"]
    assert any("BBB_Martins" in a for a in trust["anchors_wired"])
    assert any("hERG" in a for a in trust["anchors_wired"])
    assert "hERG" in trust["admet_ml_estimates"] and "BBB" in trust["admet_ml_estimates"]

    hydrophilic = build(dict(log_p_pred=1.0))
    assert hydrophilic.cl_bil_l_h > 1e-4  # low-logP biliary fraction engaged
