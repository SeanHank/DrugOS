"""Predictive-regime reliability and measured-override tests (doc/07 D25/D26).

Covers ``drugos.reliability``: regime classification on the three evidence
axes (chemistry / dose window / route), the full-true-parameters measured
override (``Measurements`` + ``apply_measurements``), the empirical-data
agreement diagnostic (``empirical_agreement`` in the trust record), regime-
aware parameter-ensemble CVs, and the machine-readable ``reliability`` block
DISCLAIMER §2 turns into an executable per-run disclosure.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pipeline import (
    EmpiricalObservations,
    Measurements,
    RunSpec,
    benchmark_data,
    run_pipeline,
)
from drugos.reliability import (
    PredictiveRegime,
    agreement_rows,
    apply_measurements,
    classify_regime,
    cv_for_spec,
    reliability_of,
)


def _novel(fast_apap: RunSpec) -> RunSpec:
    """A chemistry-identical molecule whose name has no benchmark anchor."""
    return replace(
        fast_apap,
        molecule=type(fast_apap.molecule)(
            name="novelx",
            canonical_smiles="CC(=O)Nc1ccc(O)cc1",
            log_p=0.51,
            pka_acids=[9.9],
            pka_bases=[],
        ),
        measurements=None,
    )


def test_regime_scaffold_in_range(fast_apap: RunSpec) -> None:
    scaffold = benchmark_data("acetaminophen")
    regime = classify_regime(fast_apap, scaffold=scaffold)
    assert regime == PredictiveRegime.VALIDATED_IN_RANGE_ON_LABEL
    assert cv_for_spec(fast_apap, scaffold=scaffold) == pytest.approx(0.20)


def test_regime_extrapolated_dose(fast_apap_od: RunSpec) -> None:
    scaffold = benchmark_data("acetaminophen")
    regime = classify_regime(fast_apap_od, scaffold=scaffold)
    assert regime == PredictiveRegime.VALIDATED_EXTRAPOLATED_DOSE
    assert cv_for_spec(fast_apap_od, scaffold=scaffold) == pytest.approx(0.35)


def test_regime_offlabel_route(fast_apap: RunSpec) -> None:
    off = replace(fast_apap, dose_plan=build_dose_plan("transdermal", 1000.0))
    scaffold = benchmark_data("acetaminophen")
    regime = classify_regime(off, scaffold=scaffold)
    assert regime == PredictiveRegime.VALIDATED_OFFLABEL_ROUTE
    assert cv_for_spec(off, scaffold=scaffold) == pytest.approx(0.50)


def test_regime_novel_is_least_reliable(fast_apap: RunSpec) -> None:
    novel = _novel(fast_apap)
    regime = classify_regime(novel)
    assert regime == PredictiveRegime.NOVEL_MOLECULE
    assert cv_for_spec(novel) == pytest.approx(0.60)


def test_regime_partial_evidence(fast_apap: RunSpec) -> None:
    partial = replace(fast_apap, measurements=Measurements(fup=0.75))
    regime = classify_regime(partial)
    assert regime == PredictiveRegime.PARTIAL_EVIDENCE


def test_regime_measured_full_pk_is_strongest(fast_apap: RunSpec) -> None:
    measured = replace(
        fast_apap,
        measurements=Measurements(fup=0.75, cl_hep_l_h=22.0, cl_renal_l_h=2.0),
    )
    scaffold = benchmark_data("acetaminophen")
    assert (
        classify_regime(measured, scaffold=scaffold) == PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL
    )
    assert cv_for_spec(measured, scaffold=scaffold) == pytest.approx(0.15)
    # measured chemistry outweighs the scaffold label on the molecule axis
    assert classify_regime(measured) == PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL


def test_measurements_unknown_field_rejected() -> None:
    with pytest.raises(ValueError, match="unknown measurement field"):
        Measurements.from_dict({"cl_hep_l_h": 20.0, "nonsense": 1.0})
    with pytest.raises(ValueError, match="boolean"):
        Measurements.from_dict({"fa": True})


def test_apply_measurements_overrides_and_discloses(fast_apap: RunSpec) -> None:
    spec = replace(
        fast_apap,
        measurements=Measurements(fup=0.75, cl_hep_l_h=22.0, cl_renal_l_h=2.0, qt_ic50_nm=8.0e4),
    )
    applied = apply_measurements(spec)
    assert applied.fup == pytest.approx(0.75)
    assert applied.cl_hep_l_h == pytest.approx(22.0)
    assert applied.cl_renal_l_h == pytest.approx(2.0)
    assert applied.qt_ic50_nm == pytest.approx(8.0e4)
    assert "cl_hep_l_h" in applied.estimate_basis
    assert "measured true-parameter override" in applied.estimate_basis["cl_hep_l_h"]


def test_apply_measurements_noop_without_block(fast_apap: RunSpec) -> None:
    assert apply_measurements(fast_apap) is fast_apap


def test_measured_zero_keeps_full_fidelity_engaged(fast_apap: RunSpec) -> None:
    # A measured zero for tubular secretion is a true determination ('no
    # active secretion'), not a silent opt-out: the run must not raise and the
    # term must be disclosed as measured.
    spec = replace(fast_apap, measurements=Measurements(cl_sec_l_h=0.0))
    applied = apply_measurements(spec)
    assert applied.cl_sec_l_h == 0.0
    result = run_pipeline(applied)
    trust = result.to_contract()["trust"]
    assert any("tubular secretion" in term for term in trust["mechanism_terms_engaged"])


def test_measured_run_contract_reliability(fast_apap: RunSpec) -> None:
    spec = replace(
        fast_apap,
        measurements=Measurements(fup=0.75, cl_hep_l_h=22.0, cl_renal_l_h=2.0),
    )
    result = run_pipeline(spec)
    reliability = result.to_contract()["trust"]["reliability"]
    assert reliability["regime"] == PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL.value
    assert reliability["reliability"] == "high"
    assert "measured" in reliability["basis"]
    # the run actually used the measured values, not a synthesized stand-in
    assert result.spec.cl_hep_l_h == pytest.approx(22.0)


def test_contract_reliability_present_and_deterministic(fast_apap: RunSpec) -> None:
    first = run_pipeline(fast_apap)
    second = run_pipeline(fast_apap)
    trust = first.to_contract()["trust"]
    for key in ("regime", "label", "reliability", "basis", "band_cv", "disclaimer"):
        assert key in trust["reliability"]
    assert trust["reliability"]["regime"] == PredictiveRegime.VALIDATED_IN_RANGE_ON_LABEL.value
    assert first.to_contract() == second.to_contract()


def test_empirical_agreement_disclosed(fast_apap: RunSpec) -> None:
    # observed plasma Cmax ~ realistic clinical value for paracetamol 1 g PO
    spec = replace(
        fast_apap,
        empirical=EmpiricalObservations(plasma_cmax_mg_l=12.0, peak_alt_uln=2.0),
    )
    result = run_pipeline(spec)
    agreement = result.to_contract()["trust"]["empirical_agreement"]
    assert agreement is not None
    assert "DISCLAIMER" in agreement["policy"] and "not a bug" in agreement["policy"]
    endpoints = {row["endpoint"] for row in agreement["observations"]}
    assert endpoints == {"plasma_cmax_mg_l", "peak_alt_uln"}
    cmax = next(row for row in agreement["observations"] if row["endpoint"] == "plasma_cmax_mg_l")
    assert cmax["predicted"] > 0
    assert cmax["fold_error"] == pytest.approx(12.0 / cmax["predicted"], abs=1e-3)
    assert cmax["within_2x"] is True
    # a grossly mismatched observation flips the diagnostic the other way
    spec_miss = replace(
        fast_apap,
        empirical=EmpiricalObservations(plasma_cmax_mg_l=1.0),
    )
    miss = run_pipeline(spec_miss)
    miss_row = next(
        row
        for row in miss.to_contract()["trust"]["empirical_agreement"]["observations"]
        if row["endpoint"] == "plasma_cmax_mg_l"
    )
    assert miss_row["fold_error"] < 1.0
    assert miss_row["within_2x"] is False


def test_agreement_rows_handles_missing_prediction() -> None:
    obs = EmpiricalObservations(plasma_cmax_mg_l=1.0)
    rows = agreement_rows(obs, {"plasma_cmax_mg_l": 0.0})
    assert rows[0]["fold_error"] is None
    assert rows[0]["within_2x"] is False


def test_band_cv_all_regimes_ordered(fast_apap: RunSpec) -> None:
    from drugos.pipeline import benchmark_data as bd
    from drugos.reliability import REGIMES

    scaffold = bd("acetaminophen")
    # explicit weakest->strongest construction (band CV must shrink accordingly)
    offlabel = replace(fast_apap, dose_plan=build_dose_plan("transdermal", 1000.0))
    partial = replace(fast_apap, measurements=Measurements(fa=0.9))
    extrapolated = replace(fast_apap, dose_plan=build_dose_plan("oral", 4000.0))
    measured = replace(
        fast_apap,
        measurements=Measurements(fup=0.75, cl_hep_l_h=22.0, cl_renal_l_h=2.0),
    )
    cases = [
        (PredictiveRegime.NOVEL_MOLECULE, _novel(fast_apap), None),
        (PredictiveRegime.VALIDATED_OFFLABEL_ROUTE, offlabel, scaffold),
        (PredictiveRegime.PARTIAL_EVIDENCE, partial, None),
        (PredictiveRegime.VALIDATED_EXTRAPOLATED_DOSE, extrapolated, scaffold),
        (PredictiveRegime.VALIDATED_IN_RANGE_ON_LABEL, fast_apap, scaffold),
        (PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL, measured, scaffold),
    ]
    cvs = [cv_for_spec(spec, scaffold=s) for _, spec, s in cases]
    assert cvs == sorted(cvs, reverse=True)
    for (regime, _, _), cv in zip(cases, cvs, strict=True):
        assert cv == REGIMES[regime].band_cv
    # weakest regime carries the widest band, strongest the narrowest
    assert cvs[0] == 0.60
    assert cvs[-1] == 0.15


def test_regime_aware_uncertainty_cv(fast_apap: RunSpec) -> None:
    from drugos.robustness.uncertainty import EnsembleConfig, run_uncertainty

    validated = run_uncertainty(fast_apap, EnsembleConfig(n_runs=2, seed=1))
    assert validated.config.cv == pytest.approx(0.20)
    novel = run_uncertainty(_novel(fast_apap), EnsembleConfig(n_runs=2, seed=1))
    assert novel.config.cv == pytest.approx(0.60)
    explicit = run_uncertainty(fast_apap, EnsembleConfig(n_runs=2, seed=1, cv=0.3))
    assert explicit.config.cv == pytest.approx(0.30)


def test_reliability_of_schema(fast_apap: RunSpec) -> None:
    info = reliability_of(fast_apap, scaffold=benchmark_data("acetaminophen"))
    assert info["regime"] == "validated_in_range_on_label"
    assert isinstance(info["band_cv"], float)
    first = run_pipeline(fast_apap)
    second = run_pipeline(fast_apap)
    assert (
        first.to_contract()["trust"]["reliability"] == second.to_contract()["trust"]["reliability"]
    )


def test_measurements_pk_counts() -> None:
    empty = Measurements()
    assert not empty.is_pk_complete()
    assert empty.count_pk_partial() == 0
    one = Measurements(fup=0.5)
    assert one.count_pk_partial() == 1
    full = Measurements(fup=0.5, cl_hep_l_h=10.0, cl_renal_l_h=1.0)
    assert full.is_pk_complete()


def test_measurements_validation_edges() -> None:
    for kwargs in (
        {"qt_ic50_nm": 0.0},
        {"fup": 1.5},
        {"fa": 1.1},
        {"cl_hep_l_h": -1.0},
        {"hepatic_vmax_mg_h": 10.0},
        {"hepatic_vmax_mg_h": 10.0, "hepatic_km_mg_l": -1.0},
    ):
        with pytest.raises(ValueError):
            Measurements(**kwargs)


def test_measurements_from_dict_edges() -> None:
    with pytest.raises(ValueError):
        Measurements.from_dict(["not", "a", "mapping"])
    with pytest.raises(ValueError):
        Measurements.from_dict({"fup": "abc"})
    assert Measurements.from_dict({"cl_hep_l_h": None}).cl_hep_l_h is None


def test_empirical_observations_validation_edges() -> None:
    with pytest.raises(ValueError):
        EmpiricalObservations(plasma_cmax_mg_l=-1.0)
    with pytest.raises(ValueError):
        EmpiricalObservations.from_dict(None)
    with pytest.raises(ValueError):
        EmpiricalObservations.from_dict({"nonsense_field": 1.0})
    assert EmpiricalObservations.from_dict(
        {"plasma_cmax_mg_l": 1.0}
    ).plasma_cmax_mg_l == pytest.approx(1.0)
    assert EmpiricalObservations.from_dict({"auc_last_mg_h_l": None}).auc_last_mg_h_l is None


def test_regime_meta_validation_edges() -> None:
    from drugos.reliability import RegimeMeta

    with pytest.raises(ValueError):
        RegimeMeta(label="x", reliability="high", band_cv=0.0, disclaimer="x")
    with pytest.raises(ValueError):
        RegimeMeta(label="x", reliability="bogus", band_cv=0.5, disclaimer="x")


def test_apply_measurements_all_empty_is_noop_copy(fast_apap: RunSpec) -> None:
    wrapped = replace(fast_apap, measurements=Measurements())
    out = apply_measurements(wrapped)
    assert out is not wrapped
    assert out.estimate_basis == wrapped.estimate_basis


def test_reliability_of_partial_measured_basis(fast_apap: RunSpec) -> None:
    partial = replace(fast_apap, measurements=Measurements(fa=0.9))
    info = reliability_of(partial)
    assert "1 measured PK/absorption scalar(s)" in info["basis"]
    assert "measured fields: fa" in info["basis"]
    assert info["regime"] == PredictiveRegime.PARTIAL_EVIDENCE.value


def test_reliability_of_empty_measured_basis(fast_apap: RunSpec) -> None:
    info = reliability_of(replace(fast_apap, measurements=Measurements()))
    assert "no benchmark/clinical anchor and no measured parameters" in info["basis"]


def test_empirical_agreement_rendered_in_reports(fast_apap: RunSpec) -> None:
    from drugos.report.render import render_html, render_markdown

    spec = replace(
        fast_apap,
        empirical=EmpiricalObservations(plasma_cmax_mg_l=12.0, peak_alt_uln=2.0),
    )
    result = run_pipeline(spec)
    md = render_markdown(result)
    assert "## Empirical agreement" in md
    assert "within 2x" in md
    assert "| Endpoint | Observed | Predicted | Fold error | Within 2x |" in md
    page = render_html(result)
    assert "<th>Within 2x</th>" in page
    assert "Endpoint" in page
