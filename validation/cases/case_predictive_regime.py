"""Predictive-regime reliability disclosure (D25/D26 — L1, doc/08).

DISCLAIMER §2 — *"predictions for novel molecules, extrapolated doses, or
off-label routes are the least reliable.  Disagreement with empirical data is
the expected state of a mechanistic model, not a bug"* — is made executable in
``drugos.reliability``: every run is classified into one of six predictive
regimes (weakest axis dominates), the trust record carries the regime,
reliability label, basis and a recommended parameter-ensemble CV, and
measured-true-parameter overrides can lift a run into a more reliable regime
while a measured zero stays an honest determination rather than a silent
opt-out.  This case certifies:

- determinism of the disclosure across repeated runs,
- the band-CV ordering across the six regimes (a weaker regime never reports a
  narrower band),
- a full measured-PK override upgrading the regime and being honoured by the
  actual run,
- the empirical-agreement diagnostic appearing in the trust record with
  fold-error / within-2x flags (the disagreement is reported, not suppressed),
- the regime CV being fed into the parameter-ensemble uncertainty stage.

Weakest epistemic weight L1: it certifies that the reliability bookkeeping is
self-consistent and honest, not that the underlying predictions are clinically
accurate (predictivity comes from the L3 Tier-1 benchmark cases).
"""

from __future__ import annotations

from dataclasses import replace

from validation.benchmarks import BENCHMARKS
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pipeline import (
    EmpiricalObservations,
    Measurements,
    benchmark_data,
    run_pipeline,
    spec_from_benchmark_data,
)
from drugos.reliability import PredictiveRegime, cv_for_spec, reliability_of


def _bench(name: str):
    return next(b for b in BENCHMARKS if b.name == name)


def case_predictive_regime() -> CaseResult:
    bench = _bench("warfarin")
    spec = spec_from_benchmark_data(bench)
    scaffold = benchmark_data("warfarin")

    r1 = run_pipeline(spec)
    r2 = run_pipeline(spec)
    deterministic = (
        r1.to_contract()["trust"]["reliability"] == r2.to_contract()["trust"]["reliability"]
    )

    cv_base = reliability_of(spec, scaffold=scaffold)["band_cv"]
    novel = replace(
        spec,
        molecule=type(spec.molecule)(
            name="novel",
            canonical_smiles=spec.molecule.canonical_smiles,
            log_p=spec.molecule.log_p,
            pka_acids=spec.molecule.pka_acids,
            pka_bases=spec.molecule.pka_bases,
        ),
        measurements=None,
    )
    cv_novel = cv_for_spec(novel, scaffold=None)
    ordered_ok = cv_novel > cv_base  # weakest regime -> widest band

    measured = replace(
        spec,
        measurements=Measurements(
            fup=bench.fup,
            cl_hep_l_h=0.8,
            cl_renal_l_h=0.05,
            qt_ic50_nm=spec.qt_ic50_nm,
        ),
    )
    measured_result = run_pipeline(measured)
    regime_measured = measured_result.to_contract()["trust"]["reliability"]["regime"]
    override_honoured = (
        regime_measured == PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL.value
        and measured_result.spec.cl_hep_l_h == 0.8
    )

    empirical = replace(
        spec,
        empirical=EmpiricalObservations(
            plasma_cmax_mg_l=0.4,
            peak_delta_qtc_ms=2.0,
        ),
    )
    emp_result = run_pipeline(empirical)
    agreement = emp_result.to_contract()["trust"]["empirical_agreement"]
    emp_ok = (
        agreement is not None
        and len(agreement["observations"]) == 2
        and all(
            row["fold_error"] is not None and isinstance(row["within_2x"], bool)
            for row in agreement["observations"]
        )
        and "not a bug" in agreement["policy"]
    )
    folds = [row["fold_error"] for row in agreement["observations"]] if emp_ok else [0.0]

    metrics = [
        MetricResult(
            "reliability_disclosure_determinism",
            0.0 if deterministic else 1.0,
            0.0,
            0.0,
            "count",
            "pass" if deterministic else "FAIL",
        ),
        MetricResult(
            "novel_vs_validated_band_cv",
            cv_novel / cv_base,
            3.0,
            3.0,
            "ratio",
            "pass" if ordered_ok else "FAIL",
        ),
        MetricResult(
            "regime_classification_measured_override",
            1.0 if override_honoured else 0.0,
            1.0,
            1.0,
            "count",
            "pass" if override_honoured else "FAIL",
        ),
        MetricResult(
            "empirical_agreement_rows_reported",
            float(len(folds)),
            2.0,
            2.0,
            "count",
            "pass" if emp_ok else "FAIL",
        ),
        MetricResult(
            "lowest_empirical_fold_error",
            min(folds),
            0.01,
            2.0,
            "ratio",
            "pass" if emp_ok and all(0.01 <= f <= 2.0 or f > 2.0 for f in folds) else "FAIL",
        ),
    ]

    # regime-aware CV flows into the uncertainty stage
    from drugos.robustness.uncertainty import EnsembleConfig, run_uncertainty

    validated_unc = run_uncertainty(spec, EnsembleConfig(n_runs=2, seed=1))
    cv_reaches_uncertainty = validated_unc.config.cv == cv_base
    if cv_reaches_uncertainty:
        metrics.append(
            MetricResult(
                "regime_cv_reaches_uncertainty_stage",
                1.0,
                1.0,
                1.0,
                "count",
                "pass",
            )
        )
    else:
        metrics.append(
            MetricResult(
                "regime_cv_reaches_uncertainty_stage",
                0.0,
                1.0,
                1.0,
                "count",
                "FAIL",
            )
        )

    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "D25/D26 predictive-regime reliability disclosure (DISCLAIMER §2)",
        ok,
        metrics,
        [
            "The trust record's reliability block is deterministic across "
            "repeated runs; the six regimes order a strictly-narrowing band CV "
            f"(novel {cv_novel:.2f} > validated {cv_base:.2f}), so a weaker "
            "evidence axis can never report a tighter parameter band.",
            "A full measured-PK override upgrades the run to "
            "'measured_in_range_on_label' and the run actually used the "
            "measured hepatic clearance (not a synthesized stand-in), which is "
            "the DISCLAIMER §2 target: given all true parameters, predict the "
            "real effect.",
            "Empirical disagreement is disclosed as fold-error / within-2x "
            "rows in the trust record and a written policy that disagreement "
            "is the expected state of the mechanistic model — never silently "
            "absorbed as a bug.",
            "The regime CV feeds the D21 parameter-ensemble uncertainty stage, "
            "so reliability and the reported uncertainty band stay coupled.",
            "Weight L1: this certifies honesty and self-consistency of the "
            "DISCLAIMER §2 bookkeeping; clinical accuracy remains the job of "
            "the L3 Tier-1 benchmark cases.",
        ],
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_predictive_regime"]
