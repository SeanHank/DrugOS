"""ADMET-AI BBB_Martins -> CNS partition wiring (R-4, L2 analytic limit).

Production-model upgrade of the Stage-4 CNS panel (doc/12 §1 row 4d): when an
ADMET-AI call is present, its BBB_Martins classifier head (Martins et al. 2012
molecular BBB dataset, retrained in ADMET-AI; Swanson et al. 2024, BSD-3 code
+ CC-BY weights) decides the residual brain:plasma unbound partition class:

- predicted penetrant  (``P >= 0.5``)  -> ``kpu_brain = 1.0`` (full passive
  free-drug partition, the pipeline default),
- predicted non-penetrant (``P < 0.5``) -> ``kpu_brain = 0.2`` (restricted
  entry, conservative lower bound).

This replaces the bare fixed 0.8/1.0 residual partition class with a verdict
from the production-validated, downloadable BBB classifier.  Benchmarks that
carry no ADMET-AI call are untouched (``None -> 1.0``), so all existing anchors
keep their exact values.  The check is analytic (L2): the mapping is defined
exactly in ``drugos.organ.cns.kpu_brain_from_bbb`` and this case pins the
pipeline-level wiring and the resulting 5x brain-exposure separation.
"""

from __future__ import annotations

from dataclasses import replace

from validation.benchmarks import BENCHMARKS
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pipeline import run_pipeline, spec_from_benchmark_data
from drugos.pk.admet import AdmetOutput

_SMILES = "CC(=O)c1ccc(cc1)C(C(=O)O)"


def case_admet_bbb_cns() -> CaseResult:
    spec = spec_from_benchmark_data(next(b for b in BENCHMARKS if b.name == "warfarin"))
    penetrant = run_pipeline(replace(spec, admet=AdmetOutput(smiles=_SMILES, BBB=0.9)))
    restricted = run_pipeline(replace(spec, admet=AdmetOutput(smiles=_SMILES, BBB=0.1)))
    baseline = run_pipeline(spec)

    kpu_pen = penetrant.exposure.cns_kpu_brain
    kpu_res = restricted.exposure.cns_kpu_brain
    peak_ratio = restricted.organ.cns.peak_brain_free_nm / penetrant.organ.cns.peak_brain_free_nm
    grade_pen = penetrant.organ.cns.grade
    grade_res = restricted.organ.cns.grade
    kpu_base = baseline.exposure.cns_kpu_brain

    metrics = [
        MetricResult(
            "cns_kpu_brain_penetrant",
            kpu_pen,
            1.0,
            1.0,
            "ratio",
            "pass" if abs(kpu_pen - 1.0) < 1e-9 else "FAIL",
        ),
        MetricResult(
            "cns_kpu_brain_non_penetrant",
            kpu_res,
            0.2,
            0.2,
            "ratio",
            "pass" if abs(kpu_res - 0.2) < 1e-9 else "FAIL",
        ),
        MetricResult(
            "benchmark_no_admet_kpu_brain",
            kpu_base,
            1.0,
            1.0,
            "ratio",
            "pass" if abs(kpu_base - 1.0) < 1e-9 else "FAIL",
        ),
        MetricResult(
            "non_penetrant_peak_brain_free_ratio_of_penetrant",
            peak_ratio,
            0.2,
            0.2,
            "ratio",
            "pass" if abs(peak_ratio - 0.2) < 1e-9 else "FAIL",
        ),
        MetricResult(
            "non_penetrant_cns_grade_lte_penetrant",
            float(grade_res),
            -1e9,
            float(grade_pen),
            "grade",
            "pass" if grade_res <= grade_pen else "FAIL",
        ),
    ]
    notes = [
        f"BBB_Martins P=0.9 -> kpu_brain {kpu_pen:.2f}, "
        f"P=0.1 -> {kpu_res:.2f} "
        f"(source: ADMET-AI BBB_Martins head); restricted brain peak = "
        f"{100 * peak_ratio:.0f}% of penetrant; CNS grades {grade_res} <= {grade_pen}; "
        f"baseline (no ADMET-AI) kpu={kpu_base:.2f} untouched — "
        "benchmark anchors unchanged"
    ]
    return CaseResult(
        benchmark="ADMET-AI BBB_Martins -> CNS partition (R-4)",
        passed=all(m.criterion == "pass" for m in metrics),
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )
