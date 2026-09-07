"""Validation case library (doc/08): Tier-1 benchmarks + stage analytic/CI cases.

Each case declares an ``EvidenceLevel`` (see ``base.py``) grading how much
epistemic weight its pass carries — from empirically anchored (L3) to merely
numerically self-consistent (L1).  ``run_all`` returns every case in suite
order; the runner regenerates ``validation/report.md``.
"""

from __future__ import annotations

from validation.benchmarks import BENCHMARKS
from validation.cases.base import (
    EVIDENCE_META,
    EVIDENCE_ORDER,
    CaseResult,
    EvidenceLevel,
    MetricResult,
)
from validation.cases.case_admet_bbb_cns import case_admet_bbb_cns
from validation.cases.case_benchmarks import evaluate_benchmark
from validation.cases.case_cardiac_ap_ord import case_cardiac_ap_ord
from validation.cases.case_cardiac_qtc import case_cardiac_qtc
from validation.cases.case_ckdepi_2021 import case_ckdepi_2021
from validation.cases.case_clinical_grading import case_clinical_grading
from validation.cases.case_corpus_calibration import case_corpus_calibration
from validation.cases.case_dose_proportionality import case_dose_proportionality
from validation.cases.case_kidney_gfr import case_kidney_gfr
from validation.cases.case_liver_cholestasis_pbk import case_liver_cholestasis_pbk
from validation.cases.case_liver_dose_response import case_liver_dose_response
from validation.cases.case_mass_balance import case_mass_balance
from validation.cases.case_occupancy_equilibrium import case_occupancy_equilibrium
from validation.cases.case_pathway_amplification import case_pathway_amplification
from validation.cases.case_prospective_fidelity import case_prospective_fidelity
from validation.cases.case_r_bridge import case_r_bridge
from validation.cases.case_risk_ordering import case_risk_ordering
from validation.cases.case_robustness_sanity import case_robustness_sanity
from validation.cases.case_sbml_mapk_validation import case_sbml_mapk_validation
from validation.cases.case_sc_im_depot import case_sc_im_depot
from validation.cases.case_single_pool_analytic import case_single_pool_analytic


def run_all() -> list[CaseResult]:
    results: list[CaseResult] = [evaluate_benchmark(b) for b in BENCHMARKS]
    results.append(case_mass_balance())
    results.append(case_dose_proportionality())
    results.append(case_single_pool_analytic())
    results.append(case_occupancy_equilibrium())
    results.append(case_pathway_amplification())
    results.append(case_liver_dose_response())
    results.append(case_cardiac_qtc())
    results.append(case_kidney_gfr())
    results.append(case_clinical_grading())
    results.append(case_corpus_calibration())
    results.append(case_risk_ordering())
    results.append(case_robustness_sanity())
    results.append(case_sc_im_depot())
    results.append(case_prospective_fidelity())
    results.append(case_r_bridge())
    results.append(case_cardiac_ap_ord())
    results.append(case_admet_bbb_cns())
    results.append(case_sbml_mapk_validation())
    results.append(case_ckdepi_2021())
    results.append(case_liver_cholestasis_pbk())
    return results


__all__ = [
    "BENCHMARKS",
    "CaseResult",
    "EVIDENCE_META",
    "EVIDENCE_ORDER",
    "EvidenceLevel",
    "MetricResult",
    "run_all",
]
