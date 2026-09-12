"""Validation case library (doc/08): Tier-1 benchmarks + stage analytic/CI cases.

"can be intellectually staged in a defined pipeline order, and ``run_all`` /
``run_parallel`` return every case in suite order; the runner regenerates
``validation/report.md``.
"""

from __future__ import annotations

import os
import pickle
import sys
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from typing import Any

from validation.benchmarks import BENCHMARKS
from validation.cases.base import (
    EVIDENCE_META,
    EVIDENCE_ORDER,
    CaseResult,
    EvidenceLevel,
    MetricResult,
)
from validation.cases.case_acat_multisegment_si import case_acat_multisegment_si
from validation.cases.case_admet_bbb_cns import case_admet_bbb_cns
from validation.cases.case_benchmarks import evaluate_benchmark
from validation.cases.case_bioavailability_f import case_bioavailability_f
from validation.cases.case_cardiac_ap_ord import case_cardiac_ap_ord
from validation.cases.case_cardiac_qtc import case_cardiac_qtc
from validation.cases.case_cardiac_sympathetic_suppression import (
    case_cardiac_sympathetic_suppression,
)
from validation.cases.case_cheng_prusoff_conversion import case_cheng_prusoff_conversion
from validation.cases.case_ckdepi_2021 import case_ckdepi_2021
from validation.cases.case_clearance_mechanisms import case_clearance_mechanisms
from validation.cases.case_clinical_grading import case_clinical_grading
from validation.cases.case_corpus_calibration import case_corpus_calibration
from validation.cases.case_cyp_kinetics import case_cyp_kinetics
from validation.cases.case_dili_immune_activation import case_dili_immune_activation
from validation.cases.case_dose_proportionality import case_dose_proportionality
from validation.cases.case_full_chain_admet_to_report import case_full_chain_admet_to_report
from validation.cases.case_herg_calibration import case_herg_calibration
from validation.cases.case_kidney_gfr import case_kidney_gfr
from validation.cases.case_liver_cholestasis_pbk import case_liver_cholestasis_pbk
from validation.cases.case_liver_dose_response import case_liver_dose_response
from validation.cases.case_mass_balance import case_mass_balance
from validation.cases.case_occupancy_equilibrium import case_occupancy_equilibrium
from validation.cases.case_pathway_amplification import case_pathway_amplification
from validation.cases.case_pathway_organ_coupling import case_pathway_organ_coupling
from validation.cases.case_predictive_regime import case_predictive_regime
from validation.cases.case_prospective_fidelity import case_prospective_fidelity
from validation.cases.case_r_bridge import case_r_bridge
from validation.cases.case_risk_ordering import case_risk_ordering
from validation.cases.case_robustness_sanity import case_robustness_sanity
from validation.cases.case_sbml_mapk_validation import case_sbml_mapk_validation
from validation.cases.case_sc_im_depot import case_sc_im_depot
from validation.cases.case_single_pool_analytic import case_single_pool_analytic
from validation.cases.case_tmdd_drug_disposition import case_tmdd_drug_disposition
from validation.cases.case_transdermal_multi_layer import case_transdermal_multi_layer

_ANALYTIC_CASES: tuple[Callable[[], CaseResult], ...] = (
    case_mass_balance,
    case_dose_proportionality,
    case_single_pool_analytic,
    case_occupancy_equilibrium,
    case_pathway_amplification,
    case_liver_dose_response,
    case_cardiac_qtc,
    case_kidney_gfr,
    case_clinical_grading,
    case_corpus_calibration,
    case_risk_ordering,
    case_robustness_sanity,
    case_sc_im_depot,
    case_prospective_fidelity,
    case_r_bridge,
    case_cardiac_ap_ord,
    case_admet_bbb_cns,
    case_sbml_mapk_validation,
    case_ckdepi_2021,
    case_liver_cholestasis_pbk,
    case_pathway_organ_coupling,
    case_bioavailability_f,
    case_herg_calibration,
    case_clearance_mechanisms,
    case_cheng_prusoff_conversion,
    case_cyp_kinetics,
    case_dili_immune_activation,
    case_acat_multisegment_si,
    case_full_chain_admet_to_report,
    case_tmdd_drug_disposition,
    case_transdermal_multi_layer,
    case_cardiac_sympathetic_suppression,
    case_predictive_regime,
)

_Step = tuple[Callable[..., CaseResult], tuple[Any, ...]]


def _suite_steps() -> list[_Step]:
    steps: list[_Step] = [(evaluate_benchmark, (b,)) for b in BENCHMARKS]
    for case in _ANALYTIC_CASES:
        steps.append((case, ()))
    return steps


def _run_step(step: _Step) -> CaseResult:
    fn, args = step
    return fn(*args)


def run_all() -> list[CaseResult]:
    return [_run_step(step) for step in _suite_steps()]


def run_parallel(jobs: int | None = None) -> list[CaseResult]:
    """Run the full suite, evaluating independent cases in separate processes.

    ``jobs == 0`` forces the sequential path.  Each case is a self-contained
    top-level callable (left to right in suite order), so the results come
    back in the same order as ``run_all``.  Scheduler failures fall back to a
    sequential re-run instead of silently truncating the report.
    """
    steps = _suite_steps()
    if jobs == 0:
        return [_run_step(step) for step in steps]
    n = jobs or min(max(os.cpu_count() or 1, 1), 8)
    try:
        with ProcessPoolExecutor(max_workers=n) as executor:
            return list(executor.map(_run_step, steps))
    except (OSError, BrokenProcessPool, pickle.PicklingError, RuntimeError) as exc:
        msg = f"[validation] parallel scheduler failed ({exc!r}); falling back to sequential"
        print(msg, file=sys.stderr)
        return [_run_step(step) for step in steps]


__all__ = [
    "BENCHMARKS",
    "CaseResult",
    "EVIDENCE_META",
    "EVIDENCE_ORDER",
    "EvidenceLevel",
    "MetricResult",
    "run_all",
    "run_parallel",
]
