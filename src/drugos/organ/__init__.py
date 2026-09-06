"""Stage 4: organ-level functional endpoints (doc/05 section 4).

Baseline panels convert Stage-1 PBPK free-tissue exposure and Stage-2 target
safety panel into organ functional outcomes:

* ``liver``      — DILI QST (cholestasis, mitochondrial dysfunction,
                   oxidative stress, hepatocyte death -> ALT/AST/bilirubin,
                   DILI grade, Hy's Law).
* ``cardiac``    — hERG/QT electrical axis (QTc, TdP band) + lumped
                   circulation hemodynamics (MAP, CVP, CO, SV).
* ``kidney``     — nephron injury -> GFR + serum creatinine, KDIGO AKI grade.
* ``cns``        — brain free exposure (free-drug hypothesis) + margin grade.
* ``feedback``   — organ dysfunction folded back into PK (4.5).
"""

from drugos.organ.base import auc_nm_h, free_mg_l_to_nm
from drugos.organ.cardiac import (
    CardiacParams,
    CardiacResult,
    HemodynamicsResult,
    ikr_fraction,
    predict_qtc,
    simulate_cardiac,
    simulate_hemodynamics,
    tdpr_band,
)
from drugos.organ.cns import CnsParams, CnsResult, cns_grade, simulate_cns
from drugos.organ.feedback import (
    OrganFeedback,
    apply_pk_scaling,
    clo_01,
    feedback_from_results,
    organ_feedback,
)
from drugos.organ.kidney import (
    KidneyParams,
    KidneyResult,
    aki_grade,
    gfr_trajectory,
    nephron_injury,
    scr_from_gfr,
    simulate_kidney,
)
from drugos.organ.liver import (
    LiverParams,
    LiverStress,
    LiverTrajectory,
    aten_floor_factor,
    combined_stress,
    dili_grade,
    inhibition,
    liver_params_from_panel,
    mitochondrial_block,
    redox_state,
    simulate_liver,
)

__all__ = [
    "CardiacParams",
    "CardiacResult",
    "CnsParams",
    "CnsResult",
    "HemodynamicsResult",
    "KidneyParams",
    "KidneyResult",
    "LiverParams",
    "LiverStress",
    "LiverTrajectory",
    "OrganFeedback",
    "aki_grade",
    "apply_pk_scaling",
    "aten_floor_factor",
    "auc_nm_h",
    "clo_01",
    "cns_grade",
    "combined_stress",
    "dili_grade",
    "feedback_from_results",
    "free_mg_l_to_nm",
    "gfr_trajectory",
    "ikr_fraction",
    "inhibition",
    "liver_params_from_panel",
    "mitochondrial_block",
    "nephron_injury",
    "organ_feedback",
    "predict_qtc",
    "redox_state",
    "scr_from_gfr",
    "simulate_cardiac",
    "simulate_cns",
    "simulate_hemodynamics",
    "simulate_kidney",
    "simulate_liver",
    "tdpr_band",
]
