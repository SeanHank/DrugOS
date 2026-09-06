"""Organ -> PK feedback (doc/05 4.5, monolithic direct scaling).

Organ dysfunction is fed back into the Stage-1 model as multiplicative
scaling of hepatic clearance (liver death), GFR/renal clearance (kidney
injury) and organ perfusion (cardiac output), e.g. for chronic-exposure runs
or population sensitivity sweeps.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass

from drugos.pk.pbpk_build import PBPKModel


@dataclass(frozen=True, slots=True)
class OrganFeedback:
    """Multiplicative PK scalings derived from organ state."""

    co_scale: float
    hepatic_cl_scale: float
    gfr_scale: float


def clo_01(x: float) -> float:
    """Clamp a dimensionless scale into [0, 1]."""
    return max(0.0, min(1.0, x))


def organ_feedback(
    co_factor: float,
    hepatic_factor: float,
    gfr_factor: float,
) -> OrganFeedback:
    """Build PK feedback from normalised organ function factors (0=worst)."""
    return OrganFeedback(
        co_scale=clo_01(co_factor),
        hepatic_cl_scale=clo_01(hepatic_factor),
        gfr_scale=clo_01(gfr_factor),
    )


def feedback_from_results(
    liver_dead_frac: float,
    gfr_fraction_of_base: float,
    co_fraction: float,
    hepatic_weight: float = 0.5,
) -> OrganFeedback:
    """Map Stage-4 organ outputs to PK scalings (doc/05 4.5)."""
    hepatic = 1.0 - hepatic_weight * clo_01(liver_dead_frac)
    return organ_feedback(co_fraction, hepatic, gfr_fraction_of_base)


def apply_pk_scaling(model: PBPKModel, feedback: OrganFeedback) -> PBPKModel:
    """Return a copy of ``model`` with organ dysfunction folded into PK."""
    modified = copy.copy(model)
    modified.cl_hep_l_h = model.cl_hep_l_h * feedback.hepatic_cl_scale
    modified.cl_renal_l_h = model.cl_renal_l_h * feedback.gfr_scale
    phys = copy.copy(model.physiology)
    phys.cardiac_output_l_min = model.physiology.cardiac_output_l_min * feedback.co_scale
    phys.organ_flow = {k: v * feedback.co_scale for k, v in model.physiology.organ_flow.items()}
    modified.physiology = phys
    return modified


__all__ = [
    "OrganFeedback",
    "apply_pk_scaling",
    "clo_01",
    "feedback_from_results",
    "organ_feedback",
]
