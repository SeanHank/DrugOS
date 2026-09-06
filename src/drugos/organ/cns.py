"""CNS: brain free-exposure organ panel (Phase 7 organ extension).

doc/07 Phase 7 "additional organ panels — CNS".  Rather than a full
neuropharmacodynamic model, the CNS panel applies the free-drug hypothesis:
the unbound interstitial brain concentration tracks the unbound plasma
concentration (rapid passive equilibration across an intact/passive BBB for
small molecules), scaled by a brain:plasma unbound partition ``kpu_brain``.
Inside the pipeline the driver is the PBPK brain compartment free exposure
(``simulate_cns`` receives ``unbound_tissues["brain"]`` with ``kpu_brain=1.0``)
so the brain partition used is the Rodgers–Rowland value already computed for
the brain tissue, and ``kpu_brain`` here is the residual scale for direct calls.

An exposure-ratio grade (Cmax_brain_free / IC50) is derived from CTCAE-style
margin thresholds, exactly parallel to how ``drugos.organ.kidney`` interprets
in-vitro potency against tissue exposure.  The grade feeds the CNS endpoint of
the Stage-5 composite toxicity fusion so the CNS line becomes exposure-anchored
rather than a bare class prior.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from drugos.organ.base import NDArray, free_mg_l_to_nm


@dataclass(frozen=True, slots=True)
class CnsParams:
    """CNS sub-model constants.

    ``ic50_nm`` is the in-vitro neurotoxicity potency class prior (low
    confidence, 100 uM) unless an assay- or structure-derived value overrides it
    on the ``RunSpec``; ``kpu_brain`` is the brain:plasma unbound partition
    coefficient (default 0.8, consistent with diffusible small molecules);
    ``ratio_thresholds`` are the exposure-margin ladders for grades 1..4.
    """

    ic50_nm: float = 1.0e5
    kpu_brain: float = 0.8
    ratio_thresholds: tuple[float, float, float, float] = (0.1, 0.32, 1.0, 3.2)

    def __post_init__(self) -> None:
        if self.ic50_nm <= 0:
            raise ValueError("ic50_nm must be positive")
        if self.kpu_brain <= 0:
            raise ValueError("kpu_brain must be positive")
        if len(self.ratio_thresholds) != 4 or any(
            self.ratio_thresholds[i] >= self.ratio_thresholds[i + 1]
            for i in range(len(self.ratio_thresholds) - 1)
        ):
            raise ValueError("ratio_thresholds must be 4 strictly increasing values")


@dataclass(slots=True)
class CnsResult:
    """CNS organ-panel result over the exposure horizon."""

    t_h: NDArray
    brain_free_nm: NDArray
    peak_brain_free_nm: float
    exposure_ratio: float
    grade: int


def cns_grade(ratio: float, thresholds: tuple[float, float, float, float]) -> int:
    """CTCAE-style grade from the Cmax_brain_free/IC50 exposure margin."""
    return min(4, sum(1 for t in thresholds if ratio >= t))


def simulate_cns(
    t_h: NDArray,
    plasma_free_mg_l: NDArray,
    mw: float,
    params: CnsParams = CnsParams(),
) -> CnsResult:
    """Brain free exposure from the free-drug hypothesis, graded by margin."""
    if mw <= 0:
        raise ValueError("mw must be positive")
    t = np.asarray(t_h, dtype=float)
    c_free = np.asarray(plasma_free_mg_l, dtype=float)
    if len(t) != len(c_free) or len(t) == 0:
        raise ValueError("t_h and plasma_free_mg_l must be equal, non-empty arrays")
    brain_free_nm: NDArray = free_mg_l_to_nm(c_free, mw) * params.kpu_brain
    peak = float(np.max(brain_free_nm))
    ratio = peak / params.ic50_nm
    return CnsResult(
        t_h=t,
        brain_free_nm=brain_free_nm,
        peak_brain_free_nm=peak,
        exposure_ratio=ratio,
        grade=cns_grade(ratio, params.ratio_thresholds),
    )


__all__ = ["CnsParams", "CnsResult", "cns_grade", "simulate_cns"]
