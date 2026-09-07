"""Kidney: nephron injury, GFR trajectory and AKI (KDIGO) grading.

doc/05 4.4 baseline.  A drug's free renal exposure drives nephron injury, which
depresses GFR from the physiological baseline; serum creatinine is then the
closed-form balance of production over clearance (Scr = P / GFR), and KDIGO
criteria convert the Scr rise (or GFR fall) into an acute kidney injury grade.

The **baseline** GFR is anchored to the production-validated, published
CKD-EPI 2021 race-free creatinine equation (Levey et al., NEJM 2021) when a
measured serum creatinine is available (doc/12 row 4b; R-6); otherwise the
physiology default (age/sex part) is used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from drugos.organ.base import NDArray, free_mg_l_to_nm

_KAPPA = {"M": 0.9, "F": 0.7}
_ALPHA = {"M": -0.302, "F": -0.241}
_SCR_MGDL_TO_UMOL_L = 88.4


def ckdepi_2021_egfr(
    scr_mg_dl: float, age_y: float, female: bool, bsa_m2: float | None = None
) -> float:
    """eGFR by the CKD-EPI 2021 race-free creatinine equation (mL/min/1.73 m^2).

    ``eGFR = 142 * min(Scr/k,1)^a * max(Scr/k,1)^-1.200 * 0.9938^Age``,
    with the female multiplier and sex-specific ``k``/``a`` (Levey et al.,
    N Engl J Med 2021;385:1737-49).  When ``bsa_m2`` is given the per-subject
    absolute GFR (mL/min) is returned by scaling off 1.73 m^2; otherwise the
    index value (per 1.73 m^2) is returned.
    """
    if scr_mg_dl <= 0:
        raise ValueError("scr_mg_dl must be positive")
    if not (0.0 <= age_y <= 130.0):
        raise ValueError(f"age_y out of range: {age_y}")
    key = "F" if female else "M"
    kappa: float = _KAPPA[key]
    alpha: float = _ALPHA[key]
    scr_k: float = scr_mg_dl / kappa
    egfr = (
        142.0
        * math.pow(min(scr_k, 1.0), alpha)
        * math.pow(max(scr_k, 1.0), -1.200)
        * math.pow(0.9938, age_y)
    )
    if female:
        egfr *= 1.012
    if bsa_m2 is not None:
        if bsa_m2 <= 0:
            raise ValueError("bsa_m2 must be positive")
        egfr *= bsa_m2 / 1.73
    return egfr


def scr_mg_dl_to_umol_l(scr_mg_dl: float) -> float:
    """Convert serum creatinine mg/dL to umol/L (standard clinical factor)."""
    if scr_mg_dl <= 0:
        raise ValueError("scr_mg_dl must be positive")
    return scr_mg_dl * _SCR_MGDL_TO_UMOL_L


@dataclass(frozen=True, slots=True)
class KidneyParams:
    """Kidney sub-model constants."""

    injury_ic50_nm: float = 5.0e6
    injury_hill: float = 2.0
    gfr_floor: float = 0.15


@dataclass(slots=True)
class KidneyResult:
    """Kidney QST result over the exposure horizon."""

    t_h: NDArray
    injury: NDArray
    gfr_ml_min: NDArray
    scr_umol_l: NDArray
    scr_ratio: NDArray
    peak_scr_umol_l: float
    peak_scr_ratio: float
    min_gfr_ml_min: float
    aki_grade: int

    def to_series(self) -> dict[str, NDArray | float]:
        return {
            "t_h": self.t_h,
            "injury": self.injury,
            "gfr_ml_min": self.gfr_ml_min,
            "scr_umol_L": self.scr_umol_l,
            "scr_ratio": self.scr_ratio,
        }


def nephron_injury(conc_nm: NDArray | float, ic50_nm: float, hill: float) -> NDArray:
    """Injury fraction (0..1) driven by free kidney concentration."""
    if ic50_nm <= 0 or hill <= 0:
        raise ValueError("ic50_nm and hill must be positive")
    c = np.asarray(conc_nm, dtype=float)
    c = np.maximum(c, 0.0)
    num: NDArray = c**hill
    denom = ic50_nm**hill
    out: NDArray = num / (num + denom)
    return out


def gfr_trajectory(gfr_base_ml_min: float, injury: NDArray, floor: float = 0.15) -> NDArray:
    """GFR profile at fractional injury, floored to a survival minimum."""
    if gfr_base_ml_min <= 0:
        raise ValueError("gfr_base_ml_min must be positive")
    if not (0.0 <= floor < 1.0):
        raise ValueError("floor must be in [0, 1)")
    return gfr_base_ml_min * (floor + (1.0 - floor) * (1.0 - injury))


def scr_from_gfr(gfr_ml_min: NDArray, prod_umol_min: float) -> NDArray:
    """Serum creatinine (umol/L): steady-state balance P = GFR x Scr."""
    if prod_umol_min <= 0:
        raise ValueError("creatinine production must be positive")
    g = np.asarray(gfr_ml_min, dtype=float)
    if np.any(g <= 0):
        raise ValueError("GFR must stay positive")
    return 1000.0 * prod_umol_min / g


def aki_grade(scr_ratio: float, gfr_ratio: float = 1.0) -> int:
    """KDIGO AKI stage from Scr ratio (or drop in GFR)."""
    if scr_ratio >= 3.0 or gfr_ratio <= 0.35:
        return 3
    if scr_ratio >= 2.0:
        return 2
    if scr_ratio >= 1.5 or gfr_ratio <= 0.5:
        return 1
    return 0


def simulate_kidney(
    t_h: NDArray,
    c_free_mg_l: NDArray,
    mw: float,
    gfr_base_ml_min: float,
    scr_base_umol_l: float,
    params: KidneyParams | None = None,
) -> KidneyResult:
    """Run the nephron/GFR panel on the PBPK free kidney exposure."""
    p = params or KidneyParams()
    t = np.asarray(t_h, dtype=float)
    c = free_mg_l_to_nm(np.asarray(c_free_mg_l, dtype=float), mw)
    if t.ndim != 1 or c.ndim != 1 or t.shape[0] != c.shape[0]:
        raise ValueError("t_h and c_free_mg_l must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")
    if gfr_base_ml_min <= 0 or scr_base_umol_l <= 0:
        raise ValueError("GFR and Scr bases must be positive")

    injury = nephron_injury(c, p.injury_ic50_nm, p.injury_hill)
    gfr = gfr_trajectory(gfr_base_ml_min, injury, p.gfr_floor)
    prod = scr_base_umol_l * gfr_base_ml_min / 1000.0
    scr = scr_from_gfr(gfr, prod)
    scr_ratio = scr / scr_base_umol_l
    gfr_ratio = gfr / gfr_base_ml_min
    level = int(aki_grade(float(np.max(scr_ratio)), float(np.min(gfr_ratio))))
    return KidneyResult(
        t_h=t,
        injury=injury,
        gfr_ml_min=gfr,
        scr_umol_l=scr,
        scr_ratio=scr_ratio,
        peak_scr_umol_l=float(np.max(scr)),
        peak_scr_ratio=float(np.max(scr_ratio)),
        min_gfr_ml_min=float(np.min(gfr)),
        aki_grade=level,
    )


__all__ = [
    "KidneyParams",
    "KidneyResult",
    "aki_grade",
    "ckdepi_2021_egfr",
    "gfr_trajectory",
    "nephron_injury",
    "scr_from_gfr",
    "scr_mg_dl_to_umol_l",
    "simulate_kidney",
]
