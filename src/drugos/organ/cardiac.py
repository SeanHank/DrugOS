"""Cardiovascular: hERG/QT electrical axis + lumped circulation hemodynamics.

doc/05 4.3 baseline.  The QT axis converts a Stage-2 hERG blockade fraction
(occupancy at the hERG site) into IKr current reduction via an Emax relation
and a Fridericia-corrected QTc with a torsades-de-pointes (TdP) risk band.
The hemodynamic axis is a minimal Physiome-style two-compartment Windkessel
(arterial/venous capacitances, systemic resistance) whose cardiac output is
built from heart rate and stroke volume, with optional inotropy/chronotropy
modulation feeding contractility and rate (doc/05 4.3 bullet 3).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

from drugos.organ.base import NDArray
from drugos.pk.physiology import HumanPhysiology

TdPBand = Literal["none", "low", "moderate", "high"]


@dataclass(frozen=True, slots=True)
class CardiacParams:
    """Cardiac electrical-axis and hemodynamic constants.

    ``inotropy``/``chronotropy`` scale stroke volume and heart rate from the
    stimulatory pathway tone (doc/05 4.3); ``sympathetic_tone`` (1 =
    intact) encodes the *suppressive* sympathetic branch as a fraction of
    baseline beta-adrenergic drive that remains after channel/site blockade —
    it multiplies both heart rate and contractility, so a beta-like blockade
    at IC50 exposure halves both and quarters cardiac output.
    """

    qtc_base_ms: float = 415.0
    delta_qtc_max_ms: float = 40.0
    qtc_b50: float = 0.5
    bpm: float = 70.0
    sv_ml: float | None = None
    r_sys_mmhg_min_l: float | None = None
    c_art_l_mmhg: float = 0.0018
    c_ven_l_mmhg: float = 0.09
    map_target_mmhg: float = 93.0
    cvp_ref_mmhg: float = 5.0
    inotropy: float = 1.0
    chronotropy: float = 1.0
    sympathetic_tone: float = 1.0


def sympathetic_tone_from_emax(conc_free_nm: float, ic50_nm: float) -> float:
    """Residual sympathetic tone from an Emax site-saturation axis.

    A drug that blocks a sympathetic/beta-adrenergic site suppresses the drive
    with the standard saturable relation ``suppression = C/(IC50 + C)`` (Emax
    1), so the *remaining* tone is ``1 - suppression = IC50/(IC50 + C)``.
    ``C = IC50`` halves the tone.  Applied multiplicatively to both heart rate
    and contractility (beta-1 blockade), an IC50-level exposure therefore
    produces a 4-fold cardiac-output reduction.
    """
    if not ic50_nm > 0:
        raise ValueError("ic50_nm must be positive")
    if not conc_free_nm >= 0:
        raise ValueError("conc_free_nm must be non-negative")
    return float(ic50_nm / (ic50_nm + conc_free_nm))


@dataclass(slots=True)
class CardiacResult:
    """QTc trajectory and hemodynamic state variables."""

    t_h: NDArray
    block_frac: NDArray
    ikr_frac: NDArray
    qtc_ms: NDArray
    delta_qtc_ms: NDArray
    tdpr_band: TdPBand
    tdpr_grade: int
    pa_mmhg: NDArray
    pv_mmhg: NDArray
    map_mmhg: float
    cvp_mmhg: float
    co_l_min: float
    hr_bpm: float
    sv_ml: float

    def to_series(self) -> dict[str, NDArray | str | float]:
        return {
            "t_h": self.t_h,
            "block_frac": self.block_frac,
            "ikr_frac": self.ikr_frac,
            "qtc_ms": self.qtc_ms,
            "delta_qtc_ms": self.delta_qtc_ms,
            "pa_mmHg": self.pa_mmhg,
            "pv_mmHg": self.pv_mmhg,
        }


def ikr_fraction(block_frac: NDArray | float, b50: float = 0.5) -> NDArray:
    """Fractional IKr current remaining at hERG blockade ``block_frac``."""
    if b50 <= 0:
        raise ValueError("b50 must be positive")
    b = np.asarray(block_frac, dtype=float)
    return 1.0 - b / (b + b50)


def predict_qtc(
    block_frac: NDArray,
    qtc_base_ms: float,
    delta_max_ms: float,
    b50: float = 0.5,
) -> tuple[NDArray, NDArray]:
    """QTc (ms) and delta-QTc (ms) from the Emax hERG->IKr->QTc relation."""
    if delta_max_ms < 0:
        raise ValueError("delta_max_ms must be non-negative")
    b = np.asarray(block_frac, dtype=float)
    b = np.clip(b, 0.0, 1.0)
    ikr = ikr_fraction(b, b50)
    prolong = delta_max_ms * (1.0 - ikr)
    qtc = qtc_base_ms + prolong
    return qtc, prolong


def tdpr_band(qtc_ms: NDArray) -> tuple[TdPBand, int]:
    """TdP risk band and CTCAE-like grade for a QTc trace."""
    peak = float(np.max(qtc_ms))
    if peak >= 500.0:
        return "high", 3
    if peak >= 480.0:
        return "moderate", 2
    if peak >= 450.0:
        return "low", 1
    return "none", 0


@dataclass(slots=True)
class HemodynamicsResult:
    pa_mmhg: NDArray
    pv_mmhg: NDArray
    map_mmhg: float
    cvp_mmhg: float
    co_l_min: float
    hr_bpm: float
    sv_ml: float


def _systemic_resistance(co_l_min: float, map_target_mmhg: float, cvp_ref_mmhg: float) -> float:
    if co_l_min <= 0:
        raise ValueError("cardiac output must be positive")
    return (map_target_mmhg - cvp_ref_mmhg) / co_l_min


def simulate_hemodynamics(
    physiology: HumanPhysiology,
    params: CardiacParams,
    n_eval: int = 601,
    tmax_h: float = 24.0,
    rtol: float = 1e-8,
    atol: float = 1e-9,
) -> HemodynamicsResult:
    """Two-compartment Windkessel relaxation to steady-state pressures.

    ``sympathetic_tone`` multiplies heart rate and stroke volume, so cardiac
    output scales with tone^2.  The systemic resistance is pinned to the
    *intact-tone* reference output (``co_ref``) rather than the suppressed one:
    a beta-like suppression that lowers CO with an un-compensated vasculature
    reads out as hypotension (MAP - CVP falls with tone^2), the analytic anchor
    for the doc/05 4.3 sympathetic branch.  A stimulatory ERK inotropy lifts
    co_ref itself, so it still sets resistance as before.
    """
    if not 0.0 < params.sympathetic_tone <= 1.0:
        raise ValueError("sympathetic_tone must be in (0, 1]")
    hr_ref = params.bpm * params.chronotropy
    if hr_ref <= 0:
        raise ValueError("heart rate must be positive")
    sv_ref = params.sv_ml if params.sv_ml is not None else physiology.cardiac_output_ml_min / hr_ref
    sv = sv_ref * params.inotropy
    co_ref = hr_ref * sv / 1000.0  # L/min
    hr = hr_ref * params.sympathetic_tone
    sv = sv * params.sympathetic_tone
    co = hr * sv / 1000.0  # L/min
    if co <= 0:
        raise ValueError("cardiac output must be positive")
    r_sys = params.r_sys_mmhg_min_l or _systemic_resistance(
        co_ref, params.map_target_mmhg, params.cvp_ref_mmhg
    )
    if r_sys <= 0 or params.c_art_l_mmhg <= 0 or params.c_ven_l_mmhg <= 0:
        raise ValueError("circulatory constants must be positive")

    pa0 = params.map_target_mmhg
    pv0 = params.cvp_ref_mmhg

    def rhs(sol_t: float, y: NDArray) -> NDArray:
        del sol_t
        pa, pv = float(y[0]), float(y[1])
        q = (pa - pv) / r_sys
        return np.array([(co - q) / params.c_art_l_mmhg, (q - co) / params.c_ven_l_mmhg])

    y0 = np.array([pa0, pv0], dtype=float)
    t_span = (0.0, tmax_h)
    t_eval = np.linspace(0.0, tmax_h, n_eval)
    sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, method="LSODA", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"hemodynamics solve failed: {sol.message}")
    pa = sol.y[0]
    pv = sol.y[1]
    return HemodynamicsResult(
        pa_mmhg=pa,
        pv_mmhg=pv,
        map_mmhg=float(pa[-1]),
        cvp_mmhg=float(pv[-1]),
        co_l_min=co,
        hr_bpm=hr,
        sv_ml=sv,
    )


def simulate_cardiac(
    t_h: NDArray,
    block_frac: NDArray,
    physiology: HumanPhysiology,
    params: CardiacParams | None = None,
    n_eval: int = 601,
) -> CardiacResult:
    """Full Stage-4 cardiac panel: QT axis + lumped hemodynamics."""
    p = params or CardiacParams()
    t = np.asarray(t_h, dtype=float)
    b = np.asarray(block_frac, dtype=float)
    if t.ndim != 1 or b.ndim != 1 or t.shape[0] != b.shape[0]:
        raise ValueError("t_h and block_frac must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")

    qtc, prolong = predict_qtc(
        np.interp(np.linspace(t[0], t[-1], n_eval), t, b),
        p.qtc_base_ms,
        p.delta_qtc_max_ms,
        p.qtc_b50,
    )
    band, grade = tdpr_band(qtc)
    hemo = simulate_hemodynamics(physiology, p, n_eval=n_eval)
    return CardiacResult(
        t_h=np.linspace(t[0], t[-1], n_eval),
        block_frac=np.interp(np.linspace(t[0], t[-1], n_eval), t, b),
        ikr_frac=ikr_fraction(np.interp(np.linspace(t[0], t[-1], n_eval), t, b), p.qtc_b50),
        qtc_ms=qtc,
        delta_qtc_ms=prolong,
        tdpr_band=band,
        tdpr_grade=grade,
        pa_mmhg=hemo.pa_mmhg,
        pv_mmhg=hemo.pv_mmhg,
        map_mmhg=hemo.map_mmhg,
        cvp_mmhg=hemo.cvp_mmhg,
        co_l_min=hemo.co_l_min,
        hr_bpm=hemo.hr_bpm,
        sv_ml=hemo.sv_ml,
    )


__all__ = [
    "CardiacParams",
    "CardiacResult",
    "HemodynamicsResult",
    "ikr_fraction",
    "predict_qtc",
    "simulate_cardiac",
    "simulate_hemodynamics",
    "sympathetic_tone_from_emax",
    "tdpr_band",
]
