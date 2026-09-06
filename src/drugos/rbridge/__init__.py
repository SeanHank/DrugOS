"""R bridge: literature-PK cross-check via rpy2 (required; doc/06 §1).

Every pipeline run re-computes AUC/CL/t½ with the *literature source code*
``rbridge/literature_pk.R`` (Wagner 1976; Gibaldi & Perrier 1982; Rowland &
Tozer 2010) and reports agreement with the in-repo numpy estimator in the
contract under ``r_verify``.  R is a hard runtime requirement: when rpy2 or R
is unavailable, ``verify_pk`` raises ``RuntimeError`` (no silent fallback —
G5).
"""

from __future__ import annotations

import contextlib
import io
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

_R_SOURCE = Path(__file__).with_name("literature_pk.R")
_r_loaded: bool = False
_r_imported: Any = None


@dataclass(frozen=True, slots=True)
class RVerify:
    """Result of the R literature-PK cross-check for one run."""

    r_cl_l_h: float
    py_cl_l_h: float
    r_auc_inf_mg_h_l: float
    r_t_half_beta_h: float
    agreement_frac: float
    fit: str
    verdict: str
    note: str

    def is_finite(self) -> bool:
        """True when the central (CL/AUC) R estimates are finite.

        ``r_t_half_beta_h`` is allowed to be NaN/absent: some valid profiles
        (e.g. dofetilide at 0.5 mg) have already decayed past the terminal
        log-linear measurement gate by the end of the dosing window, so no
        lambda_z exists — that is *not* a bridge failure, CL/AUC remain
        well-defined.  Only CL/AUC being non-finite is a hard bridge error.
        """
        return math.isfinite(self.r_cl_l_h) and math.isfinite(self.r_auc_inf_mg_h_l)

    def to_dict(self) -> dict[str, float | str | None]:
        """Serialize for JSON contracts; non-finite tails become ``None``.

        ``r_t_half_beta_h`` is NaN when no terminal log-linear tail exists
        (the data already decayed past the measurement gate, e.g. dofetilide
        at 0.5 mg).  NaN must never reach ``json.dumps`` (produces an invalid
        JSON ``NaN`` literal that browsers' ``JSON.parse`` rejects), so every
        non-finite float is dropped to ``None``.
        """

        def _round(v: float) -> float | None:
            return None if not math.isfinite(v) else round(float(v), 4)

        return {
            "r_cl_l_h": _round(self.r_cl_l_h),
            "py_cl_l_h": _round(self.py_cl_l_h),
            "r_auc_inf_mg_h_l": _round(self.r_auc_inf_mg_h_l),
            "r_t_half_beta_h": _round(self.r_t_half_beta_h),
            "agreement_frac": round(self.agreement_frac, 6),
            "fit": self.fit,
            "verdict": self.verdict,
            "note": self.note,
        }


def _estimator_py(
    t: NDArray[np.float64], c: NDArray[np.float64], dose: float
) -> tuple[float, float, float]:
    """numpy twin of ``verify_pk``: trapezoid AUC + tail-corrected CL.

    Mirrors the R estimator step-for-step so their agreement is an identity
    check of two independent implementations of the same literature rule.
    """
    if dose <= 0:
        raise ValueError("dose must be positive")
    keep = np.isfinite(t) & np.isfinite(c) & (t >= 0) & (c >= 0)
    t = np.asarray(t[keep], dtype=float)
    c = np.asarray(c[keep], dtype=float)
    n = t.size
    if n < 4:
        raise ValueError("at least 4 finite time-points are required")
    cmax = float(np.max(c)) if n else 0.0
    auc_t = float(np.sum(0.5 * (c[1:] + c[:-1]) * np.diff(t))) if n >= 2 else 0.0
    sel = np.where((t >= 0.75 * float(np.max(t))) & (c > 0.05 * cmax))[0]
    lam = math.nan
    auc_inf = auc_t
    tail_ok = bool(sel.size >= 3 and c[int(sel[-1])] > 0)
    if tail_ok:
        beta = np.polyfit(t[sel], np.log(np.maximum(c[sel], 1e-300)), 1)
        lam = float(-beta[0])
        auc_inf = auc_t + float(c[int(sel[-1])]) / lam
    return auc_inf, dose / auc_inf, lam


def _import_r() -> Any:
    """Import rpy2 once; raise RuntimeError with guidance when R is missing."""
    global _r_imported
    if _r_imported is not None:
        return _r_imported
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            import rpy2.robjects as ro

        _r_imported = ro
        return ro
    except Exception as exc:  # ImportError, OSError (missing libR), SystemExit
        raise RuntimeError(
            "R bridge unavailable: rpy2/R is a hard requirement (doc/06 §1); "
            "install R and 'pip install rpy2'."
        ) from exc


def _run_r_fit(t: NDArray[np.float64], c: NDArray[np.float64], dose: float) -> dict[str, float]:
    """Evaluate the literature R model and return its scalar estimates."""
    ro = _import_r()
    if not _r_loaded:
        ro.r.source(str(_R_SOURCE))
        globals()["_r_loaded"] = True
    vec = ro.r["verify_pk"](
        ro.FloatVector(t.tolist()),
        ro.FloatVector(c.tolist()),
        float(dose),
    )
    names = list(vec.names)
    return {name: float(vec.rx2(name)[0]) for name in names}


def verify_pk(t: NDArray[np.float64], c: NDArray[np.float64], dose: float) -> RVerify:
    """Run the literature R estimator and cross-check the numpy twin."""
    if dose <= 0:
        raise ValueError("dose must be positive")
    if t.size != c.size or t.size < 4:
        raise ValueError("t/c must be equal-length with at least 4 points")

    auc_py, cl_py, _lam = _estimator_py(t, c, dose)
    rdict = _run_r_fit(t, c, dose)
    r_auc = rdict["auc_inf_mg_hl"]
    r_cl = rdict["cl_l_h"]
    r_t_half = rdict["t_half_beta_h"]
    two_comp = rdict["two_comp"]
    # Both estimators can legitimately yield a non-positive CL (or a missing
    # terminal tail) when the dosing window does not span the full peak-and-
    # decay profile (e.g. t_max shorter than t1/2: the curve is still rising).
    # The bridge's contract is *agreement* between two independent
    # implementations of the same literature rule, not that the estimate be
    # physician-plausible — so only a non-finite central estimate is a hard
    # error.
    if not (math.isfinite(r_auc) and math.isfinite(r_cl) and math.isfinite(cl_py)):
        raise RuntimeError("R bridge returned non-finite central estimates")
    agreement = abs(r_cl - cl_py) / abs(cl_py) if cl_py != 0 else float("inf")
    fit = "two_comp" if two_comp > 0.5 else "loglinear"
    if agreement <= 0.02:
        verdict, note = "r:agree", "R literature estimator within 2% of the numpy twin"
    else:
        verdict, note = (
            "r:mismatch",
            (f"R literature estimator diverged from the numpy twin ({agreement:.1%})"),
        )
    return RVerify(
        r_cl_l_h=r_cl,
        py_cl_l_h=cl_py,
        r_auc_inf_mg_h_l=r_auc,
        r_t_half_beta_h=r_t_half,
        agreement_frac=agreement,
        fit=fit,
        verdict=verdict,
        note=note,
    )


__all__ = [
    "RVerify",
    "_estimator_py",
    "_import_r",
    "_run_r_fit",
    "verify_pk",
]
