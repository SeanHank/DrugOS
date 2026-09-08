"""Biomarker translation and CTCAE-style grading (Stage 5, D18).

Maps Stage-4 organ outputs onto clinical-grade biomarkers with units,
reference ranges and a 0..4 severity ladder per the CTCAE convention
(doc/05 5.1, doc/03 2.6).  Each ``BiomarkerSpec`` carries the four grade
thresholds (worse-than rows), the direction in which worse matters, and the
normal reference range.  ``grade_timeseries`` adds time-to-onset and duration
from the organ trajectory, which the clinical report consumes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from drugos.organ.base import NDArray

SEVERITY_LABELS: tuple[str, ...] = (
    "normal",
    "mild",
    "moderate",
    "severe",
    "life-threatening",
)


@dataclass(frozen=True, slots=True)
class BiomarkerSpec:
    """Grading rule for one biomarker."""

    name: str
    unit: str
    ref_lo: float
    ref_hi: float
    worse: str = "higher"  # "higher" or "lower"
    thresholds: tuple[float, float, float, float] = (1.0, 2.0, 3.0, 10.0)

    def __post_init__(self) -> None:
        if self.worse not in ("higher", "lower"):
            raise ValueError("worse must be 'higher' or 'lower'")
        if len(self.thresholds) != 4:
            raise ValueError("thresholds must contain exactly four values")


@dataclass(frozen=True, slots=True)
class BiomarkerGrade:
    """One graded biomarker outcome."""

    name: str
    value: float
    unit: str
    grade: int
    severity: str
    ref_lo: float
    ref_hi: float
    onset_h: float | None = None
    duration_h: float | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "value": round(self.value, 3),
            "unit": self.unit,
            "grade": self.grade,
            "severity": self.severity,
            "ref_lo": self.ref_lo,
            "ref_hi": self.ref_hi,
            "onset_h": self.onset_h,
            "duration_h": self.duration_h,
        }


def grade_value(value: float, spec: BiomarkerSpec) -> int:
    """CTCAE-style grade of a single value against ``spec``."""
    if spec.worse == "higher":
        return min(4, sum(1 for t in spec.thresholds if value >= t))
    return min(4, sum(1 for t in spec.thresholds if value <= t))


def peak_fraction(t_h: NDArray, values: NDArray, spec: BiomarkerSpec) -> tuple[float, float]:
    """Peak value and its time (h) for the worse direction of ``spec``."""
    if spec.worse == "higher":
        idx = int(np.argmax(values))
    else:
        idx = int(np.argmin(values))
    return float(values[idx]), float(t_h[idx])


def crossing_interval(
    t_h: NDArray, values: NDArray, spec: BiomarkerSpec, grade: int = 1
) -> tuple[float | None, float | None]:
    """Time-to-onset and duration of ``grade`` (or worse) on the trajectory."""
    thr = spec.thresholds[grade - 1]
    if spec.worse == "higher":
        active = values >= thr
    else:
        active = values <= thr
    active = np.asarray(active, dtype=bool)
    if not bool(active.any()):
        return None, None
    t = np.asarray(t_h, dtype=float)
    onset = float(t[np.argmax(active)])
    last = float(t[len(t) - 1 - int(np.argmax(active[::-1]))])
    return onset, max(0.0, last - onset)


def grade_timeseries(t_h: NDArray, values: NDArray, spec: BiomarkerSpec) -> BiomarkerGrade:
    """Grade the trajectory peak and attach onset/duration windows."""
    value, _ = peak_fraction(t_h, values, spec)
    g = grade_value(value, spec)
    onset, duration = crossing_interval(t_h, values, spec, max(g, 1))
    return BiomarkerGrade(
        name=spec.name,
        value=value,
        unit=spec.unit,
        grade=g,
        severity=SEVERITY_LABELS[g],
        ref_lo=spec.ref_lo,
        ref_hi=spec.ref_hi,
        onset_h=onset,
        duration_h=duration,
    )


def grade_absolute(value: float, spec: BiomarkerSpec) -> BiomarkerGrade:
    """Grading of a scalar (peak-anchored) value without a trajectory."""
    g = grade_value(value, spec)
    return BiomarkerGrade(
        name=spec.name,
        value=value,
        unit=spec.unit,
        grade=g,
        severity=SEVERITY_LABELS[g],
        ref_lo=spec.ref_lo,
        ref_hi=spec.ref_hi,
    )


# ---------------------------------------------------------------------------
# Baseline biomarker catalogue (units + reference ranges)
# ---------------------------------------------------------------------------
BIOMARKERS: dict[str, BiomarkerSpec] = {
    "ALT": BiomarkerSpec(
        "ALT (alanine transaminase)", "xULN", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0)
    ),
    "AST": BiomarkerSpec(
        "AST (aspartate transaminase)", "xULN", 0.0, 1.0, "higher", (1.0, 2.0, 3.0, 10.0)
    ),
    "total_bilirubin": BiomarkerSpec(
        "total bilirubin", "xULN", 0.0, 1.0, "higher", (1.0, 1.5, 2.0, 5.0)
    ),
    "QTc": BiomarkerSpec(
        "QTc (Fridericia)", "ms", 350.0, 450.0, "higher", (450.0, 480.0, 500.0, 550.0)
    ),
    "delta_QTc": BiomarkerSpec("Delta QTc", "ms", -10.0, 10.0, "higher", (20.0, 30.0, 60.0, 100.0)),
    "heart_rate": BiomarkerSpec(
        "heart rate", "bpm", 60.0, 100.0, "higher", (100.0, 120.0, 150.0, 200.0)
    ),
    "map": BiomarkerSpec(
        "mean arterial pressure", "mmHg", 70.0, 105.0, "lower", (70.0, 65.0, 60.0, 50.0)
    ),
    "gfr": BiomarkerSpec("GFR", "mL/min", 90.0, 140.0, "lower", (90.0, 60.0, 45.0, 15.0)),
    "creatinine": BiomarkerSpec(
        "serum creatinine", "xULN", 0.0, 1.0, "higher", (1.5, 2.0, 3.0, 4.0)
    ),
    "KIM_1": BiomarkerSpec("KIM-1 (urinary)", "xUNL", 0.0, 1.0, "higher", (1.5, 3.0, 5.0, 10.0)),
    "CNS_exposure": BiomarkerSpec(
        "CNS exposure ratio", "Cmax/IC50", 0.0, 0.1, "higher", (0.1, 0.32, 1.0, 3.2)
    ),
}


def summarize_liver(
    t_h: NDArray, alt_uln: NDArray, ast_uln: NDArray, tbili_uln: NDArray
) -> list[BiomarkerGrade]:
    """Grade liver trajectory rows (ALT/AST/total bilirubin)."""
    return [
        grade_timeseries(t_h, alt_uln, BIOMARKERS["ALT"]),
        grade_timeseries(t_h, ast_uln, BIOMARKERS["AST"]),
        grade_timeseries(t_h, tbili_uln, BIOMARKERS["total_bilirubin"]),
    ]


def summarize_cardiac(
    t_h: NDArray,
    qtc_ms: NDArray,
    heart_rate_bpm: float,
    map_mmhg: float,
) -> list[BiomarkerGrade]:
    """Grade cardiovascular rows (QTc trajectory + steady-state HR/MAP)."""
    qtc = grade_timeseries(t_h, qtc_ms, BIOMARKERS["QTc"])
    hr = grade_absolute(heart_rate_bpm, BIOMARKERS["heart_rate"])
    bp = grade_absolute(map_mmhg, BIOMARKERS["map"])
    return [qtc, hr, bp]


def summarize_kidney(
    t_h: NDArray,
    gfr_ml_min: NDArray,
    scr_ratio: NDArray,
    kim1_xunl: NDArray | None = None,
) -> list[BiomarkerGrade]:
    """Grade kidney rows (GFR trajectory + creatinine ratio + KIM-1)."""
    gfr = grade_timeseries(t_h, gfr_ml_min, BIOMARKERS["gfr"])
    scr = grade_timeseries(t_h, scr_ratio, BIOMARKERS["creatinine"])
    rows = [gfr, scr]
    if kim1_xunl is not None:
        rows.append(grade_timeseries(t_h, kim1_xunl, BIOMARKERS["KIM_1"]))
    return rows


__all__ = [
    "BIOMARKERS",
    "BiomarkerGrade",
    "BiomarkerSpec",
    "SEVERITY_LABELS",
    "crossing_interval",
    "grade_absolute",
    "grade_timeseries",
    "grade_value",
    "peak_fraction",
    "summarize_cardiac",
    "summarize_kidney",
    "summarize_liver",
]
