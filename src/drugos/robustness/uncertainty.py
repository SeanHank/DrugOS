"""D21 — Parameter-ensemble uncertainty propagation (doc/07 Phase 6).

A multiplicative parameter ensemble perturbs the empirically-fitted PK and
in-vitro potency inputs of a ``RunSpec`` (unbound clearance split, plasma free
fraction, hERG / DILI / CNS potencies) with independent log-normals at a fixed
coefficient of variation.  Every ensemble member is a full pipeline run; the
collected trajectories give 90% percentile bands for plasma/organ readouts and
quantiles for the endpoint risks, surfacing how much the Stage-5 verdict moves
under input uncertainty.  Deterministic via ``seed``, so results are
reproducible.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from drugos.organ.base import NDArray
from drugos.pipeline import RunResult, RunSpec, run_pipeline

DEFAULT_KEYS: tuple[str, ...] = (
    "cl_hep_l_h",
    "cl_renal_l_h",
    "fup",
    "qt_ic50_nm",
    "dili_ic50_nm",
    "cns_ic50_nm",
)

CURVE_NAMES: tuple[str, ...] = (
    "plasma_total_mg_l",
    "plasma_free_mg_l",
    "alt_u_l",
    "bilirubin_mg_dl",
    "qtc_ms",
    "delta_qtc_ms",
    "gfr_ml_min",
    "scr_ratio",
    "brain_free_nm",
)


@dataclass(frozen=True, slots=True)
class EnsembleConfig:
    """Configuration of the multiplicative parameter ensemble."""

    n_runs: int = 14
    seed: int = 7
    cv: float = 0.30
    keys: tuple[str, ...] = DEFAULT_KEYS

    def __post_init__(self) -> None:
        if self.n_runs < 1:
            raise ValueError("n_runs must be >= 1")
        if self.cv <= 0:
            raise ValueError("cv must be positive")


def sample_multipliers(config: EnsembleConfig) -> list[dict[str, float]]:
    """Deterministic log-normal multiplier sets, one per ensemble member."""
    rng = np.random.default_rng(config.seed)
    matrix: NDArray = rng.lognormal(0.0, config.cv, size=(config.n_runs, len(config.keys)))
    return [
        {name: float(matrix[i, j]) for j, name in enumerate(config.keys)}
        for i in range(config.n_runs)
    ]


def apply_multipliers(spec: RunSpec, multipliers: Mapping[str, float]) -> RunSpec:
    """Copy ``spec`` with each PK/potency scalar scaled by its multiplier."""
    import copy

    changed: dict[str, float | None] = {}
    for key, mult in multipliers.items():
        if not math.isfinite(mult) or mult <= 0:
            raise ValueError(f"multiplier for '{key}' must be finite and positive")
        value = getattr(spec, key)
        if value is None:
            continue
        scaled = value * mult
        if key == "fup":
            scaled = max(1.0e-4, min(1.0, scaled))
        else:
            scaled = max(1.0e-4, scaled)
        changed[key] = float(scaled)

    out = copy.copy(spec)
    for key, value in changed.items():
        object.__setattr__(out, key, value)
    return out


def curves_from_result(result: RunResult) -> dict[str, np.ndarray]:
    """Time-series readouts extracted from one pipeline run (common grid)."""
    return {
        "plasma_total_mg_l": np.asarray(result.pk.plasma_total, dtype=float),
        "plasma_free_mg_l": np.asarray(result.pk.plasma_free, dtype=float),
        "alt_u_l": np.asarray(result.organ.liver.alt_u_l, dtype=float),
        "bilirubin_mg_dl": np.asarray(result.organ.liver.bilirubin_mg_dl, dtype=float),
        "qtc_ms": np.asarray(result.organ.cardiac.qtc_ms, dtype=float),
        "delta_qtc_ms": np.asarray(result.organ.cardiac.delta_qtc_ms, dtype=float),
        "gfr_ml_min": np.asarray(result.organ.kidney.gfr_ml_min, dtype=float),
        "scr_ratio": np.asarray(result.organ.kidney.scr_ratio, dtype=float),
        "brain_free_nm": np.asarray(result.organ.cns.brain_free_nm, dtype=float),
    }


def points_from_result(result: RunResult) -> dict[str, float]:
    """Scalar endpoint outputs of one run, keyed by stable names."""
    points: dict[str, float] = {
        "peak_qtc_ms": float(np.max(result.organ.cardiac.qtc_ms)),
        "peak_alt_uln": float(result.organ.liver.peak_alt_uln),
        "overall_risk": result.toxicity.overall_risk(),
    }
    for risk in result.toxicity.risks:
        points[risk.endpoint.value] = risk.risk
        points[f"{risk.endpoint.value}_grade"] = risk.grade
    return points


@dataclass(slots=True)
class EnsembleResult:
    """Ensemble of pipeline runs plus percentile/determinism summaries."""

    config: EnsembleConfig
    t_h: np.ndarray
    curves: dict[str, np.ndarray]
    points: dict[str, np.ndarray]
    verdicts: list[str]

    def band(self, name: str, qs: tuple[int, int, int] = (5, 50, 95)) -> dict[str, np.ndarray]:
        """Location + 90% percentile band for curve ``name`` on the t grid."""
        raw = self.curves[name]
        return {
            "t_h": self.t_h,
            "lo": np.percentile(raw, qs[0], axis=0),
            "med": np.percentile(raw, qs[1], axis=0),
            "hi": np.percentile(raw, qs[2], axis=0),
        }

    def point_quantiles(
        self, name: str, qs: tuple[int, int, int] = (5, 50, 95)
    ) -> tuple[float, float, float]:
        """(qlo, med, qhi) of a scalar endpoint across the ensemble."""
        qlo, med, qhi = np.quantile(self.points[name], [q / 100.0 for q in qs])
        return float(qlo), float(med), float(qhi)

    def verdict_counts(self) -> dict[str, int]:
        """Verdict string -> number of ensemble members hitting it."""
        counts: dict[str, int] = {}
        for verdict in self.verdicts:
            counts[verdict] = counts.get(verdict, 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: kv[1], reverse=True))

    def to_dict(self) -> dict[str, object]:
        """Compact JSON-friendly summary for reports/playground."""
        point_q: dict[str, object] = {}
        for name in self.points:
            point_q[name] = [round(v, 4) for v in self.point_quantiles(name)]
        curve_q: dict[str, object] = {}
        for name in self.curves:
            band = self.band(name)
            step = max(1, len(self.t_h) // 60)
            curve_q[name] = {
                "t_h": [round(float(t), 3) for t in band["t_h"][::step]],
                "lo": [round(float(v), 4) for v in band["lo"][::step]],
                "med": [round(float(v), 4) for v in band["med"][::step]],
                "hi": [round(float(v), 4) for v in band["hi"][::step]],
            }
        return {
            "n_runs": self.config.n_runs,
            "cv": self.config.cv,
            "seed": self.config.seed,
            "points_q5_q50_q95": point_q,
            "band_90": curve_q,
            "verdict_counts": self.verdict_counts(),
        }


def run_uncertainty(spec: RunSpec, config: EnsembleConfig | None = None) -> EnsembleResult:
    """Run the parameter ensemble over ``spec`` and summarize the bands."""
    config = config or EnsembleConfig()
    multipliers = sample_multipliers(config)
    results = [run_pipeline(apply_multipliers(spec, m)) for m in multipliers]
    base = results[0]
    t_h = np.asarray(base.pk.t, dtype=float)
    curves = {
        name: np.stack([curves_from_result(r)[name] for r in results]) for name in CURVE_NAMES
    }
    point_rows: dict[str, list[float]] = {name: [] for name in points_from_result(base)}
    for r in results:
        for name, value in points_from_result(r).items():
            point_rows[name].append(value)
    points = {name: np.asarray(values, dtype=float) for name, values in point_rows.items()}
    verdicts = [r.verdict for r in results]
    return EnsembleResult(config=config, t_h=t_h, curves=curves, points=points, verdicts=verdicts)


__all__ = [
    "CURVE_NAMES",
    "DEFAULT_KEYS",
    "EnsembleConfig",
    "EnsembleResult",
    "apply_multipliers",
    "curves_from_result",
    "points_from_result",
    "run_uncertainty",
    "sample_multipliers",
]
