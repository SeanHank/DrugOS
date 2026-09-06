"""D23 — Sensitivity analysis and driver attribution (doc/07 Phase 6).

Two complementary tools over the pipeline's PK/potency inputs (the same
scalars the D21 ensemble perturbs):

* ``local_sensitivity`` — one-at-a-time +/- ``delta`` multiplicative
  perturbation; reports normalised log-sensitivities d ln y / d ln p and the
  ranked driver list for a chosen endpoint output.
* ``run_sobol_sensitivity`` — Saltelli first/total Sobol indices computed on
  an independent evaluator (default: full pipeline runs) with a scrambled
  LD-64 Sobol design (``scipy.stats.qmc``).  The index math is pure and
  unit-testable via ``saltelli_design`` / ``first_total_indices``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass

import numpy as np
from scipy.stats.qmc import Sobol

from drugos.pipeline import RunSpec, run_pipeline
from drugos.robustness.uncertainty import DEFAULT_KEYS, apply_multipliers, points_from_result


@dataclass(frozen=True, slots=True)
class SensitivityConfig:
    """Shared configuration for the OAT and Sobol analyses."""

    keys: tuple[str, ...] = DEFAULT_KEYS
    delta: float = 0.10
    seed: int = 3
    sobol_n: int = 8
    sobol_outputs: tuple[str, ...] = ("qt", "dili", "overall_risk")

    def __post_init__(self) -> None:
        if not (0.0 < self.delta < 1.0):
            raise ValueError("delta must be in (0, 1)")
        if self.sobol_n < 1:
            raise ValueError("sobol_n must be >= 1")
        if not self.sobol_outputs:
            raise ValueError("sobol_outputs must not be empty")


@dataclass(slots=True)
class LocalSensitivityResult:
    """Normalised one-at-a-time sensitivities for one scalar output."""

    config: SensitivityConfig
    output: str
    baseline: float
    values: dict[str, float]
    drivers: list[tuple[str, float]]


def _scalar_outputs(spec: RunSpec) -> dict[str, float]:
    return points_from_result(run_pipeline(spec))


def local_sensitivity(
    spec: RunSpec,
    output: str = "overall_risk",
    config: SensitivityConfig | None = None,
    evaluator: Callable[[RunSpec], Mapping[str, float]] = _scalar_outputs,
) -> LocalSensitivityResult:
    """OAT +/- delta sensitivities of ``output`` over the config keys."""
    config = config or SensitivityConfig()
    points = dict(evaluator(spec))
    if output not in points:
        raise KeyError(f"sensitivity output '{output}' not produced by evaluator")
    baseline = float(points[output])
    denominator = 2.0 * config.delta * (baseline if baseline != 0.0 else 1.0)
    values: dict[str, float] = {}
    for key in config.keys:
        hi = dict(evaluator(apply_multipliers(spec, {key: 1.0 + config.delta})))
        lo = dict(evaluator(apply_multipliers(spec, {key: 1.0 - config.delta})))
        values[key] = (float(hi[output]) - float(lo[output])) / denominator
    drivers = sorted(values.items(), key=lambda kv: abs(kv[1]), reverse=True)
    return LocalSensitivityResult(
        config=config, output=output, baseline=baseline, values=values, drivers=drivers
    )


def saltelli_design(n: int, d: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Independent base matrices A and B in [0,1)^(n x d)."""
    if n < 1 or d < 1:
        raise ValueError("n and d must be >= 1")
    a = Sobol(d=d, scramble=True, seed=seed).random(n)
    b = Sobol(d=d, scramble=True, seed=seed + 1).random(n)
    return np.asarray(a, dtype=float), np.asarray(b, dtype=float)


def first_total_indices(
    ya: np.ndarray, yb: np.ndarray, yc: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Saltelli first-order / total-indices estimators from Ya, Yb, Yc.

    ``yc`` is (d, n); row ``i`` is the recombination of A with column ``i`` of
    B.  Returns (s1, st) each of shape (d,) using the variance of
    ``cat(Ya, Yb)``; a degenerate (zero-variance) output yields zeros.
    """
    ya = np.asarray(ya, dtype=float)
    yb = np.asarray(yb, dtype=float)
    yc = np.asarray(yc, dtype=float)
    stacked = np.concatenate([ya, yb])
    var = float(np.var(stacked))
    if var < 1e-12 or yc.ndim != 2 or yc.shape[1] != len(ya):
        return np.zeros(yc.shape[0]), np.zeros(yc.shape[0])
    n = len(ya)
    s1 = np.zeros(yc.shape[0])
    st = np.zeros(yc.shape[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(yc.shape[0]):
            s1[i] = float(np.sum(yb * (ya - yc[i])) / (n * var))
            st[i] = float(np.sum((ya - yc[i]) ** 2) / (2.0 * n * var))
    return np.clip(s1, -1.0, 1.0), np.clip(st, 0.0, 1.0)


@dataclass(slots=True)
class SobolSensitivityResult:
    """First/total Sobol indices per scalar output, aligned to config keys."""

    config: SensitivityConfig
    design_n: int
    indices: dict[str, tuple[np.ndarray, np.ndarray]]

    def first_drivers(self, output: str) -> list[tuple[str, float]]:
        s1, _ = self.indices[output]
        ranked = sorted(
            zip(self.config.keys, s1, strict=True), key=lambda kv: abs(kv[1]), reverse=True
        )
        return [(name, float(value)) for name, value in ranked]

    def total_drivers(self, output: str) -> list[tuple[str, float]]:
        _, st = self.indices[output]
        ranked = sorted(
            zip(self.config.keys, st, strict=True), key=lambda kv: abs(kv[1]), reverse=True
        )
        return [(name, float(value)) for name, value in ranked]

    def to_dict(self) -> dict[str, object]:
        return {
            "design_n": self.design_n,
            "keys": list(self.config.keys),
            "indices": {
                output: {
                    "first": [float(v) for v in s1],
                    "total": [float(v) for v in st],
                }
                for output, (s1, st) in self.indices.items()
            },
        }


def _design_points(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """All Saltelli design rows: (A-rows, B-rows, recombined C rows)."""
    n, d = a.shape
    c = np.empty((d * n, d), dtype=float)
    block = np.empty((d, n, d), dtype=float)
    for i in range(d):
        recombined = a.copy()
        recombined[:, i] = b[:, i]
        block[i] = recombined
    c[:] = block.reshape(d * n, d)
    return a, b, c


def _default_evaluator(spec: RunSpec) -> Callable[[dict[str, float]], Mapping[str, float]]:
    return lambda m: dict(points_from_result(run_pipeline(apply_multipliers(spec, m))))


def run_sobol_sensitivity(
    spec: RunSpec,
    config: SensitivityConfig | None = None,
    evaluator: Callable[[dict[str, float]], Mapping[str, float]] | None = None,
) -> SobolSensitivityResult:
    """Sobol first/total indices over the config keys and outputs.

    The default evaluator maps a multiplier set to endpoint outputs by running
    the full pipeline; tests may inject a cheap surrogate.
    """
    config = config or SensitivityConfig()
    d = len(config.keys)
    n = config.sobol_n
    a, b = saltelli_design(n, d, config.seed)
    _, _, c = _design_points(a, b)

    if evaluator is None:
        eval_cb = _default_evaluator(spec)
    else:
        eval_cb = evaluator

    a_rows = [{config.keys[j]: float(a[i, j]) for j in range(d)} for i in range(n)]
    b_rows = [{config.keys[j]: float(b[i, j]) for j in range(d)} for i in range(n)]
    c_rows = [{config.keys[j]: float(c[i, j]) for j in range(d)} for i in range(d * n)]

    indices: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for output in config.sobol_outputs:
        ya = np.asarray([float(eval_cb(m)[output]) for m in a_rows])
        yb = np.asarray([float(eval_cb(m)[output]) for m in b_rows])
        yc = np.asarray(
            [[float(eval_cb(m)[output]) for m in c_rows[i * n : (i + 1) * n]] for i in range(d)]
        )
        indices[output] = first_total_indices(ya, yb, yc)
    return SobolSensitivityResult(config=config, design_n=n, indices=indices)


__all__ = [
    "LocalSensitivityResult",
    "SensitivityConfig",
    "SobolSensitivityResult",
    "first_total_indices",
    "local_sensitivity",
    "run_sobol_sensitivity",
    "saltelli_design",
]
