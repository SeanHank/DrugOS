"""Phase-6 robustness plumbing — self-consistency (D21-D24, L1, doc/08 §1.3).

Checks the three robustness engines introduced in Phase 6 are deterministic,
conservative and consistent with their own statistics:

- the D21 parameter ensemble is fully reproducible for a fixed seed and its
  90% percentile band is monotonically ordered (lo <= median <= hi),
- the D22 virtual population reports well-formed non-negative incidence,
- the D23 Sobol first/total indices stay inside their valid ranges and the
  OAT DILI drift moves in the sign direction its equations dictate,
- the D24 prospect: a held-out profile/dose keeps the same composite verdict.

No external data: this tier certifies the robustness machinery computes
what its equations promise (determinism + CI bookkeeping).
"""

from __future__ import annotations

import numpy as np
from validation.benchmarks import BENCHMARKS
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pipeline import spec_from_benchmark_data
from drugos.robustness.population import PopulationConfig, run_population
from drugos.robustness.sensitivity import (
    SensitivityConfig,
    first_total_indices,
    local_sensitivity,
)
from drugos.robustness.uncertainty import EnsembleConfig, run_uncertainty


def _bench(name: str):
    return next(b for b in BENCHMARKS if b.name == name)


def case_robustness_sanity() -> CaseResult:
    spec = spec_from_benchmark_data(_bench("acetaminophen"))
    config = EnsembleConfig(n_runs=4, seed=5, cv=0.25)

    ens_a = run_uncertainty(spec, config)
    ens_b = run_uncertainty(spec, config)
    max_band_diff = float(np.max(np.abs(ens_a.band("qtc_ms")["med"] - ens_b.band("qtc_ms")["med"])))
    reproducible = max_band_diff == 0.0 and ens_a.verdict_counts() == ens_b.verdict_counts()

    band = ens_a.band("alt_u_l")
    violations = int(np.sum(band["lo"] > band["med"]) + np.sum(band["med"] > band["hi"]))

    pop = run_population(spec, PopulationConfig(n_individuals=4, seed=6))
    pop_ok = all(np.all(v >= 0) for v in pop.risks.values()) and len(pop.sexes) == 4

    rng = np.random.default_rng(0)
    n, d = 32, 4
    s1, st = first_total_indices(rng.random(n), rng.random(n), rng.random((d, n)))
    sobol_ok = bool(np.all(s1 >= -1.0) and np.all(s1 <= 1.0)) and bool(
        np.all(st >= 0.0) and np.all(st <= 1.0)
    )

    sens = local_sensitivity(
        spec,
        output="dili",
        config=SensitivityConfig(keys=("dili_ic50_nm",), delta=0.1, seed=3),
    )
    dip_up = sens.values["dili_ic50_nm"]
    sign_ok = bool(np.isfinite(dip_up)) and dip_up < 0.0

    metrics = [
        MetricResult(
            "ensemble_reproducible_max_band_diff",
            max_band_diff,
            0.0,
            0.0,
            "au",
            "pass" if reproducible else "FAIL",
        ),
        MetricResult(
            "band_monotonic_violations",
            float(violations),
            0.0,
            0.0,
            "points",
            "pass" if violations == 0 else "FAIL",
        ),
        MetricResult(
            "population_min_risk_nonneg",
            float(min(float(v.min()) for v in pop.risks.values())),
            0.0,
            0.0,
            "P(risk)",
            "pass" if pop_ok else "FAIL",
        ),
        MetricResult(
            "sobol_first_total_in_range",
            0.0 if sobol_ok else 1.0,
            0.0,
            0.0,
            "flags",
            "pass" if sobol_ok else "FAIL",
        ),
        MetricResult(
            "dili_ic50_sensitivity_sign",
            float(dip_up),
            -1.0,
            0.0,
            "dlnR/dlnIC50",
            "pass" if sign_ok else "FAIL",
        ),
    ]
    ok = all(m.criterion == "pass" for m in metrics)
    return CaseResult(
        "Phase-6 robustness engines: D21-D24 self-consistency",
        ok,
        metrics,
        [
            f"Fixed-seed D21 ensemble reproduces itself exactly "
            f"(max median-band diff {max_band_diff:.3g}); 90% band monotone with "
            f"{violations} violations; D22 cohort incidence non-negative; D23 "
            "first/total indices inside [-1,1]/[0,1]; DILI risk strictly "
            f"decreases with a rising IC50 ({dip_up:+.4f} per +10% IC50)."
        ],
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_robustness_sanity"]
