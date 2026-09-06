"""D22 — Virtual-patient population simulation (doc/07 Phase 6).

Samples realistic anthropometric profiles (sex, age, height, weight via a
BMI distribution) and runs the full pipeline per individual at the same
absolute dose regimen, reporting population-level endpoint incidence and risk
quantiles.  Deterministic through ``seed``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from drugos.inputs.models import HumanProfile, Sex
from drugos.pipeline import RunSpec, run_pipeline


@dataclass(frozen=True, slots=True)
class PopulationConfig:
    """Sampling configuration for the virtual cohort."""

    n_individuals: int = 16
    seed: int = 11
    age_range: tuple[float, float] = (18.0, 80.0)
    bmi_mean: float = 26.0
    bmi_sd: float = 4.0

    def __post_init__(self) -> None:
        if self.n_individuals < 1:
            raise ValueError("n_individuals must be >= 1")
        lo, hi = self.age_range
        if not (0.0 <= lo < hi <= 100.0):
            raise ValueError("age_range must satisfy 0 <= lo < hi <= 100")
        if self.bmi_mean <= 12.0 or self.bmi_sd <= 0.0:
            raise ValueError("bmi_mean and bmi_sd must be realistic and positive")


@dataclass(slots=True)
class PopulationResult:
    """Per-individual endpoint outcomes of a virtual cohort."""

    config: PopulationConfig
    ages: list[float]
    weights: list[float]
    heights: list[float]
    sexes: list[str]
    risks: dict[str, np.ndarray]
    grades: dict[str, np.ndarray]
    overall_risk: np.ndarray
    verdicts: list[str]

    def risk_quantiles(
        self, endpoint: str, qs: tuple[int, int, int] = (5, 50, 95)
    ) -> tuple[float, float, float]:
        """Endpoint risk (qlo, med, qhi) across the cohort."""
        qlo, med, qhi = np.quantile(self.risks[endpoint], [q / 100.0 for q in qs])
        return float(qlo), float(med), float(qhi)

    def incidence(self, grade: int = 1) -> dict[str, int]:
        """Counts of individuals with endpoint grade >= ``grade``."""
        return {name: int(np.sum(self.grades[name] >= grade)) for name in self.grades}

    def verdict_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for verdict in self.verdicts:
            counts[verdict] = counts.get(verdict, 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: kv[1], reverse=True))

    def to_dict(self) -> dict[str, object]:
        return {
            "n_individuals": self.config.n_individuals,
            "seed": self.config.seed,
            "age_y_range": [round(self.config.age_range[0], 1), round(self.config.age_range[1], 1)],
            "body": {
                "age_y": [round(a, 1) for a in self.ages],
                "weight_kg": [round(w, 1) for w in self.weights],
                "height_cm": [round(h, 1) for h in self.heights],
            },
            "risk_q5_q50_q95": {
                name: [round(v, 4) for v in self.risk_quantiles(name)] for name in self.risks
            },
            "incidence_grade_ge_1": self.incidence(),
            "incidence_grade_ge_2": self.incidence(2),
            "verdict_counts": self.verdict_counts(),
        }


def sample_profiles(config: PopulationConfig) -> list[HumanProfile]:
    """Draw ``n_individuals`` anthropometrically plausible profiles."""
    rng = np.random.default_rng(config.seed)
    profiles: list[HumanProfile] = []
    for _ in range(config.n_individuals):
        is_male = bool(rng.random() < 0.5)
        sex = Sex.MALE if is_male else Sex.FEMALE
        age = rng.uniform(*config.age_range)
        height_mean = 177.5 if is_male else 163.5
        height = float(np.clip(rng.normal(height_mean, 7.0), 145.0, 205.0))
        bmi = float(np.clip(rng.normal(config.bmi_mean, config.bmi_sd), 16.0, 45.0))
        weight = float(np.clip(bmi * (height / 100.0) ** 2, 40.0, 160.0))
        profiles.append(HumanProfile(sex=sex, age_y=float(age), height_cm=height, weight_kg=weight))
    return profiles


def run_population(spec: RunSpec, config: PopulationConfig | None = None) -> PopulationResult:
    """Run the cohort: one full pipeline per sampled virtual patient."""
    config = config or PopulationConfig()
    profiles = sample_profiles(config)
    risks: dict[str, list[float]] = {}
    grades: dict[str, list[int]] = {}
    overall: list[float] = []
    verdicts: list[str] = []
    for profile in profiles:
        result = run_pipeline(replace(spec, profile=profile))
        for risk in result.toxicity.risks:
            risks.setdefault(risk.endpoint.value, []).append(risk.risk)
            grades.setdefault(risk.endpoint.value, []).append(risk.grade)
        overall.append(result.toxicity.overall_risk())
        verdicts.append(result.verdict)
    return PopulationResult(
        config=config,
        ages=[p.age_y for p in profiles],
        weights=[p.weight_kg for p in profiles],
        heights=[p.height_cm for p in profiles],
        sexes=[p.sex.value for p in profiles],
        risks={name: np.asarray(v, dtype=float) for name, v in risks.items()},
        grades={name: np.asarray(v, dtype=int) for name, v in grades.items()},
        overall_risk=np.asarray(overall, dtype=float),
        verdicts=verdicts,
    )


__all__ = [
    "PopulationConfig",
    "PopulationResult",
    "run_population",
    "sample_profiles",
]
