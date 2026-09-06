"""Phase 6 robustness & validation modules (doc/07 Phase 6, D21-D23).

* ``uncertainty`` — D21 parameter ensemble with percentile bands.
* ``population`` — D22 virtual-patient cohort with population incidence.
* ``sensitivity`` — D23 local (OAT) + global (Sobol) sensitivity and drivers.

The three analyses share the same multiplicative PK/potency scalars
(``robustness.uncertainty.apply_multipliers``) and the scalar endpoint
extractors (``points_from_result``), all deterministic by seed.
"""

from drugos.robustness.population import (
    PopulationConfig,
    PopulationResult,
    run_population,
    sample_profiles,
)
from drugos.robustness.sensitivity import (
    LocalSensitivityResult,
    SensitivityConfig,
    SobolSensitivityResult,
    first_total_indices,
    local_sensitivity,
    run_sobol_sensitivity,
    saltelli_design,
)
from drugos.robustness.uncertainty import (
    CURVE_NAMES,
    DEFAULT_KEYS,
    EnsembleConfig,
    EnsembleResult,
    apply_multipliers,
    curves_from_result,
    points_from_result,
    run_uncertainty,
    sample_multipliers,
)

__all__ = [
    "CURVE_NAMES",
    "DEFAULT_KEYS",
    "EnsembleConfig",
    "EnsembleResult",
    "LocalSensitivityResult",
    "PopulationConfig",
    "PopulationResult",
    "SensitivityConfig",
    "SobolSensitivityResult",
    "apply_multipliers",
    "curves_from_result",
    "first_total_indices",
    "local_sensitivity",
    "points_from_result",
    "run_population",
    "run_sobol_sensitivity",
    "run_uncertainty",
    "saltelli_design",
    "sample_multipliers",
    "sample_profiles",
]
