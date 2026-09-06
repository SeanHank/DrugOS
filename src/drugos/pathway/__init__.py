"""pathway package: Stage 3 — target binding -> signaling pathway (QSP).

Provides a small auditable reaction DSL (:mod:`drugos.pathway.graph`), an ODE
engine with signal-driven perturbation simulation and Emax/Hill dose-response
fitting (:mod:`drugos.pathway.simulator`), and a canonical MAPK cascade model
used for amplification validation (doc/08 Tier 2, Stage 3).
"""

from drugos.pathway.graph import (
    Activation,
    FirstOrderDegradation,
    PathwayModel,
    ReversibleBinding,
    SourceProduction,
)
from drugos.pathway.simulator import (
    CompiledPathway,
    DoseResponseFit,
    PathwayResult,
    compile_model,
    dose_response,
    mapk_cascade,
    pathway_steady_state,
    simulate_pathway,
)

__all__ = [
    "Activation",
    "CompiledPathway",
    "DoseResponseFit",
    "FirstOrderDegradation",
    "PathwayModel",
    "PathwayResult",
    "ReversibleBinding",
    "SourceProduction",
    "compile_model",
    "dose_response",
    "mapk_cascade",
    "pathway_steady_state",
    "simulate_pathway",
]
