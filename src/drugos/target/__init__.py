"""target package: Stage 2 — concentration -> target binding (doc/05 section 2).

Provides target definitions and the off-target safety panel
(:mod:`drugos.target.targets`) plus the turnover-mediated occupancy model
(:mod:`drugos.target.occupancy`) that turns Stage-1 free tissue
concentrations into per-site occupancy trajectories and the time-at-target
exposure signal consumed by Stage 3.
"""

from drugos.target.occupancy import (
    PanelEngagement,
    TargetOccupancyResult,
    free_binding_conc_mg_l_to_nm,
    simulate_occupancy,
    simulate_panel,
)
from drugos.target.targets import Target, safety_panel

__all__ = [
    "PanelEngagement",
    "Target",
    "TargetOccupancyResult",
    "free_binding_conc_mg_l_to_nm",
    "safety_panel",
    "simulate_occupancy",
    "simulate_panel",
]
