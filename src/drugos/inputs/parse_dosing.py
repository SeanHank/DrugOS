"""Dosing schedule construction and validation.

Converts route / amount / regimen descriptors into an ordered sequence of
:class:`DoseEvent` objects consumed by the PBPK simulation engine.
"""

from __future__ import annotations

from drugos.inputs.models import DoseEvent, DosePlan, Route


def build_dose_plan(
    route: str | Route,
    amount_mg: float,
    duration_h: float | None = None,
    interval_h: float | None = None,
    n_doses: int = 1,
    start_h: float = 0.0,
    food_state: str = "fasted",
) -> DosePlan:
    """Build a :class:`DosePlan` from simple descriptors.

    ``amount_mg`` is the dose per administration.  ``duration_h`` applies to
    IV infusion.  ``interval_h`` with ``n_doses > 1`` creates a repeated
    regimen (doses at 0, interval, 2*interval, ...).
    """
    route_enum = Route(route) if isinstance(route, str) else route

    if amount_mg <= 0:
        raise ValueError("amount_mg must be positive")
    if route_enum is Route.IV_INFUSION and not (duration_h and duration_h > 0):
        raise ValueError("IV infusion requires a positive duration_h")
    if n_doses < 1:
        raise ValueError("n_doses must be >= 1")

    event: DoseEvent = DoseEvent(
        time_h=start_h,
        dose_mg=amount_mg,
        route=route_enum,
        infusion_duration_h=duration_h,
        food_state=food_state,
    )
    if interval_h is None or n_doses == 1:
        plan = DosePlan(events=[event])
    else:
        plan = DosePlan(events=[event])
        plan.repeat(interval_h, n_doses)
    return plan


def dose_plan_to_events(plan: DosePlan) -> list[DoseEvent]:
    """Return the schedule as a flat, time-sorted list of events."""
    return sorted(plan.events, key=lambda e: e.time_h)


__all__ = ["DoseEvent", "DosePlan", "Route", "build_dose_plan", "dose_plan_to_events"]
