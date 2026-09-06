"""Human profile resolution into a full physiology parameter set.

Thin public wrapper over :func:`drugos.pk.physiology.build_human`.
"""

from __future__ import annotations

from typing import Any

from drugos.inputs.models import HumanProfile, Sex
from drugos.pk.physiology import HumanPhysiology, build_human


def resolve_human(profile: HumanProfile) -> HumanPhysiology:
    """Resolve a sparse :class:`HumanProfile` into a complete physiology set."""
    return build_human(profile)


def human_profile(
    sex: str | Sex,
    age_y: float = 40.0,
    height_cm: float = 170.0,
    weight_kg: float = 70.0,
    **overrides: Any,
) -> HumanProfile:
    """Construct a :class:`HumanProfile`, overriding any additional field."""
    return HumanProfile(
        sex=Sex(sex), age_y=age_y, height_cm=height_cm, weight_kg=weight_kg, **overrides
    )


__all__ = ["HumanPhysiology", "HumanProfile", "Sex", "human_profile", "resolve_human"]
