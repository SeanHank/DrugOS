"""Shared input data models for the DrugOS pipeline.

All inputs to the pipeline are representable by these models. Every field can be
overridden by the user; where a field is derived (e.g. from the chemical
structure) it carries metadata describing how it was obtained.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field

from pydantic import BaseModel, Field, field_validator


class Sex(enum.StrEnum):
    """Biological sex used for physiological scaling."""

    MALE = "male"
    FEMALE = "female"

    @classmethod
    def _missing_(cls, value: object) -> Sex:
        value = str(value).strip().lower()
        if value.startswith("m"):
            return cls.MALE
        if value.startswith("f"):
            return cls.FEMALE
        raise ValueError(f"unknown sex: {value}")


class Route(enum.StrEnum):
    """Route of administration."""

    IV_BOLUS = "iv_bolus"
    IV_INFUSION = "iv_infusion"
    ORAL = "oral"
    SUBCUTANEOUS = "subcutaneous"
    INTRAMUSCULAR = "intramuscular"
    TRANSDERMAL = "transdermal"


class IonClass(enum.StrEnum):
    """Dominant ionization state of the molecule at physiological pH."""

    NEUTRAL = "neutral"
    CATION = "cation"
    ANION = "anion"
    ZWITTERION = "zwitterion"


class RrClass(enum.StrEnum):
    """Rodgers-Rowland partition-coefficient equation class.

    ``BASE``: a moderate/strong base (basic pKa > 7) -> strong-basis equation.
    ``ACID``/``NEUTRAL``/``WEAK_BASE``: general equation; protein term uses the
    tissue albumin ratio for ionizable compounds and the lipoprotein ratio for
    neutral compounds.
    """

    STRONG_BASE = "strong_base"
    WEAK_BASE = "weak_base"
    ACID = "acid"
    ZWITTERION = "zwitterion"
    NEUTRAL = "neutral"


class Molecule(BaseModel):
    """Chemical structure together with physchem descriptors.

    ``canonical_smiles`` and the computed descriptors are filled by
    :func:`drugos.inputs.parse_structure.parse_structure`.  The user may
    override any computed value, notably ``pka_acids``/``pka_bases`` when a
    measured value supersedes the heuristic macro-pKa estimate.
    """

    name: str | None = Field(default=None, description="Free-text compound name.")
    canonical_smiles: str | None = Field(default=None, description="RDKit-canonicalized SMILES.")
    canonical_isomeric_smiles: str | None = None
    inchi: str | None = None
    inchikey: str | None = None
    molecular_formula: str | None = None

    mw: float | None = Field(default=None, description="Molecular weight (g/mol).")
    log_p: float | None = Field(default=None, description="Crippen logP (octanol/water).")
    log_d_7_4: float | None = Field(
        default=None, description="Distribution coefficient logD at pH 7.4."
    )
    tpsa: float | None = Field(default=None, description="Topological polar surface area.")
    hbd: int | None = Field(default=None, description="Lipinski H-bond donors.")
    hba: int | None = Field(default=None, description="Lipinski H-bond acceptors.")
    rotatable_bonds: int | None = None
    aromatic_rings: int | None = None
    formal_charge: int | None = None
    heavy_atoms: int | None = None

    pka_acids: list[float] = Field(
        default_factory=list, description="Heuristic macro-acid pKa values (ascending)."
    )
    pka_bases: list[float] = Field(
        default_factory=list, description="Heuristic macro-base pKa values (descending)."
    )

    ion_class: IonClass = Field(
        default=IonClass.NEUTRAL, description="Dominant ionization at pH 7.4."
    )
    rr_class: RrClass = Field(
        default=RrClass.NEUTRAL,
        description="Rodgers-Rowland partition equation class.",
    )

    @property
    def strongest_acid_pka(self) -> float | None:
        """Lowest acidic pKa (most acidic group), if any."""
        return self.pka_acids[0] if self.pka_acids else None

    @property
    def strongest_base_pka(self) -> float | None:
        """Highest basic pKa (most basic group), if any."""
        return self.pka_bases[0] if self.pka_bases else None

    @property
    def is_ionizable(self) -> bool:
        return bool(self.pka_acids or self.pka_bases)


class HumanProfile(BaseModel):
    """User-supplied sparse human profile.

    Only ``sex`` is strictly required; sparse profiles are completed against
    standard reference data by :func:`drugos.pk.physiology.build_human`.
    """

    sex: Sex = Field(..., description="Biological sex.")
    age_y: float = Field(default=40.0, ge=0.0, le=120.0)
    height_cm: float = Field(default=170.0, gt=0.0)
    weight_kg: float = Field(default=70.0, gt=0.0)

    hematocrit: float | None = Field(default=None, ge=0.2, le=0.6)
    albumin_g_l: float | None = Field(default=None, gt=0.0, le=100.0)
    agp_g_l: float | None = Field(default=None, gt=0.0, le=10.0)
    gfr_ml_min: float | None = Field(
        default=None, gt=0.0, description="Glomerular filtration rate (mL/min)."
    )
    cardiac_index_l_min_m2: float | None = Field(
        default=None, gt=0.0, le=10.0, description="Cardiac index (L/min/m^2)."
    )
    serum_creatinine_mg_dl: float | None = Field(default=None, gt=0.0, le=20.0)

    mild_renal_impairment: bool = False
    moderate_renal_impairment: bool = False
    mild_hepatic_impairment: bool = False
    moderate_hepatic_impairment: bool = False
    heart_failure: bool = False

    @field_validator("sex", mode="before")
    @classmethod
    def _coerce_sex(cls, value: object) -> Sex:
        return value if isinstance(value, Sex) else Sex(str(value))


@dataclass(slots=True)
class DoseEvent:
    """A single administration event on the dosing schedule."""

    time_h: float = 0.0
    dose_mg: float = 0.0
    route: Route = Route.ORAL
    infusion_duration_h: float | None = None
    food_state: str = "fasted"

    def __post_init__(self) -> None:
        if self.dose_mg < 0:
            raise ValueError("dose_mg must be >= 0")
        if self.route is Route.IV_INFUSION and not (
            self.infusion_duration_h and self.infusion_duration_h > 0
        ):
            raise ValueError("IV infusion requires a positive infusion_duration_h")


class DosePlan(BaseModel):
    """A dosing schedule: an ordered list of administration events."""

    events: list[DoseEvent] = Field(default_factory=list)

    @classmethod
    def iv_bolus(cls, dose_mg: float, time_h: float = 0.0) -> DosePlan:
        return cls(events=[DoseEvent(time_h=time_h, dose_mg=dose_mg, route=Route.IV_BOLUS)])

    @classmethod
    def iv_infusion(cls, dose_mg: float, duration_h: float, time_h: float = 0.0) -> DosePlan:
        return cls(
            events=[
                DoseEvent(
                    time_h=time_h,
                    dose_mg=dose_mg,
                    route=Route.IV_INFUSION,
                    infusion_duration_h=duration_h,
                )
            ]
        )

    @classmethod
    def oral(cls, dose_mg: float, time_h: float = 0.0, food_state: str = "fasted") -> DosePlan:
        return cls(
            events=[
                DoseEvent(time_h=time_h, dose_mg=dose_mg, route=Route.ORAL, food_state=food_state)
            ]
        )

    def add(self, event: DoseEvent) -> DosePlan:
        self.events.append(event)
        self.events.sort(key=lambda e: e.time_h)
        return self

    def repeat(self, interval_h: float, n: int) -> DosePlan:
        """Repeat the current schedule every ``interval_h`` hours, ``n`` times."""
        if n < 0 or interval_h <= 0:
            raise ValueError("repeat requires n >= 0 and interval_h > 0")
        base = list(self.events)
        extra = [
            DoseEvent(
                time_h=e.time_h + k * interval_h,
                dose_mg=e.dose_mg,
                route=e.route,
                infusion_duration_h=e.infusion_duration_h,
                food_state=e.food_state,
            )
            for k in range(1, n)
            for e in base
        ]
        self.events = [*base, *extra]
        self.events.sort(key=lambda e: e.time_h)
        return self

    @property
    def total_dose_mg(self) -> float:
        return sum(e.dose_mg for e in self.events)

    @property
    def routes(self) -> list[Route]:
        return list({e.route for e in self.events})

    @field_validator("events")
    @classmethod
    def _sort_events(cls, events: list[DoseEvent]) -> list[DoseEvent]:
        return sorted(events, key=lambda e: e.time_h)


@dataclass(slots=True)
class DrugParameters:
    """Physicochemical and ADME parameters required by the PBPK layer.

    ``mol`` holds structure-derived descriptors; the remaining fields are
    typically predicted by :mod:`drugos.pk.admet` or measured. Values are
    deliberately optional so that a staged pipeline degrades gracefully, until
    everything needed by the requested simulation has been supplied.
    """

    mol: Molecule | None = None
    log_p: float | None = None
    pka_acids: list[float] = field(default_factory=list)
    pka_bases: list[float] = field(default_factory=list)
    fup: float | None = None
    bp: float | None = None
    log_s: float | None = None
    cl_int_hep_l_h: float | None = None
    fa: float | None = None
    peff: float | None = None
    fu_inc: float | None = None
    t_half_h: float | None = None
    vss_prior_l: float | None = None
    permeability_class: str | None = None
    source: str = "predicted"
