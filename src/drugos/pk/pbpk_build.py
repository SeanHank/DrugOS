"""Whole-body PBPK compartment graph and mass-balance ODEs.

Model structure follows the perfusion-limited whole-body PBPK validated by Ye,
Nagar & Korzekwa (Biopharm Drug Dispos 2016;37:123-141): arterial and venous
blood, lung, and the major organs; the gut and spleen drain into the portal
vein which enters the liver together with hepatic arterial inflow.

Mass balance per tissue (amounts in mg, times in hours):

    dA_t/dt = Q_t*(C_ab - C_vbt) - CL_tissue * C_unbound,t

with C_vbt = C_t*BP/Kp_t (Kp_t = tissue-to-plasma total ratio) and
C_unbound,t = C_t/Kpu_t.  Liver clearance acts on the unbound liver
concentration; renal clearance acts on the unbound kidney concentration
(filtration of free drug).  Oral input uses an ACAT-lite first-order transit
(stomach -> small intestine -> colon); IV input is applied to the venous side;
SC/IM/transdermal input is applied to a ``depot`` compartment that feeds
venous blood by first-order absorption (defaults, or user-provided
``k_depot_absorption``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from drugos.inputs.models import DosePlan
from drugos.pk.partitions import RrPartition
from drugos.pk.physiology import PORTAL_DRAINING_TISSUES, HumanPhysiology

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]

TISSUE_LIST = [
    "lung",
    "heart",
    "liver",
    "kidney",
    "brain",
    "adipose",
    "bone_rest",
    "muscle",
    "skin",
    "gut",
    "spleen",
]

ARTERIAL_INDEX = 0
VENOUS_INDEX = 1
FIRST_TISSUE_INDEX = 2


@dataclass(slots=True)
class AbsorptionParams:
    """First-order absorption kinetics (oral ACAT-lite + SC/IM depot).

    All rate constants are per hour. ``feces_fraction`` is a model diagnostic
    (fraction of oral dose reaching the colon without absorption).  SC/IM doses
    enter a ``depot`` compartment and are absorbed first-order into venous
    blood with ``k_depot_absorption``; ``depot_bioavailability`` is the
    bioavailable fraction (the complement is lost from the depot to the
    non-absorbed sink).
    """

    k_gastric_emptying: float = 1.5
    k_stomach_absorption: float = 0.05
    k_si_absorption: float = 0.55
    k_si_transit: float = 0.25
    k_colon_absorption: float = 0.10
    k_colon_transit: float = 0.08
    k_depot_absorption: float = 0.15
    depot_bioavailability: float = 1.0

    def __post_init__(self) -> None:
        for name, value in (
            ("k_depot_absorption", self.k_depot_absorption),
            ("depot_bioavailability", self.depot_bioavailability),
        ):
            if name == "depot_bioavailability" and not (0.0 < value <= 1.0):
                raise ValueError(f"{name} must be in (0, 1]; got {value}")
            if name == "k_depot_absorption" and value <= 0.0:
                raise ValueError(f"{name} must be positive; got {value}")


def absorption_rate_from_fa(fa: float, base: float = 0.55) -> float:
    """Small-intestine absorption rate constant (h^-1) tuned to a target Fa."""
    fa = min(max(fa, 0.05), 0.98)
    return base * (fa / 0.85)


@dataclass(slots=True)
class PBPKModel:
    """The coupled whole-body compartment system."""

    physiology: HumanPhysiology
    partition: RrPartition
    bp: float
    fup: float
    cl_hep_l_h: float = 0.0
    cl_renal_l_h: float = 0.0
    dose_plan: DosePlan = field(default_factory=DosePlan)
    absorption: AbsorptionParams = field(default_factory=AbsorptionParams)

    order: list[str] = field(default_factory=lambda: [*TISSUE_LIST])

    _indices: dict[str, int] = field(init=False, default_factory=dict)
    _n_state: int = field(init=False, default=0)
    _vc_tissues: list[str] = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        self._indices: dict[str, int] = {}
        self._indices["arterial"] = ARTERIAL_INDEX
        self._indices["venous"] = VENOUS_INDEX
        idx = FIRST_TISSUE_INDEX
        for name in self.order:
            self._indices[name] = idx
            idx += 1
        self._indices["stomach"] = idx
        self._indices["si"] = idx + 1
        self._indices["colon"] = idx + 2
        self._indices["urine"] = idx + 3
        self._indices["feces"] = idx + 4
        self._indices["depot"] = idx + 5
        self._n_state = idx + 6
        self._vc_tissues = [
            t for t in self.order if t not in PORTAL_DRAINING_TISSUES and t != "lung"
        ]

    @property
    def n_state(self) -> int:
        return self._n_state

    def state_index(self, name: str) -> int:
        return self._indices[name]

    def state_total_mass(self, y: NDArray) -> float:
        """Total drug mass in every modelled compartment (mg)."""
        return float(np.sum(y))

    def _volumes(self) -> dict[str, float]:
        return self.physiology.organ_volume

    def initial_state(self) -> NDArray:
        return np.zeros(self._n_state, dtype=float)

    def _tissue_concentration(self, y: NDArray, name: str) -> float:
        vols = self._volumes()
        return float(y[self._indices[name]] / vols[name])

    def _venous_blood_out(self, y: NDArray, name: str) -> float:
        """C_vb,t: drug concentration in the blood leaving tissue ``name``."""
        c_t = self._tissue_concentration(y, name)
        return c_t * self.bp / self.partition.kp[name]

    def rhs(self, t: float, y: NDArray) -> NDArray:
        flows = self.physiology.organ_flow
        kp = self.partition.kp
        kpu = self.partition.kpu
        v_ab = self.physiology.arterial_blood_l
        v_vb = self.physiology.venous_blood_l
        qc = self.physiology.cardiac_output_l_min * 60.0  # L/h

        dydt = np.zeros_like(y)
        c_ab = y[self._indices["arterial"]] / v_ab
        c_vb = y[self._indices["venous"]] / v_vb
        c_lung = self._tissue_concentration(y, "lung")
        c_vbl = c_lung * self.bp / kp["lung"]

        # Arterial and venous blood (amount balances: mg/h).
        dydt[self._indices["arterial"]] = qc * (c_vbl - c_ab)
        mixing = -qc * c_vb
        for name in self._vc_tissues:
            mixing += (flows[name] * 60.0) * self._venous_blood_out(y, name)
        dydt[self._indices["venous"]] = mixing + self._infusion_rate(t)

        # Non-eliminating tissues (perfusion-limited exchange).
        for name in self.order:
            if name in ("lung", "liver", "kidney"):
                continue
            q = flows[name] * 60.0  # L/h
            c_vbt = self._venous_blood_out(y, name)
            dydt[self._indices[name]] = q * (c_ab - c_vbt)

        # Lung: perfused by the full cardiac output.
        dydt[self._indices["lung"]] = qc * (c_vb - c_vbl)

        # Liver: arterial + portal inflow, hepatic clearance on unbound conc.
        qh = flows["liver"] * 60.0
        qg = flows["gut"] * 60.0
        qsp = flows["spleen"] * 60.0
        qh_art = qh - qg - qsp
        c_vbg = self._venous_blood_out(y, "gut")
        c_vbsp = self._venous_blood_out(y, "spleen")
        c_vbh = self._venous_blood_out(y, "liver")
        hepatic_in = qh_art * c_ab + qg * c_vbg + qsp * c_vbsp
        c_pu_h = c_liver = self._tissue_concentration(y, "liver")
        c_pu_h = c_liver / kpu["liver"]
        dydt[self._indices["liver"]] = hepatic_in - qh * c_vbh - self.cl_hep_l_h * c_pu_h

        # Absorption compartments.
        abs_params = self.absorption
        a_st = y[self._indices["stomach"]]
        a_si = y[self._indices["si"]]
        a_col = y[self._indices["colon"]]
        st_out = abs_params.k_gastric_emptying * a_st
        si_resorb = abs_params.k_si_absorption * a_si
        col_resorb = abs_params.k_colon_absorption * a_col
        dydt[self._indices["stomach"]] = -st_out - abs_params.k_stomach_absorption * a_st
        dydt[self._indices["si"]] = (
            st_out
            + abs_params.k_stomach_absorption * a_st
            - si_resorb
            - abs_params.k_si_transit * a_si
        )
        dydt[self._indices["colon"]] = (
            abs_params.k_si_transit * a_si - col_resorb - abs_params.k_colon_transit * a_col
        )
        dydt[self._indices["feces"]] = abs_params.k_colon_transit * a_col
        dydt[self._indices["liver"]] += si_resorb + col_resorb

        # Depot (SC/IM/transdermal): first-order absorption into venous blood.
        depot = y[self._indices["depot"]]
        depot_out = abs_params.k_depot_absorption * depot
        dydt[self._indices["depot"]] = -depot_out
        dydt[self._indices["venous"]] += depot_out * abs_params.depot_bioavailability
        dydt[self._indices["feces"]] += depot_out * (1.0 - abs_params.depot_bioavailability)

        # Kidney: filtration of free drug into urine.
        c_pu_k = self._tissue_concentration(y, "kidney") / kpu["kidney"]
        urine_rate = self.cl_renal_l_h * c_pu_k
        dydt[self._indices["kidney"]] = (flows["kidney"] * 60.0) * (
            c_ab - self._venous_blood_out(y, "kidney")
        ) - urine_rate
        dydt[self._indices["urine"]] = urine_rate

        return dydt

    def _infusion_rate(self, t: float) -> float:
        """Total IV infusion rate into venous blood at time ``t`` (mg/h)."""
        rate = 0.0
        for ev in self.dose_plan.events:
            if (
                ev.route.value == "iv_infusion"
                and ev.infusion_duration_h
                and ev.time_h <= t <= ev.time_h + ev.infusion_duration_h
            ):
                rate += ev.dose_mg / ev.infusion_duration_h
        return rate

    def apply_event(self, y: NDArray, t: float) -> NDArray:
        """Apply a bolus / oral / SC / IM dose event in place; returns the state."""
        for ev in self.dose_plan.events:
            if abs(ev.time_h - t) > 1e-9:
                continue
            if ev.route.value == "iv_bolus":
                y[self._indices["venous"]] += ev.dose_mg
            elif ev.route.value == "oral":
                y[self._indices["stomach"]] += ev.dose_mg
            elif ev.route.value in ("subcutaneous", "intramuscular", "transdermal"):
                y[self._indices["depot"]] += ev.dose_mg
        return y
