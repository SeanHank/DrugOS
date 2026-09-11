"""Whole-body PBPK compartment graph and mass-balance ODEs.

Model structure follows the perfusion-limited whole-body PBPK validated by Ye,
Nagar & Korzekwa (Biopharm Drug Dispos 2016;37:123-141): arterial and venous
blood, lung, and the major organs; the gut and spleen drain into the portal
vein which enters the liver together with hepatic arterial inflow.

Mass balance per tissue (amounts in mg, times in hours):

    dA_t/dt = Q_t*(C_ab - C_vbt) - CL_tissue * C_unbound,t

with C_vbt = C_t*BP/Kp_t (Kp_t = tissue-to-plasma total ratio) and
C_unbound,t = C_t/Kpu_t.  Liver clearance acts on the unbound liver
concentration (optionally as a saturable Vmax/Km term); renal clearance acts
on the unbound kidney concentration (filtration of free drug, optionally plus
an active tubular-secretion term ``cl_sec``).  Biliary secretion
(``cl_bil``) moves unbound parent from the liver into a bile compartment that
empties (``k_bile_emptying``) into the small-intestine lumen: that mass is
then subject to the normal SI absorption/reabsorption kinetics, so a fraction
re-enters the portal vein (enterohepatic recirculation) and the rest reaches
the feaces sink.  Oral input uses an ACAT-lite first-order transit
(stomach -> small intestine -> colon); IV input is applied to the venous side;
SC/IM/transdermal input is applied to a ``depot`` compartment that feeds
venous blood by first-order absorption (defaults, or user-provided
``k_depot_absorption``).  First-pass intestinal (gut-wall) extraction
``gut_extraction_eg`` removes a fraction of the SI-absorbed flux before it
enters the portal blood.  When a solubility limit is known (mg/mL, from the
ADMET-AI ``logS`` for novel molecules) the small-intestine absorption flow is
capped at the amount that can be simultaneously in solution in the lumen
volume, so an excess dose spills forward as undissolved drug to the colon and
feces sink (solubility-limited dissolution, doc/05 1.4).  A per-CYP saturable
hepatic term (``cyp_terms``) replaces the lumped metabolic clearance with the
sum of isoform Michaelis-Menten/Hill fluxes on the unbound liver
concentration, with abundance-scaled Vmax built from the physiology
CYP-abundance table (doc/05 4.1).  A native target-mediated drug disposition
(``target_binding``) couples one binding site into a tissue mass balance so
receptor association, dissociation and internalization consume drug directly
in the ODE (the monolithic TMDD limit, doc/05 1.4/2.4).  A finite-dose
multi-layer skin-permeation membrane (``skin_layers``) replaces the generic
transdermal depot with surface -> stratum corneum -> viable epidermis ->
dermis diffusion plus dermal capillary uptake into venous blood (doc/05 1.4).
Every realistic extension (secretory/reabsorptive ``cl_sec``, MM hepatic
``hepatic_vmax``/``hepatic_km``, per-CYP ``cyp_terms``, biliary
``cl_bil``/``k_bile_emptying``, gut-wall ``gut_extraction_eg``, TMDD
``target_binding``, skin ``skin_layers``) is off by default, so baseline runs
reproduce the validated linear-clearance behavior exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from drugos.inputs.models import DosePlan
from drugos.pk.partitions import RrPartition
from drugos.pk.physiology import PORTAL_DRAINING_TISSUES, HumanPhysiology
from drugos.pk.skin import SkinLayers
from drugos.target.targets import Target

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]

_NG_PER_MG = 1e6

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
    non-absorbed sink).  ``skin_layers`` (when set) replaces the transdermal
    depot path with the finite-dose multi-layer skin-permeation membrane
    (``SkinLayers``, doc/05 1.4): surface -> stratum corneum -> viable
    epidermis -> dermis, with dermal capillary removal into venous blood
    weighted by ``depot_bioavailability`` (the complement is tallied in the
    skin unabsorbed sink).  ``solubility_mg_ml`` (when set) limits the dissolved
    amount in each small-intestine segment to ``solubility_mg_ml *
    gi_volume_ml / N`` mg where N = ``si_segments`` (default single-compartment
    SI when ``si_segments`` is None); extra dose stays undissolved and transits
    onward to the colon/feces.

    **ACAT-lite multi-segment SI (off by default).**
    ``si_segments`` splits the small intestine into N sequential sub-compartments
    of equal volume ``gi_volume_ml / N`` mL, each with its own dissolution cap
    and first-order absorption/transit; total SI transit time is preserved
    (per-segment transit = ``k_si_transit * N``).  Bile secretion enters
    segment 0 (proximal SI).  The per-segment absorbed amounts sum to the
    portal flux that enters the liver.  When ``si_segments`` is None (default),
    the model is identical to the single-SI ACAT-lite baseline.
    """

    k_gastric_emptying: float = 1.5
    k_stomach_absorption: float = 0.05
    k_si_absorption: float = 0.55
    k_si_transit: float = 0.25
    k_colon_absorption: float = 0.10
    k_colon_transit: float = 0.08
    k_depot_absorption: float = 0.15
    depot_bioavailability: float = 1.0
    solubility_mg_ml: float | None = None
    gi_volume_ml: float = 250.0
    gut_extraction_eg: float = 0.0
    skin_layers: SkinLayers | None = None
    si_segments: int | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("k_depot_absorption", self.k_depot_absorption),
            ("depot_bioavailability", self.depot_bioavailability),
        ):
            if name == "depot_bioavailability" and not (0.0 < value <= 1.0):
                raise ValueError(f"{name} must be in (0, 1]; got {value}")
            if name == "k_depot_absorption" and value <= 0.0:
                raise ValueError(f"{name} must be positive; got {value}")
        if self.solubility_mg_ml is not None and self.solubility_mg_ml <= 0.0:
            raise ValueError(f"solubility_mg_ml must be positive; got {self.solubility_mg_ml}")
        if self.gi_volume_ml <= 0.0:
            raise ValueError(f"gi_volume_ml must be positive; got {self.gi_volume_ml}")
        if not (0.0 <= self.gut_extraction_eg <= 0.9):
            raise ValueError(f"gut_extraction_eg must be in [0, 0.9]; got {self.gut_extraction_eg}")
        if self.si_segments is not None and self.si_segments < 1:
            raise ValueError(f"si_segments must be >= 1 or None; got {self.si_segments}")


def absorption_rate_from_fa(fa: float, base: float = 0.55) -> float:
    """Small-intestine absorption rate constant (h^-1) tuned to a target Fa."""
    fa = min(max(fa, 0.05), 0.98)
    return base * (fa / 0.85)


@dataclass(frozen=True, slots=True)
class CypTerm:
    """Per-isoform hepatic CYP metabolism term (doc/05 1.4, 4.1).

    Encodes one cytochrome's saturable elimination on the *unbound* liver
    concentration with Michaelis-Menten / Hill kinetics:

        rate = Vmax * c^n / (Km^n + c^n)

    Vmax is the abundance-scaled capacity in mg/h (constructed from the
    physiology CYP-abundance table via :func:`cyp_vmax_mg_h`, or given
    directly); Km is the unbound liver concentration at half-maximal velocity
    (mg/L); ``n=1`` reduces to standard Michaelis-Menten and the linear
    low-dose slope is ``Vmax/Km``.

    Off by default: ``PBPKModel`` keeps the lumped linear ``cl_hep`` hepatic
    term unless an explicit ``cyp_terms`` tuple is supplied, so every shipped
    validation run preserves the validated linear baseline.
    """

    isoform: str
    km_mg_l: float
    vmax_mg_h: float
    hill: float = 1.0

    def __post_init__(self) -> None:
        if self.km_mg_l <= 0:
            raise ValueError("km_mg_l must be positive")
        if self.vmax_mg_h < 0:
            raise ValueError("vmax_mg_h must be non-negative")
        if self.hill <= 0:
            raise ValueError("hill must be positive")

    def rate(self, c_pu_h_mg_l: float) -> float:
        """Metabolic flux (mg/h) at unbound liver concentration ``c_pu_h_mg_l``."""
        c = c_pu_h_mg_l
        if self.hill == 1.0:
            return self.vmax_mg_h * c / (self.km_mg_l + c)
        cn = c**self.hill
        return float(self.vmax_mg_h * cn / (self.km_mg_l**self.hill + cn))


def cyp_vmax_mg_h(content_nmol: float, kcat_1h: float, mw_g_per_mol: float) -> float:
    """Abundance-scaled Vmax (mg/h) of one CYP isoform.

    ``Vmax = kcat * E_total * MW`` where ``E_total`` is the per-isoform liver
    content from ``HumanPhysiology.hepatic_cyp_content_nmol`` (pmol per mg
    microsomal protein x microsomal protein per g x liver mass, Barter et al.
    2013).  Content in nmol, kcat per hour and MW in g/mol give the rate in
    mg/h, matching the compartment units of the PBPK masses.
    """
    if content_nmol <= 0:
        raise ValueError("content_nmol must be positive")
    if kcat_1h <= 0:
        raise ValueError("kcat_1h must be positive")
    if mw_g_per_mol <= 0:
        raise ValueError("mw_g_per_mol must be positive")
    return content_nmol * 1.0e-6 * kcat_1h * mw_g_per_mol


@dataclass(frozen=True, slots=True)
class TargetBinding:
    """Native target-mediated drug disposition (TMDD) coupling (doc/05 1.4/2.4).

    Couples one reversible binding site into the mass balance of a PBPK
    tissue so that receptor association, dissociation and complex
    internalization consume and ultimately clear drug directly inside the
    ODE — the monolithic TMDD limit that the sequential pipeline only
    approximates via the opt-in ``feedback_loop`` driver (doc/05 2.4-2.5).

    Binding is driven by the *unbound* tissue concentration, converted from
    mg/L to nM with the model molecular weight, and uses the same turnover
    model as ``drugos.target.occupancy.simulate_occupancy``: constant
    receptor synthesis ``ksyn = rho*R0``, degradation ``rho*R``, second-order
    association ``kon*D*R``, first-order dissociation ``koff*DR`` (with
    ``koff = kon*kd`` from the ``Target``) and irreversible internalization
    ``kint*DR`` that drains drug out of the tissue into a cleared sink.  The
    per-tissue basal abundance is ``R0 = r0_nm * V_tissue`` (nmol).

    Off by default; supply an explicit tuple on ``PBPKModel`` (together with
    ``mw_g_per_mol``) to activate.
    """

    tissue: str
    target: Target


@dataclass(slots=True)
class PBPKModel:
    """The coupled whole-body compartment system."""

    physiology: HumanPhysiology
    partition: RrPartition
    bp: float
    fup: float
    cl_hep_l_h: float = 0.0
    cl_renal_l_h: float = 0.0
    cl_sec_l_h: float = 0.0
    hepatic_vmax_mg_h: float | None = None
    hepatic_km_mg_l: float | None = None
    cyp_terms: tuple[CypTerm, ...] = ()
    cl_bil_l_h: float = 0.0
    k_bile_emptying_1h: float = 1.0
    target_binding: tuple[TargetBinding, ...] = ()
    mw_g_per_mol: float | None = None
    dose_plan: DosePlan = field(default_factory=DosePlan)
    absorption: AbsorptionParams = field(default_factory=AbsorptionParams)

    order: list[str] = field(default_factory=lambda: [*TISSUE_LIST])

    _indices: dict[str, int] = field(init=False, default_factory=dict)
    _n_state: int = field(init=False, default=0)
    _vc_tissues: list[str] = field(init=False, default_factory=list)
    _tmdd_receptor: list[int] = field(init=False, default_factory=list)
    _tmdd_complex: list[int] = field(init=False, default_factory=list)
    _tmdd_cleared: list[int] = field(init=False, default_factory=list)
    _tmdd_mg_per_nmol: float = field(init=False, default=0.0)
    _tmdd_nm_per_mg_l: float = field(init=False, default=0.0)
    _skin_indices: dict[str, int] = field(init=False, default_factory=dict)
    _skin_volumes_cm3: dict[str, float] = field(init=False, default_factory=dict)
    _si_segment_indices: list[int] = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        if self.cyp_terms and self.hepatic_vmax_mg_h is not None:
            raise ValueError(
                "cyp_terms replaces the lumped hepatic clearance; "
                "hepatic_vmax_mg_h/hepatic_km_mg_l must be left unset"
            )
        if self.target_binding:
            if self.mw_g_per_mol is None:
                raise ValueError("mw_g_per_mol is required when target_binding is used")
            if self.mw_g_per_mol <= 0:
                raise ValueError("mw_g_per_mol must be positive")
        self._indices: dict[str, int] = {}
        self._indices["arterial"] = ARTERIAL_INDEX
        self._indices["venous"] = VENOUS_INDEX
        idx = FIRST_TISSUE_INDEX
        for name in self.order:
            self._indices[name] = idx
            idx += 1
        # GI core indices: stomach, colon, urine, feces, depot, bile.
        self._indices["stomach"] = idx
        self._indices["colon"] = idx + 2  # skip si slot (populated after)
        self._indices["urine"] = idx + 3
        self._indices["feces"] = idx + 4
        self._indices["depot"] = idx + 5
        self._indices["bile"] = idx + 6
        base_n = idx + 7
        # SI segments: placed after the GI fixed block (compat).
        n_seg = self.absorption.si_segments
        if n_seg is not None and n_seg > 1:
            self._si_segment_indices = [base_n + i for i in range(n_seg)]
            base_n += n_seg
        else:
            # Single SI lump (legacy layout): re-use the old si slot.
            self._si_segment_indices = []
            self._indices["si"] = idx + 1
        self._tmdd_receptor: list[int] = []
        self._tmdd_complex: list[int] = []
        self._tmdd_cleared: list[int] = []
        for tb in self.target_binding:
            if tb.tissue not in self._indices:
                raise ValueError(f"target_binding tissue {tb.tissue!r} not in model tissues")
        for _ in self.target_binding:
            self._tmdd_receptor.append(base_n)
            self._tmdd_complex.append(base_n + 1)
            self._tmdd_cleared.append(base_n + 2)
            base_n += 3
        if self.absorption.skin_layers is not None:
            for skin_name in ("surface", "sc", "ve", "dermis", "unabsorbed"):
                self._skin_indices[f"skin_{skin_name}"] = base_n
                base_n += 1
            self._skin_volumes_cm3 = self.absorption.skin_layers.layer_volumes_cm3()
        self._n_state = base_n
        if self.target_binding:
            assert self.mw_g_per_mol is not None
            self._tmdd_mg_per_nmol = self.mw_g_per_mol * 1.0e-6
            self._tmdd_nm_per_mg_l = _NG_PER_MG / self.mw_g_per_mol
        self._vc_tissues = [
            t for t in self.order if t not in PORTAL_DRAINING_TISSUES and t != "lung"
        ]

    @property
    def n_state(self) -> int:
        return self._n_state

    def state_index(self, name: str) -> int:
        return self._indices[name]

    def state_total_mass(self, y: NDArray) -> float:
        """Total drug mass in every modelled compartment (mg).

        Bound (receptor-complexed) drug is converted from nmol to mg, and
        internalized-cleared drug is kept, so the sum stays exactly equal to
        the administered dose for models with native binding; receptor
        amounts (protein, not drug) are excluded.  Multi-layer transdermal
        skin states (surface + membranes + unabsorbed tally) are drug and
        stay inside the sum, so transdermal mass closes against the dose.
        Without ``target_binding``/``skin_layers`` this is simply ``sum(y)``.
        """
        total = float(np.sum(y))
        for r_i, c_i in zip(self._tmdd_receptor, self._tmdd_complex, strict=True):
            total += y[c_i] * (self._tmdd_mg_per_nmol - 1.0) - y[r_i]
        return total

    def _volumes(self) -> dict[str, float]:
        return self.physiology.organ_volume

    def initial_state(self) -> NDArray:
        y0 = np.zeros(self._n_state, dtype=float)
        vols = self._volumes()
        for i, tb in enumerate(self.target_binding):
            y0[self._tmdd_receptor[i]] = tb.target.r0_nm * vols[tb.tissue]
        return y0

    def tmdd_state_record(self, y: NDArray) -> tuple[dict[str, float], ...]:
        """Per-site native-TMDD snapshot of the state vector (mg/nmol).

        Returns one dict per binding site with the free-receptor amount
        (nmol), bound (receptor-complexed) amount in nmol and converted mg,
        and the cumulative internalized-cleared drug mass in mg.
        """
        rows: list[dict[str, float]] = []
        for r_i, c_i, cl_i in zip(
            self._tmdd_receptor, self._tmdd_complex, self._tmdd_cleared, strict=True
        ):
            rows.append(
                {
                    "receptor_nmol": float(y[r_i]),
                    "complex_nmol": float(y[c_i]),
                    "bound_mg": float(y[c_i] * self._tmdd_mg_per_nmol),
                    "cleared_mg": float(y[cl_i]),
                }
            )
        return tuple(rows)

    def _skin_link_fluxes(self, y: NDArray) -> dict[str, float]:
        """Per-link skin-permeation fluxes for one state vector (mg/h)."""
        skin = self.absorption.skin_layers
        assert skin is not None
        vols = self._skin_volumes_cm3
        idx = self._skin_indices
        c_surf = y[idx["skin_surface"]] / vols["surface"]
        c_sc = y[idx["skin_sc"]] / vols["sc"]
        c_ve = y[idx["skin_ve"]] / vols["ve"]
        c_der = y[idx["skin_dermis"]] / vols["dermis"]
        j1 = skin.sc_conductance_cm3_h() * (c_surf - c_sc / skin.surface_sc_partition)
        j2 = skin.ve_conductance_cm3_h() * (c_sc - c_ve / skin.sc_ve_partition)
        j3 = skin.dermis_conductance_cm3_h() * (c_ve - c_der / skin.ve_dermis_partition)
        cap = skin.k_dermal_capillary_1h * y[idx["skin_dermis"]]
        return {"surface_to_sc": j1, "sc_to_ve": j2, "ve_to_dermis": j3, "dermal_capillary": cap}

    def skin_fluxes(self, y: NDArray) -> dict[str, float]:
        """Public skin-permeation flux snapshot (mg/h); empty without layers."""
        return self._skin_link_fluxes(y) if self._skin_indices else {}

    def skin_state_record(self, y: NDArray) -> dict[str, float]:
        """Layer amounts (mg) plus cumulative dermal absorption for one state."""
        if not self._skin_indices:
            return {}
        idx = self._skin_indices
        skin = self.absorption.skin_layers
        assert skin is not None
        return {
            "surface_mg": float(y[idx["skin_surface"]]),
            "sc_mg": float(y[idx["skin_sc"]]),
            "ve_mg": float(y[idx["skin_ve"]]),
            "dermis_mg": float(y[idx["skin_dermis"]]),
            "unabsorbed_mg": float(y[idx["skin_unabsorbed"]]),
            "absorption_rate_mg_h": self.skin_absorption_rate_mg_h(y),
        }

    def skin_absorption_rate_mg_h(self, y: NDArray) -> float:
        """Dermal capillary removal into venous blood (mg/h); 0 without layers."""
        if not self._skin_indices:
            return 0.0
        skin = self.absorption.skin_layers
        assert skin is not None
        derived = skin.k_dermal_capillary_1h * y[self._skin_indices["skin_dermis"]]
        return float(derived * self.absorption.depot_bioavailability)

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
        c_pu_h = self._tissue_concentration(y, "liver") / kpu["liver"]
        if self.cyp_terms:
            hepatic_elim = sum((term.rate(c_pu_h) for term in self.cyp_terms), 0.0)
        elif self.hepatic_vmax_mg_h is not None and self.hepatic_km_mg_l is not None:
            hepatic_elim = self.hepatic_vmax_mg_h * c_pu_h / (self.hepatic_km_mg_l + c_pu_h)
        else:
            hepatic_elim = self.cl_hep_l_h * c_pu_h
        dydt[self._indices["liver"]] = (
            hepatic_in - qh * c_vbh - hepatic_elim - self.cl_bil_l_h * c_pu_h
        )

        # Absorption compartments.
        abs_params = self.absorption
        a_st = y[self._indices["stomach"]]
        a_col = y[self._indices["colon"]]
        st_out = abs_params.k_gastric_emptying * a_st
        col_resorb = abs_params.k_colon_absorption * a_col
        dydt[self._indices["stomach"]] = -st_out - abs_params.k_stomach_absorption * a_st

        # ACAT-lite multi-segment SI (doc/05 §1.6, off by default).
        n_seg = abs_params.si_segments
        if n_seg is not None and n_seg > 1 and self._si_segment_indices:
            seg_k_transit = abs_params.k_si_transit * n_seg  # preserve total SI transit time
            seg_vol_ml = abs_params.gi_volume_ml / n_seg
            dissolve_cap_seg = (
                abs_params.solubility_mg_ml * seg_vol_ml
                if abs_params.solubility_mg_ml is not None
                else None
            )
            total_portal = 0.0
            stomach_abs_out = abs_params.k_stomach_absorption * a_st
            for i, seg_idx in enumerate(self._si_segment_indices):
                a_seg = y[seg_idx]
                seg_abs = abs_params.k_si_absorption * a_seg
                if dissolve_cap_seg is not None:
                    seg_abs = abs_params.k_si_absorption * min(a_seg, dissolve_cap_seg)
                total_portal += seg_abs
                seg_in = (
                    (st_out + stomach_abs_out)
                    if i == 0
                    else seg_k_transit * y[self._si_segment_indices[i - 1]]
                )
                seg_out = seg_k_transit * a_seg
                dydt[seg_idx] = seg_in - seg_abs - seg_out
            dydt[self._indices["colon"]] = (
                seg_k_transit * y[self._si_segment_indices[-1]]
                - col_resorb
                - abs_params.k_colon_transit * a_col
            )
            portal_si = (1.0 - abs_params.gut_extraction_eg) * total_portal
        else:
            # Legacy single-SI path (unchanged).
            a_si = y[self._indices["si"]]
            si_resorb = abs_params.k_si_absorption * a_si
            if abs_params.solubility_mg_ml is not None:
                dissolve_capacity_mg = abs_params.solubility_mg_ml * abs_params.gi_volume_ml
                si_resorb = abs_params.k_si_absorption * min(a_si, dissolve_capacity_mg)
            dydt[self._indices["si"]] = (
                st_out
                + abs_params.k_stomach_absorption * a_st
                - si_resorb
                - abs_params.k_si_transit * a_si
            )
            dydt[self._indices["colon"]] = (
                abs_params.k_si_transit * a_si - col_resorb - abs_params.k_colon_transit * a_col
            )
            portal_si = (1.0 - abs_params.gut_extraction_eg) * si_resorb
        dydt[self._indices["feces"]] = abs_params.k_colon_transit * a_col
        dydt[self._indices["liver"]] += portal_si + col_resorb

        # Biliary excretion & enterohepatic recirculation: unbound parent is
        # secreted into a bile pool that empties into the SI lumen, where it is
        # then subject to the regular absorption/reabsorption kinetics (so the
        # reabsorbed fraction re-enters the portal vein and the rest transits
        # to the colon/feces).
        a_bile = y[self._indices["bile"]]
        bile_out = self.k_bile_emptying_1h * a_bile
        dydt[self._indices["bile"]] = self.cl_bil_l_h * c_pu_h - bile_out
        # Bile empties into SI segment 0 (proximal SI); when using legacy
        # single-SI, _si_segment_indices is empty and we fall back to the
        # "si" index.
        bile_target = (
            self._si_segment_indices[0] if self._si_segment_indices else self._indices["si"]
        )
        dydt[bile_target] += bile_out

        # Depot (SC/IM): first-order absorption into venous blood.  Transdermal
        # uses the same depot unless ``skin_layers`` replaces it below.
        depot = y[self._indices["depot"]]
        depot_out = abs_params.k_depot_absorption * depot
        dydt[self._indices["depot"]] = -depot_out
        dydt[self._indices["venous"]] += depot_out * abs_params.depot_bioavailability
        dydt[self._indices["feces"]] += depot_out * (1.0 - abs_params.depot_bioavailability)

        # Multi-layer transdermal skin permeation (doc/05 1.4, off by default):
        # reversible diffusion links surface -> SC -> VE -> dermis followed by
        # first-order dermal capillary removal into venous blood; the series
        # fluxes cancel inside the membrane so the skin block's net drug
        # derivative is exactly ``-cap``, with the delivered complement
        # tallied in the unabsorbed sink.
        if self._skin_indices:
            fluxes = self._skin_link_fluxes(y)
            j1, j2, j3 = fluxes["surface_to_sc"], fluxes["sc_to_ve"], fluxes["ve_to_dermis"]
            cap = fluxes["dermal_capillary"]
            dydt[self._skin_indices["skin_surface"]] = -j1
            dydt[self._skin_indices["skin_sc"]] = j1 - j2
            dydt[self._skin_indices["skin_ve"]] = j2 - j3
            dydt[self._skin_indices["skin_dermis"]] = j3 - cap
            dydt[self._indices["venous"]] += cap * abs_params.depot_bioavailability
            dydt[self._skin_indices["skin_unabsorbed"]] += cap * (
                1.0 - abs_params.depot_bioavailability
            )

        # Kidney: filtration of free drug into urine, plus optional active
        # tubular secretion (cl_sec) of the unbound kidney drug.
        c_pu_k = self._tissue_concentration(y, "kidney") / kpu["kidney"]
        urine_rate = (self.cl_renal_l_h + self.cl_sec_l_h) * c_pu_k
        dydt[self._indices["kidney"]] = (flows["kidney"] * 60.0) * (
            c_ab - self._venous_blood_out(y, "kidney")
        ) - urine_rate
        dydt[self._indices["urine"]] = urine_rate

        # Native target-mediated drug disposition (doc/05 1.4/2.4): reversible
        # binding plus internalization couples each site into its tissue mass
        # balance — the monolithic TMDD limit.  Rates follow the occupancy
        # turnover model: association kon*D*R and dissociation koff*DR in
        # nmol/h (through the volume the concentration/amount factors cancel),
        # internalization kint*DR drains the complex into a cleared sink, and
        # the *net* bound flux leaves the tissue pool.
        scale = self._tmdd_mg_per_nmol
        for i, tb in enumerate(self.target_binding):
            tgt = tb.target
            volume = self.physiology.organ_volume[tb.tissue]
            free_mg_l = self._tissue_concentration(y, tb.tissue) / kpu[tb.tissue]
            d_nm = free_mg_l * self._tmdd_nm_per_mg_l
            r_nmol = y[self._tmdd_receptor[i]]
            dr_nmol = y[self._tmdd_complex[i]]
            assoc = tgt.kon_nm_h * d_nm * r_nmol
            diss = tgt.koff_1h * dr_nmol
            ksyn = tgt.rho_h * tgt.r0_nm * volume
            dydt[self._tmdd_receptor[i]] = ksyn - tgt.rho_h * r_nmol - assoc + diss
            dydt[self._tmdd_complex[i]] = assoc - diss - tgt.kint_h * dr_nmol
            dydt[self._tmdd_cleared[i]] = tgt.kint_h * dr_nmol * scale
            dydt[self._indices[tb.tissue]] -= (assoc - diss) * scale

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
            elif ev.route.value == "transdermal":
                if self._skin_indices:
                    y[self._skin_indices["skin_surface"]] += ev.dose_mg
                else:
                    y[self._indices["depot"]] += ev.dose_mg
            elif ev.route.value in ("subcutaneous", "intramuscular"):
                y[self._indices["depot"]] += ev.dose_mg
        return y
