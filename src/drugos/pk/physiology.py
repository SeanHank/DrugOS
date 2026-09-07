"""Whole-body physiology: tissue composition, organ volumes/blood flows.

Data provenance
---------------
- Tissue phospholipid/water composition, acidic-phospholipid concentration and
  tissue-to-plasma albumin/lipoprotein ratios: Rodgers & Rowland methodology
  (J Pharm Sci 2005;94:1259/1627; Pharm Res 2006/2007), machine-readable form
  from the Metrum Research PBPK_PC repository (``tissue_comp_R&R.csv``).
- Organ volumes and blood flows for a reference 70 kg human male: Ye, Nagar &
  Korzekwa, Biopharm Drug Dispos 2016;37:123-141 (PMC4844798, Table 2).
- Allometric scaling of organ volumes and flows: Willmann et al. (2007)
  population covariate equations; cardiac output from a cardiac index applied
  to body surface area (Mosteller).
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from drugos.inputs.models import HumanProfile, Sex

REFERENCE_MALE_WEIGHT_KG = 70.0

PLASMA_PH = 7.4
INTRACELLULAR_PH = 7.0
BLOOD_CELL_PH = 7.22
HEMATOCRIT_REF_MALE = 0.45
HEMATOCRIT_REF_FEMALE = 0.40
CARDIAC_OUTPUT_REF_L_MIN = 5.69

# Plasma neutral-lipid / neutral-phospholipid volume fractions (R&R rat data).
PLASMA_F_NL = 0.0023
PLASMA_F_NP = 0.0013

# Red-blood-cell fractions (R&R).
RBC_F_NL = 0.0017
RBC_F_NP = 0.0029
RBC_F_IW = 0.603
RBC_AP_MG_G = 0.500


@dataclass(frozen=True, slots=True)
class TissueComposition:
    """Lipid/water fractions and protein ratios of a tissue (R&R methodology).

    Fractions are of whole-tissue volume; ``ap_mg_g`` is the acidic
    phospholipid concentration in mg/g tissue; ``ar``/``lr`` are the
    tissue-to-plasma ratios of albumin and lipoprotein capacities.
    """

    f_nl: float
    f_np: float
    f_ew: float
    f_iw: float
    ap_mg_g: float
    ar: float
    lr: float

    @property
    def water_fraction(self) -> float:
        return self.f_ew + self.f_iw


# Tissues of the whole-body PBPK (Ye 2016 set), with R&R composition.
_TISSUES: dict[str, tuple[float, float, float, float, float, float, float]] = {
    # name: (f_nl, f_np, f_ew, f_iw, ap_mg_g, ar, lr)
    "adipose": (0.853, 0.0016, 0.135, 0.017, 0.40, 0.049, 0.068),
    "bone_rest": (0.017, 0.0017, 0.100, 0.346, 0.67, 0.100, 0.050),
    "brain": (0.039, 0.0015, 0.162, 0.620, 0.40, 0.048, 0.041),
    "gut": (0.038, 0.0125, 0.282, 0.475, 2.41, 0.158, 0.141),
    "heart": (0.014, 0.0111, 0.320, 0.456, 2.25, 0.157, 0.160),
    "kidney": (0.012, 0.0242, 0.273, 0.483, 5.03, 0.130, 0.137),
    "liver": (0.014, 0.0240, 0.161, 0.573, 4.56, 0.086, 0.161),
    "lung": (0.022, 0.0128, 0.336, 0.446, 3.91, 0.212, 0.168),
    "muscle": (0.010, 0.0072, 0.118, 0.630, 1.53, 0.064, 0.059),
    "skin": (0.060, 0.0044, 0.382, 0.291, 1.32, 0.277, 0.096),
    "spleen": (0.0077, 0.0113, 0.207, 0.579, 3.18, 0.097, 0.207),
}

TISSUE_COMPOSITION: dict[str, TissueComposition] = {
    name: TissueComposition(f_nl, f_np, f_ew, f_iw, ap, ar, lr)
    for name, (f_nl, f_np, f_ew, f_iw, ap, ar, lr) in _TISSUES.items()
}

# Organs draining into the portal vein (their outflow reaches the liver).
PORTAL_DRAINING_TISSUES = ("gut", "spleen")

# Reference organ blood flow (L/min) and volume (L) for a 70 kg human male
# (Ye, Nagar & Korzekwa 2016, Table 2).
REF_FLOW_L_MIN: dict[str, float] = {
    "adipose": 0.296,
    "bone_rest": 0.563,
    "brain": 0.683,
    "gut": 0.967,
    "heart": 0.228,
    "kidney": 1.08,
    "liver": 1.42,
    "lung": 5.69,
    "muscle": 1.08,
    "skin": 0.330,
    "spleen": 0.114,
}

REF_VOLUME_L: dict[str, float] = {
    "adipose": 12.46,
    "bone_rest": 15.52,
    "brain": 1.40,
    "gut": 1.20,
    "heart": 0.330,
    "kidney": 0.310,
    "liver": 1.82,
    "lung": 0.530,
    "muscle": 28.0,
    "skin": 2.59,
    "spleen": 0.168,
}

REF_ARTERIAL_BLOOD_L = 0.7
REF_VENOUS_BLOOD_L = 1.4

# Allometric exponents (Willmann et al. 2007 style) per organ volume.
_VOLUME_EXPONENT: dict[str, float] = {
    "adipose": 1.0,
    "bone_rest": 1.0,
    "brain": 0.7,
    "gut": 0.85,
    "heart": 0.75,
    "kidney": 0.75,
    "liver": 0.87,
    "lung": 0.75,
    "muscle": 1.0,
    "skin": 0.78,
    "spleen": 1.0,
}

# Sex factors applied to reference volume (male 1.0, female scaling).
_SEX_FACTOR: dict[Sex, dict[str, float]] = {
    Sex.MALE: dict.fromkeys(_VOLUME_EXPONENT, 1.0),
    Sex.FEMALE: {
        "adipose": 1.44,
        "bone_rest": 0.90,
        "brain": 0.95,
        "gut": 1.0,
        "heart": 0.82,
        "kidney": 0.90,
        "liver": 0.85,
        "lung": 0.85,
        "muscle": 0.70,
        "skin": 0.85,
        "spleen": 0.95,
    },
}


@dataclass(slots=True)
class HumanPhysiology:
    """Fully-resolved physiology parameter set for one simulated individual."""

    sex: Sex
    age_y: float
    height_cm: float
    weight_kg: float
    bsa_m2: float
    cardiac_output_l_min: float
    hematocrit: float
    albumin_g_l: float
    agp_g_l: float
    gfr_l_min: float

    organ_volume: dict[str, float] = field(default_factory=dict)
    organ_flow: dict[str, float] = field(default_factory=dict)
    arterial_blood_l: float = REF_ARTERIAL_BLOOD_L
    venous_blood_l: float = REF_VENOUS_BLOOD_L

    liver_microsomal_protein_mg_g: float = 40.0
    liver_mass_g: float = 1500.0

    @property
    def cardiac_output_ml_min(self) -> float:
        return self.cardiac_output_l_min * 1000.0

    @property
    def gfr_ml_min(self) -> float:
        return self.gfr_l_min * 1000.0

    @property
    def organ_names(self) -> list[str]:
        return sorted(self.organ_volume)

    def tissue_composition(self) -> dict[str, TissueComposition]:
        return TISSUE_COMPOSITION

    @property
    def liver_flow_art_l_min(self) -> float:
        """Hepatic arterial inflow = total hepatic flow minus portal inflow."""
        return self.organ_flow["liver"] - sum(self.organ_flow[t] for t in PORTAL_DRAINING_TISSUES)

    @property
    def hepatic_portal_flow_l_min(self) -> float:
        return sum(self.organ_flow[t] for t in PORTAL_DRAINING_TISSUES)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def mosteller_bsa(height_cm: float, weight_kg: float) -> float:
    """Body surface area by the Mosteller formula (m^2)."""
    return math.sqrt(height_cm * weight_kg / 3600.0)


def _age_adjust_defaults(age_y: float, sex: Sex) -> dict[str, float]:
    shemacrit = HEMATOCRIT_REF_MALE if sex is Sex.MALE else HEMATOCRIT_REF_FEMALE
    gfr_ml_min = 125.0 if sex is Sex.MALE else 110.0
    albumin = 43.0
    agp = 0.7
    if age_y >= 60:
        albumin = 41.0
        agp = 0.9
    if age_y >= 75:
        gfr_ml_min *= 0.90
    if sex is Sex.FEMALE:
        gfr_ml_min *= 0.88
    return {
        "hematocrit": shemacrit,
        "gfr_ml_min": gfr_ml_min,
        "albumin_g_l": albumin,
        "agp_g_l": agp,
    }


def _resolve_proteins(profile: HumanProfile) -> dict[str, float]:
    base = _age_adjust_defaults(profile.age_y, profile.sex)
    albumin = profile.albumin_g_l if profile.albumin_g_l is not None else base["albumin_g_l"]
    agp = profile.agp_g_l if profile.agp_g_l is not None else base["agp_g_l"]
    return {"albumin_g_l": albumin, "agp_g_l": agp}


def _resolve_gfr(profile: HumanProfile) -> float:
    base = _age_adjust_defaults(profile.age_y, profile.sex)
    if profile.gfr_ml_min is not None:
        gfr = profile.gfr_ml_min
    elif profile.serum_creatinine_mg_dl is not None:
        # Production-validated baseline: CKD-EPI 2021 race-free eGFR from a
        # measured serum creatinine (doc/12 row 4b, R-6), scaled to the
        # per-subject absolute GFR by Mosteller BSA.
        from drugos.organ.kidney import ckdepi_2021_egfr

        bsa = mosteller_bsa(profile.height_cm, profile.weight_kg)
        gfr = ckdepi_2021_egfr(
            profile.serum_creatinine_mg_dl,
            profile.age_y,
            profile.sex is Sex.FEMALE,
            bsa_m2=bsa,
        )
    else:
        gfr = base["gfr_ml_min"]
    if profile.mild_renal_impairment:
        gfr *= 0.60
    if profile.moderate_renal_impairment:
        gfr *= 0.35
    return gfr / 1000.0  # to L/min


def _resolve_hematocrit(profile: HumanProfile) -> float:
    base = _age_adjust_defaults(profile.age_y, profile.sex)
    return profile.hematocrit if profile.hematocrit is not None else base["hematocrit"]


def _resolve_cardiac_output(profile: HumanProfile) -> float:
    bsa = mosteller_bsa(profile.height_cm, profile.weight_kg)
    if profile.cardiac_index_l_min_m2 is not None:
        ci = profile.cardiac_index_l_min_m2
    else:
        ci = 3.10 if profile.sex is Sex.MALE else 2.90
    co = ci * bsa
    if profile.heart_failure:
        co *= 0.80
    return co


def _organ_volumes(profile: HumanProfile) -> dict[str, float]:
    w = profile.weight_kg
    views: dict[str, float] = {}
    for name in REF_VOLUME_L:
        exponent = _VOLUME_EXPONENT[name]
        sex_factor = _SEX_FACTOR[profile.sex][name]
        views[name] = REF_VOLUME_L[name] * sex_factor * (w / REFERENCE_MALE_WEIGHT_KG) ** exponent
    return views


def _organ_flows(profile: HumanProfile, co: float) -> dict[str, float]:
    del profile  # flows scale with cardiac output only
    ratio = co / CARDIAC_OUTPUT_REF_L_MIN
    # Liver flow is total hepatic inflow (arterial + portal); Gut and spleen
    # drain into the portal vein and their vena-cava outflow is counted once,
    # via the liver. The sum of vena-cava-draining flows therefore equals CO.
    flows = {name: REF_FLOW_L_MIN[name] * ratio for name in REF_FLOW_L_MIN}
    return flows


def glom_filtration_clearance(gfr_ml_min: float, fup: float, fe_unchanged: float = 1.0) -> float:
    """Passive glomerular-filtration component of renal clearance (L/h).

    For a predominantly filtered drug, renal clearance CLr = fup * GFR projected
    to L/h, weighted by the fraction of the administered dose excreted unchanged
    (``fe_unchanged``).  Secretory/reabsorptive balance that is not measured is
    folded into ``fe_unchanged`` (default 1.0 = pure filtration).
    """
    if gfr_ml_min <= 0:
        raise ValueError("gfr_ml_min must be positive")
    if not (0.0 < fup <= 1.0):
        raise ValueError("fup must be in (0, 1]")
    if not (0.0 <= fe_unchanged <= 1.0):
        raise ValueError("fe_unchanged must be in [0, 1]")
    return fup * (gfr_ml_min / 1000.0) * 60.0 * fe_unchanged


def build_human(profile: HumanProfile) -> HumanPhysiology:
    """Resolve a sparse :class:`HumanProfile` into a full physiology set."""
    sex = profile.sex
    co = _resolve_cardiac_output(profile)
    gfr_l_min = _resolve_gfr(profile)

    volumes = _organ_volumes(profile)
    flows = _organ_flows(profile, co)

    if profile.mild_hepatic_impairment:
        volumes["liver"] *= 0.80
        flows["liver"] *= 0.90
    if profile.moderate_hepatic_impairment:
        volumes["liver"] *= 0.65
        flows["liver"] *= 0.85

    # Close the circulation exactly. The vena-cava return is the sum over organs
    # draining directly into it (all except lung, gut and spleen): the gut and
    # spleen flow into the portal vein and are counted once, inside the liver's
    # total flow.  By construction this sum also equals the total arterial
    # demand (direct organs + hepatic arterial + portal organ uptake), so lung,
    # arterial and venous exchange can all be driven by this one flow and the
    # PBPK mass balance is exactly conservative.
    vc_return = sum(flows[t] for t in flows if t != "lung" and t not in PORTAL_DRAINING_TISSUES)
    cardiac_output = float(vc_return)
    flows["lung"] = cardiac_output

    proteins = _resolve_proteins(profile)

    return HumanPhysiology(
        sex=sex,
        age_y=profile.age_y,
        height_cm=profile.height_cm,
        weight_kg=profile.weight_kg,
        bsa_m2=mosteller_bsa(profile.height_cm, profile.weight_kg),
        cardiac_output_l_min=cardiac_output,
        hematocrit=_resolve_hematocrit(profile),
        albumin_g_l=proteins["albumin_g_l"],
        agp_g_l=proteins["agp_g_l"],
        gfr_l_min=gfr_l_min,
        organ_volume=volumes,
        organ_flow=flows,
    )
