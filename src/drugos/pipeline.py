"""Sequential end-to-end pipeline driver (doc/03 2.4, stages 0-5).

``RunSpec`` is the single self-contained input contract for a full DrugOS
simulation: molecule + physiology profile + dose plan + empirical human PK
(clearance/fup/bp) plus the Stage-3 safety panel.  ``run_pipeline`` walks the
contract left to right — PBPK (pk), occupancy panel, pathway, organ QST
(liver/cardiac/kidney), clinical biomarker grading and composite toxicity —
and returns a ``RunResult`` whose ``to_contract`` is the serializable report
document consumed by ``drugos.report``, the CLI and the web playground.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, fields, replace
from typing import Any

import numpy as np

from drugos.clinical.biomarkers import (
    BIOMARKERS,
    BiomarkerGrade,
    grade_absolute,
    grade_timeseries,
    summarize_kidney,
)
from drugos.clinical.toxicity import Endpoint, ToxicityReport, score_toxicity
from drugos.inputs.models import DosePlan, HumanProfile, Molecule, Route, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.parse_structure import parse_structure
from drugos.inputs.resolve_human import resolve_human
from drugos.organ.base import NDArray, free_mg_l_to_nm
from drugos.organ.cardiac import (
    CardiacParams,
    CardiacResult,
    simulate_cardiac,
    sympathetic_tone_from_emax,
)
from drugos.organ.cns import CnsParams, CnsResult, kpu_brain_from_bbb, simulate_cns
from drugos.organ.feedback import OrganFeedback, apply_pk_scaling, feedback_from_results
from drugos.organ.kidney import KidneyResult, simulate_kidney
from drugos.organ.liver import LiverParams, LiverTrajectory, liver_params_from_panel, simulate_liver
from drugos.pathway.sbml_pathway import simulate_sbml_pathway
from drugos.pathway.simulator import PathwayResult
from drugos.pk.admet import AdmetOutput
from drugos.pk.partitions import partition_from_molecule
from drugos.pk.pbpk_build import (
    TISSUE_LIST,
    AbsorptionParams,
    CypTerm,
    PBPKModel,
    TargetBinding,
    absorption_rate_from_fa,
)
from drugos.pk.physiology import HumanPhysiology, build_human, glom_filtration_clearance
from drugos.pk.simulate import PBPKResult, PkMetrics, compute_pk_metrics, simulate_pbpk
from drugos.pk.skin import SkinLayers
from drugos.rbridge import RVerify
from drugos.rbridge import verify_pk as _r_verify_pk
from drugos.reliability import (
    EmpiricalObservations,
    Measurements,
    agreement_rows,
    apply_measurements,
    reliability_of,
)
from drugos.target.dti import (
    cns_ic50_nm as _dti_cns_ic50_nm,
)
from drugos.target.dti import (
    no_public_data_sites,
    resolve_offtarget_panel,
)
from drugos.target.occupancy import (
    PanelEngagement,
    TargetOccupancyResult,
    simulate_occupancy,
    simulate_panel,
)
from drugos.target.resolver import bind_site, kd_from_score
from drugos.target.targets import Target, safety_panel

_NM_PER_MG = 1e6
_DEFAULT_SCR_BASE_UMOL_L = 80.0
_ALT_ULN_U_L = 40.0
_BILI_ULN_MG_DL = 1.0
_CNS_IC50_DEFAULT_NM = 1.0e5

_REALISM_TERMS: tuple[tuple[str, Callable[[RunSpec], bool]], ...] = (
    ("tubular secretion", lambda s: s.cl_sec_l_h > 0),
    ("biliary excretion + enterohepatic recirculation", lambda s: s.cl_bil_l_h > 0),
    ("first-pass gut-wall extraction", lambda s: s.gut_extraction_eg > 0),
    (
        "saturable (Michaelis-Menten) hepatic clearance",
        lambda s: s.hepatic_vmax_mg_h is not None,
    ),
    ("per-CYP abundance-scaled kinetics", lambda s: bool(s.cyp_terms)),
    ("TMDD target binding (native mass balance)", lambda s: bool(s.target_binding)),
    ("multi-layer transdermal skin permeation", lambda s: s.skin_layers is not None),
    (
        "immune-mediated DILI axis",
        lambda s: s.dili_immune_ic50_nm is not None or s.dili_immune_weight is not None,
    ),
    (
        "sympathetic-suppression cardiac branch",
        lambda s: s.beta_block_ic50_nm is not None,
    ),
    ("organ-feedback loop (coupled clearance)", lambda s: s.feedback_loop > 0),
    ("SC/IM depot absorption", lambda s: s.sc_im_ka_per_h is not None),
    (
        "multi-segment (ACAT) small-intestine absorption",
        lambda s: s.si_segments is not None and s.si_segments >= 2,
    ),
)

#: Realism scalars a *measured* value (including a measured zero — a
#: determined absence) legitimately satisfies.  A user-supplied measured zero
#: for tubular/biliary clearance or gut-wall extraction is a true parameter,
#: never a silent opt-out (doc/07 D26): the term is disclosed as measured
#: rather than dropped from the full-fidelity contract.
_MEASURED_TERM_FIELDS: dict[str, str] = {
    "tubular secretion": "cl_sec_l_h",
    "biliary excretion + enterohepatic recirculation": "cl_bil_l_h",
    "first-pass gut-wall extraction": "gut_extraction_eg",
}

# Full-fidelity auto-anchor priors (doc/12 §7.2): values synthesized from the
# ADMET-AI evidence vector + phys-chem priors so the novel-drug path engages
# every realism term by default.  A null-effect anchor is the honest
# pharmacology of a non-blocker (IC50 several orders above any therapeutic
# free exposure), never an assertion of blockade.
_MM_KM_MG_L_FULL = 1.0
_BILIARY_CLEARANCE_FRACTION_LOG_P_HIGH = 0.30
_BILIARY_CLEARANCE_FRACTION_LOG_P_LOW = 0.15
_BILIARY_CLEARANCE_FRACTION_PRIOR = 0.20
_DILI_IMMUNE_IC50_FRACTION = 0.30
_NULL_IMMUNE_IC50_NM = 1.0e6
_NULL_SYMPATHETIC_IC50_NM = 1.0e6
_SI_SEGMENTS_FULL = 3
_DEPOT_KA_PRIOR_1H = 0.08


class RealismError(ValueError):
    """Raised when a full-fidelity run would silently degrade.

    A ``fidelity="full"`` run engages every ADME/pharmacology realism term.
    When a term cannot be engaged because the required data or an explicit
    anchor is missing, ``run_pipeline`` raises instead of quietly downgrading
    (doc/07 G5/G6, doc/12 §7.2).  The caller either supplies the missing
    values or explicitly opts into ``fidelity="baseline"``, whose every
    non-engaged term is then disclosed in the trust record.
    """


def _engaged_realism_terms(spec: RunSpec) -> list[str]:
    out = [name for name, pred in _REALISM_TERMS if pred(spec)]
    meas = spec.measurements
    if meas is not None:
        for name, field in _MEASURED_TERM_FIELDS.items():
            if name not in out and getattr(meas, field, None) is not None:
                out.append(f"{name} (measured)")
    return out


def _available_realism_terms(spec: RunSpec) -> list[str]:
    return [name for name, pred in _REALISM_TERMS if not pred(spec)]


_HEPATIC_PAIR_TERMS = (
    "saturable (Michaelis-Menten) hepatic clearance",
    "per-CYP abundance-scaled kinetics",
)


def _route_applicable(name: str, routes: set[Route]) -> bool:
    if name == "multi-layer transdermal skin permeation":
        return Route.TRANSDERMAL in routes
    if name == "SC/IM depot absorption":
        return bool(routes & {Route.SUBCUTANEOUS, Route.INTRAMUSCULAR})
    if name == "multi-segment (ACAT) small-intestine absorption":
        return Route.ORAL in routes
    if name == "first-pass gut-wall extraction":
        return Route.ORAL in routes
    return True


def _applied_off_terms(spec: RunSpec) -> list[str]:
    """Realism terms applicable to the spec's routes that are not engaged."""
    routes = {ev.route for ev in spec.dose_plan.events}
    meas = spec.measurements

    def measured(name: str) -> bool:
        field = _MEASURED_TERM_FIELDS.get(name)
        return field is not None and meas is not None and getattr(meas, field, None) is not None

    return [
        name
        for name, pred in _REALISM_TERMS
        if _route_applicable(name, routes)
        and name not in _HEPATIC_PAIR_TERMS
        and not pred(spec)
        and not measured(name)
    ]


def _missing_full_fidelity(spec: RunSpec) -> list[str]:
    """Fail-closed check: terms whose absence in a full run is a silent downgrade."""
    off = _applied_off_terms(spec)
    missing = []
    if spec.hepatic_vmax_mg_h is None and not spec.cyp_terms:
        missing.append("saturable (Michaelis-Menten) or per-CYP abundance-scaled hepatic clearance")
    missing.extend(off)
    return missing


def _degraded_realism_terms(spec: RunSpec) -> list[str]:
    """Explicit, never-silent disclosure of every non-engaged realism term.

    In a full-fidelity run the list is empty by construction (the pipeline
    raises ``RealismError`` otherwise).  In an explicitly opted-in baseline
    run every applicable non-engaged term is named as a deliberate opt-out,
    so a report can never hide what the model did not do.  Route-inapplicable
    terms (e.g. skin permeation for an oral dose) are never reported degraded.
    """
    if spec.fidelity == "baseline":
        off = _applied_off_terms(spec)
        if spec.hepatic_vmax_mg_h is None and not spec.cyp_terms:
            off.append("saturable (Michaelis-Menten) or per-CYP abundance-scaled hepatic clearance")
        return [
            f"{term}: explicit fidelity='baseline' opt-out (validated linear-baseline lane)"
            for term in off
        ]
    missing = _missing_full_fidelity(spec)
    if not missing:
        return []
    return [f"{term}: guardrail — full-fidelity run will raise RealismError" for term in missing]


def _assert_full_fidelity(spec: RunSpec) -> None:
    """Fail closed: no silent downgrade in a full-fidelity run (doc/07 G5/G6)."""
    if spec.fidelity == "baseline":
        return
    missing = _missing_full_fidelity(spec)
    if missing:
        raise RealismError(
            "full-fidelity run refused: realism terms would silently stay off -> "
            + "; ".join(missing)
            + ". Supply the missing values/estimates, or set fidelity='baseline' "
            "for an explicit linear-baseline run whose every off term is disclosed."
        )


def fidelity_provenance(
    spec: RunSpec,
    exposure: ExposureProfile,
    no_public_data_sites: tuple[str, ...] = (),
) -> dict[str, Any]:
    """Serialize the fidelity provenance (trust record) of a full-chain run.

    Every report must disclose where each admitted value came from (doc/08
    risk register, doc/07 G4/G5): the production anchors wired into organ
    readouts, the opt-in realism terms engaged vs left off, the panel sites
    still on class-typical low-confidence priors, and the ADMET-AI heads used
    as ML estimates. G5 forbids silent fallback: every term that contributes
    to a readout is engaged in this record and accounts for its own value, and
    no term is substituted by another without being named. The ``reliability``
    block classifies
    the run into its predictive regime (novel molecule / partial evidence /
    off-label route / extrapolated dose / validated in-range / fully measured)
    so the DISCLAIMER §2 reliability clause is an executable per-run disclosure
    (doc/07 D25/D26).
    """
    anchors = [
        "CKD-EPI 2021 race-free GFR baseline (R-6)",
        "de Bruijn & Rietjens 2024 GCDCA bile-acid cholestasis PBK (R-7)",
    ]
    if spec.admet is not None and spec.admet.BBB is not None:
        anchors.append("ADMET-AI BBB_Martins brain-partition head (R-4)")
    if spec.admet is not None and spec.admet.hERG is not None:
        anchors.append("corpus-calibrated hERG P->KD sieve (R-8)")
    ml_estimates = None
    if spec.admet is not None:
        ml_estimates = [
            f.name
            for f in fields(spec.admet)
            if f.name != "raw" and getattr(spec.admet, f.name) is not None
        ]
    return {
        "policy": (
            "no silent fallback (G5/G6); every admissible source and every "
            "non-engaged realism term is named in this record"
        ),
        "fidelity": spec.fidelity,
        "anchors_wired": anchors,
        "cns_grading_anchored": exposure.cns_anchored,
        "mechanism_terms_engaged": _engaged_realism_terms(spec),
        "mechanism_terms_degraded": _degraded_realism_terms(spec),
        "estimates": (
            [{"term": k, "basis": v} for k, v in sorted(spec.estimate_basis.items())] or None
        ),
        "no_public_data_sites": sorted(
            no_public_data_sites or (t.name for t in spec.panel if t.low_confidence)
        ),
        "admet_ml_estimates": ml_estimates,
        "reliability": reliability_of(spec, scaffold=benchmark_data(spec.molecule.name or "")),
    }


@dataclass(slots=True)
class RunSpec:
    """Self-contained inputs for a full pipeline run."""

    molecule: Molecule
    profile: HumanProfile
    dose_plan: DosePlan
    cl_hep_l_h: float
    cl_renal_l_h: float
    mw: float
    fup: float
    bp: float = 1.0
    panel: tuple[Target, ...] = field(default_factory=safety_panel)
    admet: AdmetOutput | None = None
    liver_params: LiverParams | None = None
    qt_ic50_nm: float | None = None
    dili_ic50_nm: float | None = None
    dili_immune_ic50_nm: float | None = None
    dili_immune_weight: float | None = None
    cns_ic50_nm: float | None = None
    beta_block_ic50_nm: float | None = None
    tmax_h: float = 48.0
    n_eval: int = 601
    include_pathway: bool = True
    feedback_loop: int = 0
    sc_im_ka_per_h: float | None = None
    fa: float | None = None
    cl_sec_l_h: float = 0.0
    cl_bil_l_h: float = 0.0
    bile_emptying_1h: float = 1.0
    hepatic_vmax_mg_h: float | None = None
    hepatic_km_mg_l: float | None = None
    cyp_terms: tuple[CypTerm, ...] = ()
    target_binding: tuple[TargetBinding, ...] = ()
    gut_extraction_eg: float = 0.0
    skin_layers: SkinLayers | None = None
    name: str = "compound"
    fidelity: str = "full"
    si_segments: int | None = None
    estimate_basis: dict[str, str] = field(default_factory=dict)
    measurements: Measurements | None = None
    empirical: EmpiricalObservations | None = None

    def __post_init__(self) -> None:
        if self.mw <= 0:
            raise ValueError("mw must be positive")
        if not (0.0 < self.fup <= 1.0):
            raise ValueError("fup must be in (0, 1]")
        if self.feedback_loop < 0:
            raise ValueError("feedback_loop must be non-negative")
        if self.sc_im_ka_per_h is not None and self.sc_im_ka_per_h <= 0:
            raise ValueError("sc_im_ka_per_h must be positive")
        if self.fa is not None and not (0.0 < self.fa <= 1.0):
            raise ValueError("fa must be in (0, 1]")
        if self.qt_ic50_nm is not None and self.qt_ic50_nm <= 0:
            raise ValueError("qt_ic50_nm must be positive")
        if self.dili_ic50_nm is not None and self.dili_ic50_nm <= 0:
            raise ValueError("dili_ic50_nm must be positive")
        if self.dili_immune_ic50_nm is not None and self.dili_immune_ic50_nm <= 0:
            raise ValueError("dili_immune_ic50_nm must be positive")
        if self.dili_immune_weight is not None and not 0.0 <= self.dili_immune_weight <= 1.0:
            raise ValueError("dili_immune_weight must be in [0, 1]")
        if self.cns_ic50_nm is not None and self.cns_ic50_nm <= 0:
            raise ValueError("cns_ic50_nm must be positive")
        if self.beta_block_ic50_nm is not None and self.beta_block_ic50_nm <= 0:
            raise ValueError("beta_block_ic50_nm must be positive")
        for label, value in (
            ("cl_sec_l_h", self.cl_sec_l_h),
            ("cl_bil_l_h", self.cl_bil_l_h),
            ("bile_emptying_1h", self.bile_emptying_1h),
        ):
            if value < 0 or value != value:  # NaN guard
                raise ValueError(f"{label} must be non-negative; got {value}")
        if self.bile_emptying_1h <= 0:
            raise ValueError("bile_emptying_1h must be positive")
        if (self.hepatic_vmax_mg_h is None) != (self.hepatic_km_mg_l is None):
            raise ValueError("hepatic_vmax_mg_h and hepatic_km_mg_l must be set together")
        if self.hepatic_km_mg_l is not None and self.hepatic_km_mg_l <= 0:
            raise ValueError("hepatic_km_mg_l must be positive")
        if self.cyp_terms and self.hepatic_vmax_mg_h is not None:
            raise ValueError(
                "cyp_terms replaces the lumped hepatic clearance; "
                "hepatic_vmax_mg_h/hepatic_km_mg_l must be left unset"
            )
        if not (0.0 <= self.gut_extraction_eg <= 0.9):
            raise ValueError("gut_extraction_eg must be in [0, 0.9]")
        if self.fidelity not in ("full", "baseline"):
            raise ValueError('fidelity must be "full" or "baseline"; got {self.fidelity!r}')
        if self.si_segments is not None and self.si_segments < 1:
            raise ValueError(f"si_segments must be >= 1 or None; got {self.si_segments}")
        for tb in self.target_binding:
            if tb.tissue not in TISSUE_LIST:
                raise ValueError(
                    f"target_binding tissue {tb.tissue!r} not in model tissues {TISSUE_LIST}"
                )


@dataclass(slots=True)
class OrganStage:
    """Stage-4 organ outcomes plus their graded biomarker rows."""

    liver: LiverTrajectory
    cardiac: CardiacResult
    kidney: KidneyResult
    cns: CnsResult
    biomarkers: dict[str, list[BiomarkerGrade]]

    def grade_by_name(self, name: str) -> BiomarkerGrade | None:
        for rows in self.biomarkers.values():
            for row in rows:
                if row.name == name:
                    return row
        return None


@dataclass(slots=True)
class ExposureProfile:
    """Free-exposure anchors used by the clinical fusion (doc/05 5.2)."""

    plasma_cmax_unbound_nm: float
    liver_cmax_free_nm: float
    kidney_cmax_free_nm: float
    qt_ic50_nm: float | None
    dili_ic50_nm: float | None
    cns_ic50_nm: float | None
    dili_immune_ic50_nm: float | None = None
    cns_anchored: bool = False
    cns_kpu_brain: float = 1.0
    beta_block_ic50_nm: float | None = None


@dataclass(slots=True)
class RunResult:
    """Full pipeline output for one compound-dose scenario."""

    name: str
    spec: RunSpec
    physiology: HumanPhysiology
    pk: PBPKResult
    metrics: PkMetrics
    panel: PanelEngagement
    primary_signal: TargetOccupancyResult | None
    pathway: PathwayResult | None
    organ: OrganStage
    exposure: ExposureProfile
    toxicity: ToxicityReport
    verdict: str
    r_verify: RVerify | None = None
    empirical_agreement: dict[str, Any] | None = None
    panel_no_data_sites: tuple[str, ...] = ()

    def to_contract(self) -> dict[str, Any]:
        """Serialize the full run to the doc/03 report contract."""
        liver = self.organ.liver
        cardiac = self.organ.cardiac
        kidney = self.organ.kidney
        cns = self.organ.cns
        biomarkers = {k: [b.to_dict() for b in rows] for k, rows in self.organ.biomarkers.items()}
        return {
            "manifest": {
                "name": self.name,
                "route": self.pk.route.value,
                "dose_mg": self.pk.dose_mg,
            },
            "pk": {
                "cmax_mg_l": self.metrics.cmax_mg_l,
                "auc_last_mg_h_l": self.metrics.auc_last_mgh_l,
                "auc_inf_mg_h_l": self.metrics.auc_inf_mgh_l,
                "tmax_h": self.metrics.tmax_h,
                "bioavailability_f": round(self.metrics.f_abs, 4),
                "plasma_total": self.pk.plasma_total.tolist(),
                "plasma_free": self.pk.plasma_free.tolist(),
                "t_h": self.pk.t.tolist(),
                "unbound_tissues": {k: v.tolist() for k, v in self.pk.unbound_tissues.items()},
            },
            "occupancy": {
                "primary": (self.primary_signal.peak_occupancy if self.primary_signal else 0.0),
                "ranked": [
                    {
                        "site": name,
                        "peak": res.peak_occupancy,
                        "time_at_target_h": res.time_at_target_h,
                    }
                    for name, res in self.panel.ranked_by_time_at_target()
                ],
            },
            "pathway": {
                "available": self.pathway is not None,
                "readout": (
                    ro
                    if self.pathway is not None and (ro := self.pathway.model.readout) is not None
                    else None
                ),
                "readout_peak": (
                    float(np.max(self.pathway.concentrations[ro]))
                    if self.pathway is not None and (ro := self.pathway.model.readout) is not None
                    else 0.0
                ),
            },
            "organ": {
                "liver": {
                    "dili_grade": liver.dili_grade,
                    "hy_law": liver.hy_law,
                    "peak_alt_uln": liver.peak_alt_uln,
                    "peak_ast_u_l": liver.peak_ast_u_l,
                    "peak_bilirubin_uln": liver.peak_bilirubin_uln,
                    "max_dead_frac": liver.max_dead_frac,
                    "regen_scale_min": round(float(np.min(liver.regen_scale)), 3),
                    "regen_scale_peak": round(float(np.max(liver.regen_scale)), 3),
                },
                "cardiac": {
                    "delta_qtc_ms": float(np.max(cardiac.delta_qtc_ms)),
                    "qtc_peak_ms": float(np.max(cardiac.qtc_ms)),
                    "tdpr_band": str(cardiac.tdpr_band),
                    "map_mmhg": cardiac.map_mmhg,
                },
                "kidney": {
                    "aki_grade": kidney.aki_grade,
                    "peak_scr_ratio": kidney.peak_scr_ratio,
                    "min_gfr_ml_min": kidney.min_gfr_ml_min,
                    "peak_kim1_xunl": round(
                        float(np.max(1.0 + _KIM1_RISE_PER_INJURY * kidney.injury)), 3
                    ),
                },
                "cns": {
                    "peak_brain_free_nm": round(float(cns.peak_brain_free_nm), 3),
                    "exposure_ratio": round(cns.exposure_ratio, 3),
                    "grade": cns.grade,
                    "kpu_brain": self.exposure.cns_kpu_brain,
                },
                "trajectories": {
                    "t_h": liver.t_h.tolist(),
                    "liver_alt_U_L": liver.alt_u_l.tolist(),
                    "liver_bilirubin_mg_dL": liver.bilirubin_mg_dl.tolist(),
                    "qtc_ms": cardiac.qtc_ms.tolist(),
                    "delta_qtc_ms": cardiac.delta_qtc_ms.tolist(),
                    "gfr_ml_min": _on_x(liver.t_h, kidney.t_h, kidney.gfr_ml_min).tolist(),
                    "scr_ratio": _on_x(liver.t_h, kidney.t_h, kidney.scr_ratio).tolist(),
                    "brain_free_nm": _on_x(liver.t_h, cns.t_h, cns.brain_free_nm).tolist(),
                    "regen_scale": liver.regen_scale.tolist(),
                },
            },
            "clinical": {
                "biomarkers": biomarkers,
                "toxicity": self.toxicity.to_dict(),
                "verdict": self.verdict,
                "exposure": {
                    "plasma_cmax_unbound_nm": round(self.exposure.plasma_cmax_unbound_nm, 3),
                    "liver_cmax_free_nm": round(self.exposure.liver_cmax_free_nm, 3),
                    "kidney_cmax_free_nm": round(self.exposure.kidney_cmax_free_nm, 3),
                    "qt_ic50_nm": self.exposure.qt_ic50_nm,
                    "dili_ic50_nm": self.exposure.dili_ic50_nm,
                    "dili_immune_ic50_nm": self.exposure.dili_immune_ic50_nm,
                    "cns_ic50_nm": self.exposure.cns_ic50_nm,
                    "beta_block_ic50_nm": self.exposure.beta_block_ic50_nm,
                },
            },
            "trust": self._trust_record(),
            "r_verify": None if self.r_verify is None else self.r_verify.to_dict(),
        }

    def _trust_record(self) -> dict[str, Any]:
        """The ``trust`` contract section (fidelity provenance + agreement)."""
        trust = fidelity_provenance(self.spec, self.exposure, self.panel_no_data_sites)
        if self.empirical_agreement is not None:
            trust["empirical_agreement"] = self.empirical_agreement
        return trust


def _mg_l_to_nm(c_mg_l: float, mw: float) -> float:
    return c_mg_l * _NM_PER_MG / mw


def _on_x(x_ref: NDArray, x_src: NDArray, y_src: NDArray) -> NDArray:
    """Resample ``y_src`` (on grid ``x_src``) onto the reference grid ``x_ref``."""
    xi = np.asarray(x_ref, dtype=float)
    xs = np.asarray(x_src, dtype=float)
    out: NDArray = np.interp(xi, xs, np.asarray(y_src, dtype=float))
    return out


def _herg_occupancy(
    t_h: NDArray,
    c_free_mg_l: NDArray,
    mw: float,
    qt_ic50_nm: float | None,
    n_eval: int = 601,
) -> NDArray:
    """hERG occupancy from the effective KD (override, ADMET-AI sieve, panel prior)."""
    if qt_ic50_nm is None:
        return np.zeros_like(t_h)
    site = Target(
        name="hERG (Kv11.1)",
        kd_nm=qt_ic50_nm,
        r0_nm=0.05,
        reference="effective hERG affinity (override/sieved/panel prior)",
    )
    occ = simulate_occupancy(t_h, c_free_mg_l, site, mw=mw, n_eval=n_eval)
    return np.interp(t_h, occ.t_h, occ.occupancy)


def _dili_ic50(spec: RunSpec) -> float | None:
    if spec.dili_ic50_nm is not None:
        return spec.dili_ic50_nm
    cands = [t.kd_nm for t in spec.panel if "BSEP" in t.name or "Mitochondrial" in t.name]
    return min(cands) if cands else None


def _effective_herg_kd_nm(spec: RunSpec) -> float | None:
    """hERG KD used by Stage-2/4 and the exposure anchors (doc/05 2.2-2.3).

    Precedence: a per-compound ``qt_ic50_nm`` override is absolute; otherwise
    a novel molecule scored by ADMET-AI maps its predicted hERG-head
    probability onto a continuous, corpus-calibrated KD that is never more
    potent than the panel prior and degrades toward the weak-kD floor as the
    blocker confidence falls (doc/12 D10, path B).  Benchmarks / no-ADMET and
    missing-head runs fall back to the panel class prior unchanged.
    """
    if spec.qt_ic50_nm is not None:
        return spec.qt_ic50_nm
    prior = next((t.kd_nm for t in spec.panel if "hERG" in t.name), None)
    if prior is None or spec.admet is None or spec.admet.hERG is None:
        return prior
    return kd_from_score(spec.admet.hERG, prior)


def _resolve_panel(spec: RunSpec) -> tuple[tuple[Target, ...], tuple[str, ...]]:
    """Resolve the safety panel for Stage-2 occupancy (doc/05 2.2-2.3, doc/12 D10).

    Site-by-site, three seams (doc/12 L13/L15, path A + path B + the hERG
    anchor):
    - path B heads first: every site whose ADMET-AI head measures the same
      interaction (CYP2D6/3A4/2C9 inhibition) is re-scored from the head,
    - path A resolver: every site with *chemotype support* in the vendored
      ChEMBL snapshot (top-match Tanimoto >= ``MIN_NEIGHBOR_TANIMOTO``) is
      re-bound to its fingerprint-kNN-predicted KD (structure-derived,
      ``low_confidence`` cleared); sites without support stay on their
      disclosed priors,
    - then the hERG site is bound to the effective KD — a measured
      ``qt_ic50_nm`` override stays absolute (a measured anchor, not a prior),
      otherwise the corpus-calibrated head sieve applies.  Keeps occupancy
      (and the pathway drive) consistent with the QT and exposure anchors.
    Returns ``(resolved_panel, no_public_data_sites)`` — computed on the final
    panel — the sites that still sit on disclosed class priors: the sites with
    no public bioactivity row in the snapshot (MRP3/MRP4, mitochondrial
    complexes II-IV, pyruvate carrier) plus any mapped site whose query has no
    chemotype support.
    """
    smiles = spec.molecule.canonical_smiles if spec.molecule is not None else None
    panel = resolve_offtarget_panel(spec.panel, smiles=smiles, admet=spec.admet)[0]
    if spec.qt_ic50_nm is not None:
        panel = bind_site(
            panel,
            "hERG (Kv11.1)",
            float(spec.qt_ic50_nm),
            "measured hERG IC50 override (absolute)",
            low_confidence=False,
        )
    else:
        eff = _effective_herg_kd_nm(spec)
        if eff is not None:
            panel = bind_site(
                panel,
                "hERG (Kv11.1)",
                eff,
                "hERG KD from ADMET-AI head, corpus-calibrated monotone P->KD (doc/12 D10)",
                low_confidence=False,
            )
    return tuple(panel), no_public_data_sites(panel)


def _resolved_cns_ic50_nm(spec: RunSpec) -> float | None:
    """CNS potency anchor for the grading seam (doc/12 L12).

    An explicit ``spec.cns_ic50_nm`` stays absolute; otherwise a real molecule
    is graded against the structure-derived worst-case prediction over the
    CNS-liability set (never more potent than the conservative default).
    ``None`` only when no structure is available.
    """
    if spec.cns_ic50_nm is not None:
        return spec.cns_ic50_nm
    if spec.molecule is not None and spec.molecule.canonical_smiles:
        return _dti_cns_ic50_nm(spec.molecule.canonical_smiles)
    return None


def _qt_ic50(spec: RunSpec) -> float | None:
    return _effective_herg_kd_nm(spec)


def _proliferation_signal(pathway: PathwayResult, t_h: NDArray) -> NDArray | None:
    """ERK readout fold-change over ``t_h`` driving hepatocyte regeneration.

    A missing readout node or a degenerate (non-positive) drug-free baseline
    yields ``None`` (no pathway coupling into the liver); otherwise the readout
    concentration over its baseline interpolated onto ``t_h``.
    """
    ro = pathway.model.readout
    if ro is None:
        return None
    base = pathway.baseline.get(ro, 0.0)
    if not base > 0.0:
        return None
    fc = pathway.concentrations[ro] / base
    return np.interp(t_h, pathway.t_h, fc)


_CARDIAC_TONE_GAIN = 0.2
_CARDIAC_TONE_MAX = 1.3


def _cardiac_tone_scale(fc: NDArray | None) -> float:
    """Mild inotropic/chronotropic tone from the ERK amplification ratio.

    doc/05 4.3: beta/catecholamine pathway effects modulate contractility and
    rate.  The shipped pathway kernel is the ERK/MAPK cascade (R-5), whose
    doubly-phosphorylated-ERK readout is the available proxy (ERK1/2 is
    downstream of beta-adrenergic E-C coupling).  The readout is a *stimulatory*
    concentration: benign drugs sit at or below the drug-free baseline (here
    ~0 fold), so only amplification above baseline (the inotropic branch) lifts
    ``inotropy=chronotropy`` via ``clamp(1 + gain*(fc-1), 1, max)``; a baseline
    or suppressed readout leaves tone exactly 1.0.  Sympathetic *suppression*
    is a separate, saturable Emax axis (`RunSpec.beta_block_ic50_nm` ->
    `CardiacParams.sympathetic_tone`, doc/05 4.3), so a benign molecule with
    no beta-block anchor never disturbs hemodynamic neutrality or the loop
    co-fraction.
    """
    if fc is None:
        return 1.0
    mean = float(np.mean(fc))
    if mean <= 1.0:
        return 1.0
    return float(
        np.clip(
            1.0 + _CARDIAC_TONE_GAIN * (mean - 1.0),
            1.0,
            _CARDIAC_TONE_MAX,
        )
    )


def _panel_model(spec: RunSpec) -> PBPKModel:
    physiology = resolve_human(spec.profile)
    part = partition_from_molecule(
        spec.molecule, fup=spec.fup, bp=spec.bp, hematocrit=physiology.hematocrit
    )
    absorption = _absorption_params(spec)
    return PBPKModel(
        physiology=physiology,
        partition=part,
        bp=spec.bp,
        fup=spec.fup,
        cl_hep_l_h=spec.cl_hep_l_h,
        cl_renal_l_h=spec.cl_renal_l_h,
        cl_sec_l_h=spec.cl_sec_l_h,
        hepatic_vmax_mg_h=spec.hepatic_vmax_mg_h,
        hepatic_km_mg_l=spec.hepatic_km_mg_l,
        cyp_terms=spec.cyp_terms,
        cl_bil_l_h=spec.cl_bil_l_h,
        k_bile_emptying_1h=spec.bile_emptying_1h,
        target_binding=spec.target_binding,
        mw_g_per_mol=spec.mw,
        dose_plan=spec.dose_plan,
        absorption=absorption,
    )


def _admet_fa(admet: AdmetOutput | None) -> float | None:
    """Fraction-absorbed target from an ADMET-AI call (HIA, else bioavailable_Ma)."""
    if admet is None:
        return None
    for attr in ("HIA", "bioavailable_Ma"):
        value = getattr(admet, attr, None)
        if value is not None and 0.0 < value <= 1.0:
            return float(value)
    return None


def _admet_solubility_mg_ml(admet: AdmetOutput | None, mw: float) -> float | None:
    """Solubility (mg/mL) from the ADMET-AI ``logS`` (log10 molar) call."""
    log_s = getattr(admet, "log_s", None) if admet is not None else None
    if log_s is None:
        return None
    log_s = float(log_s)
    if not np.isfinite(log_s):
        return None
    return float(10.0**log_s) * float(mw) / 1000.0


def _absorption_params(spec: RunSpec) -> AbsorptionParams:
    """Permeability/Fa-gated, solubility-limited first-order absorption.

    SC/IM/transdermal runs use the depot override; oral runs without an
    explicit rate gate the small-intestine absorption constant to the
    fraction-absorbed target ``fa`` (permeability-gated absorption, doc/05
    §1.4), so a low-permeability compound absorbs slower and loses more to the
    colon/feces sink while a fast one absorbs in the proximal gut.  For novel
    molecules scored by ADMET-AI the population is an explicit ``fa`` before
    an HIA-derived target (doc/05 1.2), so new chemistry never falls back to a
    fixed 0.55 h^-1 blind.  A measured/predicted ``logS`` sets a
    solubility-limited dissolution cap on the dissolved lumen pool (novel/ADMET
    runs only), so a low-solubility dose spills undissolved mass to the
    colon/feces sink.
    """
    if spec.sc_im_ka_per_h is not None or spec.skin_layers is not None:
        return AbsorptionParams(
            k_depot_absorption=spec.sc_im_ka_per_h or 0.15,
            skin_layers=spec.skin_layers,
            si_segments=spec.si_segments,
        )
    has_oral = any(ev.route is Route.ORAL for ev in spec.dose_plan.events)
    if not has_oral:
        return AbsorptionParams(
            gut_extraction_eg=spec.gut_extraction_eg,
            si_segments=spec.si_segments,
        )
    solubility = _admet_solubility_mg_ml(spec.admet, spec.mw)
    if spec.fa is not None:
        return AbsorptionParams(
            k_si_absorption=absorption_rate_from_fa(spec.fa),
            solubility_mg_ml=solubility,
            gut_extraction_eg=spec.gut_extraction_eg,
            si_segments=spec.si_segments,
        )
    admet_fa = _admet_fa(spec.admet)
    if admet_fa is not None:
        return AbsorptionParams(
            k_si_absorption=absorption_rate_from_fa(admet_fa),
            solubility_mg_ml=solubility,
            gut_extraction_eg=spec.gut_extraction_eg,
            si_segments=spec.si_segments,
        )
    return AbsorptionParams(gut_extraction_eg=spec.gut_extraction_eg, si_segments=spec.si_segments)


def _organ_state(
    pk: PBPKResult,
    model: PBPKModel,
    spec: RunSpec,
    *,
    include_pathway: bool = True,
) -> tuple[
    PanelEngagement | None,
    TargetOccupancyResult | None,
    LiverTrajectory,
    NDArray,
    CardiacResult,
    KidneyResult,
    CnsResult,
    PathwayResult | None,
    float,
    CnsParams,
    tuple[str, ...],
]:
    """Run the downstream stages (occupancy -> pathway -> organ QST) on ``pk``."""
    liver_free = pk.unbound_tissues["liver"]

    # The hERG site enters every downstream consumer (occupancy/primary signal,
    # QT drive, exposure anchors) at one effective KD — a measured override or,
    # for ADMET-AI-scored novel molecules, the predicted-hERG sieve.
    panel: PanelEngagement | None = None
    primary_signal: TargetOccupancyResult | None = None
    no_data_sites: tuple[str, ...] = ()
    if include_pathway:
        # The hERG site enters every downstream consumer (occupancy/primary
        # signal, QT drive, exposure anchors) at one effective KD — the panel
        # is first resolved per-site from ADMET-AI heads and the ChEMBL-kNN
        # resolver, then the hERG override binds everywhere at once.
        resolved_panel, no_data_sites = _resolve_panel(spec)
        panel = simulate_panel(pk.t, liver_free, resolved_panel, mw=spec.mw, n_eval=spec.n_eval)
        primary_signal = max(
            (r for r in panel.results.values()),
            key=lambda r: r.time_at_target_h,
            default=None,
        )

    liver_params = spec.liver_params or liver_params_from_panel(spec.panel)
    if spec.dili_immune_ic50_nm is not None:
        immune_weight = spec.dili_immune_weight if spec.dili_immune_weight is not None else 0.5
        liver_params = replace(
            liver_params,
            immune_ic50_nm=spec.dili_immune_ic50_nm,
            immune_weight=immune_weight,
        )

    # Stage-3 pathway before the organ QST: the ERK/MAPK readout fold-change is
    # the proliferation signal that modulates hepatocyte regeneration, closing
    # the occupancy -> pathway -> organ -> phenotype chain (doc/05 4.2).
    # Intermediate feedback-loop iterations skip it (only liver/kidney/cardiac
    # feed the scaling pass); the final organ state runs the full chain.
    pathway: PathwayResult | None = None
    proliferation: NDArray | None = None
    if include_pathway and spec.include_pathway and primary_signal is not None:
        pathway = simulate_sbml_pathway(
            primary_signal.t_h, primary_signal.occupancy, n_eval=spec.n_eval
        )
        proliferation = _proliferation_signal(pathway, pk.t)

    liver = simulate_liver(
        pk.t,
        liver_free,
        spec.mw,
        params=liver_params,
        n_eval=spec.n_eval,
        proliferation_signal=proliferation,
    )

    # hERG channel blockade is at the myocardium, so drive it with the PBPK
    # cardiac (heart) free tissue exposure, not the hepatic one.
    block = _herg_occupancy(
        pk.t,
        pk.unbound_tissues["heart"],
        spec.mw,
        _effective_herg_kd_nm(spec),
        spec.n_eval,
    )
    cardiac_params = None
    if proliferation is not None:
        tone = _cardiac_tone_scale(proliferation)
        cardiac_params = CardiacParams(inotropy=tone, chronotropy=tone)
    if spec.beta_block_ic50_nm is not None:
        heart_free_nm = np.asarray(free_mg_l_to_nm(pk.unbound_tissues["heart"], spec.mw))
        sym_tone = sympathetic_tone_from_emax(
            float(np.mean(heart_free_nm)), spec.beta_block_ic50_nm
        )
        base = cardiac_params or CardiacParams()
        cardiac_params = CardiacParams(
            inotropy=base.inotropy,
            chronotropy=base.chronotropy,
            sympathetic_tone=sym_tone,
        )
    cardiac = simulate_cardiac(
        pk.t, block, model.physiology, params=cardiac_params, n_eval=spec.n_eval
    )

    kidney_free = pk.unbound_tissues["kidney"]
    kidney = simulate_kidney(
        pk.t,
        kidney_free,
        spec.mw,
        gfr_base_ml_min=model.physiology.gfr_ml_min,
        scr_base_umol_l=_DEFAULT_SCR_BASE_UMOL_L,
    )

    # CNS is driven by the PBPK brain compartment free exposure.  The R&R brain
    # partition is already folded in, so the residual kpu_brain scale is 1.0
    # unless ADMET-AI's BBB_Martins head predicts a non-penetrant (then 0.2).
    kpu_brain = kpu_brain_from_bbb(spec.admet.BBB if spec.admet is not None else None)
    resolved_cns = _resolved_cns_ic50_nm(spec)
    cns_ic50 = resolved_cns if resolved_cns is not None else _CNS_IC50_DEFAULT_NM
    cns_params = CnsParams(ic50_nm=cns_ic50, kpu_brain=kpu_brain)
    cns = simulate_cns(
        pk.t,
        pk.unbound_tissues["brain"],
        spec.mw,
        cns_params,
    )

    return (
        panel,
        primary_signal,
        liver,
        block,
        cardiac,
        kidney,
        cns,
        pathway,
        kpu_brain,
        cns_params,
        no_data_sites,
    )


def _organ_feedback(
    model: PBPKModel,
    liver: LiverTrajectory,
    kidney: KidneyResult,
    cardiac: CardiacResult,
) -> OrganFeedback:
    """PK scalings from the current organ state (doc/05 4.5, sequential outer loop).

    Cardiac output is derived from the hemodynamic state actually simulated in
    the cardiac panel: ``co_fraction = CO_simulated / CO_reference``, so a
    reduced-cardiac-output physiology (e.g. low inotropy/chronotropy in heart
    failure) scales perfusion-limited distribution — ``apply_pk_scaling`` then
    cuts organ flows and venous return for the re-run.
    """
    gfr_fraction = kidney.min_gfr_ml_min / model.physiology.gfr_ml_min
    ref_co_l_min = model.physiology.cardiac_output_ml_min / 1000.0
    if ref_co_l_min > 0:
        co_fraction = cardiac.co_l_min / ref_co_l_min
    else:
        co_fraction = 1.0
    return feedback_from_results(liver.max_dead_frac, gfr_fraction, co_fraction)


def run_pipeline(spec: RunSpec) -> RunResult:
    """Execute the full sequential pipeline for ``spec``.

    ``spec.feedback_loop`` (>0) re-runs PK with organ-dysfunction scaling
    (doc/05 4.5): each iteration re-derives hepatic/GFR scalings from the
    current organ state and resolves the PBPK model afresh.
    Measured true parameters (``spec.measurements``) are applied first and
    replace the machine auto-anchors; user-supplied empirical observations
    (``spec.empirical``) are compared against the prediction and disclosed as
    an explicit agreement diagnostic, never absorbed.
    """
    if spec.tmax_h <= 0:
        raise ValueError("tmax_h must be positive")

    spec = apply_measurements(spec)
    _assert_full_fidelity(spec)

    model = _panel_model(spec)
    pk = simulate_pbpk(model, tmax_h=spec.tmax_h, n_eval=spec.n_eval)
    metrics = compute_pk_metrics(pk.t, pk.plasma_total, pk.dose_mg, pk.route)

    for _ in range(spec.feedback_loop):
        (
            _panel,
            _primary,
            liver,
            _block,
            _cardiac,
            kidney,
            _cns,
            _pathway,
            _kpu,
            _cns_params,
            _nds,
        ) = _organ_state(pk, model, spec, include_pathway=False)
        model = apply_pk_scaling(model, _organ_feedback(model, liver, kidney, _cardiac))
        pk = simulate_pbpk(model, tmax_h=spec.tmax_h, n_eval=spec.n_eval)

    (
        panel,
        primary_signal,
        liver,
        block,
        cardiac,
        kidney,
        cns,
        pathway,
        cns_kpu_brain,
        cns_params,
        panel_no_data_sites,
    ) = _organ_state(pk, model, spec)

    physiology = model.physiology
    liver_free = pk.unbound_tissues["liver"]
    kidney_free = pk.unbound_tissues["kidney"]

    metrics = compute_pk_metrics(
        pk.t,
        pk.plasma_total,
        pk.dose_mg,
        pk.route,
        f_abs=pk.bioavailability_f or 1.0,
    )

    organ = OrganStage(
        liver=liver,
        cardiac=cardiac,
        kidney=kidney,
        cns=cns,
        biomarkers=_grade_organ(liver, cardiac, kidney, cns, cns_params),
    )

    resolved_cns = _resolved_cns_ic50_nm(spec)
    exposure = ExposureProfile(
        plasma_cmax_unbound_nm=_mg_l_to_nm(float(np.max(pk.plasma_free)), spec.mw),
        liver_cmax_free_nm=_mg_l_to_nm(float(np.max(liver_free)), spec.mw),
        kidney_cmax_free_nm=_mg_l_to_nm(float(np.max(kidney_free)), spec.mw),
        qt_ic50_nm=_qt_ic50(spec),
        dili_ic50_nm=_dili_ic50(spec),
        dili_immune_ic50_nm=spec.dili_immune_ic50_nm,
        cns_ic50_nm=resolved_cns if resolved_cns is not None else _CNS_IC50_DEFAULT_NM,
        cns_anchored=resolved_cns is not None,
        cns_kpu_brain=cns_kpu_brain,
        beta_block_ic50_nm=spec.beta_block_ic50_nm,
    )

    toxicity = _score_clinical(organ, exposure, spec.admet)
    assert panel is not None  # final organ state always runs the full chain
    empirical_agreement = _empirical_agreement(spec, metrics, liver, cardiac)
    return RunResult(
        name=spec.name,
        spec=spec,
        physiology=physiology,
        pk=pk,
        metrics=metrics,
        panel=panel,
        primary_signal=primary_signal,
        pathway=pathway,
        organ=organ,
        exposure=exposure,
        toxicity=toxicity,
        verdict=_verdict(toxicity),
        r_verify=_r_verify_pk(pk.t, pk.plasma_total, pk.dose_mg),
        empirical_agreement=empirical_agreement,
        panel_no_data_sites=panel_no_data_sites,
    )


def _empirical_agreement(
    spec: RunSpec, metrics: PkMetrics, liver: LiverTrajectory, cardiac: CardiacResult
) -> dict[str, Any] | None:
    """Disclose observed-vs-predicted deviation when observations are given.

    Disagreement with empirical data is the expected state of a mechanistic
    model, not a bug (DISCLAIMER §2): supplied observations are never absorbed
    into the parameterization — they are folded into the trust record as an
    explicit ``empirical_agreement`` diagnostic.
    """
    if spec.empirical is None:
        return None
    predicted = {
        "plasma_cmax_mg_l": metrics.cmax_mg_l,
        "auc_last_mg_h_l": metrics.auc_last_mgh_l,
        "peak_delta_qtc_ms": float(np.max(cardiac.delta_qtc_ms)),
        "peak_alt_uln": liver.peak_alt_uln,
    }
    return {
        "observations": agreement_rows(spec.empirical, predicted),
        "policy": (
            "disagreement with empirical data is the expected state of a "
            "mechanistic model, not a bug — reported, not silenced (DISCLAIMER §2)"
        ),
    }


def _grade_organ(
    liver: LiverTrajectory,
    cardiac: CardiacResult,
    kidney: KidneyResult,
    cns: CnsResult,
    cns_params: CnsParams | None = None,
) -> dict[str, list[BiomarkerGrade]]:
    return {
        "liver": summarize_clinical_liver(liver),
        "cardiac": summarize_clinical_cardiac(cardiac),
        "kidney": summarize_clinical_kidney(kidney),
        "cns": summarize_clinical_cns(cns, cns_params)
        if cns_params is not None
        else summarize_clinical_cns(cns),
    }


def summarize_clinical_liver(liver: LiverTrajectory) -> list[BiomarkerGrade]:
    """Grade the liver trajectory (ALT/AST/bilirubin in multiples of ULN)."""
    return [
        grade_timeseries(liver.t_h, liver.alt_u_l / _ALT_ULN_U_L, BIOMARKERS["ALT"]),
        grade_timeseries(liver.t_h, liver.ast_u_l / _ALT_ULN_U_L, BIOMARKERS["AST"]),
        grade_timeseries(
            liver.t_h, liver.bilirubin_mg_dl / _BILI_ULN_MG_DL, BIOMARKERS["total_bilirubin"]
        ),
    ]


def summarize_clinical_cardiac(cardiac: CardiacResult) -> list[BiomarkerGrade]:
    """Grade the cardiovascular trajectory (QTc/ΔQTc series + HR/MAP)."""
    return [
        grade_timeseries(cardiac.t_h, cardiac.qtc_ms, BIOMARKERS["QTc"]),
        grade_timeseries(cardiac.t_h, cardiac.delta_qtc_ms, BIOMARKERS["delta_QTc"]),
        grade_absolute(cardiac.hr_bpm, BIOMARKERS["heart_rate"]),
        grade_absolute(cardiac.map_mmhg, BIOMARKERS["map"]),
    ]


_KIM1_RISE_PER_INJURY = 8.0


def summarize_clinical_kidney(kidney: KidneyResult) -> list[BiomarkerGrade]:
    """Grade the kidney trajectory (GFR + creatinine ratio + KIM-1)."""
    return summarize_kidney(
        kidney.t_h,
        kidney.gfr_ml_min,
        kidney.scr_ratio,
        kim1_xunl=1.0 + _KIM1_RISE_PER_INJURY * kidney.injury,
    )


def summarize_clinical_cns(cns: CnsResult, params: CnsParams = CnsParams()) -> list[BiomarkerGrade]:
    """Grade the CNS exposure-ratio trajectory (Cmax_brain_free / IC50)."""
    return [
        grade_timeseries(cns.t_h, cns.brain_free_nm / params.ic50_nm, BIOMARKERS["CNS_exposure"])
    ]


def _score_clinical(
    organ: OrganStage, exposure: ExposureProfile, admet: AdmetOutput | None
) -> ToxicityReport:
    qtc = organ.grade_by_name("QTc (Fridericia)")
    delta_qtc = organ.grade_by_name("Delta QTc")
    qt_grade = max(
        qtc.grade if qtc is not None else 0,
        delta_qtc.grade if delta_qtc is not None else 0,
    )
    mechanistic: dict[Endpoint, int] = {
        Endpoint.DILI: organ.liver.dili_grade,
        Endpoint.QT: qt_grade,
        Endpoint.AKI: organ.kidney.aki_grade,
    }
    ic50: dict[Endpoint, float | None] = {
        Endpoint.DILI: exposure.dili_ic50_nm,
        Endpoint.QT: exposure.qt_ic50_nm,
        Endpoint.AKI: None,
        Endpoint.CNS: exposure.cns_ic50_nm if exposure.cns_anchored else None,
    }
    if exposure.cns_anchored:
        mechanistic[Endpoint.CNS] = organ.cns.grade
    structural_mask = None if exposure.cns_anchored else frozenset({Endpoint.CNS})
    return score_toxicity(
        mechanistic,
        cmax_unbound_nm=exposure.plasma_cmax_unbound_nm,
        ic50_nm=ic50,
        admet=admet,
        structural_mask=structural_mask,
    )


def _verdict(toxicity: ToxicityReport) -> str:
    overall = toxicity.overall_risk()
    if overall >= 0.5:
        pct = f"{toxicity.overall_risk():.0%}"
        return f"High composite risk ({pct}, driver {toxicity.overall_driver()})"
    if overall >= 0.3:
        return "Elevated composite risk: monitor on the flagged endpoint(s)"
    return "No elevated composite risk detected"


@dataclass(frozen=True, slots=True)
class _BenchmarkLike:
    """Structural view of a validation ``Benchmark`` (decoupled adapter)."""

    name: str
    smiles: str
    fup: float
    bp: float
    route: str
    dose_mg: float
    tmax_h: float
    n_eval: int
    log_p: float
    pka_acids: list[float]
    pka_bases: list[float]
    published: dict[str, tuple[float, float]]


_BENCHMARKS: dict[str, _BenchmarkLike] = {
    "acetaminophen": _BenchmarkLike(
        name="acetaminophen",
        smiles="CC(=O)Nc1ccc(O)cc1",
        fup=0.8,
        bp=1.17,
        route="oral",
        dose_mg=1000.0,
        tmax_h=24.0,
        n_eval=400,
        log_p=0.51,
        pka_acids=[9.4],
        pka_bases=[],
        published={
            "cl_plasma_l_h": (15.0, 27.0),
            "t_half_h": (1.5, 3.5),
            "f_abs": (0.80, 0.98),
        },
    ),
    "warfarin": _BenchmarkLike(
        name="warfarin",
        smiles="Cc1ccc(C(=O)CC(c2ccccc2O)c2ccc(cc2)C(=O)O)cc1",
        fup=0.013,
        bp=0.55,
        route="oral",
        dose_mg=10.0,
        tmax_h=96.0,
        n_eval=800,
        log_p=2.7,
        pka_acids=[5.0],
        pka_bases=[],
        published={
            "cl_plasma_l_h": (0.10, 0.30),
            "vss_l": (4.0, 16.0),
            "t_half_h": (30.0, 45.0),
        },
    ),
    "midazolam": _BenchmarkLike(
        name="midazolam",
        smiles="Cc1ncn(-c2ccc(Cl)c(c2)C(=O)Nc2cccc(F)c2)c1C",
        fup=0.02,
        bp=0.82,
        route="iv_bolus",
        dose_mg=5.0,
        tmax_h=24.0,
        n_eval=400,
        log_p=3.94,
        pka_acids=[],
        pka_bases=[6.04],
        published={"cl_plasma_l_h": (18.0, 36.0)},
    ),
    "ciprofloxacin": _BenchmarkLike(
        name="ciprofloxacin",
        smiles="C1CN(CCN1)C1=C(F)C(=O)C(C(=O)O)=CN1C2CC2",
        fup=0.7,
        bp=0.83,
        route="oral",
        dose_mg=500.0,
        tmax_h=48.0,
        n_eval=400,
        log_p=0.28,
        pka_acids=[6.1],
        pka_bases=[8.7],
        published={
            "cl_plasma_l_h": (15.0, 27.0),
            "t_half_h": (3.0, 6.0),
            "f_abs": (0.60, 0.95),
            "urine_fraction": (0.40, 0.70),
        },
    ),
    "dofetilide": _BenchmarkLike(
        name="dofetilide",
        smiles="COc1cc(NC(=O)Nc2ccc(NS(=O)(=O)CCN(C)C)cc2)ccc1OC",
        fup=0.36,
        bp=0.86,
        route="oral",
        dose_mg=0.5,
        tmax_h=72.0,
        n_eval=400,
        log_p=2.29,
        pka_acids=[],
        pka_bases=[7.5],
        published={
            "cl_plasma_l_h": (12.0, 20.0),
            "t_half_h": (7.0, 9.0),
            "f_abs": (0.80, 1.0),
            "urine_fraction": (0.60, 0.85),
        },
    ),
}


def benchmark_names() -> list[str]:
    """Benchmark compound names available to the CLI/playground.

    Self-contained registry mirroring ``validation.benchmarks`` (kept here so
    the core package — and its wheel — never imports the validation tree).
    """
    return list(_BENCHMARKS)


def benchmark_data(name: str) -> _BenchmarkLike | None:
    """Return the canonical benchmark ``_BenchmarkLike`` for ``name`` or None."""
    return _BENCHMARKS.get(name)


def _published_band(bench: _BenchmarkLike, key: str) -> tuple[float, float]:
    lo, hi = bench.published[key]
    return lo, hi


def _clearance_from_published(bench: _BenchmarkLike) -> tuple[float, float]:
    """Unbound hepatic/renal clearances from the published CL band."""
    lo, hi = _published_band(bench, "cl_plasma_l_h")
    cl_total = (lo + hi) / 2.0
    renal = bench.published.get("urine_fraction", (0.0, 0.0))
    renal_share = (renal[0] + renal[1]) / 2.0 if renal[1] > 0 else 0.0
    cl_hep = cl_total * (1.0 - renal_share) / bench.fup
    cl_renal = cl_total * renal_share / bench.fup
    return cl_hep, cl_renal


def _other_benchmark_affinities(name: str) -> tuple[float | None, float | None]:
    """Per-compound pharmacology overrides on top of panel class priors.

    Literature anchors (doc/06): dofetilide is a high-affinity hERG blocker
    (IC50 ~2 nM), the other benchmarks (warfarin, ciprofloxacin, midazolam,
    acetaminophen) are weak/negligible hERG blockers and none is a potent
    BSEP cholestasis trigger at therapeutic exposure; high-dose acetaminophen
    is anchored on off-target mitochondrial/BSEP stress (IC50 ~0.5 mM).
    """
    strong_herg = {"dofetilide": 2.0}
    weak_herg = {"warfarin", "ciprofloxacin", "midazolam", "acetaminophen"}
    qt = strong_herg.get(name, 1e6 if name in weak_herg else None)
    dili = 5.0e5 if name == "acetaminophen" else None
    return qt, dili


def _engage_full_fidelity(spec: RunSpec) -> RunSpec:
    """Deterministically engage every realism term for novel drugs and benches.

    Full-fidelity default (doc/12 §7.2): every ADME/pharmacology realism term
    is filled from the ADMET-AI evidence vector and phys-chem priors.  Each
    synthesized value is recorded with its derivation basis in ``estimate_basis``
    so the trust record discloses the estimate and never a silent default.
    Null-effect anchors are the honest pharmacology of a non-engager (a
    non-beta-blocker and a no-immune-trigger molecule have IC50s orders above
    any therapeutic exposure), not fabricated blockade.
    """
    if spec.fidelity != "full":
        return spec
    kw: dict[str, Any] = {"estimate_basis": dict(spec.estimate_basis)}
    basis = kw["estimate_basis"]
    admet = spec.admet

    def admet_field(name: str) -> Any:
        """Read an ADMET-AI head by name, tolerating the runtime namespace."""
        return getattr(admet, name, None) if admet is not None else None

    if spec.hepatic_vmax_mg_h is None and not spec.cyp_terms:
        km = _MM_KM_MG_L_FULL
        vmax = max(spec.cl_hep_l_h * km, 1e-3)
        kw["hepatic_vmax_mg_h"] = vmax
        kw["hepatic_km_mg_l"] = km
        basis["hepatic_vmax_mg_h"] = (
            "lumped saturable MM: Vmax = CL_h*Km, Km = 1 mg/L prior; the low-dose "
            "slope Vmax/Km = CL_h preserves the linear clearance baseline"
        )

    if spec.cl_sec_l_h == 0:
        pgp = admet_field("Pgp")
        if pgp is None:
            pgp = 0.5
        cl_sec = max(spec.cl_renal_l_h * (0.5 + 0.5 * float(pgp)), 1e-4)
        kw["cl_sec_l_h"] = cl_sec
        basis["cl_sec_l_h"] = (
            f"renal tubular secretion = CL_renal x (0.5 + 0.5 x P-gp) with "
            f"P-gp substrate probability {float(pgp):.2f} (ADMET-AI Pgp_Broccatelli "
            "head, or 0.5 neutral prior when absent)"
        )

    if spec.cl_bil_l_h == 0:
        lp = admet_field("log_p_pred")
        if lp is not None and lp >= 2.0:
            frac = _BILIARY_CLEARANCE_FRACTION_LOG_P_HIGH
            lp_desc = f"logP {lp:.2f} >= 2"
        elif lp is not None:
            frac = _BILIARY_CLEARANCE_FRACTION_LOG_P_LOW
            lp_desc = f"logP {lp:.2f} < 2"
        else:
            frac = _BILIARY_CLEARANCE_FRACTION_PRIOR
            lp_desc = "logP unavailable"
        cl_bil = max(spec.cl_hep_l_h * frac, 1e-4)
        kw["cl_bil_l_h"] = cl_bil
        basis["cl_bil_l_h"] = (
            f"biliary excretion = {frac:.0%} x CL_h phys-chem prior ({lp_desc}) + "
            "enterohepatic recirculation (k_bile_emptying stays at 1.0/h)"
        )

    if spec.gut_extraction_eg == 0:
        cyp34 = admet_field("cyp3a4_inhibitor")
        if cyp34 is None:
            cyp34 = 0.5
        eg = min(max(0.10 + 0.30 * float(cyp34), 0.05), 0.60)
        kw["gut_extraction_eg"] = eg
        basis["gut_extraction_eg"] = (
            f"first-pass gut-wall extraction = 0.10 + 0.30 x P(CYP3A4 inhibition) "
            f"(head {float(cyp34):.2f}, or 0.5 neutral prior when absent)"
        )

    if spec.si_segments is None:
        kw["si_segments"] = _SI_SEGMENTS_FULL
        basis["si_segments"] = (
            f"ACAT-lite multi-segment SI = {_SI_SEGMENTS_FULL} equal sub-compartments "
            "(full-resolution default, doc/12 L10)"
        )

    if spec.feedback_loop == 0:
        kw["feedback_loop"] = 1
        basis["feedback_loop"] = "organ-feedback loop = 1 coupled re-run (doc/05 4.5)"

    routes = {ev.route for ev in spec.dose_plan.events}
    if spec.sc_im_ka_per_h is None and routes & {Route.SUBCUTANEOUS, Route.INTRAMUSCULAR}:
        kw["sc_im_ka_per_h"] = _DEPOT_KA_PRIOR_1H
        basis["sc_im_ka_per_h"] = f"SC/IM depot absorption ka = {_DEPOT_KA_PRIOR_1H}/h class prior"
    if spec.skin_layers is None and Route.TRANSDERMAL in routes:
        kw["skin_layers"] = SkinLayers()
        basis["skin_layers"] = (
            "transdermal route: default 4-layer skin membrane (SkinLayers defaults)"
        )

    if spec.dili_immune_ic50_nm is None:
        if spec.dili_ic50_nm is not None:
            immune_ic50 = spec.dili_ic50_nm * _DILI_IMMUNE_IC50_FRACTION
            immune_basis = (
                f"immune-DILI IC50 = DILI-hazard IC50 x {_DILI_IMMUNE_IC50_FRACTION} "
                f"({spec.dili_ic50_nm:g} nM x {_DILI_IMMUNE_IC50_FRACTION}) immune fraction"
            )
        else:
            immune_ic50 = _NULL_IMMUNE_IC50_NM
            immune_basis = (
                "immune-DILI null-effect anchor: no immune-mediated DILI trigger "
                f"predicted (IC50 {_NULL_IMMUNE_IC50_NM:g} nM)"
            )
        kw["dili_immune_ic50_nm"] = immune_ic50
        if spec.dili_immune_weight is None:
            kw["dili_immune_weight"] = 0.5
        basis["dili_immune_ic50_nm"] = immune_basis

    if spec.beta_block_ic50_nm is None:
        kw["beta_block_ic50_nm"] = _NULL_SYMPATHETIC_IC50_NM
        basis["beta_block_ic50_nm"] = (
            "sympathetic branch null-effect anchor: beta-blockade predicted absent "
            f"(IC50 {_NULL_SYMPATHETIC_IC50_NM:g} nM, orders above any therapeutic exposure)"
        )

    if not spec.target_binding:
        primary = min(spec.panel, key=lambda t: t.kd_nm)
        tissue = "heart" if "hERG" in primary.name else "liver"
        kw["target_binding"] = (TargetBinding(tissue=tissue, target=primary),)
        basis["target_binding"] = (
            f"native TMDD at the primary-affinity site {primary.name!r} "
            f"(kd {primary.kd_nm:g} nM, tissue {tissue}, abundance r0 = {primary.r0_nm:g} nM): "
            "the mass-balance sink scales with the site's affinity and abundance — "
            "trivially small unless a high-affinity, high-abundance site is engaged"
        )

    return replace(spec, **kw)


def spec_from_benchmark_data(
    bench: _BenchmarkLike,
    dose_override_mg: float | None = None,
    measurements: Measurements | None = None,
    empirical: EmpiricalObservations | None = None,
) -> RunSpec:
    """Build a ``RunSpec`` from a validation ``Benchmark``-like instance."""
    cl_hep, cl_renal = _clearance_from_published(bench)
    mol = parse_structure(bench.smiles, name=bench.name)
    mw = float(mol.mw or 0.0)
    if mw <= 0:
        raise ValueError(f"could not compute molecular weight for {bench.name}")
    qt_ic50, dili_ic50 = _other_benchmark_affinities(bench.name)
    liver_params = LiverParams() if bench.name == "acetaminophen" else None
    fa = None
    if bench.route.startswith("oral") and "f_abs" in bench.published:
        lo, hi = bench.published["f_abs"]
        fa = float(getattr(bench, "fa_override", None) or (lo + hi) / 2.0)
    return _engage_full_fidelity(
        RunSpec(
            name=bench.name,
            molecule=mol,
            profile=HumanProfile(sex=Sex.MALE),
            dose_plan=_benchmark_plan(bench, dose_override_mg),
            cl_hep_l_h=cl_hep,
            cl_renal_l_h=cl_renal,
            mw=mw,
            fup=bench.fup,
            bp=bench.bp,
            qt_ic50_nm=qt_ic50,
            dili_ic50_nm=dili_ic50,
            liver_params=liver_params,
            tmax_h=bench.tmax_h,
            n_eval=bench.n_eval,
            fa=fa,
            measurements=measurements,
            empirical=empirical,
        )
    )


def _benchmark_plan(bench: _BenchmarkLike, dose_override_mg: float | None = None) -> DosePlan:
    amount = bench.dose_mg if dose_override_mg is None else dose_override_mg
    return build_dose_plan(route=bench.route, amount_mg=float(amount))


def spec_from_admet(
    molecule: Molecule,
    admet: AdmetOutput,
    profile: HumanProfile,
    dose_plan: DosePlan,
    *,
    cl_hep_l_h: float | None = None,
    cl_renal_l_h: float | None = None,
    include_pathway: bool = True,
    sc_im_ka_per_h: float | None = None,
    cl_sec_l_h: float = 0.0,
    cl_bil_l_h: float = 0.0,
    bile_emptying_1h: float = 1.0,
    hepatic_vmax_mg_h: float | None = None,
    hepatic_km_mg_l: float | None = None,
    cyp_terms: tuple[CypTerm, ...] = (),
    target_binding: tuple[TargetBinding, ...] = (),
    gut_extraction_eg: float = 0.0,
    skin_layers: SkinLayers | None = None,
    beta_block_ic50_nm: float | None = None,
    dili_immune_ic50_nm: float | None = None,
    dili_immune_weight: float | None = None,
    measurements: Measurements | None = None,
    empirical: EmpiricalObservations | None = None,
) -> RunSpec:
    """Build a ``RunSpec`` for a novel molecule entirely from ADMET-AI.

    This is the shared implementation of the ADMET -> PK auto-wiring (doc/05
    1.2/4.1) so new chemistry enters the full chain without manual physchem:
    ``fup`` comes from the predicted plasma protein binding; hepatic intrinsic
    clearance (mL/min/kg) is scaled to whole-body ``cl_hep`` (L/h) by body
    weight; renal clearance defaults to filtration of the free fraction
    (``cl_renal_l_h`` overrides); the oral small-intestine absorption gate uses
    the predicted HIA as ``fa``.  Elevates when the prediction lacks the
    needed physchem signal for the chain to be meaningful.
    """
    smiles = getattr(molecule, "canonical_smiles", None) or getattr(molecule, "smiles", None) or ""
    name = getattr(molecule, "name", None) or "custom"
    parsed = parse_structure(smiles, name=name)
    mw = float(parsed.mw or getattr(molecule, "mw", None) or 0.0)
    if mw <= 0:
        raise ValueError(f"could not compute molecular weight for {name}")
    fup = admet.fup_plasma
    if fup is None or not (0.0 < fup <= 1.0):
        raise ValueError("ADMET-AI did not return a usable fup")
    physiology = build_human(profile)
    if cl_hep_l_h is None:
        cl_int = admet.cl_int_hep_ml_min_kg
        if cl_int is None:
            raise ValueError("ADMET-AI did not return hepatic intrinsic clearance")
        cl_hep_l_h = cl_int * 60.0 / 1000.0 * profile.weight_kg  # mL/min/kg -> L/h
    if cl_renal_l_h is None:
        cl_renal_l_h = glom_filtration_clearance(physiology.gfr_l_min * 1000.0, fup)
    has_oral = any(ev.route is Route.ORAL for ev in dose_plan.events)
    fa = _admet_fa(admet) if has_oral else None
    return _engage_full_fidelity(
        RunSpec(
            name=name,
            molecule=parsed,
            profile=profile,
            dose_plan=dose_plan,
            cl_hep_l_h=cl_hep_l_h,
            cl_renal_l_h=cl_renal_l_h,
            mw=mw,
            fup=fup,
            bp=1.0,
            admet=admet,
            fa=fa,
            include_pathway=include_pathway,
            sc_im_ka_per_h=sc_im_ka_per_h,
            cl_sec_l_h=cl_sec_l_h,
            cl_bil_l_h=cl_bil_l_h,
            bile_emptying_1h=bile_emptying_1h,
            hepatic_vmax_mg_h=hepatic_vmax_mg_h,
            hepatic_km_mg_l=hepatic_km_mg_l,
            cyp_terms=cyp_terms,
            target_binding=target_binding,
            gut_extraction_eg=gut_extraction_eg,
            skin_layers=skin_layers,
            beta_block_ic50_nm=beta_block_ic50_nm,
            dili_immune_ic50_nm=dili_immune_ic50_nm,
            dili_immune_weight=dili_immune_weight,
            measurements=measurements,
            empirical=empirical,
        )
    )


__all__ = [
    "EmpiricalObservations",
    "ExposureProfile",
    "Measurements",
    "OrganStage",
    "RealismError",
    "RunResult",
    "RunSpec",
    "benchmark_data",
    "benchmark_names",
    "fidelity_provenance",
    "run_pipeline",
    "spec_from_benchmark_data",
    "summarize_clinical_cardiac",
    "summarize_clinical_cns",
    "summarize_clinical_kidney",
    "summarize_clinical_liver",
]
