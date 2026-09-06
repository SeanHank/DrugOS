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

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from drugos.clinical.biomarkers import BIOMARKERS, BiomarkerGrade, grade_absolute, grade_timeseries
from drugos.clinical.toxicity import Endpoint, ToxicityReport, score_toxicity
from drugos.inputs.models import DosePlan, HumanProfile, Molecule, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.parse_structure import parse_structure
from drugos.inputs.resolve_human import resolve_human
from drugos.organ.base import NDArray
from drugos.organ.cardiac import CardiacResult, simulate_cardiac
from drugos.organ.cns import CnsParams, CnsResult, simulate_cns
from drugos.organ.feedback import OrganFeedback, apply_pk_scaling, feedback_from_results
from drugos.organ.kidney import KidneyResult, simulate_kidney
from drugos.organ.liver import LiverParams, LiverTrajectory, liver_params_from_panel, simulate_liver
from drugos.pathway.simulator import PathwayResult, mapk_cascade, simulate_pathway
from drugos.pk.admet import AdmetOutput
from drugos.pk.partitions import partition_from_molecule
from drugos.pk.pbpk_build import AbsorptionParams, PBPKModel
from drugos.pk.physiology import HumanPhysiology
from drugos.pk.simulate import PBPKResult, PkMetrics, compute_pk_metrics, simulate_pbpk
from drugos.target.occupancy import (
    PanelEngagement,
    TargetOccupancyResult,
    simulate_occupancy,
    simulate_panel,
)
from drugos.target.targets import Target, safety_panel

_NM_PER_MG = 1e6
_DEFAULT_SCR_BASE_UMOL_L = 80.0
_ALT_ULN_U_L = 40.0
_BILI_ULN_MG_DL = 1.0
_CNS_IC50_DEFAULT_NM = 1.0e5


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
    cns_ic50_nm: float | None = None
    tmax_h: float = 48.0
    n_eval: int = 601
    include_pathway: bool = True
    feedback_loop: int = 0
    sc_im_ka_per_h: float | None = None
    name: str = "compound"

    def __post_init__(self) -> None:
        if self.mw <= 0:
            raise ValueError("mw must be positive")
        if not (0.0 < self.fup <= 1.0):
            raise ValueError("fup must be in (0, 1]")
        if self.feedback_loop < 0:
            raise ValueError("feedback_loop must be non-negative")
        if self.sc_im_ka_per_h is not None and self.sc_im_ka_per_h <= 0:
            raise ValueError("sc_im_ka_per_h must be positive")
        if self.qt_ic50_nm is not None and self.qt_ic50_nm <= 0:
            raise ValueError("qt_ic50_nm must be positive")
        if self.dili_ic50_nm is not None and self.dili_ic50_nm <= 0:
            raise ValueError("dili_ic50_nm must be positive")
        if self.cns_ic50_nm is not None and self.cns_ic50_nm <= 0:
            raise ValueError("cns_ic50_nm must be positive")


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
    cns_anchored: bool = False


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
                },
                "cns": {
                    "peak_brain_free_nm": round(float(cns.peak_brain_free_nm), 3),
                    "exposure_ratio": round(cns.exposure_ratio, 3),
                    "grade": cns.grade,
                },
                "trajectories": {
                    "t_h": liver.t_h.tolist(),
                    "liver_alt_U_L": liver.alt_u_l.tolist(),
                    "liver_bilirubin_mg_dL": liver.bilirubin_mg_dl.tolist(),
                    "qtc_ms": cardiac.qtc_ms.tolist(),
                    "delta_qtc_ms": cardiac.delta_qtc_ms.tolist(),
                    "gfr_ml_min": kidney.gfr_ml_min.tolist(),
                    "scr_ratio": kidney.scr_ratio.tolist(),
                    "brain_free_nm": cns.brain_free_nm.tolist(),
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
                    "cns_ic50_nm": self.exposure.cns_ic50_nm,
                },
            },
        }


def _mg_l_to_nm(c_mg_l: float, mw: float) -> float:
    return c_mg_l * _NM_PER_MG / mw


def _herg_occupancy(
    panel: PanelEngagement,
    t_h: NDArray,
    c_free_mg_l: NDArray,
    mw: float,
    qt_ic50_nm: float | None,
    n_eval: int = 601,
) -> NDArray:
    """hERG occupancy either from an override affinity or the panel prior."""
    if qt_ic50_nm is not None:
        site = Target(
            name="hERG (Kv11.1)",
            kd_nm=qt_ic50_nm,
            r0_nm=0.05,
            reference="specified per-compound hERG affinity override",
        )
        occ = simulate_occupancy(t_h, c_free_mg_l, site, mw=mw, n_eval=n_eval)
        return np.interp(t_h, occ.t_h, occ.occupancy)
    for name, res in panel.results.items():
        if "hERG" in name:
            return np.interp(t_h, res.t_h, res.occupancy)
    return np.zeros_like(t_h)


def _dili_ic50(spec: RunSpec) -> float | None:
    if spec.dili_ic50_nm is not None:
        return spec.dili_ic50_nm
    cands = [t.kd_nm for t in spec.panel if "BSEP" in t.name or "Mitochondrial" in t.name]
    return min(cands) if cands else None


def _qt_ic50(spec: RunSpec) -> float | None:
    if spec.qt_ic50_nm is not None:
        return spec.qt_ic50_nm
    cands = [t.kd_nm for t in spec.panel if "hERG" in t.name]
    return cands[0] if cands else None


def _panel_model(spec: RunSpec) -> PBPKModel:
    physiology = resolve_human(spec.profile)
    part = partition_from_molecule(
        spec.molecule, fup=spec.fup, bp=spec.bp, hematocrit=physiology.hematocrit
    )
    absorption = (
        AbsorptionParams(k_depot_absorption=spec.sc_im_ka_per_h)
        if spec.sc_im_ka_per_h is not None
        else AbsorptionParams()
    )
    return PBPKModel(
        physiology=physiology,
        partition=part,
        bp=spec.bp,
        fup=spec.fup,
        cl_hep_l_h=spec.cl_hep_l_h,
        cl_renal_l_h=spec.cl_renal_l_h,
        dose_plan=spec.dose_plan,
        absorption=absorption,
    )


def _organ_state(
    pk: PBPKResult,
    model: PBPKModel,
    spec: RunSpec,
) -> tuple[
    PanelEngagement,
    TargetOccupancyResult | None,
    LiverTrajectory,
    NDArray,
    CardiacResult,
    KidneyResult,
    CnsResult,
    PathwayResult | None,
]:
    """Run the downstream stages (occupancy -> organ QST -> pathway) on ``pk``."""
    liver_free = pk.unbound_tissues["liver"]
    panel = simulate_panel(pk.t, liver_free, spec.panel, mw=spec.mw, n_eval=spec.n_eval)

    primary_signal = max(
        (r for r in panel.results.values()),
        key=lambda r: r.time_at_target_h,
        default=None,
    )

    liver_params = spec.liver_params or liver_params_from_panel(spec.panel)
    liver = simulate_liver(pk.t, liver_free, spec.mw, params=liver_params, n_eval=spec.n_eval)

    block = _herg_occupancy(panel, pk.t, liver_free, spec.mw, spec.qt_ic50_nm, spec.n_eval)
    cardiac = simulate_cardiac(pk.t, block, model.physiology, n_eval=spec.n_eval)

    kidney_free = pk.unbound_tissues["kidney"]
    kidney = simulate_kidney(
        pk.t,
        kidney_free,
        spec.mw,
        gfr_base_ml_min=model.physiology.gfr_ml_min,
        scr_base_umol_l=_DEFAULT_SCR_BASE_UMOL_L,
    )

    cns = simulate_cns(
        pk.t,
        pk.plasma_free,
        spec.mw,
        CnsParams(ic50_nm=spec.cns_ic50_nm or _CNS_IC50_DEFAULT_NM),
    )

    pathway: PathwayResult | None = None
    if spec.include_pathway and primary_signal is not None:
        pathway = simulate_pathway(
            mapk_cascade(), primary_signal.t_h, primary_signal.occupancy, n_eval=spec.n_eval
        )

    return panel, primary_signal, liver, block, cardiac, kidney, cns, pathway


def _organ_feedback(
    model: PBPKModel, liver: LiverTrajectory, kidney: KidneyResult
) -> OrganFeedback:
    """PK scalings from the current organ state (doc/05 4.5, sequential outer loop)."""
    gfr_fraction = kidney.min_gfr_ml_min / model.physiology.gfr_ml_min
    return feedback_from_results(liver.max_dead_frac, gfr_fraction, 1.0)


def run_pipeline(spec: RunSpec) -> RunResult:
    """Execute the full sequential pipeline for ``spec``.

    ``spec.feedback_loop`` (>0) re-runs PK with organ-dysfunction scaling
    (doc/05 4.5): each iteration re-derives hepatic/GFR scalings from the
    current organ state and resolves the PBPK model afresh.
    """
    if spec.tmax_h <= 0:
        raise ValueError("tmax_h must be positive")

    model = _panel_model(spec)
    pk = simulate_pbpk(model, tmax_h=spec.tmax_h, n_eval=spec.n_eval)
    metrics = compute_pk_metrics(pk.t, pk.plasma_total, pk.dose_mg, pk.route)

    for _ in range(spec.feedback_loop):
        (_panel, _primary, liver, _block, _cardiac, kidney, _cns, _pathway) = _organ_state(
            pk, model, spec
        )
        model = apply_pk_scaling(model, _organ_feedback(model, liver, kidney))
        pk = simulate_pbpk(model, tmax_h=spec.tmax_h, n_eval=spec.n_eval)

    (panel, primary_signal, liver, block, cardiac, kidney, cns, pathway) = _organ_state(
        pk, model, spec
    )

    physiology = model.physiology
    liver_free = pk.unbound_tissues["liver"]
    kidney_free = pk.unbound_tissues["kidney"]

    organ = OrganStage(
        liver=liver,
        cardiac=cardiac,
        kidney=kidney,
        cns=cns,
        biomarkers=_grade_organ(liver, cardiac, kidney, cns),
    )

    exposure = ExposureProfile(
        plasma_cmax_unbound_nm=_mg_l_to_nm(float(np.max(pk.plasma_free)), spec.mw),
        liver_cmax_free_nm=_mg_l_to_nm(float(np.max(liver_free)), spec.mw),
        kidney_cmax_free_nm=_mg_l_to_nm(float(np.max(kidney_free)), spec.mw),
        qt_ic50_nm=_qt_ic50(spec),
        dili_ic50_nm=_dili_ic50(spec),
        cns_ic50_nm=spec.cns_ic50_nm or _CNS_IC50_DEFAULT_NM,
        cns_anchored=spec.cns_ic50_nm is not None,
    )

    toxicity = _score_clinical(organ, exposure, spec.admet)
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
    )


def _grade_organ(
    liver: LiverTrajectory,
    cardiac: CardiacResult,
    kidney: KidneyResult,
    cns: CnsResult,
) -> dict[str, list[BiomarkerGrade]]:
    return {
        "liver": summarize_clinical_liver(liver),
        "cardiac": summarize_clinical_cardiac(cardiac),
        "kidney": summarize_clinical_kidney(kidney),
        "cns": summarize_clinical_cns(cns),
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


def summarize_clinical_kidney(kidney: KidneyResult) -> list[BiomarkerGrade]:
    """Grade the kidney trajectory (GFR + creatinine ratio)."""
    return [
        grade_timeseries(kidney.t_h, kidney.gfr_ml_min, BIOMARKERS["gfr"]),
        grade_timeseries(kidney.t_h, kidney.scr_ratio, BIOMARKERS["creatinine"]),
    ]


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


def spec_from_benchmark_data(
    bench: _BenchmarkLike, dose_override_mg: float | None = None
) -> RunSpec:
    """Build a ``RunSpec`` from a validation ``Benchmark``-like instance."""
    cl_hep, cl_renal = _clearance_from_published(bench)
    mol = parse_structure(bench.smiles, name=bench.name)
    mw = float(mol.mw or 0.0)
    if mw <= 0:
        raise ValueError(f"could not compute molecular weight for {bench.name}")
    qt_ic50, dili_ic50 = _other_benchmark_affinities(bench.name)
    liver_params = LiverParams() if bench.name == "acetaminophen" else None
    return RunSpec(
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
    )


def _benchmark_plan(bench: _BenchmarkLike, dose_override_mg: float | None = None) -> DosePlan:
    amount = bench.dose_mg if dose_override_mg is None else dose_override_mg
    return build_dose_plan(route=bench.route, amount_mg=float(amount))


__all__ = [
    "ExposureProfile",
    "OrganStage",
    "RunResult",
    "RunSpec",
    "benchmark_data",
    "benchmark_names",
    "run_pipeline",
    "spec_from_benchmark_data",
    "summarize_clinical_cardiac",
    "summarize_clinical_cns",
    "summarize_clinical_kidney",
    "summarize_clinical_liver",
]
