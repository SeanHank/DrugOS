"""Evidence ladder, result types, and shared case machinery (doc/08).

The three-level ladder grades how much *epistemic weight* each result
carries — from "reproduces clinically observed human data" down to "the
solver and bookkeeping are numerically self-consistent".  It mirrors the
doc/08 three-tier validation ladder (benchmark compounds / stage-level /
prospective + CI checks).

Levels, strongest to weakest:

- ``L3_EMPIRICAL``      — predicted vs published human clinical ranges.
- ``L2_ANALYTIC_LIMIT`` — closed-form / mechanistically forced point-match.
- ``L1_SELF_CONSISTENCY`` — ODE conservation and linearity axioms only.

Case modules call ``_build_model``, ``_simulate`` and friends from here so
that every case exercises the same single code path as the benchmark suite.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from types import SimpleNamespace

from validation.benchmarks import FOLD_ALLOWANCE, Benchmark

from drugos.inputs.models import DosePlan, HumanProfile, Molecule, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.resolve_human import resolve_human
from drugos.pk.partitions import partition_from_molecule
from drugos.pk.pbpk_build import TISSUE_LIST, AbsorptionParams, PBPKModel
from drugos.pk.simulate import PBPKResult, simulate_pbpk


class EvidenceLevel(StrEnum):
    """Epistemic weight of a validation result (doc/08 §1.1-1.4)."""

    L1_SELF_CONSISTENCY = "L1"
    L2_ANALYTIC_LIMIT = "L2"
    L3_EMPIRICAL = "L3"


@dataclass(frozen=True)
class EvidenceMeta:
    """Human-readable definition of one evidence level."""

    label: str
    basis: str
    certifies: str
    doc08: str


#: Ordered strongest -> weakest; drives report tables.
EVIDENCE_ORDER: tuple[EvidenceLevel, ...] = (
    EvidenceLevel.L3_EMPIRICAL,
    EvidenceLevel.L2_ANALYTIC_LIMIT,
    EvidenceLevel.L1_SELF_CONSISTENCY,
)

EVIDENCE_META: dict[EvidenceLevel, EvidenceMeta] = {
    EvidenceLevel.L3_EMPIRICAL: EvidenceMeta(
        label="Empirically anchored (Tier 1)",
        basis="Predicted vs published human clinical ranges "
        "(USPI / literature) under the 2x GMFE allowance.",
        certifies="The pipeline reproduces clinically observed human PK "
        "within the fold allowance — the strongest evidence in this suite.",
        doc08="doc/08 §1.1",
    ),
    EvidenceLevel.L2_ANALYTIC_LIMIT: EvidenceMeta(
        label="Analytic / mechanistic limit (Tier 2)",
        basis="Closed-form or mechanistically forced point-match derived "
        "from the stage's own equations.",
        certifies="The stage ODEs solve their intended dynamics correctly; "
        "it does not by itself certify human predictivity (needs L3).",
        doc08="doc/08 §1.2",
    ),
    EvidenceLevel.L1_SELF_CONSISTENCY: EvidenceMeta(
        label="Internal consistency / CI (Tier 3)",
        basis="ODE / compiler self-consistency: mass conservation and "
        "linearity axioms, no external data.",
        certifies="Numerical correctness of the solver and mass "
        "bookkeeping; weakest in epistemic weight — 'just computes right'.",
        doc08="doc/08 §1.3-1.4",
    ),
}

_PROFILE: HumanProfile | None = None


@dataclass(frozen=True)
class MetricResult:
    """A single quantitative validation assertion."""

    name: str
    predicted: float
    lo: float
    hi: float
    unit: str
    criterion: str


@dataclass
class CaseResult:
    benchmark: str
    passed: bool
    metrics: list[MetricResult] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    level: EvidenceLevel = EvidenceLevel.L1_SELF_CONSISTENCY


def profile() -> HumanProfile:
    global _PROFILE
    if _PROFILE is None:
        _PROFILE = resolve_human(HumanProfile(sex=Sex.MALE))
    return _PROFILE


def _clearance_from_published(b: Benchmark) -> tuple[float, float]:
    """Unbound hepatic and renal clearances from the published CL band."""
    lo, hi = b.published["cl_plasma_l_h"]
    cl_total = (lo + hi) / 2.0
    renal_frac = b.published.get("urine_fraction", (0.0, 0.0))
    renal_share = (renal_frac[0] + renal_frac[1]) / 2.0 if renal_frac[1] > 0 else 0.0
    cl_hep = cl_total * (1.0 - renal_share) / b.fup
    cl_renal = cl_total * renal_share / b.fup
    return cl_hep, cl_renal


def _build_model(b: Benchmark) -> PBPKModel:
    mol = Molecule(
        name=b.name,
        canonical_smiles=b.smiles,
        log_p=b.log_p,
        pka_acids=b.pka_acids,
        pka_bases=b.pka_bases,
    )
    h = profile()
    part = partition_from_molecule(mol, fup=b.fup, bp=b.bp, hematocrit=h.hematocrit)
    plan = build_dose_plan(route=b.route, amount_mg=b.dose_mg)
    abs_params = AbsorptionParams()
    if b.fa_override is not None:
        abs_params = AbsorptionParams(k_si_absorption=(0.55 * b.fa_override) / 0.85)
    cl_hep, cl_renal = _clearance_from_published(b)
    return PBPKModel(
        physiology=h,
        partition=part,
        bp=b.bp,
        fup=b.fup,
        cl_hep_l_h=cl_hep,
        cl_renal_l_h=cl_renal,
        dose_plan=plan,
        absorption=abs_params,
    )


def _simulate(b: Benchmark) -> PBPKResult:
    return simulate_pbpk(_build_model(b), tmax_h=b.tmax_h, n_eval=b.n_eval)


def _fraction_absorbed(b: Benchmark, res: PBPKResult) -> float:
    if res.feces_cum_mg is None:
        return 1.0
    dose = b.dose_mg
    return 1.0 - float(res.feces_cum_mg[-1]) / max(dose, 1e-12)


def _allowed_band(lo: float, hi: float, unit: str) -> tuple[float, float]:
    center = (lo * hi) ** 0.5
    hi_ok = center * FOLD_ALLOWANCE
    if unit == "fraction":
        hi_ok = min(hi_ok, 1.0)
    return center / FOLD_ALLOWANCE, hi_ok


def _single_pool(cl: float) -> PBPKModel:
    h = profile()
    names: Sequence[str] = [*TISSUE_LIST]
    kp = {n: 1.0 for n in names}
    part = SimpleNamespace(kp=kp, kpu={n: 1.0 for n in names})
    plan = DosePlan.iv_bolus(dose_mg=10.0)
    return PBPKModel(
        physiology=h,
        partition=part,
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=cl,
        cl_renal_l_h=0.0,
        dose_plan=plan,
        absorption=AbsorptionParams(),
    )


__all__ = [
    "CaseResult",
    "EVIDENCE_META",
    "EVIDENCE_ORDER",
    "EvidenceLevel",
    "MetricResult",
    "_allowed_band",
    "_build_model",
    "_clearance_from_published",
    "_fraction_absorbed",
    "_simulate",
    "_single_pool",
    "profile",
]
