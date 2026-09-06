"""Cardiac hERG/QT panel (L3 — empirically anchored, doc/08 §1.1 Tier 1).

Stage-1 plasma free exposure drives Stage-2 hERG occupancy (Kd 2 nM, the
dofetilide-class anchor) which is converted through the Emax hERG->IKr->QTc
relation (doc/05 4.3) into a Delta-QTc trajectory.  Predicted peak Delta-QTc
for a single 0.5 mg dofetilide dose is compared against the published human
QT-prolongation range; a negligible-hERG control stays flat.
"""

from __future__ import annotations

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, _simulate

from drugos.inputs.models import HumanProfile, Sex
from drugos.organ.cardiac import simulate_cardiac
from drugos.pk.physiology import build_human
from drugos.target.occupancy import simulate_occupancy
from drugos.target.targets import Target


def _bench(name: str) -> object:
    from validation.benchmarks import BENCHMARKS

    return next(b for b in BENCHMARKS if b.name == name)


def _delta_qtc_max_ms(res, mw: float, site: Target) -> float:
    occ = simulate_occupancy(res.t, res.plasma_free, site, mw=mw, n_eval=500)
    human = build_human(HumanProfile(sex=Sex.MALE, body_weight_kg=70, height_cm=175))
    out = simulate_cardiac(occ.t_h, occ.occupancy, human, n_eval=500)
    return float(out.delta_qtc_ms.max())


def case_cardiac_qtc() -> CaseResult:
    herg = Target(
        name="hERG (Kv11.1) dofetilide anchor",
        kd_nm=2.0,
        reference="low-nM hERG blocker class (e.g. dofetilide)",
    )
    weak = Target(name="hERG control", kd_nm=2.0e5, reference="negligible-hERG class")

    delta_dofetilide = _delta_qtc_max_ms(_simulate(_bench("dofetilide")), 441.537, herg)
    delta_warfarin = _delta_qtc_max_ms(_simulate(_bench("warfarin")), 308.33, weak)

    met_dof = MetricResult(
        "dofetilide_delta_qtc_ms",
        delta_dofetilide,
        5.0,
        55.0,
        "ms",
        "pass" if 5.0 <= delta_dofetilide <= 55.0 else "FAIL",
    )
    met_control = MetricResult(
        "warfarin_delta_qtc_ms_control",
        delta_warfarin,
        0.0,
        1.0,
        "ms",
        "pass" if delta_warfarin < 1.0 else "FAIL",
    )
    ok = met_dof.criterion == "pass" and met_control.criterion == "pass"
    return CaseResult(
        "cardiac QTc prolongation (dofetilide hERG)",
        ok,
        [met_dof, met_control],
        [
            f"0.5 mg dofetilide peak Delta-QTc {delta_dofetilide:.1f} ms vs "
            f"published ~20-60 ms/QTc-prolonging clinical band (low-nM hERG "
            f"block); warfarin control {delta_warfarin:.2f} ms.  "
            "Published: dofetilide (Tikosyn) USPI lists QT/QTc prolongation; "
            "peak Delta-QTc in the 0.5 mg single-dose range is around "
            "10-35 ms and TdP aggregates in QTc > 500 ms."
        ],
        level=EvidenceLevel.L3_EMPIRICAL,
    )


__all__ = ["case_cardiac_qtc"]
