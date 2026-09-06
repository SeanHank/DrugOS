"""End-to-end composite-risk ordering (L3 — empirically anchored, doc/08 §1.1).

Feeds real benchmark compounds through the *whole* pipeline contract
(``spec_from_benchmark_data`` -> ``run_pipeline`` -> ``to_contract``) and
asserts the Stage-5 composite toxicity anchors:

- dofetilide is a QT/TdP-dominant high-moderate risk (hERG IC50 ~2 nM) while
  warfarin stays low on QT; dofetilide overall risk must exceed warfarin
  overall.  The floor is calibrated on the **myocardium-exposure driver**
  (pipeline feeds cardiac hERG blockade from the PBPK heart free tissue
  exposure): a single 0.5 mg dose yields Cmax_unbound/hERG-IC50 ≈ 0.87,
  ΔQTc ≈ 19 ms and a fused QT risk ≈ 0.53 — genuinely elevated (exposure-ratio
  line 0.82) but below the legacy 0.70 floor that was calibrated on the
  over-estimating hepatic-tissue driver (doc/05 §4.3).
- acetaminophen 20 g overdose is DILI-dominant high risk; a 1 g dose stays
  low on DILI; the dose-escalation ordering must be preserved.
- unanchored CNS endpoints fall back to the 0.20 class prior exactly.

Numbers reproduce the calibrated anchors on the fixed-sim seed path.
"""

from __future__ import annotations

from dataclasses import replace

from validation.benchmarks import BENCHMARKS
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.clinical.toxicity import Endpoint
from drugos.inputs.models import HumanProfile
from drugos.pipeline import run_pipeline, spec_from_benchmark_data


def _bench(name: str):
    return next(b for b in BENCHMARKS if b.name == name)


def _endpoint_risk(result, endpoint: Endpoint) -> float:
    row = result.toxicity.by_endpoint(endpoint)
    return row.risk if row is not None else 0.0


def case_risk_ordering() -> CaseResult:
    dof = spec_from_benchmark_data(_bench("dofetilide"))
    war = spec_from_benchmark_data(_bench("warfarin"))
    apap_1g = spec_from_benchmark_data(_bench("acetaminophen"))
    apap_20g = spec_from_benchmark_data(_bench("acetaminophen"), dose_override_mg=20000.0)
    apap_f = replace(
        apap_20g,
        profile=HumanProfile(sex="female", age_y=60.0, height_cm=160.0, weight_kg=65.0),
    )

    r_dof = run_pipeline(dof)
    r_war = run_pipeline(war)
    r_ap1 = run_pipeline(apap_1g)
    r_ap20 = run_pipeline(apap_20g)
    r_apf = run_pipeline(apap_f)

    dof_qt = _endpoint_risk(r_dof, Endpoint.QT)
    war_qt = _endpoint_risk(r_war, Endpoint.QT)
    ap20_dili = _endpoint_risk(r_ap20, Endpoint.DILI)
    ap1_dili = _endpoint_risk(r_ap1, Endpoint.DILI)
    war_cns = _endpoint_risk(r_war, Endpoint.CNS)

    ordering_qt = dof_qt > war_qt
    ordering_dili = ap20_dili > ap1_dili
    driver_qt = r_dof.toxicity.overall_driver() == "qt"
    driver_dili = r_ap20.toxicity.overall_driver() == "dili"
    prior_cns = abs(war_cns - 0.20) < 1e-9
    held_out_stable = ("High composite risk" in r_ap20.verdict and "dili" in r_ap20.verdict) and (
        "High composite risk" in r_apf.verdict and "dili" in r_apf.verdict
    )

    inband = (
        0.45 <= dof_qt <= 0.95
        and 0.0 <= war_qt <= 0.35
        and 0.60 <= ap20_dili <= 1.0
        and 0.0 <= ap1_dili <= 0.40
    )
    ok = (
        inband
        and ordering_qt
        and ordering_dili
        and driver_qt
        and driver_dili
        and prior_cns
        and held_out_stable
    )

    metrics = [
        MetricResult(
            "dofetilide_qt_risk",
            dof_qt,
            0.45,
            0.95,
            "P(risk)",
            "pass" if 0.45 <= dof_qt <= 0.95 else "FAIL",
        ),
        MetricResult(
            "warfarin_qt_risk",
            war_qt,
            0.0,
            0.35,
            "P(risk)",
            "pass" if 0.0 <= war_qt <= 0.35 else "FAIL",
        ),
        MetricResult(
            "acetaminophen_20g_dili_risk",
            ap20_dili,
            0.60,
            1.0,
            "P(risk)",
            "pass" if 0.60 <= ap20_dili <= 1.0 else "FAIL",
        ),
        MetricResult(
            "acetaminophen_1g_dili_risk",
            ap1_dili,
            0.0,
            0.40,
            "P(risk)",
            "pass" if 0.0 <= ap1_dili <= 0.40 else "FAIL",
        ),
        MetricResult(
            "unanchored_cns_prior", war_cns, 0.20, 0.20, "P(risk)", "pass" if prior_cns else "FAIL"
        ),
        MetricResult(
            "held_out_dose_profile_stable",
            0.0 if held_out_stable else 1.0,
            0.0,
            0.0,
            "disagreements",
            "pass" if held_out_stable else "FAIL",
        ),
    ]
    return CaseResult(
        "Stage-5 composite risk ordering vs clinical anchors",
        ok,
        metrics,
        [
            f"dofetilide QT {dof_qt:.3f} ({r_dof.toxicity.overall_driver()}-driven) > "
            f"warfarin QT {war_qt:.3f}; APAP 20 g DILI {ap20_dili:.3f} > 1 g "
            f"DILI {ap1_dili:.3f} ({r_ap20.toxicity.overall_driver()}-driven); "
            f"unanchored CNS sits on the 0.20 class prior.  Published anchors: "
            "dofetilide (Tikosyn) is a QT-prolonging hERG blocker and is "
            "contraindicated with renal/QT risk; massive acetaminophen overdose "
            "causes centrilobular hepatic necrosis (DILI), while warfarin is not "
            "a QT liability."
        ],
        level=EvidenceLevel.L3_EMPIRICAL,
    )


__all__ = ["case_risk_ordering"]
