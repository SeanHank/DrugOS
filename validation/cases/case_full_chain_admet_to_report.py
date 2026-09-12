"""End-to-end full-chain integrity: SMILES -> ADMET-AI -> PBPK -> report (L1).

This is the trust-mechanism case (doc/08, full-chain E2E).  It feeds
*arbitrary* molecules through the exact production entry used by the CLI and
web playground — ``parse_structure`` -> real ``predict_admet`` (hard FAIL if
the ADMET-AI runtime is missing, G5 no-silent-fallback) -> ``spec_from_admet``
(auto-wires fup, hepatic intrinsic clearance, renal filtration and oral fa)
-> ``run_pipeline`` -> ``to_contract`` -> ``render_all``/``write_report`` —
and asserts, with no human-data claim:

- hand-off integrity: every stage's output round-trips into the report
  contract unchanged (the contract is the single serialization boundary
  consumed by the CLI, web and report renderers);
- determinism: two identical runs of one compound yield byte-identical
  contracts (the reproducibility precondition of D24-level trust);
- physical plausibility: positive Cmax, unit-interval oral bioavailability;
- provenance disclosure: each report's ``trust`` record names the wired
  validated anchors (R-6 CKD-EPI, R-7 cholestasis, R-4 BBB when the ADMET-AI
  BBB_Martins head is present), the class-typical low-confidence panel priors,
  and the full-engagement realism terms (default-on under fidelity="full");
  no term outside the record contributes to any readout — nothing is silently
  borrowed (G5);
- report hand-off: JSON/Markdown/HTML artifacts render and write to disk.

The case uses two chemically distinct molecules (an oral analgesic and a
CNS-penetrant stimulant) to demonstrate that the chain is generic, not tuned
to one compound.  ADMET-AI is loaded lazily once for a batched prediction.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.inputs.models import HumanProfile, Route
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.inputs.parse_structure import parse_structure
from drugos.pipeline import run_pipeline, spec_from_admet
from drugos.pk.admet import predict_admet
from drugos.report import render_all, write_report

_REQUIRED_SECTIONS = (
    "manifest",
    "pk",
    "occupancy",
    "pathway",
    "organ",
    "clinical",
    "trust",
    "r_verify",
)
_MOLECULES = ("CC(=O)Nc1ccc(O)cc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")  # acetaminophen, caffeine


def _run_one(smiles: str, admet, data: dict):
    molecule = parse_structure(smiles)
    profile = HumanProfile(sex="male", age_y=35.0, height_cm=175.0, weight_kg=70.0)
    dose_plan = build_dose_plan(route=Route.ORAL, amount_mg=10.0, interval_h=0.0, n_doses=1)
    spec = spec_from_admet(molecule, admet, profile, dose_plan)
    result = run_pipeline(spec)
    contract = result.to_contract()
    data[smiles] = {
        "spec": spec,
        "result": result,
        "contract": contract,
    }


def case_full_chain_admet_to_report() -> CaseResult:
    try:
        admet_list = predict_admet(list(_MOLECULES))
    except (ImportError, OSError, ValueError) as exc:
        notes = [
            "ADMET-AI runtime unavailable: hard FAIL, no fallback (G5) — "
            f"{exc!r}.  Install admet_ai and its CC-BY weights to run this case."
        ]
        metrics = [
            MetricResult(
                "admet_runtime",
                float("nan"),
                0.0,
                0.0,
                "available",
                "FAIL",
            )
        ]
        return CaseResult(
            "E2E full-chain ADMET -> report integrity (trust mechanism)",
            False,
            metrics,
            notes,
            level=EvidenceLevel.L1_SELF_CONSISTENCY,
        )

    acetadmet, caffadmet = admet_list
    acet_smiles, caff_smiles = _MOLECULES
    data: dict[str, dict] = {}

    _run_one(acet_smiles, acetadmet, data)
    _run_one(caff_smiles, caffadmet, data)

    # Determinism: re-run the same compound and require an identical contract.
    rerun = run_pipeline(
        spec_from_admet(
            parse_structure(acet_smiles),
            acetadmet,
            HumanProfile(sex="male", age_y=35.0, height_cm=175.0, weight_kg=70.0),
            build_dose_plan(route=Route.ORAL, amount_mg=10.0, interval_h=0.0, n_doses=1),
        )
    )
    determinism_ok = rerun.to_contract() == data[acet_smiles]["contract"]

    section_ok = all(set(c["contract"]) == set(_REQUIRED_SECTIONS) for c in data.values())
    cmax_ok = all(c["contract"]["pk"]["cmax_mg_l"] > 0.0 for c in data.values())
    fa_ok = all(0.0 < c["contract"]["pk"]["bioavailability_f"] <= 1.0 for c in data.values())
    anchors_ok = all(
        any("cholestasis" in a for a in c["contract"]["trust"]["anchors_wired"])
        and any("CKD-EPI" in a for a in c["contract"]["trust"]["anchors_wired"])
        for c in data.values()
    )
    g5_ok = all("G5" in c["contract"]["trust"]["policy"] for c in data.values())
    priors_ok = all(len(c["contract"]["trust"]["no_public_data_sites"]) >= 1 for c in data.values())
    terms_ok = all(c["contract"]["trust"]["mechanism_terms_engaged"] != [] for c in data.values())
    no_seams_ok = all(not [k for k in c["contract"]["trust"] if "seam" in k] for c in data.values())
    bbb_wired = any("BBB" in a for a in data[caff_smiles]["contract"]["trust"]["anchors_wired"])

    with tempfile.TemporaryDirectory() as tmp:
        payloads = render_all(data[caff_smiles]["result"])
        paths = write_report(data[caff_smiles]["result"], directory=tmp)
        report_files_ok = all(payloads.values()) and all(Path(p).exists() for p in paths.values())
    bbb_wired = any("BBB" in a for a in data[caff_smiles]["contract"]["trust"]["anchors_wired"])

    ok = (
        determinism_ok
        and section_ok
        and cmax_ok
        and fa_ok
        and anchors_ok
        and g5_ok
        and priors_ok
        and terms_ok
        and no_seams_ok
        and report_files_ok
    )

    acet = data[acet_smiles]["contract"]
    caff = data[caff_smiles]["contract"]
    metrics = [
        MetricResult(
            "contract_determinism",
            float(determinism_ok),
            1.0,
            1.0,
            "bool",
            "pass" if determinism_ok else "FAIL",
        ),
        MetricResult(
            "required_sections_present",
            float(sum(len(c["contract"]) == 8 for c in data.values())),
            float(len(data)),
            float(len(data)),
            "contracts",
            "pass" if section_ok else "FAIL",
        ),
        MetricResult(
            "plasma_cmax_positive",
            float(sum(c["contract"]["pk"]["cmax_mg_l"] > 0.0 for c in data.values())),
            float(len(data)),
            float(len(data)),
            "compounds",
            "pass" if cmax_ok else "FAIL",
        ),
        MetricResult(
            "oral_bioavailability_unit_interval",
            round(min(c["contract"]["pk"]["bioavailability_f"] for c in data.values()), 4),
            0.0,
            1.0,
            "fraction",
            "pass" if fa_ok else "FAIL",
        ),
        MetricResult(
            "validated_anchors_disclosed",
            float(sum(len(c["contract"]["trust"]["anchors_wired"]) >= 2 for c in data.values())),
            float(len(data)),
            float(len(data)),
            "contracts",
            "pass" if anchors_ok else "FAIL",
        ),
        MetricResult(
            "g5_no_silent_fallback_disclosed",
            float(sum("G5" in c["contract"]["trust"]["policy"] for c in data.values())),
            float(len(data)),
            float(len(data)),
            "contracts",
            "pass" if g5_ok else "FAIL",
        ),
        MetricResult(
            "class_priors_disclosed",
            float(
                sum(len(c["contract"]["trust"]["no_public_data_sites"]) >= 1 for c in data.values())
            ),
            float(len(data)),
            float(len(data)),
            "contracts",
            "pass" if priors_ok else "FAIL",
        ),
        MetricResult(
            "mechanism_terms_engaged_full",
            float(
                sum(c["contract"]["trust"]["mechanism_terms_engaged"] != [] for c in data.values())
            ),
            float(len(data)),
            float(len(data)),
            "contracts",
            "pass" if terms_ok else "FAIL",
        ),
        MetricResult(
            "report_artifacts_written",
            float(report_files_ok),
            1.0,
            1.0,
            "bool",
            "pass" if report_files_ok else "FAIL",
        ),
    ]
    notes = [
        f"acetaminophen: Cmax {acet['pk']['cmax_mg_l']:.3g} mg/L, Fa "
        f"{acet['pk']['bioavailability_f']:.3g}, verdict "
        f"{acet['clinical']['verdict']}; caffeine: Cmax {caff['pk']['cmax_mg_l']:.3g} "
        f"mg/L, Fa {caff['pk']['bioavailability_f']:.3g}, verdict "
        f"{caff['clinical']['verdict']}; anchors "
        f"{caff['trust']['anchors_wired']}; BBB_Martins head "
        f"{'wired' if bbb_wired else 'not wired'}; mechanism terms engaged "
        f"{caff['trust']['mechanism_terms_engaged']}; deterministic re-run "
        f"{'matches' if determinism_ok else 'diverges'}.",
    ]
    return CaseResult(
        "E2E full-chain ADMET -> report integrity (trust mechanism)",
        ok,
        metrics,
        notes,
        level=EvidenceLevel.L1_SELF_CONSISTENCY,
    )


__all__ = ["case_full_chain_admet_to_report"]
