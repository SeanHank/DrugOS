"""CKD-EPI 2021 race-free GFR baseline (R-6, L2 analytic limit).

Kidney-stage production anchor (doc/12 row 4b): when a measured serum
creatinine is available on the human profile, the baseline GFR is no longer
the age/sex default but the **CKD-EPI 2021 race-free creatinine equation**
(Levey et al., N Engl J Med 2021;385:1737-49), scaled to the per-subject
absolute GFR by Mosteller BSA.  This is the standard, published, openly
reproducible clinical equation for eGFR — the production-validated baseline
for the Stage-4 nephron panel (CMR nephron SBML is a planned release; see
doc/07 P7, in-house CKD-EPI + tubular transport carries the Stage-4 baseline).

The case pins: (1) reference eGFR points under the 2021 equation, (2) the
BSA scaling, (3) pipeline wiring — a profile carrying Scr drives
``gfr_ml_min`` through CKD-EPI instead of the 110/125 mL/min default, and
(4) the untouched default when no Scr is present (benchmark invariance).
"""

from __future__ import annotations

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.inputs.models import HumanProfile, Sex
from drugos.organ.kidney import ckdepi_2021_egfr
from drugos.pk.physiology import build_human, mosteller_bsa

_MALE60_SCR1 = ckdepi_2021_egfr(1.0, 60, female=False)
_BSA175_80 = mosteller_bsa(175.0, 80.0)


def _metrics() -> list[MetricResult]:
    egfr_female = ckdepi_2021_egfr(1.0, 60, female=True)
    egfr_bsa = ckdepi_2021_egfr(1.0, 60, female=False, bsa_m2=_BSA175_80)

    prof_scr = build_human(
        HumanProfile(
            sex=Sex.MALE,
            age_y=60.0,
            height_cm=175.0,
            weight_kg=80.0,
            serum_creatinine_mg_dl=1.0,
        )
    )
    expected_scr = ckdepi_2021_egfr(1.0, 60, female=False, bsa_m2=_BSA175_80)
    prof_default = build_human(HumanProfile(sex=Sex.MALE, age_y=60.0))

    return [
        MetricResult(
            "egfr_male_60_scr1_pct2_5",
            _MALE60_SCR1,
            85.0,
            92.0,
            "mL/min/1.73 m2",
            "pass" if 85.0 <= _MALE60_SCR1 <= 92.0 else "FAIL",
        ),
        MetricResult(
            "egfr_female_lower_same_scr",
            egfr_female,
            0.0,
            _MALE60_SCR1,
            "mL/min/1.73 m2",
            "pass" if egfr_female < _MALE60_SCR1 else "FAIL",
        ),
        MetricResult(
            "egfr_bsa_scaled_absolute",
            egfr_bsa,
            _MALE60_SCR1 * 1.0,
            _MALE60_SCR1 * 1.3,
            "mL/min",
            "pass" if egfr_bsa > _MALE60_SCR1 else "FAIL",
        ),
        MetricResult(
            "profile_scr_gfr_uses_ckdepi",
            prof_scr.gfr_ml_min,
            expected_scr * 0.999,
            expected_scr * 1.001,
            "mL/min",
            "pass" if abs(prof_scr.gfr_ml_min - expected_scr) < 1e-6 else "FAIL",
        ),
        MetricResult(
            "profile_without_scr_keeps_default",
            prof_default.gfr_ml_min,
            125.0 * 0.90,
            125.0,
            "mL/min",
            "pass" if 112.5 <= prof_default.gfr_ml_min <= 125.0 else "FAIL",
        ),
    ]


def case_ckdepi_2021() -> CaseResult:
    metrics = _metrics()
    male = _MALE60_SCR1
    female = ckdepi_2021_egfr(1.0, 60, female=True)
    bsa = ckdepi_2021_egfr(1.0, 60, female=False, bsa_m2=_BSA175_80)
    profile = metrics[3].predicted
    default = metrics[4].predicted
    notes = [
        f"CKD-EPI 2021 race-free: male 60 y SCR 1.0 -> {male:.1f} "
        f"mL/min/1.73 m2 (CKD-2 band); female same Scr {female:.1f}; "
        f"BSA-scaled {bsa:.1f} mL/min; Scr-carrying profile "
        f"gfr={profile:.1f} mL/min; no-Scr default untouched "
        f"({default:.1f} mL/min)."
    ]
    return CaseResult(
        benchmark="CKD-EPI 2021 race-free GFR baseline (R-6)",
        passed=all(m.criterion == "pass" for m in metrics),
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_ckdepi_2021"]
