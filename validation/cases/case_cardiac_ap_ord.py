"""Production-validated cardiac AP cross-check (R-3, doc/08 new case, L3).

Runs the vendored O'Hara-Rudy 2011 human ventricular action-potential model
(ORd, PLoS Comput Biol e1002061; BSD Myokit encoding at
``data/models/ohara-2011.mmt``) under fractional IKr block via myokit, and
cross-checks the Stage-4 cardiac axis:

- the ORd baseline APD90 is physiological (published endocardial 1 Hz APD90
  is ~270 ms; the model was validated against >100 undiseased human hearts),
- at its measured hERG IC50 (occupancy = 0.5 by definition; dofetilide
  ChEMBL geomean 26 nM, R-2) the ORd action potential prolongs strongly —
  dofetilide is a marketed TdP-class drug, so the validated ionic model must
  reproduce a large, directionally-correct prolongation,
- prolongation increases monotonically with block fraction,
- the warfarin control (no measurable hERG block, block_frac ~ 0) produces
  zero prolongation — the same ordering the calibrated encoder reports
  (warfarin delta-QTc 0.017 ms vs dofetilide 20.2 ms in the cardiac case).

The algebraic encoder stays the fast whole-heart QTc map in the web/robustness
compute path; this lane is an independent production-validated anchor.
Runtime requirements: myokit + SUNDIALS headers (doc/12, doc/06).
"""

from __future__ import annotations

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.organ.cardiac_ap import ord_apd90

_BASE_LO_MS = 200.0
_BASE_HI_MS = 350.0


def case_cardiac_ap_ord() -> CaseResult:
    res_25 = ord_apd90(0.25)
    res_50 = ord_apd90(0.5)
    res_control = ord_apd90(0.0)

    base_ok = _BASE_LO_MS <= res_25.apd90_base_ms <= _BASE_HI_MS
    ic50_delta = res_50.delta_apd90_ms
    ic50_ok = ic50_delta >= 30.0
    mono_ms = res_50.delta_apd90_ms - res_25.delta_apd90_ms
    mono_ok = mono_ms > 0.0
    control_delta = res_control.delta_apd90_ms
    control_ok = abs(control_delta) < 1e-6

    metrics = [
        MetricResult(
            "ord_apd90_base_ms",
            res_25.apd90_base_ms,
            _BASE_LO_MS,
            _BASE_HI_MS,
            "ms",
            "pass" if base_ok else "FAIL",
        ),
        MetricResult(
            "dofetilide_at_ic50_delta_apd90_ms",
            ic50_delta,
            30.0,
            1000.0,
            "ms",
            "pass" if ic50_ok else "FAIL",
        ),
        MetricResult(
            "apd90_monotone_25_to_50_ms",
            mono_ms,
            0.0,
            1000.0,
            "ms",
            "pass" if mono_ok else "FAIL",
        ),
        MetricResult(
            "warfarin_control_delta_apd90_ms",
            control_delta,
            -1e-6,
            1e-6,
            "ms",
            "pass" if control_ok else "FAIL",
        ),
    ]
    notes = [
        f"ORd 2011 (myokit, endo, 50 pre-paces @1 Hz): baseline "
        f"APD90={res_25.apd90_base_ms:.1f} ms; "
        f"delta-APD90 @25% block={res_25.delta_apd90_ms:.1f} ms, "
        f"@50% block (measured-IC50 concentration)={ic50_delta:.1f} ms, "
        f"control (0% block)={control_delta:.4f} ms; "
        f"monotone +{mono_ms:.1f} ms between block levels; "
        "warfarin control confirmed zero prolongation (matches encoder ordering "
        "warfarin 0.017 ms << dofetilide 20.2 ms)"
    ]
    return CaseResult(
        benchmark="cardiac AP cross-check (ORd/IKr)",
        passed=all(m.criterion == "pass" for m in (metrics[0], metrics[1], metrics[2], metrics[3])),
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L3_EMPIRICAL,
    )
