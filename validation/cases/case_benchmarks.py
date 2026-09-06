"""Tier-1 benchmark compounds (L3 — empirically anchored, doc/08 §1.1)."""

from __future__ import annotations

from validation.benchmarks import Benchmark
from validation.cases.base import (
    CaseResult,
    EvidenceLevel,
    MetricResult,
    _allowed_band,
    _fraction_absorbed,
    _simulate,
)


def evaluate_benchmark(b: Benchmark) -> CaseResult:
    res = _simulate(b)
    m = res.pk_metrics()
    renal = 0.0
    if res.urine_cum_mg is not None:
        renal = float(res.urine_cum_mg[-1]) / max(b.dose_mg, 1e-12)
    candidates: list[tuple[str, float, tuple[float, float], str]] = [
        ("cl_plasma_l_h", m.cl_l_h, b.published.get("cl_plasma_l_h", (0.0, 0.0)), "L/h"),
        ("vss_l", m.vss_l, b.published.get("vss_l", (0.0, 0.0)), "L"),
        ("t_half_h", m.term_half_life_h, b.published.get("t_half_h", (0.0, 0.0)), "h"),
        ("f_abs", _fraction_absorbed(b, res), b.published.get("f_abs", (0.0, 0.0)), "fraction"),
        ("urine_fraction", renal, b.published.get("urine_fraction", (0.0, 0.0)), "fraction"),
    ]
    metrics: list[MetricResult] = []
    for name, pred, band, unit in candidates:
        lo, hi = band
        if hi <= 0.0:
            continue
        allowance_lo, allowance_hi = _allowed_band(lo, hi, unit)
        ok = allowance_lo <= pred <= allowance_hi
        if ok and name not in b.assert_only:
            continue
        metrics.append(
            MetricResult(
                name=name,
                predicted=pred,
                lo=allowance_lo,
                hi=allowance_hi,
                unit=unit,
                criterion="pass" if ok else "FAIL",
            )
        )
    assert metrics, f"benchmark {b.name} gave no asserted metrics"
    passed = all(s.criterion == "pass" for s in metrics)
    reported = ", ".join(
        f"{k}={v:.3g} {u}"
        for k, v, u in (
            ("CL", m.cl_l_h, "L/h"),
            ("Vss(MRT)", m.vss_l, "L"),
            ("t1/2", m.term_half_life_h, "h"),
        )
    )
    notes = [f"reported: {reported}"]
    if "vss_l" not in b.published and "t_half_h" not in b.published:
        notes.append(
            "Vss/t1/2 reported but not asserted: lumped R&R partition "
            "overpredicts the apparent Vss of low-fup lipophilic bases "
            "(doc/08 risk #4); clearance axis is the validated output."
        )
    return CaseResult(
        benchmark=b.name,
        passed=passed,
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L3_EMPIRICAL,
    )


__all__ = ["evaluate_benchmark"]
