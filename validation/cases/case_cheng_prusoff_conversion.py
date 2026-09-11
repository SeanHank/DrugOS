"""Cheng-Prusoff IC50 -> Ki conversion (L2 — analytic, doc/08 §1.2).

doc/10 §2 (P2) demands that a measured IC50 is converted to the Kd/Ki a
pipeline binding site consumes through the Cheng & Prusoff (1973)
competitive-inhibition relation `Ki = IC50/(1 + [S]/Km)` — never by
liter-wiring IC50 onto Kd.  This case pins the recipe and its assay-context
default (`[S]/Km = 1`, the competition-binding convention, giving the
conservative Ki = IC50/2).
"""

from __future__ import annotations

import pytest
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.target import Target
from drugos.target.targets import cheng_prusoff_ki_nm


def case_cheng_prusoff_conversion() -> CaseResult:
    metrics: list[MetricResult] = []
    sites = {
        "BSEP (cholestasis)": 90.0,
        "CYP3A4 inhibition": 12.0,
        "hERG (Kv11.1)": 0.1,
    }
    ok = True
    for name, ic50_um in sites.items():
        base = Target(name=name, kd_nm=1_000_000.0)
        kd_default = base.kd_from_ic50_um(ic50_um).kd_nm
        expect_default = ic50_um * 1e3 / 2.0
        good = kd_default == pytest.approx(expect_default, rel=1e-9)
        ok &= good
        metrics.append(
            MetricResult(
                f"ki_{name}_default",
                kd_default,
                expect_default,
                expect_default,
                "nM",
                f"Ki = IC50/(1+1) = IC50/2 = {expect_default:.1f} nM"
                if good
                else f"FAIL (got {kd_default:.1f})",
            )
        )

    ic50 = 1.0
    r0 = cheng_prusoff_ki_nm(ic50, 0.0)
    r1 = cheng_prusoff_ki_nm(ic50, 1.0)
    r9 = cheng_prusoff_ki_nm(ic50, 9.0)
    monotone = r9 < r1 < r0
    ok &= monotone
    metrics.append(
        MetricResult(
            "cheng_prusoff_shape",
            r1,
            0.0,
            r0,
            "nM",
            "ratio 0 -> 1000 nM, ratio 1 -> 500 nM, ratio 9 -> 100 nM, strictly falling"
            if monotone
            else "FAIL",
        )
    )

    err = True
    try:
        cheng_prusoff_ki_nm(0.0, 1.0)
        err = False
    except ValueError:
        pass
    ok &= err
    for name, _ in (("non-positive_ic50", err),):
        metrics.append(
            MetricResult(
                name,
                float(err),
                1.0,
                1.0,
                "flag",
                "raises on non-positive IC50" if err else "FAIL",
            )
        )

    return CaseResult(
        "Cheng-Prusoff IC50->Ki conversion",
        ok,
        metrics,
        [
            "measured IC50 -> Ki via Ki = IC50/(1 + [S]/Km); "
            "doc/10 P2 'never liter-wire IC50->Kd' is now code",
            f"default assay convention [S]/Km = 1 sets Ki = IC50/2 "
            f"({sites['BSEP (cholestasis)']} uM -> {sites['BSEP (cholestasis)'] * 500:.0f} nM)",
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_cheng_prusoff_conversion"]
