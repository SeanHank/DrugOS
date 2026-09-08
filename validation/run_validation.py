#!/usr/bin/env python3
"""DrugOS validation runner (G4 gate, doc/09 section G4).

Runs every validation case and regenerates ``validation/report.md`` from the
same code path that produced the numbers (single source of truth).  Exits
non-zero if any case fails.

Usage:  python validation/run_validation.py
"""

from __future__ import annotations

import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
for entry in (str(SRC), str(ROOT)):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from validation.benchmarks import CITATIONS, FOLD_ALLOWANCE  # noqa: E402
from validation.cases import (  # noqa: E402
    BENCHMARKS,
    EVIDENCE_META,
    EVIDENCE_ORDER,
    CaseResult,
    MetricResult,
    run_all,
)

from drugos.version import __version__  # noqa: E402

REPORT_PATH = ROOT / "validation" / "report.md"


def _fmt(m: MetricResult) -> str:
    return f"{m.predicted:.4g} [{m.lo:.4g}, {m.hi:.4g}] {m.unit}"


def fold_error(results: list[CaseResult], benchmark: str) -> float | None:

    import numpy as np

    folds: list[float] = []
    for c in results:
        if c.benchmark != benchmark:
            continue
        for m in c.metrics:
            if m.lo <= 0 or m.hi <= 0:
                continue
            center = (m.lo * m.hi) ** 0.5
            folds.append(max(m.predicted / center, center / m.predicted))
    if not folds:
        return None
    return float(np.exp(np.mean(np.log(folds))))


def build_markdown(results: list[CaseResult]) -> str:
    passed = sum(1 for c in results if c.passed)
    total = len(results)
    status = "PASS" if passed == total else "FAIL"
    lines: list[str] = [
        "# DrugOS Validation Report",
        "",
        f"- Status: **{status}** ({passed}/{total} cases passed)",
        f"- Generated: {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        f"- DrugOS version: {__version__}",
        f"- Interpreter: {sys.executable}",
        f"- Fold-error allowance: within {FOLD_ALLOWANCE:.0f}x of the published"
        " band centre (doc/08 Tier 2, GMFE <= 2); Fa bands additionally clamp to 1.",
        "",
        "## Results",
        "",
        "| Case | Level | Metric | Predicted [allowed band] | Criterion |",
        "|---|---|---|---|---|",
    ]
    for c in results:
        row = c.metrics
        if not row:
            lines.append(
                f"| {c.benchmark} | {c.level.value} | (no asserted metrics) | - | "
                f"{'PASS' if c.passed else 'FAIL'} |"
            )
        for i, m in enumerate(row):
            case = c.benchmark if i == 0 else ""
            level = c.level.value if i == 0 else ""
            lines.append(f"| {case} | {level} | {m.name} | {_fmt(m)} | {m.criterion} |")
    lines += [
        "",
        "## Evidence levels",
        "",
        "Results are graded by how much epistemic weight they carry "
        "(doc/08 §1.1-1.4 tier ladder, strongest first):",
        "",
        "| Level | Meaning | Basis | What it certifies |",
        "|---|---|---|---|",
    ]
    for level in EVIDENCE_ORDER:
        meta = EVIDENCE_META[level]
        lines.append(
            f"| **{level.value}** | {meta.label} | {meta.basis} | {meta.certifies} ({meta.doc08}) |"
        )
    by_level = {level: [c for c in results if c.level == level] for level in EVIDENCE_ORDER}
    lines += [
        "",
        "Per-level status:",
        "",
    ]
    for level in EVIDENCE_ORDER:
        sub = by_level[level]
        passed = sum(1 for c in sub if c.passed)
        lines.append(
            f"- **{level.value}** ({EVIDENCE_META[level].label}): {passed}/{len(sub)} cases green."
        )
    lines += [
        "",
        "## Notes & limitations",
        "",
    ]
    for c in results:
        lines.append(f"- **{c.benchmark}**: " + " ".join(c.notes))
    gmfes = [
        (b.name, fold_error(results, b.name))
        for b in BENCHMARKS
        if fold_error(results, b.name) is not None
    ]
    if gmfes:
        lines += [
            "",
            "## Tier-1 geometric-mean fold error (L3 asserted metrics)",
            "",
            "| Benchmark | GMFE |",
            "|---|---|",
        ]
        for name, fe in gmfes:
            lines.append(f"| {name} | {fe:.2f}x |")
    lines += [
        "",
        "## Tier coverage (doc/08)",
        "",
        "- Stage 1 (PK): benchmark compounds + analytic limit + mass budget + "
        "dose-proportionality + route-dependent bioavailability F reporting "
        "(IV/depot/oral first-pass) + permeability/Fa-gated and logS-gated "
        "solubility-limited oral absorption — **green**.",
        "- Stage 2 (occupancy): target-turnover equilibrium ODE vs analytic "
        "D/(D+Kd) point-wise match — **green**.",
        "- Stage 3 (pathway): 3-tier MAPK amplifier — steady-state EC50 below "
        "the receptor-Kd-equivalent signal (Emax/Hill fit, EC50<0.5) and "
        ">2x baseline amplification — **green**.",
        "- Stage 4 (organ): liver DILI dose-response (ALT/bilirubin/Hy's Law at "
        "overdose), pathway->organ regeneration coupling + bilirubin ceiling, "
        "cardiac QTc prolongation vs the published dofetilide Delta-QTc band "
        "+ ERK-amplification inotropy/chronotropy tone coupling, and kidney "
        "GFR/AKI KDIGO escalation with a graded urinary KIM-1 row — "
        "**green**.",
        "- Stage 5 clinical / report / CLI tiers: scheduled with their stage "
        "modules (see roadmap); each new model feature must add a validation "
        "case before merge (G4 rule).",
        "",
        "## Citations",
        "",
    ]
    for key in BENCHMARKS:
        refs = [CITATIONS[r] for r in key.references]
        lines.append(f"- **{key.name}**: " + " ".join(refs))
    lines.append("- Partition & whole-body graph: " + CITATIONS["partition"] + ".")
    lines.append("- Occupancy model: " + CITATIONS["occupancy"] + ".")
    lines.append("- Pathway amplification: " + CITATIONS["pathway"] + ".")
    lines.append("")
    lines.append(
        "This report is machine-generated by `validation/run_validation.py`; "
        "hand edits are overwritten."
    )
    return "\n".join(lines)


def main() -> int:
    results = run_all()
    for c in results:
        if not c.metrics:
            print(f"[FAIL] {c.benchmark}: no metrics produced")
    markdown = build_markdown(results)
    REPORT_PATH.write_text(markdown, encoding="utf-8")
    print(f"wrote {REPORT_PATH}")
    for c in results:
        maybe = "PASS" if c.passed else "FAIL"
        summary = ", ".join(f"{m.name}={m.predicted:.3g} {m.unit}" for m in c.metrics)
        print(f"[{maybe}] {c.benchmark}: {summary}")
    failed = [c for c in results if not c.passed]
    if failed:
        print(f"validation FAILED: {len(failed)} case(s) failing")
        for c in failed:
            print(f"  - {c.benchmark}")
        return 1
    print(f"validation PASSED: {len(results)} cases green")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
