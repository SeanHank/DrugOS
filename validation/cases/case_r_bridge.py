"""R literature-PK cross-check (R-1, required-R integration, doc/06 §1).

Runs the five-benchmark corpus through the whole pipeline and asserts the
mandatory R bridge reproduced the numpy PK estimator within 2 % (the R side
implements the literature equations in ``src/drugos/rbridge/literature_pk.R``
— Wagner 1976, Gibaldi & Perrier 1982, Greenblatt & Koch-Weser 1975, Rowland
& Tozer 2010).  This certification:

- guarantees R is genuinely interrogated on every run (hard dependency),
- verifies both estimators implement the same literature rule,
- surfaces token-level agreement in validation/report.md.

The datasets feeding the run are the vendored published-PK corpus + the
measured-hERG/DILI corpora (doc/11 rows 1/8/9); see doc/10 for comparability.
"""

from __future__ import annotations

from validation.benchmarks import BENCHMARKS
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pipeline import run_pipeline, spec_from_benchmark_data

_AGREEMENT_MAX = 0.02
_COMPOUNDS = ("midazolam", "acetaminophen", "warfarin", "ciprofloxacin", "dofetilide")


def case_r_bridge() -> CaseResult:
    worst_frac = 0.0
    worst_compound = ""
    verdicts: dict[str, str] = {}
    two_comp_count = 0
    for name in _COMPOUNDS:
        bench = next(b for b in BENCHMARKS if b.name == name)
        result = run_pipeline(spec_from_benchmark_data(bench))
        verify = result.r_verify
        if verify is None:
            raise RuntimeError(f"R bridge did not run for {name}")
        verdicts[name] = verify.verdict
        if verify.fit == "two_comp":
            two_comp_count += 1
        if verify.agreement_frac > worst_frac:
            worst_frac = verify.agreement_frac
            worst_compound = name
    ok = worst_frac <= _AGREEMENT_MAX and all(v == "r:agree" for v in verdicts.values())

    metrics = [
        MetricResult(
            "r_literature_cl_agreement_max",
            worst_frac,
            0.0,
            _AGREEMENT_MAX,
            "fraction",
            "pass" if worst_frac <= _AGREEMENT_MAX else "FAIL",
        ),
        MetricResult(
            "r_verdicts_not_agree",
            sum(1 for v in verdicts.values() if v != "r:agree"),
            0.0,
            0.0,
            "count",
            "pass" if ok else "FAIL",
        ),
        MetricResult(
            "r_two_comp_fits",
            float(two_comp_count),
            0.0,
            float(len(_COMPOUNDS)),
            "count",
            "pass",
        ),
    ]
    notes = [
        "worst |CL_r - CL_py|/CL_py over "
        + " / ".join(_COMPOUNDS)
        + f": {worst_frac:.2e} ({worst_compound}); "
        + " / ".join(f"{k}={v}" for k, v in verdicts.items())
        + f"; method-of-residuals two-comp fits: {two_comp_count}/{len(_COMPOUNDS)}"
    ]
    return CaseResult(
        "R literature-PK cross-check (R-1)",
        ok,
        metrics,
        notes,
        level=EvidenceLevel.L3_EMPIRICAL,
    )


__all__ = ["case_r_bridge"]
