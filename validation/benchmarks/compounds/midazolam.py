"""Midazolam: hepatic CYP3A4-extracted benzodiazepine (doc/08 Tier 1).

Primary purpose: PBPK accuracy on the hepatic-metabolism axis.
"""

from validation.benchmarks.base import Benchmark, load_about, load_cites, load_published

BENCHMARK = Benchmark(
    name="midazolam",
    smiles="Cc1ncn(-c2ccc(Cl)c(c2)C(=O)Nc2cccc(F)c2)c1C",
    log_p=3.94,
    pka_acids=[],
    pka_bases=[6.04],
    fup=0.02,
    bp=0.82,
    about=load_about("midazolam"),
    route="iv_bolus",
    dose_mg=5.0,
    tmax_h=24.0,
    published=load_published("midazolam"),
    assert_only=("cl_plasma_l_h",),
    references=load_cites("midazolam"),
)

__all__ = ["BENCHMARK"]
