"""Ciprofloxacin: renal-elimination fluoroquinolone, zwitterion (Tier 1).

Primary purpose: renal elimination + zwitterionic partition-permeability axis.
"""

from validation.benchmarks.base import Benchmark, load_about, load_cites, load_published

BENCHMARK = Benchmark(
    name="ciprofloxacin",
    smiles="C1CN(CCN1)C1=C(F)C(=O)C(C(=O)O)=CN1C2CC2",
    log_p=0.28,
    pka_acids=[6.1],
    pka_bases=[8.7],
    fup=0.7,
    bp=0.83,
    about=load_about("ciprofloxacin"),
    route="oral",
    dose_mg=500.0,
    tmax_h=48.0,
    published=load_published("ciprofloxacin"),
    assert_only=("cl_plasma_l_h", "t_half_h", "f_abs", "urine_fraction"),
    references=load_cites("ciprofloxacin"),
)

__all__ = ["BENCHMARK"]
