"""Acetaminophen: low-binding analgesic, oral ACAT-lite absorption (Tier 1).

Primary purpose: oral-absorption + hepatic-metabolic clearance axis.
"""

from validation.benchmarks.base import Benchmark, load_about, load_cites, load_published

BENCHMARK = Benchmark(
    name="acetaminophen",
    smiles="CC(=O)Nc1ccc(O)cc1",
    log_p=0.51,
    pka_acids=[9.4],
    pka_bases=[],
    fup=0.8,
    bp=1.17,
    about=load_about("acetaminophen"),
    route="oral",
    dose_mg=1000.0,
    tmax_h=24.0,
    fa_override=0.88,
    published=load_published("acetaminophen"),
    assert_only=("cl_plasma_l_h", "t_half_h", "f_abs"),
    references=load_cites("acetaminophen"),
)

__all__ = ["BENCHMARK"]
