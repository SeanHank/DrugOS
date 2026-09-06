"""Dofetilide: hERG-active Class III anti-arrhythmic, renal elimination.

Primary purpose: renal-elimination axis; anchors Stage 4 cardiovascular QT.

"""

from validation.benchmarks.base import Benchmark, load_about, load_cites, load_published

BENCHMARK = Benchmark(
    name="dofetilide",
    smiles="COc1cc(NC(=O)Nc2ccc(NS(=O)(=O)CCN(C)C)cc2)ccc1OC",
    log_p=2.29,
    pka_acids=[],
    pka_bases=[7.5],
    fup=0.36,
    bp=0.86,
    about=load_about("dofetilide"),
    route="oral",
    dose_mg=0.5,
    tmax_h=72.0,
    published=load_published("dofetilide"),
    assert_only=("cl_plasma_l_h", "t_half_h", "f_abs", "urine_fraction"),
    references=load_cites("dofetilide"),
)

__all__ = ["BENCHMARK"]
