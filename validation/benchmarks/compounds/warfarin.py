"""Warfarin: high-protein-binding narrow-therapeutic-index anticoagulant.

Primary purpose: high fup-binding % / small apparent Vd axis (doc/08 Tier 1).
"""

from validation.benchmarks.base import Benchmark, load_about, load_cites, load_published

BENCHMARK = Benchmark(
    name="warfarin",
    smiles="Cc1ccc(C(=O)CC(c2ccccc2O)c2ccc(cc2)C(=O)O)cc1",
    log_p=2.7,
    pka_acids=[5.0],
    pka_bases=[],
    fup=0.013,
    bp=0.55,
    about=load_about("warfarin"),
    route="oral",
    dose_mg=10.0,
    tmax_h=96.0,
    n_eval=800,
    published=load_published("warfarin"),
    assert_only=("cl_plasma_l_h", "vss_l", "t_half_h"),
    references=load_cites("warfarin"),
)

__all__ = ["BENCHMARK"]
