"""Benchmark registry: the curated Tier-1 compound corpus (doc/08 §1.1)."""

from validation.benchmarks.base import Benchmark
from validation.benchmarks.compounds.acetaminophen import BENCHMARK as ACETAMINOPHEN
from validation.benchmarks.compounds.ciprofloxacin import BENCHMARK as CIPROFLOXACIN
from validation.benchmarks.compounds.dofetilide import BENCHMARK as DOFETILIDE
from validation.benchmarks.compounds.midazolam import BENCHMARK as MIDAZOLAM
from validation.benchmarks.compounds.warfarin import BENCHMARK as WARFARIN

#: Tier-1 corpus in iteration order (also the validation-suite order).
BENCHMARKS: tuple[Benchmark, ...] = (
    MIDAZOLAM,
    ACETAMINOPHEN,
    WARFARIN,
    CIPROFLOXACIN,
    DOFETILIDE,
)

__all__ = ["BENCHMARKS"]
