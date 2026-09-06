"""Benchmark corpus: compound data, citations, and the GMFE pass rule.

Split layout:

- ``base.py``        — ``Benchmark`` data class and ``FOLD_ALLOWANCE``.
- ``citations.py``   — licensing-safe literature references.
- ``compounds/``     — one module per benchmark compound.

Every entry carries *measured* input parameters (literature values) and the
*published* clinical range for the output metrics that the pipeline must
reproduce.  Referencing is licensing-safe: sources are product labelling /
standard pharmacopoeia summaries, no proprietary data.

Pass rule (doc/08, Tier 1/2): a metric passes when the prediction is within
the published range widened by the 2x geometric-mean-fold-error allowance
``center / 2 <= pred <= min(center * 2, 1)`` (upper bound clamped to 1 for
fractions such as Fa).  The report lists both the applied allowance and the
resulting GMFE per benchmark.
"""

from validation.benchmarks.base import FOLD_ALLOWANCE, Benchmark
from validation.benchmarks.citations import CITATIONS
from validation.benchmarks.compounds import BENCHMARKS

__all__ = ["BENCHMARKS", "CITATIONS", "FOLD_ALLOWANCE", "Benchmark"]
