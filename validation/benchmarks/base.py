"""Benchmark data class and the geometric-mean-fold-error pass rule.

Original experimental PK bands are **not** inlined in the compound modules:
they are loaded from the vendored, checksummed dataset
``data/benchmarks/published_pk.json`` (single source of truth). A missing or
malformed data file is an explicit error — never a silent default.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

FOLD_ALLOWANCE = 2.0

DATA_FILE = Path(__file__).resolve().parents[2] / "data" / "benchmarks" / "published_pk.json"

_CORPUS: dict[str, object] | None = None


def _corpus() -> dict[str, object]:
    global _CORPUS
    if _CORPUS is not None:
        return _CORPUS
    if not DATA_FILE.exists():
        raise FileNotFoundError(
            f"missing vendored benchmark data {DATA_FILE} "
            "(run from the repo checkout; do not delete data/)"
        )
    try:
        _CORPUS = json.loads(DATA_FILE.read_text(encoding="utf-8"))["compounds"]
    except (ValueError, KeyError, TypeError) as exc:
        raise RuntimeError(f"corrupt benchmark data file {DATA_FILE}: {exc}") from exc
    return _CORPUS


def load_published(name: str) -> dict[str, tuple[float, float]]:
    """Original experimental PK pass-bands (lo, hi) for a benchmark compound."""
    entry = _corpus().get(name)
    if entry is None or not isinstance(entry, Mapping):
        raise KeyError(f"benchmark '{name}' missing from {DATA_FILE}")
    published_raw = entry.get("published")
    if not isinstance(published_raw, Mapping):
        raise ValueError(f"benchmark '{name}': missing 'published' bands in {DATA_FILE}")
    return {metric: (float(lo), float(hi)) for metric, (lo, hi) in published_raw.items()}


def load_cites(name: str) -> tuple[str, ...]:
    """Citation keys (into citations.py) for a benchmark compound."""
    entry = _corpus().get(name)
    if entry is None or not isinstance(entry, Mapping):
        raise KeyError(f"benchmark '{name}' missing from {DATA_FILE}")
    cites = entry.get("cites", ())
    return tuple(cites) if isinstance(cites, (list, tuple)) else ()


def load_about(name: str) -> str:
    """One-line mechanism note for a benchmark compound."""
    entry = _corpus().get(name)
    if entry is None or not isinstance(entry, Mapping):
        raise KeyError(f"benchmark '{name}' missing from {DATA_FILE}")
    about = entry.get("about", "")
    return str(about)


@dataclass(frozen=True)
class Benchmark:
    """A benchmark compound with literature inputs and published output bands.

    Published metric bands are encoded as (lo, hi) on the *clinical* value;
    the validation pass rule applies the 2x fold-error allowance defined
    above, so ``hi`` here is the raw literature band, not the allowance.
    """

    name: str
    smiles: str
    log_p: float
    pka_acids: list[float]
    pka_bases: list[float]
    fup: float
    bp: float
    about: str
    route: str
    dose_mg: float
    tmax_h: float
    n_eval: int = 400
    fa_override: float | None = None
    published: dict[str, tuple[float, float]] = field(default_factory=dict)
    assert_only: tuple[str, ...] = ()
    references: tuple[str, ...] = ()


__all__ = ["FOLD_ALLOWANCE", "Benchmark", "load_published", "load_cites", "load_about"]
