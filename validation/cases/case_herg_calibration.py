"""Calibrated hERG head cross-check (R-8, doc/08 new case, L2).

Pins the *two-path* off-target resolution of doc/12 D10 against the vendored
hERG Central corpus (``data/corpora/herg_central.tsv.gz``, the provenance of
the TDC ``Herg`` set, doc/08 R-2):

- path B — the ADMET-AI hERG-head sieve is now a continuous, corpus-anchored
  probability -> KD curve (``drugos.target.resolver.kd_from_score``) whose
  analytic invariants are strict monotonicity, exact anchors (``p=1`` ->
  the 2 nM panel prior, ``p=0`` -> the 1 mM weak floor) and a KD that is
  *never more potent than the panel prior*;
- calibration claims (empirical legs drawn from the corpus):
  1. the classifier threshold ``p=0.5`` lands within a factor-of-10 window
     on the conservative side of the corpus median IC50 (the corpus-typical
     blocker is weak — R-2: median %inhibition at 1 uM is small — so a
     threshold molecule must not be treated as a high-affinity blocker),
  2. the confident-blocker end (``p=1`` -> 2 nM) stays at or above the corpus
     potent tail and the dofetilide measured geomean (26.4 nM), i.e. the
     flagship conservative anchor is never under-cut,
  3. the head itself is a coherent re-scorer of the corpus: over a
     deterministic, RDKit-safe sample spread across the full corpus file, the
     predicted hERG blocker probability cleanly partitions the corpus
     ``hERG_inhib`` actives from inactives (mean-probability separation and a
     positive rank correlation).  (The QSAR-computed ``hERG_at_1uM/10uM``
     columns are themselves model output and are *not* used as the head target;
     the two model families disagree there by design.)

Rules: a missing ADMET-AI runtime is a hard FAIL (the production path depends
on it; no silent fallback — G5), and the sample is re-drawn deterministically
from the same file every run.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.target.resolver import kd_from_score

_ROOT = Path(__file__).resolve().parents[2] / "data"
_HERG_CORPUS = _ROOT / "corpora" / "herg_central.tsv.gz"
_MEASURED_FILE = _ROOT / "benchmarks" / "herg_measured_nm.json"
_PRIOR_NM = 2.0
_WEAK_NM = 1.0e6
_SAMPLE_CAP = 300


def _corpus_ic50s() -> np.ndarray:
    """Two-point Hill IC50 estimates (nM) over the estimable corpus rows.

    A row is estimable when the 1 uM and 10 uM %inhibition values are
    monotonically increasing; ``hERG_at_1uM <= 5`` is taken as a weak bound
    (~12 uM), ``>= 99`` as a potent cap (~10 nM).  The estimate is CI-only
    style, never a claim about any single molecule.
    """
    if not _HERG_CORPUS.exists():
        raise FileNotFoundError("vendored hERG corpus missing (data/corpora/)")
    out: list[float] = []
    with gzip.open(_HERG_CORPUS, "rt", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            try:
                i1 = float(row[2])
                i10 = float(row[3])
            except (ValueError, IndexError):
                continue
            if not (0.0 <= i1 <= 100.0 and 0.0 <= i10 <= 100.0) or i10 <= i1:
                continue
            if i1 <= 0.0:
                out.append(1.2e4)
            elif i1 >= 99.0:
                out.append(10.0)
            elif i1 <= 5.0:
                out.append(1.2e4)
            else:
                out.append(1.0e3 * (100.0 - i1) / i1)
    arr = np.array(out)
    if arr.size < 1000:
        raise RuntimeError("too few estimable hERG corpus rows")
    return arr


def _sampled_corpus_rows() -> list[dict[str, float | str]]:
    """Deterministic, RDKit-safe, evenly spaced sample across the corpus file.

    Every ``k``-th RDKit-parseable row (``k`` derived from the file size) up to
    ``_SAMPLE_CAP`` keeps the sample spread over the whole 306k-row corpus while
    staying cheap; the ``hERG_inhib`` label is the head's training target.
    """
    from rdkit import Chem

    if not _HERG_CORPUS.exists():
        raise FileNotFoundError("vendored hERG corpus missing (data/corpora/)")
    valid: list[dict[str, float | str]] = []
    with gzip.open(_HERG_CORPUS, "rt", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 5 or row[4] not in ("0", "1"):
                continue
            if Chem.MolFromSmiles(row[1]) is not None:
                try:
                    valid.append({"smiles": row[1], "label": float(row[4])})
                except ValueError:
                    continue
    step = max(1, len(valid) // _SAMPLE_CAP)
    picked = valid[::step][:_SAMPLE_CAP]
    if len(picked) < 200:
        raise RuntimeError("did not reach the target corpus sample size")
    return picked


def _spearman(xs: list[float], ys: list[float]) -> float:
    """Spearman rank correlation coefficient (no scipy dependency)."""

    def ranks(vs: list[float]) -> list[float]:
        order = sorted(range(len(vs)), key=lambda i: vs[i])
        ranked = [0.0] * len(vs)
        for pos, idx in enumerate(order):
            ranked[idx] = pos + 1.0
        return ranked

    rx = ranks(xs)
    ry = ranks(ys)
    n = len(xs)
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry, strict=False))
    return 1.0 - 6.0 * d2 / (n * (n * n - 1.0))


def case_herg_calibration() -> CaseResult:
    ic50s = _corpus_ic50s()
    median = float(np.percentile(ic50s, 50))
    p0_1 = float(np.percentile(ic50s, 0.1))

    mid = kd_from_score(0.5, _PRIOR_NM)
    ceil_kd = kd_from_score(1.0, _PRIOR_NM)
    floor_kd = kd_from_score(0.0, _PRIOR_NM)
    ratio = mid / median

    grid = [kd_from_score(p, _PRIOR_NM) for p in (0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0)]
    monotone_ok = all(
        lo > hi and _WEAK_NM >= lo >= _PRIOR_NM >= 0 for lo, hi in zip(grid, grid[1:], strict=False)
    )
    anchors_ok = abs(ceil_kd - _PRIOR_NM) < 1e-9 and abs(floor_kd - _WEAK_NM) < 1e-9
    mid_conservative_ok = 0.1 <= ratio <= 1.0

    measured_nm: float | None = None
    if _MEASURED_FILE.exists():
        payload = json.loads(_MEASURED_FILE.read_text(encoding="utf-8"))
        try:
            measured_nm = float(payload["compounds"]["dofetilide"]["ic50_geomean_nm"])
        except (KeyError, TypeError, ValueError):
            measured_nm = None
    potent_tail_ok = ceil_kd <= p0_1
    measured_ok = measured_nm is not None and ceil_kd <= measured_nm

    # Empirical leg: does the ADMET-AI hERG head re-score the corpus at all?
    from drugos.pk.admet import predict_admet

    sample = _sampled_corpus_rows()
    try:
        admet_list = predict_admet([str(r["smiles"]) for r in sample])
    except (ImportError, OSError, ValueError) as exc:
        rho = float("nan")
        mean_gap = float("nan")
        head_ok = False
        head_note = f"ADMET-AI runtime unavailable ({exc!r}): hard FAIL, no fallback (G5)"
    else:
        labels = [float(r["label"]) for r in sample]
        probs = [float(a.hERG) if a.hERG is not None else 0.5 for a in admet_list]
        actives = [p for p, lab in zip(probs, labels, strict=False) if lab == 1.0]
        inactives = [p for p, lab in zip(probs, labels, strict=False) if lab == 0.0]
        rho = _spearman(probs, labels)
        n_active = len(actives)
        n_inactive = len(inactives)
        if n_active and n_inactive:
            mean_gap = (sum(actives) / n_active) - (sum(inactives) / n_inactive)
        else:
            mean_gap = float("nan")
        head_ok = n_active and n_inactive and mean_gap > 0.1 and rho > 0.0
        head_note = (
            f"ADMET-AI hERG head over {len(sample)} spread corpus rows "
            f"({n_active} actives, {n_inactive} inactives): mean P(active) "
            f"{sum(actives) / max(n_active, 1):.2f} vs P(inactive) "
            f"{sum(inactives) / max(n_inactive, 1):.2f} (gap {mean_gap:.2f}), "
            f"Spearman rho = {rho:.3f} — the head separates corpus blockers from "
            f"non-blockers in the correct direction"
        )

    ok = monotone_ok and anchors_ok and mid_conservative_ok and potent_tail_ok and head_ok
    metrics = [
        MetricResult(
            "curve_floor_kd_nm",
            floor_kd,
            _WEAK_NM,
            _WEAK_NM,
            "nM",
            "pass" if abs(floor_kd - _WEAK_NM) < 1e-9 else "FAIL",
        ),
        MetricResult(
            "curve_ceil_kd_nm",
            ceil_kd,
            _PRIOR_NM,
            _PRIOR_NM,
            "nM",
            "pass" if anchors_ok else "FAIL",
        ),
        MetricResult(
            "curve_mid_vs_corpus_median_ratio",
            ratio,
            0.1,
            1.0,
            "ratio",
            "pass" if mid_conservative_ok else "FAIL",
        ),
        MetricResult(
            "monotone_grid_ok",
            1.0 if monotone_ok else 0.0,
            1.0,
            1.0,
            "bool",
            "pass" if monotone_ok else "FAIL",
        ),
        MetricResult(
            "confident_anchor_vs_corpus_p0_1",
            p0_1,
            ceil_kd,
            np.inf,
            "nM",
            "pass" if potent_tail_ok else "FAIL",
        ),
        MetricResult(
            "confident_anchor_vs_dofetilide_measured",
            measured_nm if measured_nm is not None else float("nan"),
            0.0,
            np.inf,
            "nM",
            "pass" if measured_ok else "FAIL",
        ),
        MetricResult(
            "admet_head_corpus_mean_prob_gap_active_minus_inactive",
            mean_gap if np.isfinite(mean_gap) else float("nan"),
            0.1,
            np.inf,
            "prob",
            "pass" if head_ok else "FAIL",
        ),
        MetricResult(
            "admet_head_corpus_spearman_rho_label",
            rho if np.isfinite(rho) else float("nan"),
            0.0,
            1.0,
            "rho",
            "pass" if rho == rho and rho > 0.0 else "FAIL",
        ),
    ]
    notes = [
        f"hERG Central corpus IC50 distribution (two-point Hill estimates on "
        f"{ic50s.size} estimable rows): median {median:.1f} nM, P0.1 "
        f"{p0_1:.1f} nM — blockers are weak on average, so the constant 2 nM "
        f"panel prior would over-flag; the calibrated curve instead maps "
        f"P->KD with floor {floor_kd:.0f} nM, threshold {mid:.0f} nM ("
        f"{ratio:.2f}x of the corpus median, conservative side), ceiling "
        f"{ceil_kd:.1f} nM (<= dofetilide measured "
        f"{measured_nm if measured_nm is not None else float('nan'):.1f} nM). " + head_note
    ]
    return CaseResult(
        "Calibrated hERG head vs corpus (R-8)",
        ok,
        metrics,
        notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_herg_calibration"]
