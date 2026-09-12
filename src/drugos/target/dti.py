"""General off-target DTI resolver: ChEMBL fingerprint kNN (doc/12 L15, path A).

For a novel molecule carrying no measured affinity, every panel site backed by
public record is now re-scored from *structure* instead of a class median:

- the training set is the vendored, size-capped bioactivity snapshot of
  ``data/chembl/offtarget_snapshot.json`` (measured ``IC50``/``Ki`` rows,
  geomean-central per canonical SMILES, potency-stratified for review),
- the query molecule and every snapshot molecule are featurized as Morgan
  fingerprints (RDKit, radius 2, 2048 bits),
- the predicted pChemBL is a Tanimoto-similarity-weighted mean over the ``k``
  nearest snapshot neighbours for the *same target* (a nearest-estimator of
  the molecule's affinity to that site — never extrapolated outside the
  observed potency band of the target),
- ``resolve_offtarget_panel`` then re-binds every mapped panel site whose
  query has *chemotype support* in the snapshot to the predicted KD
  (``low_confidence=False``, reference naming the kNN target) and returns the
  remaining ``no_public_data_sites``: the panel entries that stay on their
  disclosed class priors — the sites with no public bioactivity row in the
  snapshot (MRP3/MRP4, mitochondrial complexes II-IV and the pyruvate carrier)
  plus any mapped site whose nearest-neighbour Tanimoto similarity to the
  query falls below ``MIN_NEIGHBOR_TANIMOTO`` (a ligand-similarity estimator
  sitting in a chemotype desert is noise, not evidence, so it is disclosed on
  the prior instead of silently driving occupancy).

The CNS grading anchor (doc/12 L12) is the same resolver's *worst-case* result
over the CNS-liability set — the most potent predicted binding to the curated
brain-exposed off-target panel (dopamine D2, muscarinic M1, alpha-1A
adrenergic, mu-opioid and histamine H1 receptors) — but only from targets
whose snapshot has chemotype support for the query, capped at the conservative
class-prior default.  Conservatism: a genuinely CNS-active molecule whose
chemotype is represented in the record (e.g. haloperidol, DRD2 top-match
Tanimoto 0.67) predicts a potent (low nM) value that is honored; a clean
molecule sitting in a chemotype desert (e.g. warfarin, top-match 0.25) fails
the support gate and the seam discloses the bare 0.20 prior instead of a
structure-derived false potency.

Admission rule (G4): a predictor must ship a corpus-calibration + equivalence
case before it is wired; that is ``validation/cases/case_dti_resolver_calibration``,
which leaves molecules out of the snapshot, recomputes the kNN on the held-in
set and pins the leave-out correlation and the within-2x accuracy against the
measured rows.
"""

from __future__ import annotations

import json
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import cast

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

rdkit.RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]

from drugos.target.resolver import resolve_admet_panel  # noqa: E402
from drugos.target.targets import Target  # noqa: E402

K_DEFAULT = 7
_RADIUS = 2
_NBITS = 2048
_CNS_DEFAULT_IC50_NM = 1.0e5  # same conservative default the pipeline discloses

# Chemotype-support gate for the evidence seam (L12 CNS anchor, L13 panel
# re-bind).  The kNN is a nearest estimator: when the query's best Tanimoto
# match in a target's snapshot rows is negligible the "prediction" is single
# weak-similarity-neighbour noise, not evidence.  Below this threshold the
# prediction is refused and the site/endpoint stays on its disclosed prior.
# 0.35 separates the flagships honestly — haloperidol/DRD2 0.67 (anchor),
# diazepam/OPRM1 0.43 (anchor), dofetilide/warfarin/acetaminophen 0.15-0.30
# (no anchor).
MIN_NEIGHBOR_TANIMOTO = 0.35

# Panel site -> snapshot target key.  Only sites with public bioactivity rows
# in the vendored snapshot are mapped; every other site stays on its prior and
# is listed in ``no_public_data_sites`` by ``resolve_offtarget_panel`` (doc/12
# L13/L15).
SITE_TO_TARGET: dict[str, str] = {
    "hERG (Kv11.1)": "KCNH2",
    "CYP3A4 inhibition": "CYP3A4",
    "CYP2D6 inhibition": "CYP2D6",
    "CYP2C9 inhibition": "CYP2C9",
    "BSEP (cholestasis)": "ABCB11",
    "OATP1B1": "SLCO1B1",
    "P-gp (MDR1)": "ABCB1",
    "Glucocorticoid receptor": "NR3C1",
    "Estrogen receptor": "ESR1",
    "Androgen receptor": "AR",
    "Mitochondrial complex I": "NADH_DEHYDROGENASE",
}

# CNS-liability off-target set for a structure-derived CNS potency anchor
# (doc/12 L12).  A molecule's CNS IC50 is the most potent prediction across
# these brain-exposed receptors.
CNS_TARGETS: tuple[str, ...] = ("DRD2", "CHRM1", "ADRA1A", "OPRM1", "HRH1")


@dataclass(frozen=True)
class _TargetRows:
    smiles: tuple[str, ...]
    pchembl: tuple[float, ...]
    pchembl_min: float
    pchembl_max: float


@dataclass(frozen=True)
class Snapshot:
    """The vendored measured-bioactivity training set, keyed by target."""

    targets: dict[str, _TargetRows] = field(default_factory=dict)

    def __contains__(self, target: str) -> bool:
        return target in self.targets


@lru_cache(maxsize=1)
def load_snapshot(path: Path | None = None) -> Snapshot:
    """Load and cache the vendored off-target snapshot."""
    if path is None:
        path = Path(__file__).resolve().parents[3] / "data" / "chembl" / "offtarget_snapshot.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    targets: dict[str, _TargetRows] = {}
    for name, info in raw["targets"].items():
        rows = info["rows"]
        targets[name] = _TargetRows(
            smiles=tuple(str(r["smiles"]) for r in rows),
            pchembl=tuple(float(r["pchembl"]) for r in rows),
            pchembl_min=float(info["pchembl_min"]),
            pchembl_max=float(info["pchembl_max"]),
        )
    return Snapshot(targets=targets)


@lru_cache(maxsize=4096)
def _fingerprint(smiles: str) -> DataStructs.ExplicitBitVect | None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return cast(
            DataStructs.ExplicitBitVect,
            AllChem.GetMorganFingerprintAsBitVect(  # type: ignore[attr-defined]
                mol, radius=_RADIUS, nBits=_NBITS
            ),
        )


def _tanimoto(a: DataStructs.ExplicitBitVect, b: DataStructs.ExplicitBitVect) -> float:
    return float(DataStructs.TanimotoSimilarity(a, b))


def predict_pchembl(target: str, smiles: str, k: int = K_DEFAULT) -> float | None:
    """Predicted pChemBL (-log10 IC50/M) of ``smiles`` against ``target``.

    Tanimoto-similarity-weighted mean over the ``k`` nearest snapshot
    neighbours of the same target — always a weighted mean of *in-band*
    measured values, so a prediction can never extrapolate outside the
    target's observed potency band.  Returns ``None`` when the target has no
    snapshot rows, the SMILES does not parse, or all nearest neighbours are
    orthogonal to the query.
    """
    rows = load_snapshot().targets.get(target)
    if rows is None or not rows.smiles:
        return None
    fp = _fingerprint(smiles)
    if fp is None:
        return None
    neighbours: list[tuple[float, float]] = []
    for smi, pchembl in zip(rows.smiles, rows.pchembl, strict=True):
        other = _fingerprint(smi)
        if other is None:
            continue
        sim = _tanimoto(fp, other)
        if sim > 0.0:
            neighbours.append((sim, pchembl))
    if not neighbours:
        return None
    neighbours.sort(key=lambda pair: pair[0], reverse=True)
    top = neighbours[:k]
    weight = sum(sim for sim, _ in top)
    # A similarity-weighted mean of in-band neighbour values is in-band by
    # construction, so the prediction can never extrapolate outside the
    # target's observed potency band.
    return sum(sim * pchembl for sim, pchembl in top) / weight


def predict_kd_nm(target: str, smiles: str, k: int = K_DEFAULT) -> float | None:
    """Predicted KD (nM) of ``smiles`` against ``target`` (IC50-derived)."""
    pchembl = predict_pchembl(target, smiles, k)
    if pchembl is None:
        return None
    return float(10.0 ** (9.0 - pchembl))


def top_tanimoto(target: str, smiles: str) -> float | None:
    """Best Tanimoto match of ``smiles`` in ``target``'s snapshot rows.

    The evidence metric for the support gate: the maximum similarity against
    the whole vendored row set, ``None`` when the target/SMILES is absent or
    the query is orthogonal (no positive match).  A low top match means the
    kNN prediction would interpolate across empty chemistry — refuse it.
    """
    rows = load_snapshot().targets.get(target)
    if rows is None or not rows.smiles:
        return None
    fp = _fingerprint(smiles)
    if fp is None:
        return None
    best = 0.0
    for smi in rows.smiles:
        other = _fingerprint(smi)
        if other is None:
            continue
        sim = _tanimoto(fp, other)
        if sim > best:
            best = sim
    return best or None


def cns_ic50_nm(
    smiles: str,
    k: int = K_DEFAULT,
    max_nm: float = _CNS_DEFAULT_IC50_NM,
    min_tanimoto: float = MIN_NEIGHBOR_TANIMOTO,
) -> float | None:
    """Structure-derived CNS potency anchor (nM) for the grading seam (L12).

    Worst case over the CNS-liability set — the most potent predicted binding —
    capped at ``max_nm`` so a molecule whose predictions are all weak is bounded
    by (not more potent than) the disclosed conservative default.  Only targets
    with chemotype support for the query contribute (``top_tanimoto >=
    min_tanimoto``); when *no* CNS target has a supported prediction the anchor
    is ``None`` and the seam falls back to the disclosed class prior instead of
    a structure-derived value with no evidence behind it.  ``None`` also for an
    unparseable SMILES.
    """
    predicted: list[float] = []
    for target in CNS_TARGETS:
        best = top_tanimoto(target, smiles)
        if best is None or best < min_tanimoto:
            continue
        kd = predict_kd_nm(target, smiles, k)
        if kd is not None:
            predicted.append(kd)
    if not predicted:
        return None
    return min(max_nm, min(predicted))


def no_public_data_sites(resolved: Sequence[Target]) -> tuple[str, ...]:
    """Panel sites that remain on class priors after head + kNN resolution."""
    return tuple(sorted(t.name for t in resolved if t.low_confidence))


def resolve_offtarget_panel(
    panel: Sequence[Target],
    smiles: str | None = None,
    admet: object = None,
    k: int = K_DEFAULT,
    min_tanimoto: float = MIN_NEIGHBOR_TANIMOTO,
) -> tuple[tuple[Target, ...], tuple[str, ...]]:
    """Resolve the safety panel site-by-site (path B heads, then path A kNN).

    Every site whose ADMET-AI head measures the same interaction is re-scored
    first (``resolve_admet_panel``); then every mapped site with *chemotype
    support* for the query (``top_tanimoto >= min_tanimoto``) is re-bound to
    its kNN-predicted KD (``low_confidence=False``).  Sites without support
    keep their head re-score or class prior and stay disclosed in
    ``no_public_data_sites``.
    Returns ``(resolved_panel, no_public_data_sites)``.
    """
    out = list(resolve_admet_panel(panel, admet))
    if smiles is not None:
        for i, t in enumerate(out):
            target = SITE_TO_TARGET.get(t.name)
            if target is None:
                continue
            best = top_tanimoto(target, smiles)
            if best is None or best < min_tanimoto:
                continue
            kd = predict_kd_nm(target, smiles, k)
            if kd is None:
                continue
            out[i] = Target(
                name=t.name,
                kd_nm=kd,
                kon_nm_h=t.kon_nm_h,
                r0_nm=t.r0_nm,
                rho_h=t.rho_h,
                kint_h=t.kint_h,
                low_confidence=False,
                reference=(f"ChEMBL kNN {target} snapshot (fingerprint resolver, doc/12 L15)"),
            )
    resolved = tuple(out)
    return resolved, no_public_data_sites(resolved)


__all__ = [
    "CNS_TARGETS",
    "MIN_NEIGHBOR_TANIMOTO",
    "SITE_TO_TARGET",
    "Snapshot",
    "cns_ic50_nm",
    "load_snapshot",
    "no_public_data_sites",
    "predict_kd_nm",
    "predict_pchembl",
    "resolve_offtarget_panel",
    "top_tanimoto",
]
