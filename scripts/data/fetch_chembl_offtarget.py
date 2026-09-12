"""Regenerate data/chembl/offtarget_snapshot.json from the live ChEMBL API.

Purpose
-------
The general off-target DTI resolver (``drugos.target.dti``) needs a bounded,
vendored, reproducibly-pulled bioactivity snapshot covering the safety panel
sites *and* the CNS-liability set the resolver is calibrated on (doc/12 L13,
L15, L12).  Call this script any time a re-issue should be pulled::

    /opt/anaconda3/envs/drug_os/bin/python scripts/data/fetch_chembl_offtarget.py

For each curated ChEMBL target it queries measured (``standard_relation='='``)
IC50/Ki activities expressed in nM with a resolved molecule structure, keeps
the canonical-SMILES central pChemBL per unique molecule, and stores a
potency-stratified, deterministic, size-capped slice (so the vendored file is
human-reviewable and never a bulk leak).  Targets without public data are
simply omitted — absence is recorded, not guessed.

The output is deliberately bounded (<= ``MAX_ROWS_PER_TARGET`` molecules per
target) and potency-stratified: pChemBL is binned in half-log steps and at most
``ROWS_PER_BIN`` molecules per bin survive, preserving the potency spread a
nearest-neighbour regressor needs instead of dumping one crowded section.
"""

from __future__ import annotations

import json
import math
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
UA = "drugos-data/2026.9"

MAX_ROWS_PER_TARGET = 200
ROWS_PER_BIN = 25

# Curated ChEMBL target set: panel sites + CNS-liability sites (doc/05 2.2,
# doc/10 P2).  KEYS are consumed by drugos/target/dti.py site->target mapping.
TARGETS: dict[str, str] = {
    "KCNH2": "CHEMBL240",  # hERG (Kv11.1)
    "CYP2D6": "CHEMBL289",
    "CYP3A4": "CHEMBL340",
    "CYP2C9": "CHEMBL3397",
    "ABCB11": "CHEMBL335",  # BSEP / canalicular efflux
    "SLCO1B1": "CHEMBL1697668",  # OATP1B1
    "ABCB1": "CHEMBL4302",  # P-glycoprotein / MDR1
    "NR3C1": "CHEMBL2034",  # glucocorticoid receptor
    "ESR1": "CHEMBL206",  # estrogen receptor
    "AR": "CHEMBL1871",  # androgen receptor
    "NADH_DEHYDROGENASE": "CHEMBL4552",  # mitochondrial complex I (rotenone-site)
    "DRD2": "CHEMBL217",  # CNS liability: dopamine D2
    "CHRM1": "CHEMBL226",  # CNS: muscarinic M1
    "ADRA1A": "CHEMBL229",  # CNS: alpha-1A adrenergic
    "OPRM1": "CHEMBL233",  # CNS: mu-opioid
    "HRH1": "CHEMBL231",  # CNS: histamine H1
}

OUT = Path(__file__).resolve().parents[2] / "data" / "chembl" / "offtarget_snapshot.json"

_TYPES = ("IC50", "Ki")


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def _n_m_from_row(act: dict) -> float | None:
    sv, su = act.get("standard_value"), act.get("standard_units")
    if sv is None or su is None:
        return None
    try:
        value = float(sv)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value) or value <= 0:
        return None
    if su.lower() in ("nm", "nm/l", "nm/u"):
        return value
    if su.lower() in ("um", "µm", "microm", "um/l"):
        return value * 1e3
    if su.lower() in ("mm", "mmol/l") and value * 1e6 >= 1.0:
        return value * 1e6
    return None


def _fetch_target(chembl_id: str) -> list[dict[str, object]]:
    url = (
        f"{API}?target_chembl_id={chembl_id}&pchembl_value__isnull=false"
        f"&_fields=canonical_smiles,standard_value,standard_units,standard_relation,"
        f"standard_type,pchembl_value&limit=1000"
    )
    payload = _get(url)
    per_molecule: dict[str, list[float]] = {}
    for act in payload.get("activities", []):
        if act.get("standard_relation") not in (None, "="):
            continue
        if act.get("standard_type") not in _TYPES:
            continue
        smiles = act.get("canonical_smiles")
        if not smiles:
            continue
        nM = _n_m_from_row(act)
        if nM is None:
            continue
        per_molecule.setdefault(smiles, []).append(nM)
    rows: list[dict[str, object]] = []
    for smiles, vals in per_molecule.items():
        geomean = math.exp(sum(math.log(v) for v in vals) / len(vals))
        rows.append(
            {
                "smiles": smiles,
                "ic50_nm": round(geomean, 3),
                "pchembl": round(9.0 - math.log10(geomean), 3),  # -log10(IC50 / M)
                "n_rows": len(vals),
            }
        )
    rows.sort(key=lambda r: float(r["pchembl"]), reverse=True)
    # Half-log potency-stratified cap (deterministic, spread-preserving).
    binned: dict[int, list[dict[str, object]]] = {}
    for row in rows:
        bucket = int(math.floor(float(row["pchembl"]) * 2.0))
        binned.setdefault(bucket, []).append(row)
    kept: list[dict[str, object]] = []
    for bucket in sorted(binned):
        for row in binned[bucket][:ROWS_PER_BIN]:
            kept.append(row)
        if len(kept) >= MAX_ROWS_PER_TARGET:
            break
    return kept[:MAX_ROWS_PER_TARGET]


def main() -> None:
    targets: dict[str, object] = {}
    for name, chembl_id in TARGETS.items():
        try:
            rows = _fetch_target(chembl_id)
        except Exception as exc:  # network / API failure: record, don't invent
            print(f"[warn] {name} ({chembl_id}) fetch failed: {exc!r}")
            continue
        pc = [float(r["pchembl"]) for r in rows]
        targets[name] = {
            "target_chembl_id": chembl_id,
            "n_molecules": len(rows),
            "pchembl_min": round(min(pc), 3) if pc else None,
            "pchembl_max": round(max(pc), 3) if pc else None,
            "rows": rows,
        }
        print(f"{name:20s} {chembl_id} {len(rows)} molecules")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        json.dumps(
            {
                "about": (
                    "ChEMBL bioactivity snapshot for the general DTI resolver "
                    "(drugos/target/dti.py). Measured (relation '=') IC50/Ki rows "
                    "with resolved structures, geomean-central per canonical SMILES, "
                    "potency-stratified and size-capped for review. Absent targets "
                    "or molecules are absent on purpose; nothing is guessed."
                ),
                "retrieved_utc": datetime.now(UTC).isoformat(timespec="seconds"),
                "targets": targets,
            },
            indent=1,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {OUT}  ({len(targets)} targets)")


if __name__ == "__main__":
    main()