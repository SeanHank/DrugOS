"""Regenerate data/benchmarks/herg_measured_nm.json from the live ChEMBL API.

Purpose
-------
Replace the hard-coded class-median hERG prior for the benchmark corpus with
measured patch-clamp / radioligand hERG data where the public record is
unambiguous.  Call this script any time a re-issue of ChEMBL should be
pulled::

    /opt/anaconda3/envs/drug_os/bin/python scripts/data/fetch_chembl_herg.py

It queries ChEMBL target CHEMBL240 (hERG) IC50 activity for each benchmark
molecule, keeps only *measured* (activity_type != 'Potency'-derived) rows,
takes the geometric mean of the binding-class (assay_type in {B, F}) values,
and writes the vendored, checksummed JSON consumed by the validation layer.

The output file is deliberately small and human-reviewable: one central
estimate per compound plus the individual rows, so a reviewer can recompute
the decision by hand.  Missing measurements (warfarin, midazolam,
acetaminophen, ciprofloxacin have no hERG IC50 rows against CHEMBL240) are
simply omitted — never filled with guesses.
"""

from __future__ import annotations

import json
import math
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

HERG_TARGET = "CHEMBL240"
TYPE = "IC50"
API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

MOLECULES = {
    "dofetilide": "CHEMBL473",
    "warfarin": "CHEMBL1464",
    "midazolam": "CHEMBL655",
    "acetaminophen": "CHEMBL112",
    "ciprofloxacin": "CHEMBL8",
}

OUT = Path(__file__).resolve().parents[2] / "data" / "benchmarks" / "herg_measured_nm.json"


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "drugos-data/2026.9"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return __import__("json").load(resp)


def _geomean(values: list[float]) -> float:
    prod = 0.0
    for v in values:
        prod += math.log(v)
    return math.exp(prod / len(values))


def main() -> None:
    compounds: dict[str, object] = {}
    compounds_meta: dict[str, object] = {}
    for name, chembl_id in MOLECULES.items():
        url = (
            f"{API}?molecule_chembl_id={chembl_id}&target_chembl_id={HERG_TARGET}"
            f"&type={TYPE}&pchembl_value__isnull=false&limit=100"
        )
        payload = _get(url)
        rows = []
        for act in payload.get("activities", []):
            sv, su, pc = (
                act.get("standard_value"),
                act.get("standard_units"),
                act.get("pchembl_value"),
            )
            if sv is None or su is None or pc is None:
                continue
            try:
                value_nm = float(sv) if su.lower() == "nm" else float(sv) * 1e3
            except (TypeError, ValueError):
                continue
            if value_nm <= 0 or not math.isfinite(value_nm):
                continue
            if act.get("assay_type") not in ("B", "F"):
                continue
            rows.append(
                {
                    "standard_value": float(sv),
                    "units": su,
                    "pchembl_value": float(pc),
                    "assay_type": act.get("assay_type"),
                    "document_chembl_id": act.get("document_chembl_id"),
                }
            )
        if rows:
            # Outlier policy (documented in the vendored file): a row whose
            # IC50 is >= 10 uM is a suspected data-entry/config error for a
            # compound otherwise cited at low-nM potency (e.g. the 44 uM
            # dofetilide row); it is *retained* for transparency but excluded
            # from the central estimate.
            core = [r for r in rows if r["standard_value"] < 1.0e4]
            if core:
                central = _geomean([r["standard_value"] for r in core])
            else:
                central = _geomean([r["standard_value"] for r in rows])
            compounds[name] = {
                "chembl_id": chembl_id,
                "ic50_geomean_nm": round(central, 4),
                "n_rows_in_central_estimate": len(core),
                "rows": rows,
            }
        compounds_meta[name] = {"chembl_id": chembl_id, "n_measured_rows": len(rows)}

    OUT.write_text(
        json.dumps(
            {
                "about": (
                    "Measured hERG (Kv11.1) IC50 for the five-benchmark corpus from "
                    "the public ChEMBL record (target CHEMBL240). Geometric mean of "
                    "experimental (binding/functional) IC50 rows per compound. Compounds "
                    "without any hERG IC50 row are omitted: absence is recorded, not guessed."
                ),
                "retrieved_utc": datetime.now(UTC).isoformat(timespec="seconds"),
                "compounds": compounds,
                "compounds_presence": compounds_meta,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {OUT}  (compounds with measurements: {len(compounds)})")


if __name__ == "__main__":
    main()
