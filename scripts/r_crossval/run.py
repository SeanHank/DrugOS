#!/usr/bin/env python3
"""Optional R cross-validation loader (off-gate).

Runs a DrugOS benchmark in Python, dumps the simulated plasma curve + published
bands (read straight from ``data/benchmarks/published_pk.json``) to TSV for the
R harness ``run.r``, and optionally executes Rscript when an R interpreter is
available.  See ``scripts/r_crossval/README.md``.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from drugos.cli import spec_from_cli
from drugos.pipeline import run_pipeline

_DATA = Path(__file__).resolve().parents[2] / "data" / "benchmarks" / "published_pk.json"
_BAND_KEYS = ("cl_plasma_l_h", "t_half_h")


def _bands(benchmark: str) -> dict[str, list[float]]:
    corpus = json.loads(_DATA.read_text())
    if benchmark not in corpus["compounds"]:
        raise ValueError(f"benchmark '{benchmark}' has no published bands")
    published = corpus["compounds"][benchmark]["published"]
    bands = {k: list(published[k]) for k in _BAND_KEYS if k in published}
    if not bands:
        raise ValueError(f"benchmark '{benchmark}' has none of {_BAND_KEYS}")
    return bands


def _dump(benchmark: str, dose: float) -> tuple[Path, Path, Path]:
    spec = spec_from_cli(benchmark=benchmark, dose=dose, route="oral")
    result = run_pipeline(spec)
    contract = result.to_contract()
    pk = contract["pk"]

    out = Path(tempfile.mkdtemp(prefix="r_crossval_")) / benchmark
    out.mkdir(parents=True, exist_ok=True)

    sim_f = out / "sim.tsv"
    with sim_f.open("w") as fh:
        fh.write("time_h\tplasma_total_mg_l\n")
        for t, c in zip(pk["t_h"], pk["plasma_total"], strict=True):
            fh.write(f"{t}\t{c}\n")

    bands_f = out / "bands.tsv"
    bands = _bands(benchmark)
    with bands_f.open("w") as fh:
        fh.write("quantity\tlo\thi\n")
        for k, (lo, hi) in bands.items():
            fh.write(f"{k}\t{lo}\t{hi}\n")

    meta_f = out / "meta.json"
    meta_f.write_text(
        json.dumps(
            {
                "benchmark": benchmark,
                "dose_mg": dose,
                "route": contract["manifest"]["route"],
                "bands": bands,
            },
            indent=2,
        )
    )
    return sim_f, bands_f, meta_f


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Dump a DrugOS run for R cross-validation.")
    parser.add_argument(
        "--benchmark",
        default="warfarin",
        help="benchmark name with published bands in published_pk.json",
    )
    parser.add_argument("--dose", type=float, default=5.0)
    parser.add_argument("--with-r", action="store_true", help="run Rscript after dumping")
    args = parser.parse_args(argv)

    sim_f, bands_f, _ = _dump(args.benchmark, args.dose)
    print(f"sim.tsv  -> {sim_f}")
    print(f"bands.tsv-> {bands_f}")

    if not args.with_r:
        print("Run the R verdict with:")
        print(
            "  Rscript scripts/r_crossval/run.r"
            f" --sim {sim_f} --bands {bands_f} --out {sim_f.parent}/report.md"
        )
        return 0

    rscript = shutil.which("Rscript")
    if rscript is None:
        print("Rscript not found on PATH; skipping R execution (see scripts/r_crossval/README.md).")
        return 1
    report = sim_f.parent / "report.md"
    subprocess.check_call(
        [
            rscript,
            os.path.join(os.path.dirname(__file__), "run.r"),
            "--sim",
            str(sim_f),
            "--bands",
            str(bands_f),
            "--out",
            str(report),
        ]
    )
    print(f"R report -> {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
