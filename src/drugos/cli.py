"""Command-line interface for DrugOS (``python -m drugos`` / ``drugos``).

Subcommands:

- ``run``: execute the full pipeline for a benchmark compound (or a SMILES
  structure with predicted physiology) and emit the JSON / markdown / HTML
  report to stdout or a directory.
- ``benchmarks``: list the built-in benchmark compounds.
- ``serve``: start the web playground (Flask, deep-purple dark theme).
- ``version``: print the package version.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from drugos.inputs.models import HumanProfile, Sex
from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pipeline import (
    RunSpec,
    benchmark_data,
    benchmark_names,
    run_pipeline,
    spec_from_benchmark_data,
)
from drugos.report import render_html, render_json, render_markdown, write_report
from drugos.version import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="drugos",
        description=(
            "DrugOS: multiscale mechanism-based model of drug response in the human body "
            "(drug -> concentration -> target -> pathway -> organ -> phenotype)."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser("run", help="run the full pipeline and emit a report")
    run_p.add_argument("--benchmark", default="acetaminophen", help="benchmark compound name")
    run_p.add_argument("--smiles", help="SMILES structure (RDKit descriptors + ADMET-AI)")
    run_p.add_argument("--dose", type=float, help="total dose in mg (overrides benchmark)")
    run_p.add_argument(
        "--route",
        choices=["iv_bolus", "iv_infusion", "oral", "subcutaneous", "intramuscular", "transdermal"],
        default="oral",
        help="route of administration",
    )
    run_p.add_argument(
        "--sc-im-ka",
        type=float,
        help="SC/IM first-order depot absorption rate (1/h); SC/IM only",
    )
    run_p.add_argument("--sex", choices=["male", "female"], default="male")
    run_p.add_argument("--age", type=float, default=40.0)
    run_p.add_argument("--height", type=float, default=170.0)
    run_p.add_argument("--weight", type=float, default=70.0)
    run_p.add_argument(
        "--format",
        choices=["json", "markdown", "html"],
        default="markdown",
        help="stdout format",
    )
    run_p.add_argument("--out", help="directory to write report.json/md/html (creates it)")
    run_p.add_argument(
        "--no-pathway", action="store_true", help="skip the Stage-3 pathway simulation"
    )

    sub.add_parser("benchmarks", help="list benchmark compounds")
    sub.add_parser("version", help="print the package version")
    sv = sub.add_parser("serve", help="start the web playground")
    sv.add_argument("--host", default="127.0.0.1")
    sv.add_argument("--port", type=int, default=8080)

    st = sub.add_parser("study", help="Phase-6 study (uncertainty + population + sensitivity)")
    st.add_argument("--benchmark", default="acetaminophen", help="benchmark compound name")
    st.add_argument("--dose", type=float, help="total dose in mg (overrides benchmark)")
    st.add_argument(
        "--route",
        choices=["oral", "iv_bolus", "iv_infusion"],
        help="route of administration",
    )
    st.add_argument("--sex", choices=["male", "female"], default="male")
    st.add_argument("--age", type=float, default=40.0)
    st.add_argument("--height", type=float, default=170.0)
    st.add_argument("--weight", type=float, default=70.0)
    st.add_argument("--no-pathway", action="store_true", help="skip the Stage-3 pathway simulation")
    st.add_argument("--n-unc", type=int, default=17, help="uncertainty ensemble size")
    st.add_argument("--n-pop", type=int, default=16, help="virtual population size")
    st.add_argument("--n-sobol", type=int, default=8, help="Sobol base-design size")
    st.add_argument("--seed", type=int, default=7, help="shared random seed")
    st.add_argument("--format", choices=["json", "markdown"], default="json", help="stdout format")
    st.add_argument("--out", help="directory to write study.json/md (creates it)")
    return parser


def _benchmark_lookup(name: str) -> Any:
    """Return the canonical benchmark adapter from the self-contained registry."""
    return benchmark_data(name)


def _smiles_spec(args: _RunInputs) -> RunSpec:
    """Spec for a custom SMILES using structure + ADMET-AI predicted PK."""
    from drugos.inputs.parse_structure import parse_structure
    from drugos.pk.admet import predict_admet
    from drugos.pk.physiology import build_human, glom_filtration_clearance

    assert args.smiles is not None
    mol = parse_structure(args.smiles, name="custom")
    mw = float(mol.mw or 0.0)
    if mw <= 0:
        raise ValueError("could not compute molecular weight")
    pred = predict_admet(args.smiles)
    pred_obj = pred[0] if isinstance(pred, list) else pred
    fup = pred_obj.fup_plasma
    if fup is None or fup <= 0:
        raise ValueError("ADMET-AI did not return a usable fup")
    cl_int = pred_obj.cl_int_hep_ml_min_kg
    if cl_int is None:
        raise ValueError("ADMET-AI did not return hepatic intrinsic clearance")
    weight_kg = args.weight
    cl_hep = cl_int * 60.0 / 1000.0 * weight_kg  # mL/min/kg -> L/h (gross scaling)
    profile = HumanProfile(
        sex=Sex(args.sex), age_y=args.age, height_cm=args.height, weight_kg=weight_kg
    )
    physiology = build_human(profile)
    cl_renal = glom_filtration_clearance(physiology.gfr_l_min * 1000.0, fup)
    return RunSpec(
        name=mol.name or "custom",
        molecule=mol,
        profile=profile,
        dose_plan=build_dose_plan(args.route or "oral", args.dose or 10.0),
        cl_hep_l_h=cl_hep,
        cl_renal_l_h=cl_renal,
        mw=mw,
        fup=fup,
        bp=1.0,
        admet=pred_obj,
        include_pathway=not args.no_pathway,
    )


def spec_from_cli(
    *,
    benchmark: str | None = None,
    smiles: str | None = None,
    dose: float | None = None,
    route: str | None = None,
    sex: str = "male",
    age: float = 40.0,
    height: float = 170.0,
    weight: float = 70.0,
    no_pathway: bool = False,
    sc_im_ka: float | None = None,
) -> RunSpec:
    """Build a ``RunSpec`` from playground/CLI string inputs."""
    if smiles:
        spec = _smiles_spec(
            _RunInputs(
                smiles=smiles,
                dose=dose,
                route=route,
                sex=sex,
                age=age,
                height=height,
                weight=weight,
                no_pathway=no_pathway,
                sc_im_ka=sc_im_ka,
            )
        )
    else:
        bench = _benchmark_lookup(benchmark or "acetaminophen")
        if bench is None:
            raise ValueError(f"unknown benchmark '{benchmark}'")
        spec = spec_from_benchmark_data(bench, dose_override_mg=dose)
        spec.profile = HumanProfile(sex=Sex(sex), age_y=age, height_cm=height, weight_kg=weight)
        if route is not None:
            spec.dose_plan = build_dose_plan(route, spec.dose_plan.total_dose_mg)
    if no_pathway:
        spec.include_pathway = False
    if sc_im_ka is not None:
        spec.sc_im_ka_per_h = sc_im_ka
    return spec


@dataclass(frozen=True, slots=True)
class _RunInputs:
    """Typed projection of the CLI/playground run inputs (defaults overridable)."""

    smiles: str | None = None
    dose: float | None = None
    route: str | None = None
    sex: str = "male"
    age: float = 40.0
    height: float = 170.0
    weight: float = 70.0
    no_pathway: bool = False
    sc_im_ka: float | None = None


def cmd_run(args: argparse.Namespace) -> int:
    try:
        spec = spec_from_cli(
            benchmark=args.benchmark,
            smiles=args.smiles,
            dose=args.dose,
            route=args.route,
            sex=args.sex,
            age=args.age,
            height=args.height,
            weight=args.weight,
            no_pathway=args.no_pathway,
            sc_im_ka=args.sc_im_ka,
        )
        result = run_pipeline(spec)
    except (ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.out:
        paths = write_report_directory(args.out, result)
        for kind, path in paths.items():
            print(f"wrote {kind} report -> {path}")
    else:
        text = {
            "json": lambda: render_json(result),
            "markdown": lambda: render_markdown(result),
            "html": lambda: render_html(result),
        }[args.format]()
        print(text)
    return 0


def write_report_directory(directory: str, result: Any) -> dict[str, str]:
    """Write the three report artifacts into ``directory`` (created if needed)."""
    import os

    os.makedirs(directory, exist_ok=True)
    return write_report(result, directory)


def cmd_benchmarks(args: argparse.Namespace | None = None) -> int:
    for name in benchmark_names():
        print(name)
    return 0


def cmd_version(args: argparse.Namespace | None = None) -> int:
    print(f"drugos {__version__}")
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    module = importlib.import_module("drugos.web.app")
    app = module.create_app(host=args.host, port=args.port)
    print(f"DrugOS playground running at http://{args.host}:{args.port} (Ctrl-C to stop)")
    app.run(host=args.host, port=args.port, debug=False, use_reloader=False)
    return 0


def _study_objects(spec: RunSpec, args: argparse.Namespace) -> tuple[Any, Any, Any, Any]:
    """Compute the D21-D23 study artifacts for ``spec``."""
    from drugos.robustness import (
        EnsembleConfig,
        PopulationConfig,
        SensitivityConfig,
        run_population,
        run_sobol_sensitivity,
        run_uncertainty,
    )
    from drugos.robustness.sensitivity import local_sensitivity

    seed = args.seed
    ensemble = run_uncertainty(spec, EnsembleConfig(n_runs=args.n_unc, seed=seed))
    population = run_population(spec, PopulationConfig(n_individuals=args.n_pop, seed=seed + 1))
    config = SensitivityConfig(seed=seed + 2, sobol_n=args.n_sobol)
    local = local_sensitivity(spec, "overall_risk", config)
    sobol = run_sobol_sensitivity(spec, config)
    return ensemble, population, local, sobol


def _study_contract(spec: RunSpec, args: argparse.Namespace) -> dict[str, object]:
    ensemble, population, local, sobol = _study_objects(spec, args)
    sens = {
        "keys": list(local.config.keys),
        "local_overall_risk": local.values,
        "drivers_overall_risk": local.drivers,
        "sobol": sobol.to_dict(),
    }
    return {
        "compound": spec.name,
        "dose_mg": spec.dose_plan.total_dose_mg,
        "uncertainty": ensemble.to_dict(),
        "population": population.to_dict(),
        "sensitivity": sens,
    }


def _study_markdown(spec: RunSpec, args: argparse.Namespace) -> str:
    ensemble, population, local, sobol = _study_objects(spec, args)
    del sobol, args
    unc = ensemble.to_dict()
    pop = population.to_dict()
    counts = dict(unc["verdict_counts"])
    incidence = dict(pop["incidence_grade_ge_1"])
    quantiles = dict(pop["risk_q5_q50_q95"])
    ranked = sorted(local.drivers, key=lambda kv: abs(kv[1]), reverse=True)
    lines = [
        f"# DrugOS robustness study — {spec.name}",
        "",
        f"- Dose {spec.dose_plan.total_dose_mg:g} mg",
        "",
        "## D21 Uncertainty (parameter ensemble)",
        "",
        f"- {unc['n_runs']} ensemble runs, CV={unc['cv']:g}, seed={unc['seed']}.",
        "- Verdict distribution: " + ", ".join(f"{k} x{v}" for k, v in counts.items()),
        "",
        "## D22 Virtual population",
        "",
        f"- {pop['n_individuals']} virtual patients, age {pop['age_y_range'][0]}"
        f"-{pop['age_y_range'][1]} y.",
        "- Endpoint incidence (grade >= 1): " + ", ".join(f"{k}={v}" for k, v in incidence.items()),
        "- Endpoint risk q5/q50/q95: " + ", ".join(f"{k}={v}" for k, v in quantiles.items()),
        "",
        "## D23 Sensitivity",
        "",
        "- Top local drivers of overall risk: "
        + ", ".join(f"{k} ({v:+.3f})" for k, v in ranked[:4]),
        "",
        "> Research-grade model output; not for clinical decision-making (doc/08).",
    ]
    return "\n".join(lines)


def cmd_study(args: argparse.Namespace) -> int:
    try:
        spec = spec_from_cli(
            benchmark=args.benchmark,
            dose=args.dose,
            route=args.route,
            sex=args.sex,
            age=args.age,
            height=args.height,
            weight=args.weight,
            no_pathway=args.no_pathway,
        )
        study = _study_contract(spec, args)
    except (ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.out:
        import os

        os.makedirs(args.out, exist_ok=True)
        paths = {
            "json": f"{args.out}/study.json",
            "markdown": f"{args.out}/study.md",
        }
        with open(paths["json"], "w", encoding="utf-8") as fh:
            json.dump(study, fh, indent=2)
        with open(paths["markdown"], "w", encoding="utf-8") as fh:
            fh.write(_study_markdown(spec, args))
        for kind, path in paths.items():
            print(f"wrote {kind} study -> {path}")
    else:
        text = json.dumps(study, indent=2) if args.format == "json" else _study_markdown(spec, args)
        print(text)
    return 0


def _dispatch(
    command: str,
) -> Callable[[argparse.Namespace], int] | None:
    """Map a subcommand name to its handler (for main and tests)."""
    handlers: dict[str, Callable[[argparse.Namespace], int]] = {
        "run": cmd_run,
        "benchmarks": cmd_benchmarks,
        "serve": cmd_serve,
        "study": cmd_study,
        "version": cmd_version,
    }
    return handlers.get(command)


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handler = _dispatch(args.command)
    if handler is None:
        return 0
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser",
    "cmd_benchmarks",
    "cmd_run",
    "cmd_serve",
    "cmd_study",
    "cmd_version",
    "main",
    "spec_from_cli",
    "write_report_directory",
]
