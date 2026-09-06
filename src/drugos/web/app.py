"""DrugOS web playground (deep-purple dark theme).

A small Flask application exposing two endpoints on top of the pipeline:

- ``GET /``: the single-page playground UI.
- ``GET /api/benchmarks``: benchmark names for the form.
- ``POST /api/run``: run the pipeline for a payload and return the doc/03
  contract (JSON), which the browser renders with canvas charts.

"""

from __future__ import annotations

import os
from typing import Any

from flask import Flask, jsonify, render_template_string, request

from drugos.cli import spec_from_cli
from drugos.pipeline import benchmark_names, run_pipeline
from drugos.robustness import (
    EnsembleConfig,
    PopulationConfig,
    SensitivityConfig,
    run_population,
    run_uncertainty,
)
from drugos.robustness.sensitivity import local_sensitivity, run_sobol_sensitivity


def create_app(host: str = "127.0.0.1", port: int = 8080) -> Flask:
    del host, port
    app = Flask(__name__)
    template_dir = os.path.join(os.path.dirname(__file__), "templates")
    with open(os.path.join(template_dir, "index.html"), encoding="utf-8") as fh:
        index_html = fh.read()

    @app.get("/")
    def index() -> str:
        return render_template_string(index_html)

    @app.get("/api/benchmarks")
    def benchmarks() -> Any:
        return jsonify({"benchmarks": benchmark_names()})

    @app.post("/api/run")
    def api_run() -> Any:
        payload = request.get_json(force=True, silent=True)
        payload = payload if isinstance(payload, dict) else {}
        try:
            spec = spec_from_cli(
                benchmark=payload.get("benchmark"),
                smiles=payload.get("smiles"),
                dose=_float(payload.get("dose")),
                route=payload.get("route") or None,
                sex=payload.get("sex") or "male",
                age=_float(payload.get("age")) or 40.0,
                height=_float(payload.get("height")) or 170.0,
                weight=_float(payload.get("weight")) or 70.0,
                no_pathway=bool(payload.get("no_pathway")),
            )
            result = run_pipeline(spec)
        except (ValueError, RuntimeError) as exc:
            return jsonify({"error": str(exc)}), 400
        contract = result.to_contract()
        if payload.get("with_uncertainty"):
            seed = _int(payload.get("seed"), 7)
            ensemble = run_uncertainty(
                spec, EnsembleConfig(n_runs=_int(payload.get("n_unc"), 8), seed=seed)
            )
            contract["uncertainty"] = ensemble.to_dict()
        if payload.get("with_population"):
            seed = _int(payload.get("seed"), 7)
            population = run_population(
                spec, PopulationConfig(n_individuals=_int(payload.get("n_pop"), 8), seed=seed + 1)
            )
            contract["population"] = population.to_dict()
        if payload.get("with_sensitivity"):
            seed = _int(payload.get("seed"), 7)
            config = SensitivityConfig(seed=seed + 2, sobol_n=_int(payload.get("n_sobol"), 4))
            local = local_sensitivity(spec, "overall_risk", config)
            contract["sensitivity"] = {
                "keys": list(config.keys),
                "local_overall_risk": local.values,
                "drivers_overall_risk": [list(t) for t in local.drivers],
                "sobol": run_sobol_sensitivity(spec, config).to_dict(),
            }
        return jsonify(contract)

    return app


def _float(value: Any) -> float | None:
    """Parse a float field; ``None`` passes through, malformed input raises.

    The caller (``/api/run``) catches ``ValueError`` and reports a 400 — an
    explicit error. There is no silent default for bad user input.
    """
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"expected number, got boolean {value!r}")
    if not isinstance(value, (int, float, str)):
        raise ValueError(f"expected number, got {type(value).__name__}")
    return float(value)


def _int(value: Any, default: int) -> int:
    """Parse a positive int field; ``None`` -> ``default``, malformed input raises.

    A client-supplied bad value must fail loudly (400), never silently
    substitute a default.
    """
    if value is None:
        return default
    if isinstance(value, bool):
        raise ValueError(f"expected integer, got boolean {value!r}")
    parsed = int(value)
    if parsed < 1:
        raise ValueError(f"expected positive integer, got {parsed}")
    return parsed


if __name__ == "__main__":
    create_app().run()
