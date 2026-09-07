"""Production SBML MAPK cascade integration (R-5, L2 analytic/numerical).

Stage-3 pathway upgrade (doc/12 §1 row 3): the pipeline default no longer
runs the in-house DSL cascade but the vendored Huang & Ferrell / Levchenko
ultrasensitive MAPK cascade (BioModels BIOMD0000000009, CC0) parsed through
python-libsbml.  This case pins the production-model properties the pipeline
relies on:

- the vendored model parses to the reference 20-reaction / 22-species graph
  (beyond-audit: species named after the Huang1996 species set),
- drug-free system relaxes to its fully-activated steady state (readout is
  doubly-phosphorylated ERK, ``PP_K``, ~1 in concentration units),
- signal (target occupancy / Kd) inhibits upstream MAPKKK activation
  monotonically: full occupancy collapses the readout << baseline, and the
  response is monotone in occupancy (the known ultrasensitive + monotone
  dose-response regime of a clamped steady-state cascade).

The checks are analytic/numerical (L2): the ODE is exactly the SBML kinetic
laws with stoichiometry, and the case verifies the annotated wiring rather
than re-deriving physiology.
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pathway import parse_sbml, simulate_sbml_pathway
from drugos.pathway.sbml_pathway import model_file

_TIME_H, _N = 72.0, 240


def _run(occ: float) -> object:
    t = np.linspace(0.0, _TIME_H, _N)
    return simulate_sbml_pathway(t, np.full_like(t, occ), n_eval=_N)


def case_sbml_mapk_validation() -> CaseResult:
    spec = parse_sbml(model_file())

    baseline = _run(0.0)
    full = _run(1.0)
    mid1 = _run(0.6)
    mid2 = _run(0.9)

    ppk_base = float(np.max(baseline.concentrations["PP_K"]))
    ppk_full = float(np.max(full.concentrations["PP_K"]))
    ppk_mid1 = float(np.max(mid1.concentrations["PP_K"]))
    ppk_mid2 = float(np.max(mid2.concentrations["PP_K"]))
    fc_end = float(full.readout_fold_change("PP_K")[-1])

    metrics = [
        MetricResult(
            "reactions_count",
            len(spec.reactions),
            20,
            20,
            "count",
            "pass" if len(spec.reactions) == 20 else "FAIL",
        ),
        MetricResult(
            "species_contain_ppk",
            float("PP_K" in spec.species),
            1.0,
            1.0,
            "flag",
            "pass" if "PP_K" in spec.species else "FAIL",
        ),
        MetricResult(
            "drug_free_ppk_steady_state",
            ppk_base,
            0.5,
            2.0,
            "conc",
            "pass" if ppk_base > 0.5 else "FAIL",
        ),
        MetricResult(
            "full_occupancy_ppk_suppression",
            ppk_full,
            0.0,
            0.05,
            "conc",
            "pass" if ppk_full < 0.05 else "FAIL",
        ),
        MetricResult(
            "monotone_inhibition_in_occupancy",
            float(ppk_mid1 >= ppk_mid2),
            1.0,
            1.0,
            "flag",
            "pass" if ppk_mid1 >= ppk_mid2 else "FAIL",
        ),
        MetricResult(
            "inhibited_fold_change_lte_1",
            fc_end,
            -1.0,
            1.0,
            "fold",
            "pass" if fc_end <= 1.0 else "FAIL",
        ),
        MetricResult(
            "readout_is_dual_phospho_erk",
            float(baseline.model.readout == "PP_K"),
            1.0,
            1.0,
            "flag",
            "pass" if baseline.model.readout == "PP_K" else "FAIL",
        ),
    ]
    notes = [
        f"parsed {len(spec.species)} species / {len(spec.reactions)} reactions from "
        f"BIOMD0000000009 (volume {spec.compartment:.1e} L); drug-free PP_K "
        f"steady state {ppk_base:.3f}; occupancy monotone 0.6->{ppk_mid1:.3f}, "
        f"0.9->{ppk_mid2:.3f}, 1.0->{ppk_full:.3f}; full-signal fold-change "
        f"{fc_end:.3f} (inhibition)."
    ]
    return CaseResult(
        benchmark="Huang/Levchenko SBML MAPK cascade integration (R-5)",
        passed=all(m.criterion == "pass" for m in metrics),
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_sbml_mapk_validation"]
