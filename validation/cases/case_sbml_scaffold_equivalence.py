"""Second SBML scaffold: non-MAPK pathway parse + steady-state equivalence (L18).

doc/12 L18 called for "Physiome/Reactome pathway SBML scaffolds ... parse +
steady-state equivalence per scaffold, gated like R-5".  The SBML lane is now
generic (:func:`simulate_sbml_pathway` takes ``input_species`` / ``readout``),
and a second openly-licensed scaffold is vendored: the **Rohwer 2000 E. coli
PEP:glucose phosphotransferase system** (BioModels BIOMD0000000038, CC0) — a
non-MAPK, metabolism-domain pathway model, driving the phosphotransfer of the
energy donor phosphoenolpyruvate (PEP) into glucose and pyruvate.

The equivalence gate applies the same three checks the R-5 MAPK gate uses:

- **faithful parse** — every reaction in the artifact compiled into its ODE
  kinetic contribution (reaction counts agree with a direct libsbml read, so
  nothing is silently dropped), and all 17 species are resolved;
- **steady-state equivalence** — the drug-free system is at a true fixed
  point: a 200 h and a 400 h horizon land on the same pyruvate level (the
  trajectory tail is time-translation-invariant) and both match the computed
  per-scaffold baseline;
- **monotone substrate-gate** — scaling the PEP donor down (the honest
  energy-supply input: glucose is consumed at equilibrium, so Glc is not the
  probe) suppresses pyruvate non-increasingly, with full donor depletion
  well below the quarter-supply level.

R-5 is re-gated in the same file to prove the generalization did not disturb
the wired MAPK default.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.pathway import parse_sbml, simulate_sbml_pathway
from drugos.pathway.sbml_pathway import _import_libsbml, model_file

_DATA = Path(__file__).resolve().parents[2] / "data" / "models"
_PTS_FILE = _DATA / "rohwer2000-pts.xml"

_N = 120


def _sweep(
    path: Path,
    occs: tuple[float, ...],
    input_species: str,
    readout: str,
    horizon_h: float,
) -> list[float]:
    t = np.linspace(0.0, horizon_h, 240)
    folds: list[float] = []
    for occ in occs:
        r = simulate_sbml_pathway(
            t,
            np.full_like(t, occ),
            path=path,
            input_species=input_species,
            readout=readout,
            n_eval=_N,
        )
        folds.append(float(r.readout_fold_change(readout)[-1]))
    return folds


def case_sbml_scaffold_equivalence() -> CaseResult:
    metrics: list[MetricResult] = []
    notes: list[str] = []
    ok_flags: list[bool] = []

    def pin(name: str, value: float, lo: float, hi: float, unit: str) -> bool:
        good = lo <= value <= hi
        metrics.append(MetricResult(name, value, lo, hi, unit, "pass" if good else "FAIL"))
        ok_flags.append(good)
        return good

    # 1. Faithful parse: direct libsbml counts agree with the compiled spec.
    spec = parse_sbml(_PTS_FILE)
    lib = _import_libsbml()
    doc = lib.readSBML(str(_PTS_FILE))
    model = doc.getModel()
    n_sp, n_rx = model.getNumSpecies(), model.getNumReactions()
    pin(
        "pts_species_full_resolution",
        float(len(spec.species)),
        float(n_sp),
        float(n_sp),
        "species (libsbml == spec)",
    )
    pin(
        "pts_reactions_all_compiled",
        float(len(spec.reactions)),
        float(n_rx),
        float(n_rx),
        "reactions (no silent drop)",
    )
    pin("pts_distinct_from_mapk", float("PP_K" not in spec.species), 1.0, 1.0, "flag")

    # 2. Steady-state equivalence across horizons (time-translation invariant).
    t_lo = np.linspace(0.0, 200.0, 240)
    t_hi = np.linspace(0.0, 400.0, 241)
    r_lo = simulate_sbml_pathway(
        t_lo, np.zeros_like(t_lo), path=_PTS_FILE, input_species="PEP", readout="Pyr", n_eval=_N
    )
    r_hi = simulate_sbml_pathway(
        t_hi, np.zeros_like(t_hi), path=_PTS_FILE, input_species="PEP", readout="Pyr", n_eval=_N
    )
    pyr_lo = float(r_lo.concentrations["Pyr"][-1])
    pyr_hi = float(r_hi.concentrations["Pyr"][-1])
    base = r_lo.baseline["Pyr"]
    drift = abs(pyr_hi - pyr_lo) / max(pyr_lo, 1e-12)
    pin("pts_fixed_point_drift", drift, 0.0, 2e-2, "rel. drift (200h vs 400h)")
    pin("pts_baseline_matches_equilated", abs(pyr_lo / base - 1.0), 0.0, 2e-2, "fold")

    # 3. Monotone substrate-gate on the PEP donor.
    folds = _sweep(_PTS_FILE, (0.25, 0.5, 0.75, 1.0), "PEP", "Pyr", 200.0)
    mono = all(folds[i] >= folds[i + 1] for i in range(len(folds) - 1))
    pin("pts_monotone_in_donor_supply", float(mono), 1.0, 1.0, "flag")
    pin("pts_full_donor_depletion", folds[-1], 0.0, 0.85 * folds[0], "fold @ s=1")

    # 4. R-5 still gated (generalized defaults preserved the MAPK lane).
    fc_mapk = _sweep(model_file(), (0.0, 0.6, 0.9, 1.0), "E1", "PP_K", 72.0)
    mapk_ok = fc_mapk[-1] < min(0.05, fc_mapk[1]) and fc_mapk[1] >= fc_mapk[2] >= fc_mapk[3]
    pin("mapk_default_lane_still_gated", float(mapk_ok), 1.0, 1.0, "flag")

    notes.append(
        f"Rohwer PTS: {len(spec.species)} species / {len(spec.reactions)} reactions "
        f"(libsbml {n_sp}/{n_rx}); pyruvate drift 200h->400h {drift:.1e}; "
        f"donor sweep {[round(x, 4) for x in folds]}; MAPK sweep "
        f"{[round(x, 4) for x in fc_mapk]}"
    )

    return CaseResult(
        "Second SBML scaffold: Rohwer 2000 PTS equivalence gate (L18, non-MAPK)",
        all(ok_flags),
        metrics,
        notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_sbml_scaffold_equivalence"]
