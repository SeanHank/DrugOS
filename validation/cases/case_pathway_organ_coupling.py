"""Pathway -> organ regeneration coupling (L2 — mechanistic limit, doc/08 §1.2).

Closes the occupancy -> pathway -> organ -> phenotype chain at the liver:
the Stage-3 ERK/MAPK readout fold-change (relative to the drug-free baseline)
scales hepatocyte regeneration through the damped, bounded
``regeneration_scale`` mapping (doc/05 4.2).  At a fixed liver exposure a
suppressed proliferative readout must slow regeneration and *increase* the
death fraction; a stimulated readout must *decrease* it; a neutral signal
must be idempotent with the no-pathway run; and the reported total-bilirubin
rise must honour its ``bile_rise_max_fold`` ceiling while the underlying
dose-response keeps climbing.  All checks are against the stage's own
equations, so this is analytic, not external-data evidence.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, _simulate

from drugos.organ.liver import LiverParams, simulate_liver
from drugos.pk.simulate import PBPKResult

ACETAMINOPHEN_MW = 151.2


def _run(res: PBPKResult, dose_mg: float, **kwargs: object) -> object:
    drift = dose_mg / 1000.0
    c = res.unbound_tissues["liver"] * drift
    return simulate_liver(res.t, c, ACETAMINOPHEN_MW, n_eval=500, **kwargs)


def case_pathway_organ_coupling() -> CaseResult:
    res = _simulate(_bench())
    base = _run(res, 20000.0)
    neutral = _run(res, 20000.0, proliferation_signal=np.ones_like(res.t))
    suppressed = _run(res, 20000.0, proliferation_signal=np.zeros_like(res.t))
    stimulated = _run(res, 20000.0, proliferation_signal=np.full_like(res.t, 3.0))

    base_dead = max(base.max_dead_frac, 1e-9)
    met_suppression = MetricResult(
        "suppressed_regen_death_fold",
        suppressed.max_dead_frac / base_dead,
        1.001,
        1e3,
        "fold",
        "pass" if suppressed.max_dead_frac > base.max_dead_frac else "FAIL",
    )
    met_stimulation = MetricResult(
        "stimulated_regen_death_fold",
        stimulated.max_dead_frac / base_dead,
        1e-6,
        0.999,
        "fold",
        "pass" if stimulated.max_dead_frac < base.max_dead_frac else "FAIL",
    )
    met_monotone = MetricResult(
        "suppression_to_stimulation_fold",
        suppressed.max_dead_frac / max(stimulated.max_dead_frac, 1e-9),
        1.01,
        1e3,
        "fold",
        "pass" if suppressed.max_dead_frac > stimulated.max_dead_frac else "FAIL",
    )
    met_neutral = MetricResult(
        "neutral_signal_dead_absdiff",
        abs(neutral.max_dead_frac - base.max_dead_frac),
        0.0,
        1e-6,
        "fraction",
        "pass" if abs(neutral.max_dead_frac - base.max_dead_frac) <= 1e-6 else "FAIL",
    )
    met_floor = MetricResult(
        "regen_scale_floor",
        float(np.min(suppressed.regen_scale)),
        0.499,
        0.501,
        "scale",
        "pass" if 0.499 <= float(np.min(suppressed.regen_scale)) <= 0.501 else "FAIL",
    )
    met_ceiling = MetricResult(
        "regen_scale_ceiling",
        float(np.max(stimulated.regen_scale)),
        1.499,
        1.501,
        "scale",
        "pass" if 1.499 <= float(np.max(stimulated.regen_scale)) <= 1.501 else "FAIL",
    )

    params_hi = LiverParams(bili_rise_per_bsep=4.0, bili_rise_per_dead=8.0, bile_rise_max_fold=2.0)
    capped = _run(res, 20000.0, params=params_hi)
    uncapped = _run(res, 20000.0, params=replace(params_hi, bile_rise_max_fold=3.0))
    met_bili_capped = MetricResult(
        "bilirubin_capped_xULN",
        capped.peak_bilirubin_uln,
        1.999,
        2.001,
        "xULN",
        "pass" if abs(capped.peak_bilirubin_uln - 2.0) <= 1e-3 else "FAIL",
    )
    met_bili_uncapped = MetricResult(
        "bilirubin_uncapped_rise_xULN",
        uncapped.peak_bilirubin_uln,
        2.01,
        3.0,
        "xULN",
        "pass" if 2.01 <= uncapped.peak_bilirubin_uln <= 3.0 else "FAIL",
    )

    metrics = (
        met_suppression,
        met_stimulation,
        met_monotone,
        met_neutral,
        met_floor,
        met_ceiling,
        met_bili_capped,
        met_bili_uncapped,
    )
    ok = sum(1 for m in metrics if m.criterion == "pass") == len(metrics)
    return CaseResult(
        "pathway->organ regeneration coupling + bilirubin ceiling",
        ok,
        list(metrics),
        [
            f"at {suppressed.max_dead_frac:.4f} suppressed vs "
            f"{base.max_dead_frac:.4f} baseline vs {stimulated.max_dead_frac:.4f} "
            f"stimulated max dead fraction with regen_scale clamped to "
            f"[{float(np.min(suppressed.regen_scale)):.2f}, "
            f"{float(np.max(stimulated.regen_scale)):.2f}]; bilirubin capped at "
            f"{capped.peak_bilirubin_uln:.2f}xULN while uncapped hits "
            f"{uncapped.peak_bilirubin_uln:.2f}xULN. "
            "The blocked ERK/proliferation readout attenuates (never ablates) "
            "hepatocyte regeneration, tipping the same direct stress into more "
            "cell death (occupancy -> pathway -> organ -> phenotype)."
        ],
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


def _bench() -> object:
    """The acetaminophen benchmark entry used to drive the PBPK exposure."""
    from validation.benchmarks import BENCHMARKS

    return next(b for b in BENCHMARKS if b.name == "acetaminophen")


__all__ = ["case_pathway_organ_coupling"]
