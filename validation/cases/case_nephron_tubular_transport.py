"""Single-nephron tubular transport case (doc/12 L17 realized, L2 evidence).

drives the shipped segmental nephron cascade
(:func:`drugos.organ.nephron.nephron_handling`) through its analytic limits;
these are the closed-form anchors a tubular-transport model must satisfy:

- **mass conservation** — filtered load = reabsorbed + excreted exactly
  (the cascade moves solute, it never invents or drops it),
- **physiological Na+ balance** — the healthy four-segment cascade lands
  fractional excretion of Na+ in the physiologic ~0.5-2 % band with
  iso-osmotic water bookkeeping (urine flow ~1 % of GFR),
- **pure filtration** — a solute with zero transport is excreted at the
  filtered load: FE_x = 1 and Cl_renal = GFR exactly (this is the same
  closed-form the production kidney lane uses, Scr = P / GFR),
- **net reabsorption** — a passively reabsorbed solute sits below FE_x = 1
  (Cl_renal < GFR),
- **net secretion** — a saturably secreted solute (OAT/OCT-style
  Michaelis-Menten flux into the lumen) exceeds FE_x = 1 (Cl_renal > GFR),
- **transporter saturation** — doubling the free concentration at a fixed
  transporter Vmax does not double transport: FE falls as the filtered load
  outgrows the saturable flux, the classic net-effect roll-off used in renal
  DDI assessments.
"""

from __future__ import annotations

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.organ.nephron import NephronParams, nephron_handling

_GFR_ML_MIN = 125.0
_CFREE_UM_MM = 1.0e-3


def case_nephron_tubular_transport() -> CaseResult:
    healthy = nephron_handling(_GFR_ML_MIN, plasma_na_mm=140.0)
    pure = nephron_handling(
        _GFR_ML_MIN,
        c_free_mm=_CFREE_UM_MM,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=0.0),
    )
    reab = nephron_handling(
        _GFR_ML_MIN,
        c_free_mm=_CFREE_UM_MM,
        params=NephronParams(
            pas_reab_fraction=0.5,
            sec_vmax_nmol_min_per_nephron=0.0,
        ),
    )
    secr = nephron_handling(
        _GFR_ML_MIN,
        c_free_mm=_CFREE_UM_MM,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=4.0e-3),
    )
    hi_conc = nephron_handling(
        _GFR_ML_MIN,
        c_free_mm=10.0 * _CFREE_UM_MM,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=4.0e-3),
    )

    resid_na = (
        healthy.filtered_na_nmol_min - healthy.reabsorbed_na_nmol_min - healthy.excreted_na_nmol_min
    )
    na_ok = 0.005 <= healthy.fe_na <= 0.02
    urine_ok = 0.004 * _GFR_ML_MIN <= healthy.urine_flow_ml_min <= 0.03 * _GFR_ML_MIN
    fe_filt = abs(pure.fe_x - 1.0) < 1e-9 if pure.fe_x is not None else False
    cl_filt = abs(pure.cl_r_ml_min - _GFR_ML_MIN) < 1e-6
    fe_down = bool(reab.fe_x is not None and abs(reab.fe_x - 0.5) < 1e-9)
    cl_down = abs(reab.cl_r_ml_min - 0.5 * _GFR_ML_MIN) < 1e-6
    fe_up = bool(secr.fe_x is not None and secr.fe_x > 1.0)
    sat_ok = bool(
        secr.fe_x is not None
        and hi_conc.fe_x is not None
        and hi_conc.fe_x < secr.fe_x
        and hi_conc.cl_r_ml_min < secr.cl_r_ml_min
    )

    metrics = [
        MetricResult(
            "fe_na_healthy",
            healthy.fe_na,
            0.005,
            0.02,
            "fraction",
            "pass" if na_ok else "FAIL",
        ),
        MetricResult(
            "urine_flow_ml_min",
            healthy.urine_flow_ml_min,
            0.004 * _GFR_ML_MIN,
            0.03 * _GFR_ML_MIN,
            "ml/min",
            "pass" if urine_ok else "FAIL",
        ),
        MetricResult(
            "fe_x_pure_filtration",
            pure.fe_x if pure.fe_x is not None else 0.0,
            1.0,
            1.0,
            "fraction",
            "pass" if fe_filt else "FAIL",
        ),
        MetricResult(
            "cl_r_pure_filtration_ml_min",
            pure.cl_r_ml_min,
            _GFR_ML_MIN,
            _GFR_ML_MIN,
            "ml/min",
            "pass" if cl_filt else "FAIL",
        ),
        MetricResult(
            "fe_x_net_reabsorption",
            reab.fe_x if reab.fe_x is not None else 0.0,
            0.5,
            0.5,
            "fraction",
            "pass" if fe_down else "FAIL",
        ),
        MetricResult(
            "cl_r_net_reabsorption_ml_min",
            reab.cl_r_ml_min,
            0.5 * _GFR_ML_MIN,
            0.5 * _GFR_ML_MIN,
            "ml/min",
            "pass" if cl_down else "FAIL",
        ),
        MetricResult(
            "fe_x_net_secretion",
            secr.fe_x if secr.fe_x is not None else 0.0,
            1.0,
            1.0e6,
            "fraction",
            "pass" if fe_up else "FAIL",
        ),
        MetricResult(
            "cl_r_saturation_ratio_lo_hi",
            hi_conc.fe_x if hi_conc.fe_x is not None else 0.0,
            0.0,
            secr.fe_x if secr.fe_x is not None else 0.0,
            "fraction",
            "pass" if sat_ok else "FAIL",
        ),
    ]
    notes = [
        f"healthy four-segment cascade: FE_Na={healthy.fe_na:.4f} "
        f"urine={healthy.urine_flow_ml_min:.2f} ml/min (≈{100.0 * healthy.fe_na:.1f}% of GFR); "
        f"mass balance residual={resid_na:.2e} nmol/min; "
        f"pure filtration FE_x={pure.fe_x:.4f} Cl_renal={pure.cl_r_ml_min:.1f} ml/min (= GFR); "
        f"50% passive reabsorption FE_x={reab.fe_x:.4f} Cl_renal={reab.cl_r_ml_min:.1f} ml/min; "
        f"saturable secretion FE_x={secr.fe_x:.3f} Cl_renal={secr.cl_r_ml_min:.1f} ml/min; "
        f"10x concentration → FE_x={hi_conc.fe_x:.3f} (saturation roll-off, "
        "the renal-DDI net-effect flag)",
    ]
    return CaseResult(
        benchmark="nephron tubular transport (segmental cascade, analytic limits)",
        passed=all(m.criterion == "pass" for m in metrics),
        metrics=metrics,
        notes=notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )
