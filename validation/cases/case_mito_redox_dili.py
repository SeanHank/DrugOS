"""L2 analytic-Limit DILI evidence: mitochondrial + redox + ATP-floor axes.

The Mito/redox/hepatocyte-death axes of doc/05 4.3 are real shipped QST
machinery in ``src/drugos/organ/liver.py`` (the documented calibration
model; DILIsym-equivalent closiness is proprietary and is not claimed).
This case pins the *closed-form* behavior the axes must exhibit so the
end-to-end DILI lane is anchored on physics, not on a fitted constant:

- **mitochondrial_block** — worst-case ETC inhibition across the measured
  complexes (:func:`mitochondrial_block`); at C == IC50 a single hit is
  exactly 0.5 (Michaelis), and several hits reduce to their maximum
  (no averaging: the failing complex binds the electron transport);
- **redox_state** — ROS is normalized c/(c+IC50) and GSH fraction is
  1 - ROS with the declared 0.15 viability floor (below which the
  analytical curve would predict sanitizing the cell of glutathione);
  at C == IC50 exactly (ROS, GSH, product) = (0.5, 0.5, 0.25);
- **aten_floor_factor** — the adaptive (mitogenic) ATP-recovery term
  re-scales the *unblocked* supply only; a fully blocked transport can
  never be rescued below the declared floor, so ATP = atp_floor exactly
  at total block and 1.0 at zero block;
- **combined_stress** — default-determined composition (cholestasis,
  ATP, GSH = 0.5/0.3/0.2 weights), so a lone cholestasis hit is 0.5 and
  an all-axes hit is exactly 1.0, with the immune axis inert at its
  default weight 0;
- **integrated exposure differential** — engaging the axis through
  realistic panel IC50s (ETC complex and redox 1 µM) at a saturating
  free-hepatic exposure must raise end-of-horizon hepatocyte death and
  bury ATP / GSH below the inert controls, monotone in exposure.
"""

from __future__ import annotations

import numpy as np
from validation.cases.base import CaseResult, EvidenceLevel, MetricResult

from drugos.organ.liver import (
    LiverParams,
    aten_floor_factor,
    combined_stress,
    mitochondrial_block,
    redox_state,
    simulate_liver,
)


def case_mito_redox_dili() -> CaseResult:
    metrics: list[MetricResult] = []
    notes: list[str] = []
    ok_flags: list[bool] = []

    def pin(name: str, value: float, lo: float, hi: float, unit: str) -> bool:
        good = lo <= value <= hi
        metrics.append(MetricResult(name, value, lo, hi, unit, "pass" if good else "FAIL"))
        ok_flags.append(good)
        return good

    # 1. mitochondrial_block closed forms.
    pin("mito_empty_tuple_inert", mitochondrial_block(1.0, ()), 0.0, 0.0, "block")
    pin("mito_zero_exposure_inert", mitochondrial_block(0.0, (10.0, 100.0)), 0.0, 0.0, "block")
    pin(
        "mito_single_hit_at_ic50",
        mitochondrial_block(50.0, (50.0,)),
        0.5,
        0.5,
        "block (C==IC50 -> 0.5 Michaelis)",
    )
    worst = mitochondrial_block(50.0, (10.0, 100.0))
    single = max(50.0 / 60.0, 50.0 / 150.0)
    pin(
        "mito_worst_case_max",
        worst,
        single - 1e-12,
        single + 1e-12,
        "block (max across complexes)",
    )

    # 2. redox_state closed forms.
    ros0, gsh0, net0 = redox_state(0.0, 6.0e5)
    pin("redox_zero_ros", ros0, 0.0, 0.0, "ROS")
    pin("redox_zero_gsh", gsh0, 1.0, 1.0, "GSH frac")
    pin("redox_zero_product", net0, 0.0, 0.0, "ROS*GSH")
    ros_m, gsh_m, net_m = redox_state(100.0, 100.0)
    pin("redox_at_ic50_ros", ros_m, 0.5, 0.5, "ROS")
    pin("redox_at_ic50_gsh", gsh_m, 0.5, 0.5, "GSH frac")
    pin("redox_at_ic50_net", net_m, 0.25, 0.25, "ROS*GSH")
    _, gsh_sat, _ = redox_state(1.0e4, 100.0)
    pin("redox_floor_binds", gsh_sat, 0.15, 0.15, "GSH frac (floor)")

    # 3. aten_floor_factor closed forms.
    pin("atp_full_block_at_floor", aten_floor_factor(1.0, 0.2, 0.0), 0.2, 0.2, "ATP")
    pin("atp_zero_block_supplied", aten_floor_factor(0.0, 0.2, 0.0), 1.0, 1.0, "ATP")
    pin(
        "atp_half_block_closed_form",
        aten_floor_factor(0.4, 0.2, 0.0),
        0.68 - 1e-9,
        0.68 + 1e-9,
        "ATP (0.2 + 0.8*0.6)",
    )
    pin(
        "atp_adaptive_boost",
        aten_floor_factor(0.4, 0.2, 0.5),
        0.92 - 1e-9,
        0.92 + 1e-9,
        "ATP (0.2 + 0.8*0.6*1.5)",
    )
    pin(
        "atp_adaptive_inert_at_full_block",
        aten_floor_factor(1.0, 0.2, 99.0),
        0.2,
        0.2,
        "ATP (shortfall zero -> adaptive inert)",
    )

    # 4. combined_stress default composition.
    pin("stress_lone_cholestasis", combined_stress(1.0, 1.0, 1.0), 0.5, 0.5, "stress")
    pin("stress_all_axes", combined_stress(1.0, 0.0, 0.0), 1.0, 1.0, "stress")
    pin(
        "stress_immune_inert_at_default",
        combined_stress(0.5, 1.0, 1.0, immune=1.0),
        combined_stress(0.5, 1.0, 1.0),
        combined_stress(0.5, 1.0, 1.0),
        "stress (immune weight 0)",
    )

    # 5. Degenerate inputs raise.
    raises = 0
    for call in (
        lambda: mitochondrial_block(50.0, (0.0,)),
        lambda: redox_state(50.0, 0.0),
        lambda: aten_floor_factor(0.5, 1.5, 0.0),
    ):
        try:
            call()
        except ValueError:
            raises += 1
    pin("degenerate_mito_redox_rejected", float(raises), 3.0, 3.0, "count")

    # 6. Integrated exposure differential through simulate_liver.
    t = np.linspace(0.0, 72.0, 241)
    c = np.full_like(t, 10.0)  # 10 mg/L, MW 400 -> 25 uM free hepatic
    inert = LiverParams(mito_ic50_nm=(1.0e12,), redox_ic50_nm=1.0e12)
    mito = LiverParams(mito_ic50_nm=(1000.0,), redox_ic50_nm=1.0e12)
    redox = LiverParams(mito_ic50_nm=(1.0e12,), redox_ic50_nm=1000.0)

    r_inert = simulate_liver(t, c, 400.0, inert)
    r_mito = simulate_liver(t, c, 400.0, mito)
    r_redox = simulate_liver(t, c, 400.0, redox)
    r_half = simulate_liver(t, 0.5 * c, 400.0, mito)

    pin(
        "mito_axis_buried_atp",
        float(np.min(r_mito.atp_frac)),
        0.0,
        float(np.min(r_inert.atp_frac)) - 0.05,
        "min ATP frac (inert >= mito + 0.05)",
    )
    pin(
        "mito_axis_increases_dead",
        float(r_mito.dead_frac[-1]),
        r_inert.dead_frac[-1] + 1e-6,
        1.0,
        "dead_frac @72h (mito engaged)",
    )
    pin(
        "redox_axis_increases_dead",
        float(r_redox.dead_frac[-1]),
        r_inert.dead_frac[-1] + 1e-6,
        1.0,
        "dead_frac @72h (redox engaged)",
    )
    pin(
        "mito_death_monotone_in_exposure",
        float(r_half.dead_frac[-1]),
        0.0,
        r_mito.dead_frac[-1] + 1e-6,
        "dead_frac @72h (half exposure)",
    )

    notes.append(
        f"worst-case block {worst:.4f} == max of single-complex hits {single:.4f}; "
        f"redox (ROS,GSH,net)@IC50 = ({ros_m:.2f},{gsh_m:.2f},{net_m:.2f}); "
        f"GSH floor binds at {gsh_sat:.2f}; ATP closed-forms {aten_floor_factor(0.4, 0.2, 0.0):.2f}"
        f"/{aten_floor_factor(0.4, 0.2, 0.5):.2f}; "
        f"engaged dead_frac @72h: mito={r_mito.dead_frac[-1]:.4f}, "
        f"redox={r_redox.dead_frac[-1]:.4f}, inert={r_inert.dead_frac[-1]:.4f}; "
        f"mint ATP mito={np.min(r_mito.atp_frac):.3f} vs inert={np.min(r_inert.atp_frac):.3f}"
    )

    return CaseResult(
        "Mito/redox/ATP-floor DILI axes (closed-form + end-to-end differential)",
        all(ok_flags),
        metrics,
        notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_mito_redox_dili"]
