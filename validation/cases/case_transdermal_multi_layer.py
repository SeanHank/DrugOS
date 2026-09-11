"""Finite-dose multi-layer skin permeation for transdermal (L2, doc/05 1.4).

Transdermal was a generic first-order ``depot``; the deferred item asked for a
real multi-layer skin-permeation membrane.  That is now
``PBPKModel.absorption.skin_layers`` (``SkinLayers``, off by default): four
reversible, diffusion-limited compartments surface -> stratum corneum (SC) ->
viable epidermis (VE) -> dermis, with first-order dermal capillary removal
into venous blood:

    J_k = D_k * area / L_k * (C_up - C_down / K_k)

The checks pin, on a kp=1 single pool:

- transdermal drug mass stays exactly at the dose (no metabolic/renal
  clearance) across surface + membranes + systemic + unabsorbed;
- the serial 3-resistance membrane reaches the Fick steady state — the three
  links carry one common flux equal to ``P_eff * A * C_surface`` with
  ``P_eff = 1/(L_sc/D_sc + L_ve/(D_ve*K_sc) + L_der/(D_der*K_ve*K_sc))``;
- zero-flux equilibrium returns the partition coefficients at every interface
  (C_down/C_up = K at long time with a sealed dermis);
- the membrane is a genuine barrier: a 10x thicker stratum corneum retains a
  large fraction of a 6 h finite dose, and a 10x higher SC diffusivity faster
  systemic absorption at 1 h;
- the extension is off by default (transdermal then uses the depot) and a
  degenerate (zero-area) membrane is rejected.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

from validation.cases.base import CaseResult, EvidenceLevel, MetricResult, profile

from drugos.inputs.parse_dosing import build_dose_plan
from drugos.pk.pbpk_build import AbsorptionParams, PBPKModel
from drugos.pk.simulate import simulate_pbpk
from drugos.pk.skin import SkinLayers


def _flat_partition() -> SimpleNamespace:
    names = list(profile().organ_volume) + ["arterial", "venous"]
    return SimpleNamespace(kp={n: 1.0 for n in names}, kpu={n: 1.0 for n in names})


def _transdermal(dose_mg: float, **skin_kw: float) -> PBPKModel:
    return PBPKModel(
        physiology=profile(),
        partition=_flat_partition(),
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        dose_plan=build_dose_plan("transdermal", dose_mg),
        absorption=AbsorptionParams(skin_layers=SkinLayers(**skin_kw)),
    )


def case_transdermal_multi_layer() -> CaseResult:
    metrics: list[MetricResult] = []
    dose = 20.0
    notes: list[str] = []

    # 1. Transdermal mass closes exactly against the dose (zero clearance).
    m_mass = _transdermal(dose)
    r_mass = simulate_pbpk(m_mass, tmax_h=72.0, n_eval=400)
    assert r_mass.final_state is not None
    mass_closed = float(m_mass.state_total_mass(r_mass.final_state))
    ok1 = math.isclose(mass_closed, r_mass.dose_mg, rel_tol=1e-9)
    metrics.append(
        MetricResult(
            "mass_conserved_with_skin_layers",
            mass_closed,
            r_mass.dose_mg * 0.999999999,
            r_mass.dose_mg * 1.000000001,
            "mg",
            "pass" if ok1 else "FAIL",
        )
    )

    # 2. Fick steady state: composite-permeability flux against a large,
    #    slowly-depleting surface reservoir with the dermis a strong sink.
    m_flux = _transdermal(100.0, surface_thickness_um=5000.0, k_dermal_capillary_1h=50.0)
    r_flux = simulate_pbpk(m_flux, tmax_h=6.0, n_eval=600)
    assert r_flux.final_state is not None
    skin = m_flux.absorption.skin_layers
    assert skin is not None
    y = r_flux.final_state
    j = m_flux.skin_fluxes(y)
    j_obs = j["sc_to_ve"]
    c_surf = y[m_flux._skin_indices["skin_surface"]] / skin.layer_volumes_cm3()["surface"]
    j_theory = skin.composite_permeability_cm_h() * c_surf * skin.area_cm2
    ratio = j_obs / j_theory
    linked = math.isclose(j["surface_to_sc"], j_obs, rel_tol=1e-2) and math.isclose(
        j["ve_to_dermis"], j_obs, rel_tol=1e-2
    )
    ok2 = 0.98 <= ratio <= 1.02 and linked
    metrics.append(
        MetricResult(
            "steady_flux_matches_composite_permeability",
            ratio,
            0.98,
            1.02,
            "J_obs/J_Fick",
            "pass" if ok2 else "FAIL",
        )
    )

    # 3. Zero-flux partition equilibrium returns K at every interface (sealed
    #    dermis, surface kept near-constant by a large vehicle).
    m_eq = _transdermal(
        10.0,
        k_dermal_capillary_1h=1e-4,
        surface_sc_partition=2.0,
        sc_ve_partition=3.0,
        ve_dermis_partition=4.0,
        surface_thickness_um=2000.0,
    )
    r_eq = simulate_pbpk(m_eq, tmax_h=400.0, n_eval=800)
    assert r_eq.final_state is not None
    wy = r_eq.final_state
    vols = m_eq.absorption.skin_layers.layer_volumes_cm3()
    idx = m_eq._skin_indices
    c_s = wy[idx["skin_surface"]] / vols["surface"]
    c_sc = wy[idx["skin_sc"]] / vols["sc"]
    c_ve = wy[idx["skin_ve"]] / vols["ve"]
    c_de = wy[idx["skin_dermis"]] / vols["dermis"]
    deviations = [
        (c_sc / c_s, 2.0),
        (c_ve / c_sc, 3.0),
        (c_de / c_ve, 4.0),
    ]
    max_dev = max(abs(math.log2(actual / target)) for actual, target in deviations)
    ok3 = max_dev < 0.02
    metrics.append(
        MetricResult(
            "partition_equilibrium_recovers_k",
            max_dev,
            0.0,
            0.02,
            "max log2 deviation",
            "pass" if ok3 else "FAIL",
        )
    )

    # 4. The stratum corneum is a rate-limiting barrier: a 10x thicker SC
    #    retains most of a 6 h finite dose versus the thin membrane.
    def absorbed_at_tmax_h(tmax: float, n_eval: int, **kw: float) -> float:
        res = simulate_pbpk(_transdermal(10.0, **kw), tmax_h=tmax, n_eval=n_eval)
        assert res.skin is not None
        return float(res.skin["absorbed_mg"][-1])

    absorbed_thin = absorbed_at_tmax_h(6.0, 400, sc_thickness_um=20.0)
    absorbed_thick = absorbed_at_tmax_h(6.0, 400, sc_thickness_um=200.0)
    ok4 = absorbed_thin > 8.0 and absorbed_thick < 0.6 * absorbed_thin
    metrics.append(
        MetricResult(
            "sc_barrier_retains_finite_dose",
            absorbed_thick,
            0.0,
            0.6 * absorbed_thin,
            "mg absorbed @6h (thick SC)",
            "pass" if ok4 else "FAIL",
        )
    )

    # 5. Higher SC diffusivity delivers more systemically at an early time.
    absorbed_slow = absorbed_at_tmax_h(1.0, 300, sc_diffusivity_cm2_h=1.0e-5)
    absorbed_fast = absorbed_at_tmax_h(1.0, 300, sc_diffusivity_cm2_h=1.0e-4)
    ok5 = 0.0 < absorbed_slow < absorbed_fast and absorbed_fast > 1.5 * absorbed_slow
    metrics.append(
        MetricResult(
            "diffusivity_speeds_systemic_absorption",
            absorbed_fast,
            1.5 * absorbed_slow,
            10.0,
            "mg absorbed @1h (fast SC)",
            "pass" if ok5 else "FAIL",
        )
    )

    # 6. Off by default: plain transdermal keeps the depot (no skin states).
    plain = PBPKModel(
        physiology=profile(),
        partition=_flat_partition(),
        bp=1.0,
        fup=1.0,
        cl_hep_l_h=0.0,
        cl_renal_l_h=0.0,
        dose_plan=build_dose_plan("transdermal", dose),
    )
    expected_n = 2 + len(profile().organ_volume) + 7
    base_n = m_mass.n_state - 5
    ok6 = base_n == expected_n and plain.n_state == base_n
    metrics.append(
        MetricResult(
            "off_by_default_state_count",
            base_n,
            expected_n,
            expected_n,
            "state dim without skin layers",
            "pass" if ok6 else "FAIL",
        )
    )

    # 7. Degenerate membranes raise instead of corrupting the ODE.
    raised = 0.0
    try:
        SkinLayers(area_cm2=0.0)
    except ValueError:
        raised = 1.0
    metrics.append(
        MetricResult(
            "degenerate_skin_rejected",
            raised,
            1.0,
            1.0,
            "flag",
            "pass" if raised == 1.0 else "FAIL",
        )
    )

    ok = all(m.criterion == "pass" for m in metrics)
    notes.append(
        f"mass={mass_closed:.3f}/{dose:.0f} mg; J_obs/J_Fick={ratio:.3f}; "
        f"max partition deviation={max_dev:.3f} log2; "
        f"absorbed@6h thin={absorbed_thin:.2f} vs thick={absorbed_thick:.2f} mg; "
        f"absorbed@1h slow={absorbed_slow:.2f} vs fast={absorbed_fast:.2f} mg"
    )
    return CaseResult(
        "Multi-layer transdermal skin permeation (finite-dose membrane)",
        ok,
        metrics,
        notes,
        level=EvidenceLevel.L2_ANALYTIC_LIMIT,
    )


__all__ = ["case_transdermal_multi_layer"]
