"""Finite-dose multi-layer skin permeation for transdermal delivery (doc/05 1.4).

The transdermal route defaults to the generic first-order ``depot``
compartment; ``SkinLayers`` (opt-in, off by default) replaces it with a
four-layer membrane: vehicle/surface reservoir -> stratum corneum (SC) ->
viable epidermis (VE) -> dermis, with first-order removal from the dermis into
the dermal capillary bed (systemic venous blood).

Each inter-layer link is a reversible, diffusion-limited flux in the classic
serial-resistance compartment form

    J_k = D_k * area / L_k * (C_up - C_down / K_k)

where ``D_k`` is the *receiving* layer's effective diffusivity, ``L_k`` its
thickness, and ``K_k`` the equilibrium partition ``C_down/C_up`` at the
boundary.  At steady state with the dermis sink held at zero the three links
carry one identical flux, which for a surface concentration ``C0`` yields the
serial 3-resistance permeability (cm/h)

    P_eff = 1 / (L_sc/D_sc + L_ve/(D_ve*K_sc) + L_der/(D_der*K_ve*K_sc))

so ``J = P_eff * A * C0`` — the Fick steady-state anchor exercised by the L2
case (``case_transdermal_multi_layer``).  The dermal capillary removal
``k_dermal_capillary_1h * A_dermis`` is weighted by ``depot_bioavailability``
(the complement is tallied in the skin ``unabsorbed`` sink), matching the
existing depot semantics.  All layer amounts stay inside the PBPK mass tally,
so administered transdermal drug is conserved across surface + membranes +
systemic + unabsorbed.
"""

from __future__ import annotations

from dataclasses import dataclass

_UM_PER_CM = 1e-4


@dataclass(frozen=True, slots=True)
class SkinLayers:
    """Four-layer skin membrane parameters for a transdermal finite dose.

    Thicknesses are in micrometres, diffusivities in cm^2/h (effective
    membrane diffusivities), areas in cm^2.  Partitions are the equilibrium
    concentration ratios ``C_down/C_up`` at each interface (``>1`` means the
    receiving layer concentrates the drug).
    """

    area_cm2: float = 40.0
    surface_thickness_um: float = 50.0
    sc_thickness_um: float = 20.0
    ve_thickness_um: float = 100.0
    dermis_thickness_um: float = 150.0
    sc_diffusivity_cm2_h: float = 1.0e-5
    ve_diffusivity_cm2_h: float = 2.0e-3
    dermis_diffusivity_cm2_h: float = 5.0e-3
    surface_sc_partition: float = 1.0
    sc_ve_partition: float = 1.0
    ve_dermis_partition: float = 1.0
    k_dermal_capillary_1h: float = 5.0

    def __post_init__(self) -> None:
        values: list[tuple[float, str]] = [
            (self.area_cm2, "area_cm2"),
            (self.surface_thickness_um, "surface_thickness_um"),
            (self.sc_thickness_um, "sc_thickness_um"),
            (self.ve_thickness_um, "ve_thickness_um"),
            (self.dermis_thickness_um, "dermis_thickness_um"),
            (self.sc_diffusivity_cm2_h, "sc_diffusivity_cm2_h"),
            (self.ve_diffusivity_cm2_h, "ve_diffusivity_cm2_h"),
            (self.dermis_diffusivity_cm2_h, "dermis_diffusivity_cm2_h"),
            (self.surface_sc_partition, "surface_sc_partition"),
            (self.sc_ve_partition, "sc_ve_partition"),
            (self.ve_dermis_partition, "ve_dermis_partition"),
            (self.k_dermal_capillary_1h, "k_dermal_capillary_1h"),
        ]
        for value, name in values:
            if not value > 0.0:
                raise ValueError(f"{name} must be positive; got {value}")

    def thickness_cm(self, thickness_um: float) -> float:
        return thickness_um * _UM_PER_CM

    def layer_volumes_cm3(self) -> dict[str, float]:
        """Physical layer volumes = area * thickness (cm^3)."""
        return {
            "surface": self.area_cm2 * self.thickness_cm(self.surface_thickness_um),
            "sc": self.area_cm2 * self.thickness_cm(self.sc_thickness_um),
            "ve": self.area_cm2 * self.thickness_cm(self.ve_thickness_um),
            "dermis": self.area_cm2 * self.thickness_cm(self.dermis_thickness_um),
        }

    def transfer_cm3_h(self, diffusivity_cm2_h: float, thickness_um: float) -> float:
        """Diffusion-link conductance ``D * area / L`` for one interface (cm^3/h)."""
        return diffusivity_cm2_h * self.area_cm2 / self.thickness_cm(thickness_um)

    def sc_conductance_cm3_h(self) -> float:
        return self.transfer_cm3_h(self.sc_diffusivity_cm2_h, self.sc_thickness_um)

    def ve_conductance_cm3_h(self) -> float:
        return self.transfer_cm3_h(self.ve_diffusivity_cm2_h, self.ve_thickness_um)

    def dermis_conductance_cm3_h(self) -> float:
        return self.transfer_cm3_h(self.dermis_diffusivity_cm2_h, self.dermis_thickness_um)

    def composite_permeability_cm_h(self) -> float:
        """Serial 3-resistance steady-state permeability (Fick, doc/05 1.4).

        ``J_ss = P_eff * A * C_surface`` with the dermis held at zero by the
        capillary sink.  Partition corrections enter each downstream
        resistance so a concentrating barrier (K>1) slows the flux.
        """
        l_sc = self.thickness_cm(self.sc_thickness_um)
        l_ve = self.thickness_cm(self.ve_thickness_um)
        l_der = self.thickness_cm(self.dermis_thickness_um)
        resistance = (
            l_sc / self.sc_diffusivity_cm2_h
            + l_ve / (self.ve_diffusivity_cm2_h * self.surface_sc_partition)
            + l_der
            / (self.dermis_diffusivity_cm2_h * self.sc_ve_partition * self.surface_sc_partition)
        )
        return 1.0 / resistance
