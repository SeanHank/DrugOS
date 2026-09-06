"""Pathway graph DSL for Stage 3 (doc/05 section 3).

Small, auditable reaction vocabulary that compiles to a sparse ODE system:

* :class:`Activation` — enzyme- or signal-driven activation of a substrate
  into a product (saturating Michaelis-Menten term on the substrate and a
  Hill-coefficient driver term; ``input_drive`` reactions are driven directly
  by the Stage-2 occupancy signal, 0..1).  Optional competitive inhibition by
  another species.
* :class:`ReversibleBinding` — mass-action complex formation A + B <-> AB.
* :class:`FirstOrderDegradation` — linear loss of a species.
* :class:`SourceProduction` — constant basal production into a species.

Graphs are pruned to the tractable size of doc/05 3.2 (tens of reactions) and
are exported to the shared ODE contract by the simulator.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True)
class Activation:
    """Signal-driven activation: substrate -> product.

    rate = drive * Vmax * S^n_s/(Km^n_s + S^n_s) * (1/(1 + I/Ki) if inhibited)

    ``drive`` equals the occupancy signal ``s(t)`` when ``input_drive`` (the
    drug-derived perturbation from Stage 2) or the (saturating) Hill term of
    the driver species otherwise.
    """

    substrate: str
    product: str
    driver: str | None = None
    vmax: float = 1.0
    km: float = 1.0
    substrate_hill: float = 1.0
    driver_half: float = 1.0
    driver_hill: float = 1.0
    input_drive: bool = False
    inhibitor: str | None = None
    ki: float = 1.0

    def __post_init__(self) -> None:
        if self.vmax <= 0 or self.km <= 0:
            raise ValueError("vmax and km must be positive")
        if self.input_drive and self.driver is not None:
            raise ValueError("input_drive reactions must not name a driver species")
        if not self.input_drive and self.driver is None:
            raise ValueError("indirectly driven reactions must name a driver")


@dataclass(frozen=True, slots=True)
class ReversibleBinding:
    """Mass-action reversible complex formation A + B <-> AB."""

    a: str
    b: str
    complex: str
    kon: float = 1.0
    koff: float = 0.1

    def __post_init__(self) -> None:
        if self.kon <= 0 or self.koff <= 0:
            raise ValueError("kon and koff must be positive")


@dataclass(frozen=True, slots=True)
class FirstOrderDegradation:
    """First-order loss of a species (degradation / decay)."""

    species: str
    rate: float = 0.1

    def __post_init__(self) -> None:
        if self.rate <= 0:
            raise ValueError("rate must be positive")


@dataclass(frozen=True, slots=True)
class SourceProduction:
    """Constant basal production into a species."""

    species: str
    rate: float = 0.1

    def __post_init__(self) -> None:
        if self.rate <= 0:
            raise ValueError("rate must be positive")


@dataclass(slots=True)
class PathwayModel:
    """A pathway: species initial conditions plus a reaction list."""

    name: str
    species: dict[str, float] = field(default_factory=dict)
    reactions: list[Reaction] = field(default_factory=list)
    input_node: str | None = None
    readout: str | None = None

    def add_activation(self, activation: Activation) -> None:
        self.reactions.append(activation)

    def add_reaction(self, reaction: Reaction) -> None:
        self.reactions.append(reaction)

    def node_names(self) -> list[str]:
        return sorted(self.species)


#: A pathway reaction: any of the primitive kinetic steps defined above.
Reaction = Activation | ReversibleBinding | FirstOrderDegradation | SourceProduction


__all__ = [
    "Activation",
    "FirstOrderDegradation",
    "PathwayModel",
    "ReversibleBinding",
    "SourceProduction",
]
