"""Target parameterization: target definitions and the off-target safety panel.

Implements doc/05 section 2.2-2.3: for known drugs, primary + off-target sites
are parameterized with affinities from literature; for novel molecules the
safety-critical panel below supplies pharmacology-informed class priors
(doc/05 2.3 priority 4), each flagged ``low_confidence`` so downstream
reporting can temper its weight.  Units are documented on every field; binding
is expressed against the *unbound* (free) drug concentration in nanomolar,
consistent with the occupancy model in ``occupancy.py``.

The panel covers the baseline safety set of doc/05 2.2: hERG, CYP
inhibition sites, bile-acid / transporter axis and mitochondrial ETC
complexes, plus endocrine nuclear receptors.
"""

from __future__ import annotations

from dataclasses import dataclass

_NM_PER_MICRO_M = 1e3


@dataclass(frozen=True, slots=True)
class Target:
    """A drug-binding site with equilibrium and turnover parameters.

    Args:
        name: canonical site label (e.g. "hERG (Kv11.1)").
        kd_nm: dissociation constant against free drug, nM.
        kon_nm_h: second-order association rate, 1/(nM*h).  Defaults to an
            equilibrium-consistent effective rate (0.5); ``koff = kon*kd``.
        r0_nm: basal target abundance (total receptor per tissue volume), nM.
        rho_h: basal target turnover rate (degradation), 1/h.  ``ksyn = rho*R0``.
        kint_h: drug-receptor complex internalization rate, 1/h.
        low_confidence: True when the value is a class-typical prior, not a
            measured affinity (flagging for report).
        reference: traceable source for the affinity value.
    """

    name: str
    kd_nm: float
    kon_nm_h: float = 0.5
    r0_nm: float = 1.0
    rho_h: float = 0.05
    kint_h: float = 0.05
    low_confidence: bool = True
    reference: str = ""

    def __post_init__(self) -> None:
        if self.kd_nm <= 0:
            raise ValueError("kd_nm must be positive")
        if self.kon_nm_h <= 0:
            raise ValueError("kon_nm_h must be positive")
        if self.r0_nm <= 0:
            raise ValueError("r0_nm must be positive")

    @property
    def koff_1h(self) -> float:
        return self.kon_nm_h * self.kd_nm

    def kd_from_ic50_um(self, ic50_um: float) -> Target:
        """Return a copy whose affinity derives from a measured IC50 (µM)."""
        kd = ic50_um * _NM_PER_MICRO_M
        return Target(
            name=self.name,
            kd_nm=kd,
            kon_nm_h=self.kon_nm_h,
            r0_nm=self.r0_nm,
            rho_h=self.rho_h,
            kint_h=self.kint_h,
            low_confidence=False,
            reference=self.reference,
        )


def safety_panel() -> tuple[Target, ...]:
    """The baseline off-target safety panel (doc/05 2.2).

    Affinities are class-typical pharmacology priors (IC50 median proxies
    converted via ``kd_from_ic50_um``); abundance values are rough tissue
    expression priors.  All panel entries are ``low_confidence=True`` until a
    measured value is supplied via ``Target.kd_from_ic50_um`` or a re-bound
    ``Target``.
    """

    def site(name: str, ic50_um: float, r0_nm: float = 1.0) -> Target:
        return Target(
            name=name,
            kd_nm=ic50_um * _NM_PER_MICRO_M,
            kon_nm_h=0.5,
            r0_nm=r0_nm,
            rho_h=0.05,
            kint_h=0.05,
            low_confidence=True,
            reference="class-median pharmacology prior (flag: low confidence)",
        )

    return (
        Target(
            name="hERG (Kv11.1)",
            kd_nm=2.0,
            r0_nm=0.05,
            reference="typical low-nM hERG blocker affinity (e.g. dofetilide)",
        ),
        site("CYP3A4 inhibition", 12.0, r0_nm=20.0),
        site("CYP2D6 inhibition", 15.0, r0_nm=10.0),
        site("CYP2C9 inhibition", 15.0, r0_nm=10.0),
        site("BSEP (cholestasis)", 90.0, r0_nm=2.0),
        site("OATP1B1", 8.0, r0_nm=2.0),
        site("P-gp (MDR1)", 25.0, r0_nm=1.0),
        site("MRP3", 60.0, r0_nm=1.0),
        site("MRP4", 60.0, r0_nm=1.0),
        site("Mitochondrial complex I", 2.0, r0_nm=10.0),
        site("Mitochondrial complex II", 12.0, r0_nm=10.0),
        site("Mitochondrial complex III", 6.0, r0_nm=10.0),
        site("Mitochondrial complex IV", 0.5, r0_nm=10.0),
        site("Mitochondrial pyruvate carrier", 100.0, r0_nm=10.0),
        Target(
            name="Glucocorticoid receptor",
            kd_nm=30.0,
            r0_nm=0.1,
            reference="endocrine nuclear-receptor class prior",
        ),
        Target(
            name="Estrogen receptor",
            kd_nm=20.0,
            r0_nm=0.05,
            reference="endocrine nuclear-receptor class prior",
        ),
        Target(
            name="Androgen receptor",
            kd_nm=20.0,
            r0_nm=0.05,
            reference="endocrine nuclear-receptor class prior",
        ),
    )


__all__ = ["Target", "safety_panel"]
