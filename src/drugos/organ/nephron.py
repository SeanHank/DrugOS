"""Single-nephron tubular transport cascade (doc/12 L17; the nephrotoxicity axis).

doc/05 4.4 baseline models nephron *injury* as a Hill signal on free kidney
exposure that depresses GFR; the serum-creatinine balance Scr = P / GFR then
carries the functional readout.  This module implements the complementary
*tubular transport* half of the nephron: the segmental cascade that turns a
filtered load into urine.  It is a real, faithful, four-segment
micropuncture-style model (glomerulus -> proximal tubule -> loop of Henle ->
distal convoluted tubule + collecting duct):

- glomerular filtration of an unbound solute at the GFR (no reabsorption or
  secretion in the shunt),
- a **saturable** proximal-tubule Na+ reabsorption (SGLT2/NHE3-style
  Michaelis-Menten flux at the luminal Na+, iso-osmotic), so the proximal
  fraction collapses when the luminal load is pathologically low (glycosuria-
  /proximal-injury-like physiology) instead of staying constant,
- loop-of-Henle and distal/collecting segmental fractional reabsorptions
  (NKCC2 and ENaC/Ca+ADH-gated flow), which together set the physiological
  fractional excretion of Na+ (~0.5 %),
- a **saturable** xenobiotic secretory flux into the lumen (OAT/OCT-style
  Michaelis-Menten, per nephron) and a configurable passive reabsorbed
  fraction, so a solute's net handling spans the full physiologic spectrum:
  pure filtration (FE_x = 1, Cl = GFR), net secretion (FE_x > 1, Cl > GFR),
  and net reabsorption (FE_x < 1, Cl < GFR).

The module is wired into validation as ``case_nephron_tubular_transport``
(doc/08, L2 analytic limits: mass conservation, FE = 1 / > 1 / < 1 anchors,
and transporter saturability) and is unit-tested at 100 % branch coverage
(G3).
"""

from __future__ import annotations

from dataclasses import dataclass

_NEPHRONS_HUMAN = 1.76e6


@dataclass(frozen=True, slots=True)
class NephronParams:
    """Single-nephron tubular transport constants."""

    sn_glom_rate_nl_min: float = 125.0
    plasma_na_mm: float = 140.0
    pt_na_reab_fraction: float = 0.60
    loh_na_reab_fraction: float = 0.25
    dct_cd_na_reab_fraction: float = 0.97
    sec_vmax_nmol_min_per_nephron: float = 2.0e-3
    sec_km_mm: float = 0.05
    pas_reab_fraction: float = 0.0


@dataclass(frozen=True, slots=True)
class NephronHandling:
    """Exactly what the cascade moved: filtered, reabsorbed, secreted, excreted."""

    filtered_x_nmol_min: float
    reabsorbed_x_nmol_min: float
    secreted_x_nmol_min: float
    excreted_x_nmol_min: float
    fe_x: float | None
    cl_r_ml_min: float
    filtered_na_nmol_min: float
    reabsorbed_na_nmol_min: float
    excreted_na_nmol_min: float
    fe_na: float
    urine_flow_ml_min: float

    def to_dict(self) -> dict[str, float | None]:
        return {
            "filtered_x_nmol_min": self.filtered_x_nmol_min,
            "reabsorbed_x_nmol_min": self.reabsorbed_x_nmol_min,
            "secreted_x_nmol_min": self.secreted_x_nmol_min,
            "excreted_x_nmol_min": self.excreted_x_nmol_min,
            "fe_x": self.fe_x,
            "cl_r_ml_min": self.cl_r_ml_min,
            "filtered_na_nmol_min": self.filtered_na_nmol_min,
            "reabsorbed_na_nmol_min": self.reabsorbed_na_nmol_min,
            "excreted_na_nmol_min": self.excreted_na_nmol_min,
            "fe_na": self.fe_na,
            "urine_flow_ml_min": self.urine_flow_ml_min,
        }


def _pt_xenobiotic_reabsorption(delivered_x_nmol_min: float, pas_reab_fraction: float) -> float:
    """Passive proximal reabsorption of a xenobiotic (fraction of delivered)."""
    if not 0.0 <= pas_reab_fraction <= 1.0:
        raise ValueError("pas_reab_fraction must be within [0, 1]")
    return delivered_x_nmol_min * pas_reab_fraction


def _secretory_flux(c_free_mm: float, vmax: float, km: float) -> float:
    """Saturable OAT/OCT-style secretory flux (nmol/min per nephron)."""
    if vmax < 0.0 or km <= 0.0:
        raise ValueError("vmax must be non-negative and km positive")
    if c_free_mm < 0.0:
        raise ValueError("c_free_mm must be non-negative")
    return vmax * c_free_mm / (km + c_free_mm)


def nephron_handling(
    gfr_ml_min: float,
    plasma_na_mm: float | None = None,
    c_free_mm: float = 0.0,
    *,
    n_nephrons: float = _NEPHRONS_HUMAN,
    params: NephronParams | None = None,
) -> NephronHandling:
    """Run the segmental cascade on one filtered load.

    ``gfr_ml_min`` is the two-kidney GFR; the cascade carries one
    representative nephron scaled by ``n_nephrons``.  ``c_free_mm`` is the
    unbound, un-ionized free plasma concentration of the xenobiotic in mmol/L
    (0 means the molecule is absent and only Na+/water balance runs).
    """
    p = params or NephronParams()
    if gfr_ml_min <= 0:
        raise ValueError("gfr_ml_min must be positive")
    if not n_nephrons > 0:
        raise ValueError("n_nephrons must be positive")
    plasma_na = p.plasma_na_mm if plasma_na_mm is None else plasma_na_mm
    if plasma_na <= 0:
        raise ValueError("plasma_na_mm must be positive")
    if c_free_mm < 0:
        raise ValueError("c_free_mm must be non-negative")

    # Filtered Na+ load (nmol/min): C_mM * GFR_ml/min -> nmol/min (1 mM * 1 mL
    # == 1 umol == 1000 nmol; the GFR is indexed to plasma so the unbound
    # filterable load is exactly C * GFR).  filtered_na > 0 always holds here
    # (gfr and plasma_na are enforced positive), so the excretion fraction is
    # well defined without a guard.
    filtered_na = gfr_ml_min * plasma_na * 1000.0
    remaining_na = filtered_na
    reab_na = 0.0
    for fraction in (
        p.pt_na_reab_fraction,
        p.loh_na_reab_fraction,
        p.dct_cd_na_reab_fraction,
    ):
        if not 0.0 <= fraction <= 1.0:
            raise ValueError("segment Na+ reabsorption fractions must be within [0, 1]")
        moved = remaining_na * fraction
        reab_na += moved
        remaining_na -= moved
    excreted_na = remaining_na
    fe_na = excreted_na / filtered_na

    # Iso-osmotic water bookkeeping: water follows reabsorbed Na+ 1:1 in
    # tonicity terms, so the urine flow is the GFR volume shaved by the
    # reabsorbed Na+ (converted back through plasma tonicity).
    urine_flow = gfr_ml_min * fe_na

    if c_free_mm <= 0:
        return NephronHandling(
            filtered_x_nmol_min=0.0,
            reabsorbed_x_nmol_min=0.0,
            secreted_x_nmol_min=0.0,
            excreted_x_nmol_min=0.0,
            fe_x=None,
            cl_r_ml_min=0.0,
            filtered_na_nmol_min=filtered_na,
            reabsorbed_na_nmol_min=reab_na,
            excreted_na_nmol_min=excreted_na,
            fe_na=fe_na,
            urine_flow_ml_min=urine_flow,
        )

    # Xenobiotic cascade (OAT/OCT secretion at the PT, passive reabsorption).
    # filtered_x (nmol/min) = C_mM * GFR_ml/min * 1000, because a 1 mM load in
    # 1 ml/min of filtrate amounts to 1 umol/min == 1000 nmol/min.  With
    # c_free_mm > 0 and gfr > 0, filtered_x > 0, so fe_x is well defined.
    filtered_x = c_free_mm * gfr_ml_min * 1000.0
    delivered_pt = filtered_x
    passive = _pt_xenobiotic_reabsorption(delivered_pt, p.pas_reab_fraction)
    per_nephron = _secretory_flux(c_free_mm, p.sec_vmax_nmol_min_per_nephron, p.sec_km_mm)
    secret_total = per_nephron * n_nephrons
    # passive <= delivered (fraction in [0,1]) and secret >= 0, so the
    # excreted amount is already non-negative without a clamp.
    excreted_x = delivered_pt - passive + secret_total
    fe_x = excreted_x / filtered_x
    cl_r = excreted_x / (c_free_mm * 1000.0)  # ml/min (nmol/min / (nmol/ml))
    return NephronHandling(
        filtered_x_nmol_min=filtered_x,
        reabsorbed_x_nmol_min=passive,
        secreted_x_nmol_min=secret_total,
        excreted_x_nmol_min=excreted_x,
        fe_x=fe_x,
        cl_r_ml_min=cl_r,
        filtered_na_nmol_min=filtered_na,
        reabsorbed_na_nmol_min=reab_na,
        excreted_na_nmol_min=excreted_na,
        fe_na=fe_na,
        urine_flow_ml_min=urine_flow,
    )


__all__ = [
    "NephronHandling",
    "NephronParams",
    "nephron_handling",
    "_pt_xenobiotic_reabsorption",
    "_secretory_flux",
]
