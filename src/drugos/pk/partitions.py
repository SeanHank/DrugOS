"""Rodgers & Rowland tissue-to-plasma partition coefficients.

Implements the equations as published (Rodgers, Leahy & Rowland, J Pharm Sci
2005;94:1259-76 [bases] and 94:1627-41 [acids/neutrals/zwitterions]; corrected
in the 2007 erratum, J Pharm Sci 96:3151-52). Plain-text renderings of the
same equations appear in Korzekwa et al., DMD 2019;47:152-64 and Coutinho et
al. 2024.

Notation follows the 2007 erratum: X, Y are pure 10^raises of ionization
ratios (X = 10^(pKa - pH_IW), Y = 10^(pKa - pH_P) for bases; reversed for
acids). A rendering in which the protein-binding term of the acid/neutral
equation is multiplied by X/(1+Y) or vanishes for neutral drugs reflects a
transcription divergence found in some secondary implementations (mrgsolve,
Utsey et al. 2020) and is NOT used here.

Constants: pH(plasma)=7.4, pH(intracellular water)=7.0, pH(blood cell)=7.22;
plasma neutral-lipid fraction 0.0023, neutral-phospholipid 0.0013. For adipose
tissue the neutral-lipid and neutral-phospholipid terms use the vegetable
oil/water partition coefficient: logP_vo = 1.115*logP - 1.35.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from drugos.inputs.models import Molecule, RrClass
from drugos.pk.physiology import (
    BLOOD_CELL_PH,
    INTRACELLULAR_PH,
    PLASMA_F_NL,
    PLASMA_F_NP,
    PLASMA_PH,
    RBC_AP_MG_G,
    RBC_F_IW,
    RBC_F_NL,
    RBC_F_NP,
    TISSUE_COMPOSITION,
    TissueComposition,
)

_PH_IW = INTRACELLULAR_PH
_PH_P = PLASMA_PH
_PH_BC = BLOOD_CELL_PH


@dataclass(slots=True)
class RrPartition:
    """Computed partition coefficients for all tissues.

    ``kpu`` is the tissue-to-unbound-plasma concentration ratio; ``kp`` is the
    tissue-to-plasma (total) ratio: ``kp = kpu * fup``.
    """

    kpu: dict[str, float]
    kp: dict[str, float]
    ka_ap: float | None = None
    kpu_bc: float | None = None
    ka_pr: float | None = None
    rr_class: RrClass = RrClass.NEUTRAL
    extras: dict[str, float] = field(default_factory=dict)


def _ionization(base_pka: float, ph: float) -> float:
    """Ionization ratio 10^(pKa-pH) for a base group."""
    return math.pow(10.0, base_pka - ph)


def _acid_ionization(acid_pka: float, ph: float) -> float:
    """Ionization ratio 10^(pH-pKa) for an acid group."""
    return math.pow(10.0, ph - acid_pka)


def _xy(pka_acids: list[float], pka_bases: list[float], ph: float) -> tuple[float, float]:
    """Summed X (at ``ph``) and Y (at plasma pH) ionization ratios."""
    x = sum(_acid_ionization(pka, ph) for pka in pka_acids) + sum(
        _ionization(pka, ph) for pka in pka_bases
    )
    y = sum(_acid_ionization(pka, _PH_P) for pka in pka_acids) + sum(
        _ionization(pka, _PH_P) for pka in pka_bases
    )
    return x, y


def kpu_blood_cell(bp: float, fup: float, hematocrit: float) -> float:
    """Unbound blood-cell-to-plasma-water ratio (Korzekwa et al. eq. 7)."""
    return (bp + hematocrit - 1.0) / (hematocrit * fup)


def ka_ap_strong_base(
    log_p: float,
    pka_base: float,
    fup: float,
    bp: float,
    hematocrit: float,
) -> float:
    """Calibrate the acidic-phospholipid association constant from RBC data."""
    if not 0.0 < fup <= 1.0:
        raise ValueError("fup must be in (0, 1]")
    if bp <= 0.0:
        raise ValueError("blood-to-plasma ratio must be positive")

    p_oct = math.pow(10.0, log_p)
    y = _ionization(pka_base, _PH_P)
    z = _ionization(pka_base, _PH_BC)
    kpu_bc = kpu_blood_cell(bp, fup, hematocrit)

    water_term = (1.0 + z) / (1.0 + y) * RBC_F_IW
    lipid_term = (p_oct * RBC_F_NL + (0.3 * p_oct + 0.7) * RBC_F_NP) / (1.0 + y)
    ka_ap = (kpu_bc - water_term - lipid_term) * (1.0 + y) / (RBC_AP_MG_G * z)
    return ka_ap


def _plasma_binding_excess(fup: float, p_oct: float, y: float) -> float:
    """Plasma binding capacity in excess of the tissue lipid/water phase."""
    plasma_lipids = (p_oct * PLASMA_F_NL + (0.3 * p_oct + 0.7) * PLASMA_F_NP) / (1.0 + y)
    return (1.0 / fup - 1.0) - plasma_lipids


def _pr_ratio(
    comp: TissueComposition,
    mode: str,
    ionizable: bool,
) -> float:
    """Tissue:plasma protein-capacity ratio ([PR]_t / [PR]_p)."""
    if mode == "albumin":
        return comp.ar
    if mode == "lipoprotein":
        return comp.lr
    return comp.ar if ionizable else comp.lr


def rodgers_rowland_partition(
    fup: float,
    log_p: float,
    pka_acids: list[float] | None = None,
    pka_bases: list[float] | None = None,
    *,
    bp: float | None = None,
    hematocrit: float = 0.45,
    rr_class: RrClass | None = None,
    protein_ratio: str = "auto",
) -> RrPartition:
    """Compute Rodgers-Rowland partition coefficients for a compound.

    Parameters
    ----------
    fup : fraction unbound in plasma, in (0, 1].
    log_p : octanol/water logP of the neutral species.
    pka_acids, pka_bases : macro-pKa lists (may be empty).
    bp : blood-to-plasma ratio; REQUIRED for strong bases (calibrates Ka_AP).
    hematocrit : blood cell fraction (default 0.45).
    protein_ratio : "auto" (albumin if the compound is ionizable, otherwise
        lipoprotein), "albumin", or "lipoprotein". Used by the acid/neutral
        equation only.
    """
    pka_acids = sorted(pka_acids or [])
    pka_bases = sorted(pka_bases or [], reverse=True)
    if not 0.0 < fup <= 1.0:
        raise ValueError(f"fup must be in (0, 1]: got {fup}")

    if rr_class is None:
        from drugos.inputs.parse_structure import _rr_class

        rr_class = _rr_class(pka_acids, pka_bases)

    p_oct = math.pow(10.0, log_p)
    p_vo = math.pow(10.0, 1.115 * log_p - 1.35)
    is_strong_base = rr_class is RrClass.STRONG_BASE

    if is_strong_base:
        if bp is None:
            raise ValueError("blood-to-plasma ratio (bp) is required for strong bases")
        pka_ref = pka_bases[0]
        x, y = _xy([], pka_bases, _PH_IW)
        ka_ap: float | None = ka_ap_strong_base(log_p, pka_ref, fup, bp, hematocrit)
        kpu_bc: float | None = kpu_blood_cell(bp, fup, hematocrit)
        ka_pr: float | None = None
    else:
        x, y = _xy(pka_acids, pka_bases, _PH_IW)
        ka_ap = None
        kpu_bc = None
        ka_pr = _plasma_binding_excess(fup, p_oct, y)

    ionizable = bool(pka_acids or pka_bases)
    kpu: dict[str, float] = {}
    for name, comp in TISSUE_COMPOSITION.items():
        p = p_vo if name == "adipose" else p_oct
        value = comp.f_ew + (1.0 + x) / (1.0 + y) * comp.f_iw
        value += (p * comp.f_nl + (0.3 * p + 0.7) * comp.f_np) / (1.0 + y)
        if is_strong_base and ka_ap is not None:
            value += (ka_ap * comp.ap_mg_g * x) / (1.0 + y)
        else:
            assert ka_pr is not None
            value += ka_pr * _pr_ratio(comp, protein_ratio, ionizable)
        kpu[name] = value

    kp = {name: kpu[name] * fup for name in kpu}
    return RrPartition(
        kpu=kpu,
        kp=kp,
        ka_ap=ka_ap,
        kpu_bc=kpu_bc,
        ka_pr=ka_pr,
        rr_class=rr_class,
    )


def partition_from_molecule(
    mol: Molecule,
    fup: float,
    *,
    bp: float | None = None,
    hematocrit: float = 0.45,
    protein_ratio: str = "auto",
) -> RrPartition:
    """Convenience wrapper accepting a :class:`Molecule` directly."""
    if mol.log_p is None:
        raise ValueError("Molecule.log_p is required for partition coefficient computation")
    return rodgers_rowland_partition(
        fup=fup,
        log_p=mol.log_p,
        pka_acids=mol.pka_acids,
        pka_bases=mol.pka_bases,
        bp=bp,
        hematocrit=hematocrit,
        rr_class=mol.rr_class,
        protein_ratio=protein_ratio,
    )
