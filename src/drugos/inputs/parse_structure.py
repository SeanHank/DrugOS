"""Chemical structure parsing and physicochemical descriptor computation.

The molecule may be provided as SMILES, InChI or an SDF file path. RDKit
computes classical descriptors; macro-pKa values are estimated with a curated
SMARTS-based heuristic (no external predictor dependency).  Any value may be
overridden by the caller (e.g. measured pKa).

The heuristic macro-pKa table is intentionally small and documented; it is a
first-pass approximation to be superseded by measured or ML-based pKa where
available.
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable
from typing import Any

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from drugos.inputs.models import IonClass, Molecule, RrClass


# Typed shims over upstream RDKit APIs that carry no (or partial) type
# information.  These keep the analysis strictly typed per
# doc/09-quality-gate.md (G2) without error suppressions.
def _rdkit_fn(module: object, name: str) -> Any:
    return getattr(module, name)


_MolFromInchi: Callable[[str], Chem.Mol | None] = _rdkit_fn(Chem, "MolFromInchi")
_MolToInchi: Callable[[Chem.Mol], str] = _rdkit_fn(Chem, "MolToInchi")
_InchiToInchiKey: Callable[[str], str] = _rdkit_fn(Chem, "InchiToInchiKey")
_MolLogP: Callable[[Chem.Mol], float] = _rdkit_fn(Crippen, "MolLogP")
_MolWt: Callable[[Chem.Mol], float] = _rdkit_fn(Descriptors, "MolWt")
_NumHDonors: Callable[[Chem.Mol], int] = _rdkit_fn(Lipinski, "NumHDonors")
_NumHAcceptors: Callable[[Chem.Mol], int] = _rdkit_fn(Lipinski, "NumHAcceptors")
_NumRotatableBonds: Callable[[Chem.Mol], int] = _rdkit_fn(Lipinski, "NumRotatableBonds")

ACIDIC_SMARTS: list[tuple[str, float]] = [
    ("[CX3](=O)[OX2H1]", 4.2),  # carboxylic acid
    ("[SX2H1]", 10.5),  # thiol
    ("[cX3][OX2H1]", 10.0),  # phenol / aromatic OH
    ("[NX3;H1][SX4](=O)(=O)", 10.1),  # sulfonamide N-H
]

BASIC_SMARTS: list[tuple[str, float]] = [
    ("[NX3][CX3](=[NX2])[NX3]", 13.5),  # guanidine
    ("[NX3][CX3]=[NX2]", 12.0),  # amidine
    ("[CX4][NX3;H2,H1,H0;!$(NC=O);!$(N#C)]", 10.3),  # aliphatic amine
    ("[cX3][NX3;H1,H2;!$(NC=O)]", 4.6),  # aniline
    ("[nX3;H0]", 5.2),  # pyridine-like ring N
]


class StructureParseError(ValueError):
    """Raised when a chemical structure cannot be interpreted."""


def _as_mol(structure: str) -> Chem.Mol:
    text = structure.strip()
    low = text.lower()
    mol: Chem.Mol | None
    if low.endswith(".sdf"):
        path = os.path.abspath(text)
        if not os.path.exists(path):
            raise StructureParseError(f"SDF file not found: {path}")
        supplier = Chem.SDMolSupplier(path, removeHs=False)
        mol = supplier[0]
    elif text.startswith("InChI="):
        mol = _MolFromInchi(text)
    else:
        mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise StructureParseError(f"Could not parse structure: {structure!r}")
    return mol


def _macro_pka(mol: Chem.Mol, table: list[tuple[str, float]]) -> list[float]:
    values: list[float] = []
    for smarts, pka in table:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            continue
        if mol.HasSubstructMatch(pat):
            values.append(pka)
    return sorted(values)


def _pka_acids(mol: Chem.Mol) -> list[float]:
    """Macro-acid pKas in ascending order (most acidic first)."""
    return _macro_pka(mol, ACIDIC_SMARTS)


def _pka_bases(mol: Chem.Mol) -> list[float]:
    """Macro-base pKas in descending order (most basic first)."""
    return sorted(_macro_pka(mol, BASIC_SMARTS), reverse=True)


def _ion_class(pka_acids: list[float], pka_bases: list[float]) -> IonClass:
    strong_acid = min(pka_acids) if pka_acids else None
    strong_base = max(pka_bases) if pka_bases else None
    f_anion = 1.0 / (1.0 + 10 ** (strong_acid - 7.4)) if strong_acid is not None else 0.0
    f_cation = 1.0 / (1.0 + 10 ** (7.4 - strong_base)) if strong_base is not None else 0.0
    if f_anion > 0.1 and f_cation > 0.1:
        return IonClass.ZWITTERION
    if f_cation > f_anion and f_cation > 0.1:
        return IonClass.CATION
    if f_anion > 0.1:
        return IonClass.ANION
    return IonClass.NEUTRAL


def _rr_class(pka_acids: list[float], pka_bases: list[float]) -> RrClass:
    strong_base = max(pka_bases) if pka_bases else None
    if strong_base is not None and strong_base > 7.0:
        return RrClass.STRONG_BASE
    if pka_acids and pka_bases:
        return RrClass.ZWITTERION
    if pka_acids:
        return RrClass.ACID
    if pka_bases:
        return RrClass.WEAK_BASE
    return RrClass.NEUTRAL


def _fraction_neutral(pka_acids: list[float], pka_bases: list[float], ph: float) -> float:
    """Product of per-group neutral fractions (independent-group approximation)."""
    f_acid = math.prod(1.0 - 1.0 / (1.0 + 10 ** (pka - ph)) for pka in pka_acids)
    f_base = math.prod(1.0 - 1.0 / (1.0 + 10 ** (ph - pka)) for pka in pka_bases)
    return f_acid * f_base


def parse_structure(structure: str, name: str | None = None) -> Molecule:
    """Parse a structure (SMILES / InChI / SDF path) into a :class:`Molecule`."""
    raw = _as_mol(structure)
    mol = Chem.Mol(raw)  # explicit copy before hydrogen handling
    mol = Chem.AddHs(mol)
    mol = Chem.RemoveHs(mol)

    canonical_noniso = Chem.MolToSmiles(mol, isomericSmiles=False)
    canonical_iso = Chem.MolToSmiles(mol, isomericSmiles=True)

    pka_acids = _pka_acids(mol)
    pka_bases = _pka_bases(mol)
    log_p = _MolLogP(mol)
    f_neut = _fraction_neutral(pka_acids, pka_bases, 7.4)
    log_d = log_p + (math.log10(f_neut) if f_neut > 0 else 0.0)

    heavy = mol.GetNumHeavyAtoms()
    formula = rdMolDescriptors.CalcMolFormula(mol)
    inchi = _MolToInchi(mol)
    inchikey = _InchiToInchiKey(inchi) if inchi else None

    return Molecule(
        name=name,
        canonical_smiles=canonical_noniso,
        canonical_isomeric_smiles=canonical_iso,
        inchi=inchi,
        inchikey=inchikey,
        molecular_formula=formula,
        mw=_MolWt(mol),
        log_p=log_p,
        log_d_7_4=log_d,
        tpsa=rdMolDescriptors.CalcTPSA(mol),
        hbd=_NumHDonors(mol),
        hba=_NumHAcceptors(mol),
        rotatable_bonds=_NumRotatableBonds(mol),
        aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
        formal_charge=Chem.GetFormalCharge(mol),
        heavy_atoms=heavy,
        pka_acids=pka_acids,
        pka_bases=pka_bases,
        ion_class=_ion_class(pka_acids, pka_bases),
        rr_class=_rr_class(pka_acids, pka_bases),
    )


def fraction_neutral(mol: Molecule, ph: float) -> float:
    """Fraction of the molecule in the neutral form at ``ph``."""
    return _fraction_neutral(mol.pka_acids, mol.pka_bases, ph)
