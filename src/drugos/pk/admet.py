"""ADMET-AI (GNN) predictor wrapper.

Wraps the admet-ai package (Swanson et al. 2020; Livermore) whose output
endpoints are mapped onto the fields the DrugOS pipeline consumes. The heavy
torch model is loaded lazily and kept for reuse.
"""

from __future__ import annotations

import io
import logging
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass, field
from dataclasses import fields as dataclass_fields
from typing import TYPE_CHECKING, Any, Protocol

if TYPE_CHECKING:
    from admet_ai import ADMETModel


class _Loc(Protocol):
    def __getitem__(self, i: int) -> Any: ...


class _Series(Protocol):
    @property
    def iloc(self) -> _Loc: ...


class _Df(Protocol):
    def __getitem__(self, key: str) -> _Series: ...


logging.getLogger("pytorch_lightning").setLevel(logging.ERROR)
logging.getLogger("lightning").setLevel(logging.ERROR)


@dataclass(slots=True)
class AdmetOutput:
    """Canonical ADMET predictions for one molecule.

    ``fup_plasma`` is the fraction unbound in plasma (0..1) derived from the
    predicted plasma protein binding; ``log_s`` aqueous solubility (log mol/L);
    ``cl_int_hep`` / ``cl_int_mic`` intrinsic clearance (mL/min/kg);
    ``half_life_h`` predicted half-life (h); ``vdss_log_l_kg`` log10 steady-state
    volume (L/kg).  Probability endpoints (0..1): hERG, DILI, AMES, BBB
    penetration, P-glycoprotein substrate and the CYP inhibition panel.
    """

    smiles: str
    fup_plasma: float | None = None
    log_s: float | None = None
    cl_int_hep_ml_min_kg: float | None = None
    cl_int_mic_ml_min_kg: float | None = None
    half_life_h: float | None = None
    vdss_log_l_kg: float | None = None
    hERG: float | None = None
    DILI: float | None = None
    AMES: float | None = None
    BBB: float | None = None
    Pgp: float | None = None
    HIA: float | None = None
    bioavailable_Ma: float | None = None
    caco2_log: float | None = None
    log_p_pred: float | None = None
    cyp2d6_inhibitor: float | None = None
    cyp3a4_inhibitor: float | None = None
    cyp2c9_inhibitor: float | None = None
    cyp1a2_inhibitor: float | None = None
    cyp2c19_inhibitor: float | None = None
    raw: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, float | None]:
        return {f.name: getattr(self, f.name) for f in dataclass_fields(self) if f.name != "raw"}


class ADMETPredictor:
    """Lazy-wrapping predictor around ``admet_ai.ADMETModel``."""

    def __init__(self) -> None:
        self._model: ADMETModel | None = None

    @property
    def model(self) -> ADMETModel:
        if self._model is None:
            from admet_ai import ADMETModel

            self._model = ADMETModel()
        return self._model

    def predict(self, smiles: str | list[str]) -> AdmetOutput | list[AdmetOutput]:
        if isinstance(smiles, str):
            return _map_row(self.model, [smiles])[0]
        return _map_row(self.model, list(smiles))


def predict_admet(
    smiles: str | list[str], predictor: ADMETPredictor | None = None
) -> AdmetOutput | list[AdmetOutput]:
    """One-shot ADMET-AI prediction for one or more SMILES strings."""
    own = predictor is None
    pred = predictor or ADMETPredictor()
    try:
        return pred.predict(smiles)
    finally:
        if own:
            pred._model = None


def _map_row(model: ADMETModel, smiles: list[str]) -> list[AdmetOutput]:
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        df = model.predict(smiles)

    outputs: list[AdmetOutput] = []
    for i, smi in enumerate(smiles):
        ppbr = _num(df, i, "PPBR_AZ")
        fup = 1.0 - ppbr / 100.0 if ppbr is not None else None
        outputs.append(
            AdmetOutput(
                smiles=smi,
                fup_plasma=fup,
                log_s=_num(df, i, "Solubility_AqSolDB"),
                cl_int_hep_ml_min_kg=_positive(_num(df, i, "Clearance_Hepatocyte_AZ")),
                cl_int_mic_ml_min_kg=_positive(_num(df, i, "Clearance_Microsome_AZ")),
                half_life_h=_positive(_num(df, i, "Half_Life_Obach")),
                vdss_log_l_kg=_num(df, i, "VDss_Lombardo"),
                hERG=_prob(_num(df, i, "hERG")),
                DILI=_prob(_num(df, i, "DILI")),
                AMES=_prob(_num(df, i, "AMES")),
                BBB=_prob(_num(df, i, "BBB_Martins")),
                Pgp=_prob(_num(df, i, "Pgp_Broccatelli")),
                HIA=_prob(_num(df, i, "HIA_Hou")),
                bioavailable_Ma=_prob(_num(df, i, "Bioavailability_Ma")),
                caco2_log=_num(df, i, "Caco2_Wang"),
                log_p_pred=_num(df, i, "Lipophilicity_AstraZeneca"),
                cyp2d6_inhibitor=_prob(_num(df, i, "CYP2D6_Veith")),
                cyp3a4_inhibitor=_prob(_num(df, i, "CYP3A4_Veith")),
                cyp2c9_inhibitor=_prob(_num(df, i, "CYP2C9_Veith")),
                cyp1a2_inhibitor=_prob(_num(df, i, "CYP1A2_Veith")),
                cyp2c19_inhibitor=_prob(_num(df, i, "CYP2C19_Veith")),
            )
        )
    return outputs


def _num(df: _Df, i: int, key: str) -> float | None:
    """Read a numeric cell; NaN/None -> None (prediction absent).

    A missing *column* (``KeyError``) and a non-numeric value (``TypeError``)
    propagate as explicit errors: the ADMET-AI schema is pinned and a schema
    drift must never silently degrade the input vector.
    """
    value = df[key].iloc[i]
    if value is None or value != value:  # NaN is the only value != itself
        return None
    return float(value)


def _positive(value: float | None) -> float | None:
    if value is None or not (value > 0):
        return None
    return value


def _prob(value: float | None) -> float | None:
    if value is None:
        return None
    return min(max(value, 0.0), 1.0)
