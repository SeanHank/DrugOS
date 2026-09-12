"""Predictive-regime classification and reliability disclosure (DISCLAIMER §2).

The disclaimer's honesty clause — *"predictions for novel molecules,
extrapolated doses, or off-label routes are the least reliable; disagreement
with empirical data is the expected state of a mechanistic model, not a bug"* —
is implemented here as an executable per-run contract, not prose (doc/07
D25/D26).

Every run is classified into a predictive regime from three evidence axes:

- **chemistry** — fully measured PK (unbound fraction + hepatic + renal
  clearance), a validated benchmark/clinical scaffold, partially measured
  PK/absorption, or machine-predicted only (novel);
- **dose window** — inside or outside the validated reference window
  (one geometric fold each way of the scaffold reference dose);
- **route** — on-label (the scaffold's validated route) or off-label.

The weakest axis sets the regime, so a run never reports more confidence than
its least supported axis.  Each regime carries a recommended
parameter-ensemble CV (``band_cv``): novel-molecule runs sweep parameters far
wider than measured, on-label runs, so reported uncertainty widens exactly
where confidence is lowest.

``Measurements`` is the full-true-parameters path: any user-supplied measured
scalar (permeability/absorption, clearances, saturable Michaelis–Menten
kinetics, potencies) is applied by ``apply_measurements`` as a true-parameter
override that replaces whatever auto-anchor the pipeline synthesized — it is
never overwritten back by ``_engage_full_fidelity``.  A run carrying the full
PK set is classified into the most reliable regime and its trust record shows
the measured basis.  The goal is: *given all true parameters of a drug and a
human body, predict the effects they actually produce* — this module is the
contract that says a run used those measured parameters instead of its
synthesized stand-in.

``EmpiricalObservations`` + ``agreement_rows`` make the second half of the
clause operational: the user may attach observed human PK/PD values, and the
run reports back an explicit fold-error diagnostic in the trust record
(``empirical_agreement``) instead of silently absorbing or hiding the
deviation.  Disagreement is the expected state of a mechanistic model; it is
surfaced as data, never assigned blame and never suppressed (doc/08 risk
register, doc/12 §7).

This module deliberately imports nothing from ``drugos.pipeline`` (avoiding a
cycle); callers supply benchmark-scaffold context through the ``scaffold``
argument.
"""

from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from drugos.inputs.models import Route

_DISAGREEMENT = (
    "disagreement with empirical data is the expected state of a mechanistic "
    "model, not a bug (DISCLAIMER §2)"
)


class PredictiveRegime(StrEnum):
    """Evidence regime of a prediction, weakest to strongest."""

    NOVEL_MOLECULE = "novel_molecule"
    PARTIAL_EVIDENCE = "partial_evidence"
    VALIDATED_OFFLABEL_ROUTE = "validated_offlabel_route"
    VALIDATED_EXTRAPOLATED_DOSE = "validated_extrapolated_dose"
    VALIDATED_IN_RANGE_ON_LABEL = "validated_in_range_on_label"
    MEASURED_IN_RANGE_ON_LABEL = "measured_in_range_on_label"


@dataclass(frozen=True, slots=True)
class RegimeMeta:
    """Human-readable definition of one predictive regime."""

    label: str
    reliability: str
    band_cv: float
    disclaimer: str

    def __post_init__(self) -> None:
        if not (0.0 < self.band_cv <= 1.0):
            raise ValueError("band_cv must be in (0, 1]")
        if self.reliability not in ("highest", "high", "medium", "medium-low", "low", "lowest"):
            raise ValueError(f"unknown reliability label {self.reliability!r}")


REGIMES: dict[PredictiveRegime, RegimeMeta] = {
    PredictiveRegime.NOVEL_MOLECULE: RegimeMeta(
        label="novel molecule — chemistry, PK and pharmacology all machine-predicted",
        reliability="lowest",
        band_cv=0.60,
        disclaimer=(
            "no benchmark/clinical anchor and no measured parameter: this run is the "
            f"least reliable prediction the platform produces; {_DISAGREEMENT}."
        ),
    ),
    PredictiveRegime.PARTIAL_EVIDENCE: RegimeMeta(
        label="partially characterized molecule (some measured PK/absorption)",
        reliability="low",
        band_cv=0.45,
        disclaimer=(
            "only part of the PK/enabled-absorption picture is measured; the rest is "
            f"machine-anchored — predictions remain fragile; {_DISAGREEMENT}."
        ),
    ),
    PredictiveRegime.VALIDATED_OFFLABEL_ROUTE: RegimeMeta(
        label="known chemistry, off-label route",
        reliability="low",
        band_cv=0.50,
        disclaimer=(
            "a route is predicted that the scaffold's validated label does not cover; "
            "absorption/elimination mechanisms are extrapolated into new tissue paths, so "
            f"this is among the least reliable predictions; {_DISAGREEMENT}."
        ),
    ),
    PredictiveRegime.VALIDATED_EXTRAPOLATED_DOSE: RegimeMeta(
        label="known chemistry, dose extrapolated outside the validated window",
        reliability="medium-low",
        band_cv=0.35,
        disclaimer=(
            "the dose lies outside the scaffold's validated reference window; "
            "saturable/clearance mechanisms are extrapolated beyond their anchored range, "
            f"so disagreement with empirical data is expected; {_DISAGREEMENT}."
        ),
    ),
    PredictiveRegime.VALIDATED_IN_RANGE_ON_LABEL: RegimeMeta(
        label="validated molecule, on-label route, in-range dose",
        reliability="medium",
        band_cv=0.20,
        disclaimer=(
            "the run sits inside the scaffold's empirical anchor window (benchmark corpus, "
            f"2x fold-error GMFE allowance); disagreement within the fold allowance is still "
            f"expected, not a bug; {_DISAGREEMENT}."
        ),
    ),
    PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL: RegimeMeta(
        label="measured molecule (full PK), on-label route, in-range dose",
        reliability="high",
        band_cv=0.15,
        disclaimer=(
            "the run is driven by a full measured parameterization of the drug: the closest "
            "the platform comes to 'give all true parameters, predict the true effects'; "
            f"residual disagreement is still disclosed, never blamed; {_DISAGREEMENT}."
        ),
    ),
}

_WINDOW_LO_FRAC = 0.5
_WINDOW_HI_FRAC = 2.0

APPLIED_FIELDS: tuple[str, ...] = (
    "fup",
    "cl_hep_l_h",
    "cl_renal_l_h",
    "cl_sec_l_h",
    "cl_bil_l_h",
    "fa",
    "hepatic_vmax_mg_h",
    "hepatic_km_mg_l",
    "qt_ic50_nm",
    "dili_ic50_nm",
    "dili_immune_ic50_nm",
    "cns_ic50_nm",
    "beta_block_ic50_nm",
)

EMPIRICAL_FIELDS: tuple[str, ...] = (
    "plasma_cmax_mg_l",
    "auc_last_mg_h_l",
    "peak_delta_qtc_ms",
    "peak_alt_uln",
)


class Measurements:
    """User-supplied true-parameter overrides (the measured-data path).

    Any non-``None`` field replaces the pipeline's auto-anchor or prior for
    that scalar (``apply_measurements``).
    """

    __slots__ = APPLIED_FIELDS

    def __init__(
        self,
        *,
        fup: float | None = None,
        cl_hep_l_h: float | None = None,
        cl_renal_l_h: float | None = None,
        cl_sec_l_h: float | None = None,
        cl_bil_l_h: float | None = None,
        fa: float | None = None,
        hepatic_vmax_mg_h: float | None = None,
        hepatic_km_mg_l: float | None = None,
        qt_ic50_nm: float | None = None,
        dili_ic50_nm: float | None = None,
        dili_immune_ic50_nm: float | None = None,
        cns_ic50_nm: float | None = None,
        beta_block_ic50_nm: float | None = None,
    ) -> None:
        """Validate and store measured scalars (units match ``RunSpec``)."""
        for name in (
            "qt_ic50_nm",
            "dili_ic50_nm",
            "dili_immune_ic50_nm",
            "cns_ic50_nm",
            "beta_block_ic50_nm",
        ):
            value = locals()[name]
            if value is not None and value <= 0:
                raise ValueError(f"{name} must be positive; got {value}")
        if fup is not None and not (0.0 < fup <= 1.0):
            raise ValueError(f"fup must be in (0, 1]; got {fup}")
        if fa is not None and not (0.0 < fa <= 1.0):
            raise ValueError(f"fa must be in (0, 1]; got {fa}")
        for name in ("cl_hep_l_h", "cl_renal_l_h", "cl_sec_l_h", "cl_bil_l_h"):
            value = locals()[name]
            if value is not None and (value < 0 or not math.isfinite(value)):
                raise ValueError(f"{name} must be non-negative and finite; got {value}")
        if (hepatic_vmax_mg_h is None) != (hepatic_km_mg_l is None):
            raise ValueError("hepatic_vmax_mg_h and hepatic_km_mg_l must be set together")
        if hepatic_km_mg_l is not None and hepatic_km_mg_l <= 0:
            raise ValueError(f"hepatic_km_mg_l must be positive; got {hepatic_km_mg_l}")
        for name in APPLIED_FIELDS:
            object.__setattr__(self, name, locals()[name])

    def is_pk_complete(self) -> bool:
        """True when the full PK picture (fup + hepatic + renal) is measured."""
        return all(
            getattr(self, name) is not None for name in ("fup", "cl_hep_l_h", "cl_renal_l_h")
        )

    def count_pk_partial(self) -> int:
        """Number of measured PK/absorption scalars supplied."""
        return sum(
            getattr(self, name) is not None for name in ("fup", "cl_hep_l_h", "cl_renal_l_h", "fa")
        )

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Measurements:
        if not isinstance(data, Mapping):
            raise ValueError("measurements must be a JSON object")
        unknown = set(data) - set(APPLIED_FIELDS)
        if unknown:
            raise ValueError(
                "unknown measurement field(s): "
                + ", ".join(sorted(unknown))
                + "; valid fields: "
                + ", ".join(APPLIED_FIELDS)
            )
        return cls(**{name: _coerce(value) for name, value in data.items()})


@dataclass(frozen=True, slots=True)
class EmpiricalObservations:
    """Observed human PK/PD values a user may pin a run against (optional).

    Supplied values are never absorbed into the parameterization — they are
    only compared with the prediction and reported back as the
    ``empirical_agreement`` trust-record diagnostic.
    """

    plasma_cmax_mg_l: float | None = None
    auc_last_mg_h_l: float | None = None
    peak_delta_qtc_ms: float | None = None
    peak_alt_uln: float | None = None

    def __post_init__(self) -> None:
        for name in EMPIRICAL_FIELDS:
            value = getattr(self, name)
            if value is not None and value < 0:
                raise ValueError(f"{name} must be non-negative; got {value}")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EmpiricalObservations:
        if not isinstance(data, Mapping):
            raise ValueError("empirical must be a JSON object")
        unknown = set(data) - set(EMPIRICAL_FIELDS)
        if unknown:
            raise ValueError(
                "unknown empirical observation field(s): "
                + ", ".join(sorted(unknown))
                + "; valid fields: "
                + ", ".join(EMPIRICAL_FIELDS)
            )
        return cls(**{name: _coerce(value) for name, value in data.items()})


def _coerce(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError("field values must be numbers, got boolean")
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"field value must be a number; got {value!r}") from exc


def apply_measurements(spec: Any) -> Any:
    """Return a copy of ``spec`` with every supplied measurement applied.

    Measured scalars are true-parameter overrides: they replace the pipeline's
    auto-anchor or prior, the run's invariants are re-validated, and the trust
    record's ``estimates`` gains a "measured true-parameter override" basis
    line that names each value.  The measurement block is preserved on the
    returned spec so the provenance stays visible.
    """
    meas = getattr(spec, "measurements", None)
    if meas is None:
        return spec
    out = copy.copy(spec)
    basis = dict(getattr(spec, "estimate_basis", None) or {})
    changed = False
    for name in APPLIED_FIELDS:
        value = getattr(meas, name, None)
        if value is None:
            continue
        setattr(out, name, value)
        basis[name] = (
            f"measured true-parameter override: {name} = {value:g} "
            "(user-supplied value replaces the auto-anchor/prior)"
        )
        changed = True
    if changed:
        out.estimate_basis = basis
        out.__post_init__()
    return out


EMPIRICAL_ENDPOINT_LABELS: dict[str, str] = {
    "plasma_cmax_mg_l": "plasma Cmax (mg/L)",
    "auc_last_mg_h_l": "AUC0-t (mg·h/L)",
    "peak_delta_qtc_ms": "peak ΔQTc (ms)",
    "peak_alt_uln": "peak ALT (xULN)",
}


def agreement_rows(
    observations: EmpiricalObservations, predicted: Mapping[str, float]
) -> list[dict[str, Any]]:
    """Fold-error report of observed vs predicted endpoints (DISCLAIMER §2)."""
    rows: list[dict[str, Any]] = []
    for name in EMPIRICAL_FIELDS:
        observed = getattr(observations, name)
        if observed is None:
            continue
        pred = float(predicted.get(name) or float("nan"))
        if math.isfinite(pred) and pred > 0.0:
            fold = observed / pred
            within = 0.5 <= fold <= 2.0
        else:
            fold = None
            within = False
        rows.append(
            {
                "endpoint": name,
                "label": EMPIRICAL_ENDPOINT_LABELS[name],
                "observed": round(float(observed), 5),
                "predicted": round(pred, 5),
                "fold_error": None if fold is None else round(fold, 4),
                "within_2x": within,
            }
        )
    return rows


def classify_regime(spec: Any, *, scaffold: Any | None = None) -> PredictiveRegime:
    """Classify ``spec`` into its predictive regime (weakest axis dominates).

    ``scaffold`` is a benchmark ``_BenchmarkLike`` (reference dose + validated
    route) when the molecule has a clinical anchor; ``None`` when chemistry is
    measured-only or novel.
    """
    meas = getattr(spec, "measurements", None)
    pk_full = meas is not None and meas.is_pk_complete()
    pk_partial = meas is not None and meas.count_pk_partial() >= 1

    if scaffold is not None:
        molecule_axis = 3 if pk_full else 2
    elif pk_full:
        molecule_axis = 3
    elif pk_partial:
        molecule_axis = 1
    else:
        return PredictiveRegime.NOVEL_MOLECULE

    routes = {ev.route for ev in spec.dose_plan.events}
    route_off_label = False
    dose_extrapolated = False
    if scaffold is not None:
        label_route = Route(scaffold.route)
        route_off_label = any(r is not label_route for r in routes)
        total = float(spec.dose_plan.total_dose_mg)
        ref = float(scaffold.dose_mg)
        dose_extrapolated = not (_WINDOW_LO_FRAC * ref <= total <= _WINDOW_HI_FRAC * ref)

    if route_off_label:
        return PredictiveRegime.VALIDATED_OFFLABEL_ROUTE
    if dose_extrapolated:
        return PredictiveRegime.VALIDATED_EXTRAPOLATED_DOSE
    if molecule_axis == 1:
        return PredictiveRegime.PARTIAL_EVIDENCE
    if molecule_axis == 3:
        return PredictiveRegime.MEASURED_IN_RANGE_ON_LABEL
    return PredictiveRegime.VALIDATED_IN_RANGE_ON_LABEL


def _basis(spec: Any, scaffold: Any | None, regime: PredictiveRegime) -> str:
    del regime
    meas = getattr(spec, "measurements", None)
    routes = ", ".join(sorted({ev.route.value for ev in spec.dose_plan.events}))
    name = getattr(spec.molecule, "name", None) or spec.name
    parts: list[str] = []
    if scaffold is not None:
        parts.append(
            f"scaffold {scaffold.name!r} anchors chemistry (reference dose "
            f"{scaffold.dose_mg:g} mg, validated route {scaffold.route})"
        )
    if meas is not None:
        if meas.is_pk_complete():
            parts.append("full measured PK (fup + hepatic + renal clearance)")
        elif meas.count_pk_partial() >= 1:
            parts.append(f"{meas.count_pk_partial()} measured PK/absorption scalar(s)")
        measured_cov = [f for f in APPLIED_FIELDS if getattr(meas, f) is not None]
        if measured_cov:
            parts.append("measured fields: " + ", ".join(measured_cov))
    if scaffold is None and (meas is None or meas.count_pk_partial() == 0):
        parts.append("no benchmark/clinical anchor and no measured parameters")
    if scaffold is not None:
        total = float(spec.dose_plan.total_dose_mg)
        ref = float(scaffold.dose_mg)
        parts.append(
            "dose within validated window"
            if _WINDOW_LO_FRAC * ref <= total <= _WINDOW_HI_FRAC * ref
            else "dose outside validated reference window"
        )
        label = Route(scaffold.route)
        off = sorted({ev.route.value for ev in spec.dose_plan.events if ev.route is not label})
        status = (
            "route on-label (" + scaffold.route + ")"
            if not off
            else "route off-label: " + ", ".join(off)
        )
        parts.append(status)
    else:
        parts.append("no validated dose window declared; route not judged against a label")
    parts.append(f"run routes: {routes}; compound {name!r}")
    return "; ".join(parts)


def reliability_of(spec: Any, *, scaffold: Any | None = None) -> dict[str, Any]:
    """Serialize the reliability disclosure for one run (trust record)."""
    regime = classify_regime(spec, scaffold=scaffold)
    meta = REGIMES[regime]
    return {
        "regime": regime.value,
        "label": meta.label,
        "reliability": meta.reliability,
        "basis": _basis(spec, scaffold, regime),
        "band_cv": meta.band_cv,
        "disclaimer": meta.disclaimer,
    }


def cv_for_spec(spec: Any, *, scaffold: Any | None = None) -> float:
    """Recommended parameter-ensemble CV for ``spec``'s predictive regime."""
    return REGIMES[classify_regime(spec, scaffold=scaffold)].band_cv


__all__ = [
    "APPLIED_FIELDS",
    "EMPIRICAL_ENDPOINT_LABELS",
    "EMPIRICAL_FIELDS",
    "EmpiricalObservations",
    "Measurements",
    "PredictiveRegime",
    "REGIMES",
    "RegimeMeta",
    "agreement_rows",
    "apply_measurements",
    "classify_regime",
    "cv_for_spec",
    "reliability_of",
]
