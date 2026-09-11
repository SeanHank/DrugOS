"""PBPK simulation: ODE solve, PK metrics and result container."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

from drugos.inputs.models import Route
from drugos.pk.pbpk_build import PBPKModel, TargetBinding

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]
_METHOD = Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"]
_DEFAULT_N_EVAL = 601


def _site_label(binding: TargetBinding) -> str:
    return f"{binding.tissue}::{binding.target.name}"


_DEPOT_ROUTES = frozenset({Route.SUBCUTANEOUS, Route.INTRAMUSCULAR, Route.TRANSDERMAL})
_IV_ROUTES = frozenset({Route.IV_BOLUS, Route.IV_INFUSION})


def bioavailable_fraction(
    model: PBPKModel,
    route: Route,
    dose_mg: float,
    feces_cum_mg: NDArray | None,
) -> float:
    """Model-predicted systemic bioavailability ``F`` (doc/05 §2.1, 0..1).

    IV doses are 100% bioavailable; depot routes (SC/IM/transdermal) deliver
    the depot availability fraction directly (they bypass the portal vein).
    For oral dosing the absorbed fraction (dose minus the colon-transit feces
    sink) reaches the liver through the portal blood and suffers a
    well-stirred first-pass hepatic extraction ``Eh = CL_h/(Q_h + CL_h)`` and
    (when set) a first-pass intestinal extraction ``Eg`` from the absorbed
    flux, so ``F = Fa * (1 - Eh) * (1 - Eg)``.  The PBPK model already routes
    absorbed oral drug into the liver compartment, so this is the reported
    first-pass-corrected bioavailability.
    """
    if route in _IV_ROUTES:
        return 1.0
    if route in _DEPOT_ROUTES:
        return min(1.0, max(0.0, float(model.absorption.depot_bioavailability)))
    absorbed = 1.0
    if feces_cum_mg is not None and dose_mg > 0:
        absorbed = 1.0 - float(feces_cum_mg[-1]) / dose_mg
    qh = model.physiology.organ_flow["liver"] * 60.0  # L/h
    eh = model.cl_hep_l_h / (qh + model.cl_hep_l_h) if (qh + model.cl_hep_l_h) > 0 else 0.0
    eg = float(model.absorption.gut_extraction_eg)
    value = max(0.0, min(1.0, absorbed * (1.0 - eh) * (1.0 - eg)))
    return float(value)


@dataclass(slots=True)
class PkMetrics:
    """Non-compartmental PK metrics computed from the plasma profile."""

    cmax_mg_l: float
    tmax_h: float
    auc_last_mgh_l: float
    auc_inf_mgh_l: float
    term_half_life_h: float
    lambda_z_1h: float
    cl_l_h: float
    mrt_h: float
    vss_l: float
    c_last_mg_l: float
    dose_mg: float
    route: Route
    f_abs: float = 1.0

    def as_dict(self) -> dict[str, float | bool | int | None]:
        return {
            "cmax_mg_l": self.cmax_mg_l,
            "tmax_h": self.tmax_h,
            "auc_last_mg_h_l": self.auc_last_mgh_l,
            "auc_inf_mg_h_l": self.auc_inf_mgh_l,
            "t_half_h": self.term_half_life_h,
            "lambda_z_1_h": self.lambda_z_1h,
            "cl_l_h": self.cl_l_h,
            "mrt_h": self.mrt_h,
            "vss_l": self.vss_l,
            "c_last_mg_l": self.c_last_mg_l,
            "bioavailability_f": self.f_abs,
        }


@dataclass(slots=True)
class PBPKResult:
    """Concentration-time output of a PBPK simulation.

    Concentrations are in mg/L, times in hours.  Tissue entries are total
    (free + bound) tissue concentrations; ``unbound_tissues`` are the free
    concentrations that drive pharmacological and toxicological responses.
    ``tmdd`` carries the native-TMDD bound and internalized-cleared
    trajectories per binding site when ``model.target_binding`` is used.
    ``skin`` carries the multi-layer transdermal layer amounts plus the
    cumulative absorbed mass when ``model.absorption.skin_layers`` is set.
    """

    t: NDArray
    plasma_total: NDArray
    plasma_free: NDArray
    venous_plasma: NDArray
    tissues: dict[str, NDArray] = field(default_factory=dict)
    unbound_tissues: dict[str, NDArray] = field(default_factory=dict)
    urine_cum_mg: NDArray | None = None
    feces_cum_mg: NDArray | None = None
    dose_mg: float = 0.0
    route: Route = Route.ORAL
    n_eval: int = _DEFAULT_N_EVAL
    solver: str = "LSODA"
    final_state: NDArray | None = None
    bioavailability_f: float | None = None
    tmdd: dict[str, dict[str, NDArray]] | None = None
    skin: dict[str, NDArray] | None = None

    def pk_metrics(self) -> PkMetrics:
        return compute_pk_metrics(
            self.t,
            self.plasma_total,
            self.dose_mg,
            self.route,
            f_abs=self.bioavailability_f or 1.0,
        )

    def to_data_contract(self) -> dict[str, object]:
        """Stage-1 slice of the pipeline data contract (doc/03, section 3)."""
        return {
            "time": self.t.tolist(),
            "plasma_total": self.plasma_total.tolist(),
            "plasma_free": self.plasma_free.tolist(),
            "venous_plasma": self.venous_plasma.tolist(),
            "tissues": {k: v.tolist() for k, v in self.tissues.items()},
            "unbound_tissues": {k: v.tolist() for k, v in self.unbound_tissues.items()},
            "urine_cum_mg": None if self.urine_cum_mg is None else self.urine_cum_mg.tolist(),
            "feces_cum_mg": None if self.feces_cum_mg is None else self.feces_cum_mg.tolist(),
            "pk_metrics": self.pk_metrics().as_dict(),
            "tmdd": None
            if self.tmdd is None
            else {k: {kk: vv.tolist() for kk, vv in v.items()} for k, v in self.tmdd.items()},
            "skin": None if self.skin is None else {k: v.tolist() for k, v in self.skin.items()},
        }


def simulate_pbpk(
    model: PBPKModel,
    tmax_h: float = 48.0,
    n_eval: int = _DEFAULT_N_EVAL,
    rtol: float = 1e-8,
    atol: float = 1e-9,
    method: _METHOD = "LSODA",
) -> PBPKResult:
    """Integrate the PBPK system over ``tmax_h`` hours.

    Events (bolus/oral doses, infusion starts/ends) are handled piecewise:
    the state is integrated segment by segment and bolus/oral doses are applied
    instantaneously at event boundaries.
    """
    if tmax_h <= 0:
        raise ValueError("tmax_h must be positive")

    events = sorted(model.dose_plan.events, key=lambda e: e.time_h)
    event_times = [e.time_h for e in events if e.time_h <= tmax_h]
    # 0.0 is always present below, so ``boundaries`` is non-empty and starts
    # at 0 by construction; no extra guard is needed.
    boundaries = sorted(set([0.0, tmax_h, *event_times]))

    y = model.initial_state()
    t_parts: list[NDArray] = []
    y_tissues: dict[str, list[NDArray]] = {name: [] for name in model.order}
    y_plasma: list[NDArray] = []
    y_plasma_free: list[NDArray] = []
    y_venous: list[NDArray] = []
    y_urine: list[NDArray] = []
    y_feces: list[NDArray] = []
    y_skin: dict[str, list[float]] = {key: [] for key in _SKIN_KEYS}
    y_skin_rate: list[float] = []

    def record(t: NDArray, y: np.ndarray[tuple[int, int], np.dtype[np.float64]]) -> None:
        c_ab = y[model.state_index("arterial")] / model.physiology.arterial_blood_l
        c_vb = y[model.state_index("venous")] / model.physiology.venous_blood_l
        y_plasma.append(c_ab / model.bp)
        y_plasma_free.append(model.fup * c_ab / model.bp)
        y_venous.append(c_vb / model.bp)
        y_urine.append(y[model.state_index("urine")])
        y_feces.append(y[model.state_index("feces")])
        for name in model.order:
            conc = y[model.state_index(name)] / model.physiology.organ_volume[name]
            y_tissues[name].append(conc)
        if model._skin_indices:
            for key in _SKIN_KEYS:
                y_skin[key].extend(y[model._skin_indices[f"skin_{key}"], :].tolist())
            y_skin_rate.extend(
                float(model.skin_absorption_rate_mg_h(y[:, j])) for j in range(y.shape[1])
            )
        if model.target_binding:
            for j in range(y.shape[1]):
                for site, row in zip(
                    model.target_binding,
                    model.tmdd_state_record(y[:, j]),
                    strict=True,
                ):
                    tmdd_sites[_site_label(site)].append(row)

    tmdd_sites: dict[str, list[dict[str, float]]] = {
        _site_label(site): [] for site in model.target_binding
    }

    for i in range(len(boundaries) - 1):
        a, b = boundaries[i], boundaries[i + 1]
        model.apply_event(y, a)
        sub = np.linspace(a, b, max(2, int(n_eval * (b - a) / tmax_h) + 1))
        sol = solve_ivp(
            model.rhs,
            (a, b),
            y,
            method=method,
            t_eval=sub,
            rtol=rtol,
            atol=atol,
        )
        if not sol.success:
            raise RuntimeError(f"solver failed on segment [{a}, {b}]: {sol.message}")
        y = sol.y[:, -1]
        t_parts.append(sol.t)
        record(sol.t, sol.y)

    # ``boundaries`` always contains at least [0, tmax_h] (tmax_h > 0 is
    # enforced above), so at least one segment is integrated and ``t_parts``
    # is never empty here.
    t_all = np.concatenate(t_parts)
    plasma = np.concatenate(y_plasma)
    plasma_free = np.concatenate(y_plasma_free)
    venous = np.concatenate(y_venous)
    urine = np.concatenate(y_urine)
    feces = np.concatenate(y_feces)
    tissues = {name: np.concatenate(ys) for name, ys in y_tissues.items()}
    unbound_tissues = {name: tissues[name] / model.partition.kpu[name] for name in tissues}

    tmdd = None
    if model.target_binding:
        tmdd = {
            site: {
                key: np.asarray([row[key] for row in rows], dtype=float)
                for key in ("receptor_nmol", "complex_nmol", "bound_mg", "cleared_mg")
            }
            for site, rows in tmdd_sites.items()
        }

    return PBPKResult(
        t=t_all,
        plasma_total=plasma,
        plasma_free=plasma_free,
        venous_plasma=venous,
        tissues=tissues,
        unbound_tissues=unbound_tissues,
        urine_cum_mg=urine,
        feces_cum_mg=feces,
        dose_mg=model.dose_plan.total_dose_mg,
        route=_route_of(model),
        n_eval=n_eval,
        solver=method,
        final_state=y,
        bioavailability_f=bioavailable_fraction(
            model, _route_of(model), model.dose_plan.total_dose_mg, feces
        ),
        tmdd=tmdd,
        skin=_skin_result(model, t_all, y_skin, y_skin_rate),
    )


_SKIN_KEYS = ("surface", "sc", "ve", "dermis", "unabsorbed")


def _skin_result(
    model: PBPKModel,
    t_all: NDArray,
    y_skin: dict[str, list[float]],
    y_skin_rate: list[float],
) -> dict[str, NDArray] | None:
    """Assemble the per-layer skin trajectories and cumulative absorption."""
    if not model._skin_indices or not y_skin:
        return None
    rate = np.asarray(y_skin_rate, dtype=float)
    dt = np.diff(np.asarray(t_all, dtype=float))
    cum = np.concatenate(([0.0], np.cumsum(0.5 * dt * (rate[:-1] + rate[1:]))))
    return {
        **{f"{k}_mg": np.asarray(v, dtype=float) for k, v in y_skin.items()},
        "absorbed_mg": cum,
    }


def _route_of(model: PBPKModel) -> Route:
    routes = model.dose_plan.routes
    return routes[0] if routes else Route.ORAL


def compute_pk_metrics(
    t: NDArray,
    concentration: NDArray,
    dose_mg: float,
    route: Route,
    *,
    f_abs: float = 1.0,
) -> PkMetrics:
    """Non-compartmental PK metrics (linear trapezoid + log-linear terminal).

    ``f_abs`` is the model-predicted systemic bioavailability (default 1.0,
    e.g. IV); callers with access to the PBPK model (``PBPKResult``) pass the
    orally-extracted first-pass value from :func:`bioavailable_fraction`.
    """
    if not (0.0 <= f_abs <= 1.0):
        raise ValueError("f_abs must be in [0, 1]")
    if t.size < 3:
        raise ValueError("at least 3 time points are required for PK metrics")
    conc = np.asarray(concentration, dtype=float)
    t = np.asarray(t, dtype=float)

    i_max = int(np.argmax(conc))
    cmax = float(conc[i_max])
    tmax = float(t[i_max])

    auc_last = float(np.trapezoid(conc, t))

    # Terminal half-life: log-linear fit on the tail after Cmax.
    tail = conc[i_max + 1 :]
    t_tail = t[i_max + 1 :]
    if tail.size >= 3 and tail[-1] > 0 and np.all(tail > 0):
        floor = cmax * 0.05
        keep = tail > floor
        if keep.sum() >= 3:
            t_sel, c_sel = t_tail[keep], tail[keep]
        else:
            t_sel, c_sel = t_tail[-3:], tail[-3:]
        (slope, _) = np.polyfit(t_sel, np.log(c_sel), 1)
        lambda_z = -float(slope)
    else:
        lambda_z = 0.0

    if lambda_z > 0:
        c_last = float(conc[-1])
        auc_inf = auc_last + c_last / lambda_z
        term_half_life = np.log(2.0) / lambda_z
    else:
        auc_inf = auc_last
        c_last = float(conc[-1])
        term_half_life = np.inf

    if auc_inf > 0 and dose_mg > 0:
        cl = dose_mg / auc_inf
        aumc = float(np.trapezoid(conc * t, t))
        mrt = aumc / auc_last if auc_last > 0 else 0.0
        vss = cl * mrt
    else:
        cl = 0.0
        mrt = 0.0
        vss = 0.0

    return PkMetrics(
        cmax_mg_l=cmax,
        tmax_h=tmax,
        auc_last_mgh_l=auc_last,
        auc_inf_mgh_l=auc_inf,
        term_half_life_h=term_half_life,
        lambda_z_1h=lambda_z,
        cl_l_h=cl,
        mrt_h=mrt,
        vss_l=vss,
        c_last_mg_l=float(conc[-1]),
        dose_mg=dose_mg,
        route=route,
        f_abs=f_abs,
    )
