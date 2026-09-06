"""Pathway ODE engine and perturbation simulation (doc/05 section 3.3-3.4).

Compiles a :class:`PathwayModel` into a sparse ODE system, reaches the
drug-free steady state as baseline, then drives the ``input_drive`` reactions
with the Stage-2 occupancy signal s(t) (0..1).  Dose-response utility fits an
Emax/Hill model to the steady-state readout and reports EC50 — used to verify
that signal amplification pushes EC50 below the receptor Kd (doc/05 3.5).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

from drugos.pathway.graph import (
    Activation,
    FirstOrderDegradation,
    PathwayModel,
    ReversibleBinding,
    SourceProduction,
)

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]
_Y2D = np.ndarray[tuple[int, int], np.dtype[np.float64]]

_SIGNAL_SUPPORT = (0.0, 1.0)


@dataclass(slots=True)
class PathwayResult:
    """Time-resolved pathway response to an occupancy/drug signal."""

    model: PathwayModel
    t_h: NDArray
    concentrations: dict[str, NDArray]
    signal: NDArray
    baseline: dict[str, float] = field(default_factory=dict)

    def readout_fold_change(self, readout: str | None = None) -> NDArray:
        node = readout or self.model.readout
        if node is None:
            raise ValueError("model has no readout node")
        base = self.baseline.get(node, np.nan)
        if not base > 0:
            raise ValueError(f"readout {node}: invalid baseline {base}")
        scaled: NDArray = self.concentrations[node] * np.float64(1.0 / base)
        return scaled


@dataclass(slots=True)
class DoseResponseFit:
    """Emax/Hill fit parameters of the steady-state readout vs drive signal."""

    r0: float
    emax: float
    ec50: float
    hill: float

    def predict(self, s: NDArray) -> NDArray:
        num: NDArray = (self.emax - self.r0) * (s**self.hill)
        den: NDArray = self.ec50**self.hill + s**self.hill
        out: NDArray = self.r0 + num / den
        return out

    def drug_ec50_nm(self, kd_nm: float) -> float:
        """Drug concentration (nM) at half-maximum steady-state readout."""
        if not 0.0 < self.ec50 < 1.0:
            raise ValueError("ec50 (fractional occupancy) must be in (0, 1)")
        return kd_nm * self.ec50 / (1.0 - self.ec50)


class _SignalContext:
    """Mutable signal time course consumed by the compiled RHS at solve time."""

    def __init__(self) -> None:
        self.t_h: NDArray = np.array([0.0, 1.0])
        self.signal: NDArray = np.array([0.0, 0.0])

    def at(self, t: float) -> float:
        value = float(np.interp(t, self.t_h, self.signal))
        return min(max(value, _SIGNAL_SUPPORT[0]), _SIGNAL_SUPPORT[1])


@dataclass(slots=True)
class CompiledPathway:
    """Species ordering and a compiled RHS bound to a signal context."""

    names: list[str]
    idx: dict[str, int]
    context: _SignalContext
    n_species: int
    rhs: Callable[[float, NDArray], NDArray]

    def state_names(self) -> list[str]:
        return self.names

    def initial_condition(self, model: PathwayModel) -> NDArray:
        return np.array([model.species[n] for n in self.names], dtype=float)


def compile_model(model: PathwayModel) -> CompiledPathway:
    names = model.node_names()
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    context = _SignalContext()

    def rhs(t: float, y: np.ndarray[tuple[int], np.dtype[np.float64]]) -> NDArray:
        dydt = np.zeros_like(y)
        dydt_t: NDArray = dydt
        mods = {n: float(y[i]) for n, i in idx.items()}
        activations = [r for r in model.reactions if isinstance(r, Activation)]
        bindings = [r for r in model.reactions if isinstance(r, ReversibleBinding)]
        degradations = [r for r in model.reactions if isinstance(r, FirstOrderDegradation)]
        productions = [r for r in model.reactions if isinstance(r, SourceProduction)]
        for reaction in activations:
            if reaction.input_drive:
                drive = context.at(t)
            else:
                driver = reaction.driver or ""
                drive = mods[driver] ** reaction.driver_hill / (
                    reaction.driver_half**reaction.driver_hill
                    + mods[driver] ** reaction.driver_hill
                )
            sub_hill = reaction.substrate_hill
            sub = mods[reaction.substrate]
            michaelis = sub**sub_hill / (reaction.km**sub_hill + sub**sub_hill)
            rate = drive * reaction.vmax * michaelis
            if reaction.inhibitor is not None:
                rate *= 1.0 / (1.0 + mods[reaction.inhibitor] / reaction.ki)
            dydt[idx[reaction.substrate]] -= rate
            dydt[idx[reaction.product]] += rate
        for binding in bindings:
            on = binding.kon * mods[binding.a] * mods[binding.b]
            off = binding.koff * mods[binding.complex]
            dydt[idx[binding.a]] += off - on
            dydt[idx[binding.b]] += off - on
            dydt[idx[binding.complex]] += on - off
        for degradation in degradations:
            dydt[idx[degradation.species]] -= degradation.rate * mods[degradation.species]
        for production in productions:
            dydt[idx[production.species]] += production.rate
        return dydt_t

    return CompiledPathway(names=names, idx=idx, context=context, n_species=n, rhs=rhs)


def _integrate(
    compiled: CompiledPathway,
    model: PathwayModel,
    signal: NDArray,
    t_h: NDArray,
    n_eval: int,
    rtol: float,
    atol: float,
) -> _Y2D:
    signal_arr = np.asarray(signal, dtype=float)
    t_arr = np.asarray(t_h, dtype=float)
    if signal_arr.ndim != 1 or t_arr.shape[0] != signal_arr.shape[0]:
        raise ValueError("signal and t_h must be equal-length 1-D arrays")
    compiled.context.t_h = t_arr
    compiled.context.signal = signal_arr
    y0 = compiled.initial_condition(model)
    t_eval = np.linspace(t_arr[0], t_arr[-1], n_eval)
    sol = solve_ivp(
        compiled.rhs,
        (t_arr[0], t_arr[-1]),
        y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"pathway solve failed: {sol.message}")
    ymat: _Y2D = sol.y  # species x time
    return ymat


def pathway_steady_state(
    model: PathwayModel,
    horizon_h: float = 400.0,
    n_eval: int = 200,
    rtol: float = 1e-8,
    atol: float = 1e-9,
) -> dict[str, float]:
    """Drug-free steady-state concentrations of every species."""
    t = np.linspace(0.0, horizon_h, n_eval)
    signal = np.zeros_like(t)
    compiled = compile_model(model)
    ymat = _integrate(compiled, model, signal, t, n_eval, rtol, atol)
    return {name: float(ymat[i, -1]) for i, name in enumerate(compiled.names)}


def simulate_pathway(
    model: PathwayModel,
    t_h: NDArray,
    signal: NDArray,
    n_eval: int = 601,
    rtol: float = 1e-8,
    atol: float = 1e-9,
) -> PathwayResult:
    """Drive the pathway with a Stage-2 occupancy signal and sample readouts."""
    signal_arr = np.asarray(signal, dtype=float)
    t_arr = np.asarray(t_h, dtype=float)
    compiled = compile_model(model)
    baseline = pathway_steady_state(model, horizon_h=float(t_arr[-1]), n_eval=100)
    ymat = _integrate(compiled, model, signal_arr, t_arr, n_eval, rtol, atol)
    t_eval = np.linspace(t_arr[0], t_arr[-1], n_eval)
    concentrations = {name: ymat[i] for i, name in enumerate(compiled.names)}
    return PathwayResult(
        model=model,
        t_h=t_eval,
        concentrations=concentrations,
        signal=np.clip(np.interp(t_eval, t_arr, signal_arr), *_SIGNAL_SUPPORT),
        baseline=baseline,
    )


def _emax_model(s: NDArray, r0: float, emax: float, ec50: float, n: float) -> NDArray:
    num: NDArray = (emax - r0) * (s**n)
    den: NDArray = ec50**n + s**n
    out: NDArray = r0 + num / den
    return out


def dose_response(
    model: PathwayModel,
    signals: NDArray,
    readout: str | None = None,
    horizon_h: float = 400.0,
    n_eval: int = 200,
) -> DoseResponseFit:
    """Fit an Emax/Hill curve to steady-state readout across signal levels.

    ``signals`` is a grid of constant occupancy-drive values in [0, 1]; the
    steady-state readout concentration at each level is fit to
    R(s) = r0 + (emax - r0) * s^n/(ec50^n + s^n).  By construction
    ``drug_ec50_nm(kd)`` reports the equivalent free-drug EC50.
    """
    node = readout or model.readout
    if node is None:
        raise ValueError("model has no readout node")
    sig = np.asarray(signals, dtype=float)
    if sig.ndim != 1:
        raise ValueError("signals must be a 1-D array")
    levels: list[float] = []
    responses: list[float] = []
    for s in sig:
        t = np.linspace(0.0, horizon_h, n_eval)
        signal = np.full_like(t, s)
        res = simulate_pathway(model, t, signal, n_eval=n_eval)
        levels.append(float(s))
        responses.append(float(res.concentrations[node][-1]))
    x = np.array(levels, dtype=float)
    y = np.array(responses, dtype=float)
    p0 = [float(y[0]), float(y[-1]), 0.5, 1.0]
    bounds = ([0.0, 0.0, 1e-6, 0.2], [np.inf, np.inf, 1.0, 8.0])
    popt, _ = curve_fit(_emax_model, x, y, p0=p0, bounds=bounds, maxfev=5000)
    r0, emax, ec50, n = popt
    return DoseResponseFit(r0=float(r0), emax=float(emax), ec50=float(ec50), hill=float(n))


def mapk_cascade() -> PathwayModel:
    """Canonical 3-tier MAPK cascade with amplification (doc/05 3.2 example).

    Each tier has a weak *basal* arm (keeps a small drug-free readout) and a
    signal-driven arm fed by the Stage-2 occupancy signal; downstream tiers are
    driven with Hill ultrasensitivity, so the steady-state readout EC50 falls
    below the signal that half-saturates the first step (EC50 < receptor Kd
    equivalence; doc/05 3.5).  Sources replenish inactive pools and active
    species degrade, keeping a bounded steady state.
    """

    def arm(substrate: str, product: str, vmax: float, driver_half: float) -> Activation:
        return Activation(
            substrate=substrate,
            product=product,
            driver=substrate,
            vmax=vmax,
            km=1.0,
            driver_half=driver_half,
        )

    model = PathwayModel(
        name="MAPK",
        species={
            "raf": 10.0,
            "raf_active": 0.0,
            "mek": 100.0,
            "mek_active": 0.0,
            "erk": 1000.0,
            "erk_active": 0.0,
        },
        input_node="signal",
        readout="erk_active",
    )
    model.add_activation(arm("raf", "raf_active", vmax=0.005, driver_half=20.0))
    model.add_activation(
        Activation(
            substrate="raf",
            product="raf_active",
            vmax=2.0,
            km=1.0,
            input_drive=True,
        )
    )
    model.add_activation(
        Activation(
            substrate="mek",
            product="mek_active",
            driver="raf_active",
            vmax=2.0,
            km=2.0,
            driver_half=1.0,
            driver_hill=2.0,
        )
    )
    model.add_activation(
        Activation(
            substrate="erk",
            product="erk_active",
            driver="mek_active",
            vmax=4.0,
            km=4.0,
            driver_half=2.0,
            driver_hill=2.0,
        )
    )
    for species in ("raf", "mek", "erk"):
        model.add_reaction(SourceProduction(species=species, rate=0.5))
        model.add_reaction(FirstOrderDegradation(species=species, rate=0.05))
    for species in ("raf_active", "mek_active", "erk_active"):
        model.add_reaction(FirstOrderDegradation(species=species, rate=0.05))
    return model


__all__ = [
    "CompiledPathway",
    "DoseResponseFit",
    "PathwayResult",
    "compile_model",
    "dose_response",
    "mapk_cascade",
    "pathway_steady_state",
    "simulate_pathway",
]
