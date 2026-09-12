"""Production-validated SBML pathway lane (doc/12 row Stage-3).

Loads a vendored, openly-licensed SBML model with `libsbml`, turns every
reaction into its kinetic-law ODE contribution (reactant/product
stoichiometry + kinetic-law formula over species and local parameters), and
integrates with scipy.  The canonical wired model is the ultrasensitive
**three-tier MAPK cascade** of Huang & Ferrell (*PNAS* 1996; Kholodenko &
wider RAS/MAPK lineage), from BioModels `BIOMD0000000009` (CC0), vendored at
``data/models/huang1996-mapk-cascade.xml``.

This gives the Stage-3 signaling pathway a **parsed-from-artifact production
model** instead of the hand-built toy cascade (`mapk_cascade`), so the pathway
constants are the published model's.  No silent fallback (G5): a missing
vendored file or libsbml raises an explicit error.

Signal coupling (biologically grounded): the drug-occupancy signal s(t) (0..1)
is interpreted as inhibition of the pathway stimulus.  In the Huang model the
pathway's input is the fixed MAPKKK activator **E1** (Ras/activator density).
A target-bound drug that inhibits the upstream signal lowers the effective
stimulus: ``E1_eff = E1 * (1 - s)``.  Because the cascade has intrinsic
ultrasensitivity, even a modest occupancy produces a sharp, monotone drop in
the pathway readout (doubly-phosphorylated ERK, ``PP_K``) — the
amplification-to-inhibition behaviour the in-house design intended.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.integrate import solve_ivp

from drugos.pathway.simulator import PathwayResult

if TYPE_CHECKING:
    from libsbml import ListOfSpeciesReferences

NDArray = np.ndarray[tuple[int], np.dtype[np.float64]]


def _import_libsbml() -> Any:
    """Lazy-import python-libsbml so the native binding is loaded on demand.

    The distribution ships no type annotations; the subset DrugOS uses is
    typed by ``stubs/libsbml/__init__.pyi`` (doc/06-technology-stack.md, G2).
    """
    import libsbml

    return libsbml


def model_file() -> Path:
    """Path to the vendored, sha256-pinned Huang/Levchenko MAPK SBML."""
    return Path(__file__).resolve().parents[3] / "data" / "models" / "huang1996-mapk-cascade.xml"


_SAFE_MATH = frozenset({"pow", "exp", "sqrt", "log", "sin", "cos", "tan", "log10"})
_INPUT_SPECIES = "E1"
_READOUT_SPECIES = "PP_K"


def _rate(expr: str, spec: SbmlSpec, mods: dict[str, float]) -> float:
    """Evaluate a sanitized per-reaction term against the current state."""
    return float(
        eval(expr, {"__builtins__": {}}, {"np": np, "compartment": spec.compartment, **mods})
    )


def _safe_formula(formula: str, allowed: set[str], local_params: set[str]) -> str:
    """Sanitize an SBML math string to a restricted, attribute-safe expression.

    Only known species/parameter identifiers and a closed math-function
    allowlist are permitted; ``compartment`` is kept as a dimensionless factor
    (the model uses it as a 1.0 volume scale).  Anything unexpected raises.
    """
    tokens = set(re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", formula))
    bad = tokens - allowed - local_params - {"compartment"} - _SAFE_MATH
    if bad:
        raise ValueError(f"SBML formula references unknown identifiers: {sorted(bad)}")
    body = formula.replace("^", "**")
    body = re.sub(
        r"\b(" + "|".join(re.escape(w) for w in sorted(_SAFE_MATH)) + r")\b",
        r"np.\1",
        body,
    )
    return body


@dataclass(frozen=True, slots=True)
class SbmlSpec:
    """A parsed SBML model: species initials + compiled reactions."""

    name: str
    species: dict[str, float]
    reactions: tuple[_Reaction, ...]
    compartment: float = 1.0


@dataclass(frozen=True, slots=True)
class _Reaction:
    """One SBML reaction reduced to the numbers an integrator needs."""

    id: str
    formula: str
    reactants: dict[str, float]
    products: dict[str, float]
    local_params: dict[str, float]


def parse_sbml(path: Path) -> SbmlSpec:
    """Parse an SBML file into a :class:`SbmlSpec` (species, reactions, params)."""
    lib = _import_libsbml()
    doc = lib.readSBML(str(path))
    model = doc.getModel() if doc is not None else None
    if model is None:
        raise ValueError(f"failed to parse SBML: {path}")

    compartment = 1.0
    if model.getNumCompartments() > 0 and model.getCompartment(0).isSetSize():
        compartment = float(model.getCompartment(0).getSize())

    species: dict[str, float] = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        if sp.isSetInitialConcentration():
            species[sp.getId()] = float(sp.getInitialConcentration())

    def _stoich(term_list: ListOfSpeciesReferences) -> dict[str, float]:
        out: dict[str, float] = {}
        for t in term_list:
            out[t.getSpecies()] = out.get(t.getSpecies(), 0.0) + float(t.getStoichiometry())
        return out

    parsed: list[_Reaction] = []
    allowed = set(species) | {"compartment"}
    for i in range(model.getNumReactions()):
        rx = model.getReaction(i)
        kl = rx.getKineticLaw()
        if kl is None or not kl.isSetFormula():
            continue
        params = {
            p.getId(): float(p.getValue()) for p in kl.getListOfParameters() if p.isSetValue()
        }
        formula = _safe_formula(kl.getFormula(), allowed, set(params))
        parsed.append(
            _Reaction(
                id=rx.getId(),
                formula=formula,
                reactants=_stoich(rx.getListOfReactants()),
                products=_stoich(rx.getListOfProducts()),
                local_params=params,
            )
        )
    return SbmlSpec(
        name=model.getName() or model.getId(),
        species=species,
        reactions=tuple(parsed),
        compartment=compartment,
    )


def _build_term_lists(
    spec: SbmlSpec, names: list[str], input_species: str
) -> tuple[dict[str, list[str]], dict[str, set[str]]]:
    """Per-species algebraic terms and which terms reference the input species.

    Each reaction's kinetic law (constant-folded local params) contributes
    ``+stoich*rate`` to each product and ``-stoich*rate`` to each reactant.
    Terms that mention the input species are separated so the signal can scale
    them independently of the (constant) species values in the term string.
    """
    terms: dict[str, list[str]] = {n: [] for n in names}
    input_refs: dict[str, set[str]] = {n: set() for n in names}
    for r in spec.reactions:
        params = {k: f"{v:.12g}" for k, v in r.local_params.items()}
        expr = r.formula
        for k, v in params.items():
            expr = re.sub(rf"\b{re.escape(k)}\b", v, expr)
        # SBML semantics: each kinetic law carries an explicit `compartment`
        # volume factor, and the BioModels concentration ODE multiplies it back
        # out with 1/compartment (the net ODE is sum over reactions of
        # +/- stoich * rate/volume).  We therefore fold `compartment` into a
        # per-term 1/compartment scale so the rates are in concentration/time.
        scale = f"({1.0 / spec.compartment:.12g})"
        for sp, coef in r.products.items():
            terms[sp].append(f"({scale})*({coef:g})*({expr})")
            if input_species in set(re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", expr)):
                input_refs[sp].add(r.id)
        for sp, coef in r.reactants.items():
            terms[sp].append(f"({scale})*(-{coef:g})*({expr})")
            if input_species in set(re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", expr)):
                input_refs[sp].add(r.id)
    return terms, input_refs


def _steady_state(
    names: list[str], spec: SbmlSpec, terms: dict[str, list[str]]
) -> dict[str, float]:
    """Drug-free steady state: integrate the un-driven system to equilibrium."""
    y0 = np.array([spec.species[n] for n in names], dtype=float)

    def rhs(t: float, y: np.ndarray[tuple[int], np.dtype[np.float64]]) -> NDArray:
        del t
        mods = {n: float(y[i]) for i, n in enumerate(names)}
        dydt = np.zeros_like(y)
        for n, tlist in terms.items():
            total = sum(_rate(e, spec, mods) for e in tlist)
            dydt[names.index(n)] += total
        return dydt

    sol = solve_ivp(rhs, (0.0, 1e6), y0, method="LSODA", rtol=1e-9, atol=1e-14)
    if not sol.success:
        raise RuntimeError(f"SBML steady-state solve failed: {sol.message}")
    return {n: float(sol.y[i, -1]) for i, n in enumerate(names)}


def simulate_sbml_pathway(
    t_h: NDArray,
    signal: NDArray,
    path: Path | None = None,
    n_eval: int = 601,
    input_species: str = _INPUT_SPECIES,
    readout: str = _READOUT_SPECIES,
) -> PathwayResult:
    """Drive a vendored SBML scaffold and return a :class:`PathwayResult`.

    ``signal`` (per-occupancy 0..1 over ``t_h``) scales the pathway stimulus
    down: effective input ``<input_species>`` becomes ``X_0 * (1 - s(t))``,
    an upstream inhibition.  ``readout`` selects which species backs the
    canonical output; the default scaffold's is doubly-phosphorylated ERK
    (``PP_K``).  Baseline is the drug-free steady state, so
    ``PathwayResult.readout_fold_change`` reports signal gain (inhibition shows
    as ``< 1`` for increasing signal).
    """
    p = path or model_file()
    if not p.is_file():
        raise FileNotFoundError(f"vendored SBML model missing: {p}")
    spec = parse_sbml(p)
    if input_species not in spec.species:
        raise KeyError(f"SBML model has no input species {input_species!r}")
    if readout not in spec.species:
        raise KeyError(f"SBML model has no readout species {readout!r}")

    t = np.asarray(t_h, dtype=float)
    s = np.asarray(signal, dtype=float)
    if t.ndim != 1 or s.ndim != 1 or t.shape[0] != s.shape[0]:
        raise ValueError("t_h and signal must be equal-length 1-D arrays")

    names = sorted(spec.species)
    x_0 = spec.species[input_species]
    terms, _input_refs = _build_term_lists(spec, names, input_species)

    baseline = _steady_state(names, spec, terms)

    # Signal-dependent term builder: substitute input -> X_0*(1-s) at run time.
    def rhs_for(s_frac: float) -> Callable[[float, NDArray], NDArray]:
        x_eff = x_0 * (1.0 - s_frac)

        def rhs(t: float, y: np.ndarray[tuple[int], np.dtype[np.float64]]) -> NDArray:
            del t
            mods = {n: float(y[i]) for i, n in enumerate(names)}
            mods[input_species] = x_eff
            dydt = np.zeros_like(y)
            for n, tlist in terms.items():
                total = sum(_rate(e, spec, mods) for e in tlist)
                dydt[names.index(n)] += total
            return dydt

        return rhs

    y0 = np.array([spec.species[n] for n in names], dtype=float)
    t_eval = np.linspace(t[0], t[-1], n_eval)
    sig_arr = np.clip(s, 0.0, 1.0)

    def f(tt: float, yy: np.ndarray[tuple[int], np.dtype[np.float64]]) -> NDArray:
        return rhs_for(float(sig_arr[min(int(np.searchsorted(t, tt)), len(t) - 1)]))(tt, yy)

    sol = solve_ivp(f, (t[0], t[-1]), y0, t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"SBML pathway solve failed: {sol.message}")
    concentrations = {n: sol.y[i] for i, n in enumerate(names)}

    from drugos.pathway.graph import PathwayModel

    model = PathwayModel(
        name=spec.name, species=dict(spec.species), input_node="signal", readout=readout
    )
    return PathwayResult(
        model=model,
        t_h=t_eval,
        concentrations=concentrations,
        signal=np.interp(t_eval, t, sig_arr),
        baseline=baseline,
    )


__all__ = [
    "SbmlSpec",
    "model_file",
    "parse_sbml",
    "simulate_sbml_pathway",
]
