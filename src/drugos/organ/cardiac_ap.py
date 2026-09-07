"""Production-validated cardiac APD90 cross-check via O'Hara-Rudy 2011 (ORd).

The Stage-4 QTc encoder (:mod:`drugos.organ.cardiac`) is a calibrated,
algebraic Emax relation from hERG blockade to QTc prolongation.  This module
additionally runs the *production-validated* O'Hara-Rudy 2011 human
ventricular action-potential model (ORd, PLoS Comput Biol e1002061; BSD
Myokit encoding vendored at ``data/models/ohara-2011.mmt``) under fractional
IKr block and returns the resulting APD90 prolongation.

It backs validation case R-3 (doc/08): the algebraic encoder and the
mechanistic ionic model must agree in *direction* and in *ordering* for the
benchmark hERG binders, giving an independent, open-source, regression-tested
anchor for the cardiac axis without replacing the fast encoder in the
web/robustness compute path (see doc/12 for the integration decision record).

Runtime requirements: ``myokit>=1.39`` (pip) plus the SUNDIALS headers
(``conda install -n <env> -c conda-forge sundials`` or ``apt/brew``) that its
C code-generator needs.  No silent fallback exists: missing either surfaces as
a hard error (doc/09-quality-gate.md, G5).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from drugos.organ.base import NDArray

_MODEL_FILE = Path(__file__).resolve().parents[3] / "data" / "models" / "ohara-2011.mmt"

Runner = Callable[[float], tuple[NDArray, NDArray]]


def _import_myokit() -> Any:
    """Lazy-import myokit so the plain import chain (scipy.special --> numpy.fft)

    is not pulled at package-import time (it clashes with the coverage tracer on
    numpy 2.5.x).  The only entry point is the runner below; a missing myokit
    surfaces as an ImportError just like any other missing dependency."""
    import myokit

    return myokit


@dataclass(frozen=True, slots=True)
class Apd90Result:
    """APD90 outcome of a fractional IKr block in the ORd ventricular model."""

    block_frac: float
    apd90_base_ms: float
    apd90_block_ms: float
    delta_apd90_ms: float

    def to_dict(self) -> dict[str, float]:
        return {
            "block_frac": self.block_frac,
            "apd90_base_ms": self.apd90_base_ms,
            "apd90_block_ms": self.apd90_block_ms,
            "delta_apd90_ms": self.delta_apd90_ms,
        }


def apd90_from_trace(time_ms: NDArray, v_mv: NDArray) -> float:
    """APD90 (ms) of one paced beat, from a voltage trace and its time base.

    Definition (standard ORd post-processing): the interval between the
    upstroke (maximum dV/dt) and the first later crossing of the 90 %
    repolarisation level, computed as 10 % of the plateau-to-rest swing from
    the resting level.  Returns ``nan`` when no repolarisation crossing is
    observed in the window.
    """
    t = np.asarray(time_ms, dtype=float)
    v = np.asarray(v_mv, dtype=float)
    if t.ndim != 1 or v.ndim != 1 or t.shape[0] != v.shape[0]:
        raise ValueError("time_ms and v_mv must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("time_ms needs at least two samples")
    v_rest = float(np.percentile(v, 5.0))
    v_peak = float(np.max(v))
    level = v_rest + 0.10 * (v_peak - v_rest)
    upstroke = int(np.argmax(np.diff(v)))
    crossings = np.where(v < level)[0]
    crossings = crossings[crossings > upstroke]
    if crossings.shape[0] == 0:
        return float("nan")
    return float(t[crossings[0]] - t[upstroke])


def _load_and_run(
    model_path: Path,
    block_frac: float,
    cell_mode: int,
    pre_paces: int,
    run_ms: float,
    log_interval_ms: float,
) -> tuple[NDArray, NDArray]:
    """Run one logged, paced beat of ORd at a fractional IKr block."""
    mkm = _import_myokit()
    model, protocol, _ = mkm.load(str(model_path))
    sim = mkm.Simulation(model, protocol)
    sim.set_constant("cell.mode", cell_mode)
    g_base = float(model.get("ikr.gKr").value())
    if block_frac > 0:
        sim.set_constant("ikr.gKr", g_base * (1.0 - block_frac))
    sim.pre(float(pre_paces) * 1000.0)
    data = sim.run(run_ms, log=["engine.time", "membrane.V"], log_interval=log_interval_ms)
    t = np.asarray(data["engine.time"], dtype=float)
    v = np.asarray(data["membrane.V"], dtype=float)
    return t, v


def ord_apd90(
    block_frac: float,
    *,
    runner: Runner | None = None,
    cell_mode: int = 0,
    pre_paces: int = 50,
    run_ms: float = 1000.0,
    log_interval_ms: float = 0.1,
) -> Apd90Result:
    """APD90 at ``block_frac`` IKr block vs the un-blocked ORd baseline."""
    if not 0.0 <= block_frac <= 1.0:
        raise ValueError("block_frac must be within [0, 1]")
    run = (
        runner
        if runner is not None
        else (
            lambda bf: _load_and_run(_MODEL_FILE, bf, cell_mode, pre_paces, run_ms, log_interval_ms)
        )
    )
    t0, v0 = run(0.0)
    tb, vb = run(block_frac)
    base = apd90_from_trace(t0, v0)
    blocked = apd90_from_trace(tb, vb)
    return Apd90Result(
        block_frac=block_frac,
        apd90_base_ms=base,
        apd90_block_ms=blocked,
        delta_apd90_ms=blocked - base,
    )


__all__ = [
    "Apd90Result",
    "Runner",
    "apd90_from_trace",
    "ord_apd90",
]
