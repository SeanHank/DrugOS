"""Liver QST (DILI): bile-acid, mitochondrial, redox and hepatocyte-death axes.

DILIsym-style baseline (doc/05 4.2): PBPK liver free exposure drives four
mechanistic axes —

1. **Cholestasis** — canalicular BSEP inhibition raises the intrahepatic
   bile-acid pool and, with it, serum bile-acid / bilirubin load.
2. **Mitochondrial dysfunction** — ETC complex inhibition cuts ATP production
   to a floor, with an adaptive (mitogenic) ATP recovery term.
3. **Oxidative stress** — drug-driven ROS generation vs the glutathione
   buffer depletes GSH.
4. **Hepatocyte death** — combined stress kills hepatocytes (sigmoid) with a
   regeneration term; ALT/AST release and bilirubin rise track the dead
   fraction.  DILI grading follows ALT > 3x ULN and Hy's Law criteria.

In-vitro IC50s (BSEP, ETC complexes, redox) are explicit inputs so each
compound's measured ChEMBL/consortium values wire in directly; the defaults
are class-typical calibration constants (low confidence, doc/06).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from drugos.organ.base import NDArray, free_mg_l_to_nm
from drugos.target.targets import Target


@dataclass(frozen=True, slots=True)
class LiverParams:
    """DILI sub-model constants (IC50s against free liver nM).

    Defaults are class-typical calibration constants (low confidence, doc/06);
    supply compound-specific measured/dose-response values where available.
    """

    bsep_ic50_nm: float = 9.0e4
    mito_ic50_nm: tuple[float, ...] = (3.0e5,)
    redox_ic50_nm: float = 6.0e5
    atp_floor: float = 0.20
    gsh_floor: float = 0.15
    kill_max_1h: float = 0.03
    kill_ec50: float = 0.25
    kill_hill: float = 3.0
    regeneration_1h: float = 0.008
    alt_uln_u_l: float = 40.0
    ast_uln_u_l: float = 40.0
    bilirubin_uln_mg_dl: float = 1.0
    alt_release_per_dead: float = 12.0
    ast_release_per_dead: float = 9.0
    bile_rise_max_fold: float = 3.0
    bili_rise_per_bsep: float = 1.5
    bili_rise_per_dead: float = 1.5
    kexpz_atp: float = 0.05  # adaptive mitogenesis recovery per unit ATP shortfall


@dataclass(frozen=True, slots=True)
class LiverStress:
    """Instantaneous mechanistic state at one time point."""

    free_nm: float
    cholestasis: float
    mito_block: float
    atp_frac: float
    ros: float
    gsh_frac: float
    stress: float
    kill_rate_1h: float


@dataclass(slots=True)
class LiverTrajectory:
    """Liver QST result over the exposure horizon."""

    t_h: NDArray
    free_nm: NDArray
    cholestasis: NDArray
    atp_frac: NDArray
    gsh_frac: NDArray
    stress: NDArray
    dead_frac: NDArray
    alt_u_l: NDArray
    ast_u_l: NDArray
    bilirubin_mg_dl: NDArray
    peak_alt_uln: float
    peak_alt_u_l: float
    peak_ast_u_l: float
    peak_bilirubin_mg_dl: float
    peak_bilirubin_uln: float
    max_dead_frac: float
    dili_grade: int
    hy_law: bool

    def to_series(self) -> dict[str, NDArray]:
        return {
            "t_h": self.t_h,
            "free_nM": self.free_nm,
            "cholestasis": self.cholestasis,
            "atp_frac": self.atp_frac,
            "gsh_frac": self.gsh_frac,
            "stress": self.stress,
            "dead_frac": self.dead_frac,
            "alt_U_L": self.alt_u_l,
            "ast_U_L": self.ast_u_l,
            "bilirubin_mg_dL": self.bilirubin_mg_dl,
        }


def inhibition(c: float, ic50_nm: float) -> float:
    """Fractional inhibition at concentration ``c`` (Michaelis saturating)."""
    if ic50_nm <= 0:
        raise ValueError("ic50_nm must be positive")
    if c <= 0:
        return 0.0
    return c / (c + ic50_nm)


def mitochondrial_block(c: float, ic50s: tuple[float, ...]) -> float:
    """Worst-case ETC inhibition across the supplied complexes."""
    if not ic50s:
        return 0.0
    return max(inhibition(c, ic50) for ic50 in ic50s)


def redox_state(c: float, ic50_nm: float) -> tuple[float, float, float]:
    """ROS generation and resulting GSH fraction at free nM ``c``."""
    if ic50_nm <= 0:
        raise ValueError("ic50_nm must be positive")
    ros = 0.0 if c <= 0 else c / (c + ic50_nm)
    return ros, max(0.15, 1.0 - ros), ros * max(0.0, 1.0 - ros)


def combined_stress(
    cholestasis: float,
    atp_frac: float,
    gsh_frac: float,
    weights: tuple[float, float, float] = (0.5, 0.3, 0.2),
) -> float:
    """Normalized 0..1 hepatocyte stress from the three toxicity axes."""
    return max(
        0.0,
        weights[0] * cholestasis
        + weights[1] * max(0.0, 1.0 - atp_frac)
        + weights[2] * max(0.0, 1.0 - gsh_frac),
    )


def _kill_rate(stress: float, kill_max_1h: float, ec50: float, hill: float) -> float:
    if ec50 <= 0 or hill <= 0 or kill_max_1h <= 0:
        raise ValueError("kill constants must be positive")
    if stress <= 0:
        return 0.0
    s_n = stress**hill
    life_n = ec50**hill
    return float(kill_max_1h * s_n / (s_n + life_n))


def _death_rhs(dead: float, stress: float, p: LiverParams) -> float:
    kill = _kill_rate(stress, p.kill_max_1h, p.kill_ec50, p.kill_hill)
    return kill * (1.0 - dead) - p.regeneration_1h * dead


def aten_floor_factor(mito_block: float, atp_floor: float, adaptive: float) -> float:
    """ATP supply fraction with an adaptive (mitogenic) recovery term."""
    if atp_floor < 0 or atp_floor > 1:
        raise ValueError("atp_floor must be in [0, 1]")
    shortfall = 1.0 - mito_block
    return min(1.0, atp_floor + (1.0 - atp_floor) * shortfall * (1.0 + adaptive))


def dili_grade(alt_uln: float, tbili_uln: float) -> tuple[int, bool]:
    """DILI severity per the doc/05 cut-offs (ALT>3x ULN; Hy's Law)."""
    hy = alt_uln >= 3.0 and tbili_uln >= 2.0
    if hy:
        return 3, True
    if alt_uln >= 3.0:
        return 2, False
    if alt_uln >= 2.0:
        return 1, False
    return 0, False


def liver_params_from_panel(panel: tuple[Target, ...]) -> LiverParams:
    """Wire Stage-2 safety-panel binding data into the liver QST constants.

    Reads the canalicular efflux (BSEP) and mitochondrial ETC complex
    affinities delivered by :func:`drugos.target.targets.safety_panel`
    (kd_nm, already converted from assay IC50).
    """

    def kd_of(name_prefix: str) -> float | None:
        for target in panel:
            if target.name.startswith(name_prefix):
                return target.kd_nm
        return None

    bsep = kd_of("BSEP")
    mito = []
    for prefix in (
        "Mitochondrial complex I",
        "Mitochondrial complex II",
        "Mitochondrial complex III",
        "Mitochondrial complex IV",
    ):
        kd = kd_of(prefix)
        if kd is not None:
            mito.append(kd)
    return LiverParams(
        bsep_ic50_nm=bsep if bsep is not None else 1.0e6,
        mito_ic50_nm=tuple(mito) if mito else (1.0e6,),
    )


def simulate_liver(
    t_h: NDArray,
    c_free_mg_l: NDArray,
    mw: float,
    params: LiverParams | None = None,
    n_eval: int = 601,
    rtol: float = 1e-8,
    atol: float = 1e-9,
) -> LiverTrajectory:
    """Integrate the four-axis liver QST model on the PBPK exposure grid."""
    p = params or LiverParams()
    t = np.asarray(t_h, dtype=float)
    c = free_mg_l_to_nm(np.asarray(c_free_mg_l, dtype=float), mw)
    if t.ndim != 1 or c.ndim != 1 or t.shape[0] != c.shape[0]:
        raise ValueError("t_h and c_free_mg_l must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")

    def stress_at(sol_t: float) -> float:
        cc = float(np.interp(sol_t, t, c))
        chol = inhibition(cc, p.bsep_ic50_nm)
        mito = mitochondrial_block(cc, p.mito_ic50_nm)
        atp_frac = aten_floor_factor(mito, p.atp_floor, p.kexpz_atp * mito)
        gsh = redox_state(cc, p.redox_ic50_nm)[1]
        return combined_stress(chol, atp_frac, gsh)

    def rhs(sol_t: float, y: NDArray) -> NDArray:
        dead = float(y[0])
        return np.array([_death_rhs(dead, stress_at(sol_t), p)], dtype=float)

    y0 = np.array([0.0], dtype=float)
    t_eval = np.linspace(t[0], t[-1], n_eval)
    sol = solve_ivp(rhs, (t[0], t[-1]), y0, t_eval=t_eval, method="LSODA", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"liver solve failed: {sol.message}")
    dead = sol.y[0]

    chol_g = np.array([inhibition(float(cc), p.bsep_ic50_nm) for cc in c], dtype=float)
    mito = np.array([mitochondrial_block(float(cc), p.mito_ic50_nm) for cc in c], dtype=float)
    atp_g = np.array(
        [aten_floor_factor(float(m), p.atp_floor, p.kexpz_atp * float(m)) for m in mito],
        dtype=float,
    )
    gsh_g = np.array([redox_state(float(cc), p.redox_ic50_nm)[1] for cc in c], dtype=float)
    stress_g = np.array(
        [
            combined_stress(float(ch), float(a), float(g))
            for ch, a, g in zip(chol_g, atp_g, gsh_g, strict=True)
        ],
        dtype=float,
    )
    chol: NDArray = np.interp(t_eval, t, chol_g)
    atp: NDArray = np.interp(t_eval, t, atp_g)
    gsh: NDArray = np.interp(t_eval, t, gsh_g)
    stress: NDArray = np.interp(t_eval, t, stress_g)
    free: NDArray = np.interp(t_eval, t, c)

    alt: NDArray = p.alt_uln_u_l * (1.0 + p.alt_release_per_dead * dead)
    ast: NDArray = p.ast_uln_u_l * (1.0 + p.ast_release_per_dead * dead)
    bili: NDArray = p.bilirubin_uln_mg_dl * (
        1.0 + p.bili_rise_per_bsep * chol + p.bili_rise_per_dead * dead
    )
    peak_alt_u_l = float(np.max(alt))
    peak_ast_u_l = float(np.max(ast))
    peak_bili = float(np.max(bili))
    alt_uln = peak_alt_u_l / p.alt_uln_u_l
    bili_uln = peak_bili / p.bilirubin_uln_mg_dl
    grade, hy = dili_grade(alt_uln, bili_uln)
    return LiverTrajectory(
        t_h=t_eval,
        free_nm=free,
        cholestasis=chol,
        atp_frac=atp,
        gsh_frac=gsh,
        stress=stress,
        dead_frac=dead,
        alt_u_l=alt,
        ast_u_l=ast,
        bilirubin_mg_dl=bili,
        peak_alt_uln=alt_uln,
        peak_alt_u_l=peak_alt_u_l,
        peak_ast_u_l=peak_ast_u_l,
        peak_bilirubin_mg_dl=peak_bili,
        peak_bilirubin_uln=bili_uln,
        max_dead_frac=float(np.max(dead)),
        dili_grade=grade,
        hy_law=hy,
    )


__all__ = [
    "LiverParams",
    "LiverStress",
    "LiverTrajectory",
    "aten_floor_factor",
    "combined_stress",
    "dili_grade",
    "inhibition",
    "liver_params_from_panel",
    "mitochondrial_block",
    "redox_state",
    "simulate_liver",
]
