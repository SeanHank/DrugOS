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
are class-typical calibration constants (low confidence, doc/06).  The
**cholestasis axis is production-validated**: it is anchored to the open,
clinically-validated GCDCA bile-acid PBK of de Bruijn & Rietjens
(Arch. Toxicol. 2024, CC BY 4.0, R-7), where free-hepatic drug competitively
inhibits BSEP efflux and bile-acid pool accumulation above a 1.5x risk
threshold drives cholestatic stress (see :func:`simulate_gcdca_pbk`).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from drugos.organ.base import NDArray, free_mg_l_to_nm
from drugos.target.targets import Target

_KM_FN = Callable[[float], float]


@dataclass(frozen=True, slots=True)
class LiverParams:
    """DILI sub-model constants (IC50s against free liver nM).

    The BSEP IC50 (driving the validated bile-acid PBK cholestasis axis, R-7)
    defaults to a **benign** low-affinity value (Ki ~ 150 uM, an order of
    magnitude weaker than the cholestatic reference dataset); the measured
    efflux-inhibition affinities from Stage-2 replace it per-compound where
    available.  Mito/redox/kill constants are class-typical calibration
    constants (low confidence, doc/06).

    ``regen_pathway_*`` couple the Stage-3 pathway readout into hepatocyte
    regeneration (the occupancy -> pathway -> organ -> phenotype chain): the
    ERK/MAPK proliferative fold-change scales ``regeneration_1h``
    (``regeneration_scale``), with a damped gain and a floor so an off-target
    pathway perturbation can attenuate but never fully ablate regeneration.
    """

    bsep_ic50_nm: float = 3.0e5
    mito_ic50_nm: tuple[float, ...] = (3.0e5,)
    redox_ic50_nm: float = 6.0e5
    atp_floor: float = 0.20
    gsh_floor: float = 0.15
    kill_max_1h: float = 0.03
    kill_ec50: float = 0.25
    kill_hill: float = 3.0
    regeneration_1h: float = 0.008
    regen_pathway_gain: float = 0.5
    regen_pathway_min: float = 0.5
    regen_pathway_max: float = 1.5
    alt_uln_u_l: float = 40.0
    ast_uln_u_l: float = 40.0
    bilirubin_uln_mg_dl: float = 1.0
    alt_release_per_dead: float = 12.0
    ast_release_per_dead: float = 9.0
    bile_rise_max_fold: float = 3.0
    bili_rise_per_bsep: float = 1.5
    bili_rise_per_dead: float = 1.5
    kexpz_atp: float = 0.05  # adaptive mitogenesis recovery per unit ATP shortfall
    immune_weight: float = 0.0  # immune-mediated DILI axis; inert at 0 (default)
    immune_ic50_nm: float | None = None  # immune-hazard IC50; None disables the axis
    immune_hill: float = 1.0
    immune_recruit_1h: float = 0.02  # adaptive immune-cell recruitment rate
    immune_decay_1h: float = 0.05  # immune-response decay rate


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
    immune: float = 0.0


@dataclass(slots=True)
class LiverTrajectory:
    """Liver QST result over the exposure horizon."""

    t_h: NDArray
    free_nm: NDArray
    cholestasis: NDArray
    atp_frac: NDArray
    gsh_frac: NDArray
    stress: NDArray
    immune: NDArray
    dead_frac: NDArray
    regen_scale: NDArray
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
            "immune": self.immune,
            "dead_frac": self.dead_frac,
            "regen_scale": self.regen_scale,
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


def immune_hazard(c: float, ic50_nm: float, hill: float = 1.0) -> float:
    """Immune-mediated DILI hazard from free hepatic nM ``c`` (Hill sigmoid).

    The hazard represents haptenized reactive-intermediate / danger-signal
    formation driving adaptive immune-cell recruitment (doc/05 4.2 seam, now
    engaged).  It is inert at zero exposure and saturable at high
    exposure; the *dynamic* immune response derives from it via
    ``immune_recruit_1h``/``immune_decay_1h`` in :func:`simulate_liver`.
    """
    if ic50_nm <= 0:
        raise ValueError("ic50_nm must be positive")
    if hill <= 0:
        raise ValueError("hill must be positive")
    if c <= 0:
        return 0.0
    cn = c**hill
    return float(cn / (cn + ic50_nm**hill))


# ---------------------------------------------------------------------------
# Production-validated cholestasis anchor (R-7): GCDCA bile-acid PBK
#
# The cholestasis axis is driven by the open, clinically-validated bile-acid
# kinetic model of de Bruijn & Rietjens (Arch. Toxicol. 98:3077-3095, 2024,
# doi:10.1007/s00204-024-03775-6, CC BY 4.0).  The published equations for
# glycochenodeoxycholic acid (GCDCA) enterohepatic circulation are ported
# here as a PBK: NTCP/ASBT uptake, BSEP-mediated canalicular efflux and de
# novo synthesis == faecal loss.  Drug-induced cholestasis enters as
# **competitive inhibition of BSEP efflux**, Km_BSEP_app = Km_BSEP*(1 + C/Ki),
# where C is the free hepatic drug concentration and Ki the BSEP-efflux
# inhibition constant.  A >1.5-fold increase of the intrahepatic bile-acid
# pool is the authors' validated cholestasis risk threshold.
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BileAcidParams:
    """GCDCA enterohepatic-circulation constants (de Bruijn & Rietjens 2024).

    Reference individual (70 kg); values are the published kinetic constants,
    including the Sf9-vesicle BSEP Vmax/Km and the hepatocyte NTCP/ASBT
    parameters scaled to the whole liver / ileum.
    """

    bw: float = 70.0
    vil: float = 5.9
    qc: float = 389.988
    fq_ad: float = 0.05
    fq_bo: float = 0.05
    fq_br: float = 0.12
    fq_gu: float = 0.146462
    fq_he: float = 0.04
    fq_ki: float = 0.19
    fq_h: float = 0.215385
    fq_mu: float = 0.17
    fq_sk: float = 0.05
    fq_sp: float = 0.017231
    fq_re: float = 0.103855
    fv_ad: float = 0.21
    fv_bo: float = 0.085629
    fv_br: float = 0.02
    fv_gu: float = 0.0171
    fv_li: float = 0.021
    fv_liw: float = 0.017535
    fv_lew: float = 0.003465
    fv_mu: float = 0.4
    fv_sk: float = 0.0371
    fv_sp: float = 0.0026
    fv_te: float = 0.01
    fv_ve: float = 0.0514
    fv_pv: float = 0.0064
    fv_ar: float = 0.0257
    fv_re: float = 0.0997711
    gdose: float = 3020.0
    ka_cba: float = 0.09
    ktj_fed: float = 2.155
    ktj_fasted: float = 2.44
    kti_fed: float = 1.2
    kti_fasted: float = 2.76
    ge: float = 1.2
    fec: float = 0.05
    qgb: float = 0.5
    a_bsep: float = 0.839
    mw_bsep: float = 140000.0
    hep: float = 99.0
    wl: float = 1400.0
    vmax_bsep_c: float = 5.848
    km_bsep: float = 4.3
    vmax_ntcp_c: float = 510.0
    km_ntcp: float = 5.284
    sf_oatp: float = 1.25
    vmax_asbt_c: float = 203.2
    km_asbt: float = 0.662
    sf_asbt: float = 4712.0 * 2.8e-6
    bp_ba: float = 0.55
    kp_ad_cba: float = 0.05
    kp_gu_cba: float = 0.18
    kp_slp_cba: float = 0.19
    kp_ra_cba: float = 0.136
    kp_lew_cba: float = 0.31
    cb_fs_cba: float = 0.45


def bsep_ki_from_ic50_nm(ic50_nm: float) -> float:
    """BSEP-efflux inhibition constant Ki (umol/L) from an IC50 (nM).

    The reference model treats drug transport-inhibition as competitive
    (Ki = IC50/2); the dataset IC50s were measured in suspension-cultured
    human hepatocytes (SHH) as a worst-case estimate.
    """
    if ic50_nm <= 0:
        raise ValueError("ic50_nm must be positive")
    return ic50_nm / 2000.0


def bile_acid_stress(fold_ratio: float, risk_fold: float = 1.5) -> float:
    """Map an intrahepatic bile-acid fold-ratio onto a 0..1 cholestasis stress.

    ``fold_ratio`` is the drug-driven GCDCA pool increase over baseline.  Below
    unity there is no cholestasis; the stress rises sigmoidally to 1 beyond
    the validated ``risk_fold`` threshold (1.5-fold, de Bruijn & Rietjens
    2024).
    """
    if risk_fold <= 1.0:
        raise ValueError("risk_fold must be > 1")
    s = max(0.0, fold_ratio - 1.0)
    denom = s * s + (risk_fold - 1.0) * (risk_fold - 1.0)
    return s * s / denom


def simulate_gcdca_pbk(
    t_h: NDArray,
    c_free_umol_l: NDArray,
    ki_umol_l: float,
    params: BileAcidParams | None = None,
    rtol: float = 1e-6,
    atol: float = 1e-8,
) -> NDArray:
    """Intrahepatic GCDCA fold-ratio over time under BSEP inhibition.

    Ports the reference GDCCA PBK (NTCP uptake, BSEP canalicular efflux with
    competitive drug inhibition, ASBT ileal reabsorption, gallbladder and
    intestinal compartments) driven by the supplied *free hepatic* drug
    exposure.  Returns the intrahepatic intracellular-water GCDCA
    concentration relative to the no-drug baseline peak.
    """
    t = np.asarray(t_h, dtype=float)
    c = np.asarray(c_free_umol_l, dtype=float)
    if t.ndim != 1 or c.ndim != 1 or t.shape[0] != c.shape[0]:
        raise ValueError("t_h and c_free_umol_l must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")
    if ki_umol_l <= 0:
        raise ValueError("ki_umol_l must be positive")
    if t[-1] <= t[0]:
        raise ValueError("t_h must be increasing")
    p = params or BileAcidParams()

    bw = p.bw
    vad = bw * p.fv_ad
    vgu = bw * p.fv_gu
    vliw = bw * p.fv_liw
    vlew = bw * p.fv_lew
    vpv = bw * p.fv_pv
    vbl = bw * (p.fv_ve + p.fv_ar - p.fv_pv)
    vslp = bw * (p.fv_bo + p.fv_sk + p.fv_re)
    vra = bw * (1.0 - p.fv_ad - p.fv_gu - p.fv_li - p.fv_ve - p.fv_ar - p.fv_bo - p.fv_sk - p.fv_re)

    qad = p.qc * p.fq_ad
    qgu = p.qc * p.fq_gu
    qh = p.qc * p.fq_h
    qsp = p.qc * p.fq_sp
    qpv = qgu + qsp
    qha = qh - qpv
    qslp = p.qc * (p.fq_bo + p.fq_sk + p.fq_re)
    qrp = p.qc * (1.0 - p.fq_ad - p.fq_gu - p.fq_h - p.fq_bo - p.fq_sk - p.fq_re)

    vmax_bsep = p.vmax_bsep_c * (p.a_bsep * p.mw_bsep * p.hep * p.wl * 1e-9) * 60.0
    vmax_ntcp = p.vmax_ntcp_c * (p.hep * p.wl * 1e-6) * p.sf_oatp * 60.0
    vmax_asbt = p.vmax_asbt_c * p.sf_asbt * 60.0
    qib = 1.0 - p.qgb

    iw, ew = 2, 3

    def efflux_km(sol_t: float, drug: bool) -> float:
        if not drug:
            return p.km_bsep
        drive = float(np.interp(sol_t, t, c))
        return p.km_bsep * (1.0 + drive / ki_umol_l)

    def rhs(sol_t: float, y: NDArray, km_fn: _KM_FN) -> NDArray:
        cad = y[0] / vad
        cgu = y[1] / vgu
        ciw = y[iw] / vliw
        cew = y[ew] / vlew
        cportal = y[4] / vpv
        cslp = y[5] / vslp
        cra = y[6] / vra
        cileum = y[9] / p.vil
        cblood = y[11] / vbl

        lewfrac = cew * p.bp_ba / p.kp_lew_cba
        uptake_ntcp = vmax_ntcp * lewfrac / (p.km_ntcp + lewfrac)
        efflux = vmax_bsep * ciw / (km_fn(sol_t) + ciw)
        uptake_asbt = vmax_asbt * cileum / (p.km_asbt + cileum)
        empties = p.ge * (1.0 if sol_t % 24.0 < 1.5 else 0.0) * y[7]
        ktj = p.ktj_fasted
        kti = p.kti_fasted

        d = np.empty(12, dtype=float)
        d[0] = qad * (cblood - cad / p.kp_ad_cba * p.bp_ba)
        d[1] = qgu * (cblood - cgu / p.kp_gu_cba * p.bp_ba)
        d[2] = uptake_ntcp + p.fec * y[10] - efflux
        d[3] = qha * (cblood - lewfrac) + qpv * (cportal - lewfrac) - uptake_ntcp
        d[4] = p.ka_cba * y[8] + uptake_asbt + p.ka_cba * y[10] - qpv * cportal
        d[5] = qslp * (cblood - cslp / p.kp_slp_cba * p.bp_ba)
        d[6] = qrp * (cblood - cra / p.kp_ra_cba * p.bp_ba)
        d[7] = -empties + efflux * p.qgb
        d[8] = empties + efflux * qib - p.ka_cba * y[8] - ktj * y[8]
        d[9] = ktj * y[8] - kti * y[9] - uptake_asbt
        d[10] = kti * y[9] - p.ka_cba * y[10] - p.fec * y[10]
        d[11] = (
            qad * cad / p.kp_ad_cba * p.bp_ba
            + qgu * cgu / p.kp_gu_cba * p.bp_ba
            + qh * cew / p.kp_lew_cba * p.bp_ba
            + qslp * cslp / p.kp_slp_cba * p.bp_ba
            + qrp * cra / p.kp_ra_cba * p.bp_ba
            - (qad + qha + qslp + qrp + qgu) * cblood
        )
        return d

    y0 = np.zeros(12, dtype=float)
    y0[7] = p.gdose

    def integrate(
        y: NDArray, t_start: float, t_end: float, km_fn: _KM_FN
    ) -> tuple[NDArray, NDArray]:
        times: list[list[float]] = []
        rows: list[NDArray] = []
        day = t_start
        while day < t_end:
            seg_end = min(day + 24.0, t_end)
            t_seg = np.arange(day, seg_end + 1e-9, 0.1)
            t_seg = t_seg[t_seg <= seg_end]
            if t_seg[-1] < seg_end:
                t_seg = np.concatenate((t_seg, np.asarray([seg_end], dtype=float)))
            sol = solve_ivp(
                rhs,
                (day, seg_end),
                y,
                args=(km_fn,),
                t_eval=t_seg,
                method="LSODA",
                rtol=rtol,
                atol=atol,
            )
            if not sol.success:
                raise RuntimeError(f"bile-acid PBK solve failed: {sol.message}")
            y = sol.y[:, -1]
            y[7] = p.gdose
            times.append(list(t_seg))
            rows.append(sol.y[:, :])
            day += 24.0
        time_arr = np.asarray([v for seg in times for v in seg], dtype=float)
        state: NDArray = np.concatenate(rows, axis=1)
        return time_arr, state

    def ctrl_km(sol_t: float) -> float:
        return efflux_km(sol_t, drug=False)

    def drug_km(sol_t: float) -> float:
        return efflux_km(sol_t, drug=True)

    # Warm up the enterohepatic loop to its converged no-drug baseline.
    _, y_base = integrate(y0, 0.0, 72.0, ctrl_km)
    baseline = y_base[:, -1].copy()

    time_arr, state_ctrl = integrate(baseline, t[0], t[-1], ctrl_km)
    _, state_drug = integrate(baseline, t[0], t[-1], drug_km)
    baseline_peak = float(np.max(state_ctrl[iw] / vliw))
    if baseline_peak <= 0:
        raise RuntimeError("bile-acid PBK produced a non-positive baseline pool")
    ratio = (state_drug[iw] / vliw) / baseline_peak
    interp = np.interp(t, time_arr, ratio)
    return np.asarray(interp, dtype=float)


def combined_stress(
    cholestasis: float,
    atp_frac: float,
    gsh_frac: float,
    immune: float = 0.0,
    weights: tuple[float, float, float] = (0.5, 0.3, 0.2),
    immune_weight: float = 0.0,
) -> float:
    """Normalized 0..1 hepatocyte stress from the toxicity axes.

    The fourth (immune) axis is wired: the ``immune`` hazard input
    is weighted by ``immune_weight``, which defaults to 0 so every shipped
    default run keeps the validated cholestasis/ATP/GSH composition exactly.
    Full runs engage the axis through the ``dili_immune_*`` anchors
    (doc/12 D14; validation case ``case_dili_immune_activation``).
    """
    w_immune = max(0.0, immune_weight)
    w_sum = sum(max(0.0, w) for w in weights) + w_immune
    if w_sum <= 0:
        raise ValueError("at least one stress weight must be positive")
    return max(
        0.0,
        (
            max(0.0, weights[0]) * cholestasis
            + max(0.0, weights[1]) * max(0.0, 1.0 - atp_frac)
            + max(0.0, weights[2]) * max(0.0, 1.0 - gsh_frac)
            + w_immune * max(0.0, immune)
        )
        / w_sum,
    )


def _kill_rate(stress: float, kill_max_1h: float, ec50: float, hill: float) -> float:
    if ec50 <= 0 or hill <= 0 or kill_max_1h <= 0:
        raise ValueError("kill constants must be positive")
    if stress <= 0:
        return 0.0
    s_n = stress**hill
    life_n = ec50**hill
    return float(kill_max_1h * s_n / (s_n + life_n))


def regeneration_scale(
    fold_change: NDArray | float,
    gain: float = 0.5,
    floor: float = 0.5,
    ceiling: float = 1.5,
) -> NDArray:
    """Map an ERK proliferative fold-change onto a regeneration-rate scale.

    ``fold_change`` is the Stage-3 pathway readout relative to its drug-free
    baseline (1.0 = unperturbed signal).  Higher readout (more MAPK/ERK
    mitogenic drive) accelerates regeneration; a drug-suppressed readout
    attenuates it.  The linear relation is damped by ``gain`` and clamped to
    ``[floor, ceiling]`` so a pathway perturbation can never fully ablate (or
    infinitely potentiate) hepatocyte regeneration.
    """
    if gain <= 0:
        raise ValueError("gain must be positive")
    if floor < 0.0 or ceiling < floor:
        raise ValueError("floor must be non-negative and ceiling >= floor")
    fc = np.asarray(fold_change, dtype=float)
    scale = 1.0 + gain * (fc - 1.0)
    return np.clip(scale, floor, ceiling)


def _death_rhs(dead: float, stress: float, p: LiverParams, regen_scale: float = 1.0) -> float:
    kill = _kill_rate(stress, p.kill_max_1h, p.kill_ec50, p.kill_hill)
    return kill * (1.0 - dead) - p.regeneration_1h * regen_scale * dead


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
    proliferation_signal: NDArray | None = None,
) -> LiverTrajectory:
    """Integrate the four-axis liver QST model on the PBPK exposure grid.

    The cholestasis axis is anchored to the production-validated GCDCA
    bile-acid PBK of de Bruijn & Rietjens (2024): drug free-hepatic exposure
    competitively inhibits BSEP efflux and the intrahepatic bile-acid pool
    fold-ratio above the validated 1.5x risk threshold drives ``cholestasis``
    (doc/12 row 4c, R-7).  The Mito/redox/hepatocyte-death axes remain the
    documented calibration model (DILIsym-equivalent closiness is proprietary);
    ALT/AST release follows cell death.

    ``proliferation_signal`` (optional, on ``t_h``) is the Stage-3 pathway
    readout fold-change over time; it scales the hepatocyte regeneration rate
    via :func:`regeneration_scale`, closing the occupancy -> pathway -> organ
    chain (a suppressed ERK/MAPK signal slows regeneration and lets the same
    direct stress accumulate to more cell death).  The total-bilirubin rise is
    capped at ``bile_rise_max_fold`` x ULN (the declared cholestasis ceiling;
    the un-clamped rise above it stays in ``stress``).

    **Immune-mediated axis** (model layer default-off; auto-anchored in full
    fidelity — ``IC50_immune = 0.30 * DILI IC50`` or a disclosed null effect):
    when ``LiverParams.immune_ic50_nm``
    is set, the hapten/danger hazard :func:`immune_hazard` drives an adaptive
    immune-response state ``imm`` (ODEs ``dimm = k_recruit*hazard*(1-imm) -
    k_decay*imm``) whose level loads the fourth ``combined_stress`` axis via
    ``immune_weight``.  With ``immune_weight`` or ``immune_ic50_nm`` at their
    defaults the axis is inert and the trajectory is the validated
    cholestasis/ATP/GSH composition exactly.
    """
    p = params or LiverParams()
    t = np.asarray(t_h, dtype=float)
    c = free_mg_l_to_nm(np.asarray(c_free_mg_l, dtype=float), mw)
    if t.ndim != 1 or c.ndim != 1 or t.shape[0] != c.shape[0]:
        raise ValueError("t_h and c_free_mg_l must be equal-length 1-D arrays")
    if t.shape[0] < 2:
        raise ValueError("t_h needs at least two time points")
    if proliferation_signal is not None:
        s = np.asarray(proliferation_signal, dtype=float)
        if s.ndim != 1 or s.shape[0] != t.shape[0]:
            raise ValueError("proliferation_signal must be equal-length to t_h")
    if p.immune_ic50_nm is not None and p.immune_ic50_nm <= 0:
        raise ValueError("immune_ic50_nm must be positive")
    if p.immune_hill <= 0 or p.immune_recruit_1h < 0 or p.immune_decay_1h < 0:
        raise ValueError("immune constants must be positive/non-negative")

    chol_raw = simulate_gcdca_pbk(
        t,
        np.asarray(c, dtype=float) / 1000.0,
        bsep_ki_from_ic50_nm(p.bsep_ic50_nm),
    )
    chol_g = np.array([bile_acid_stress(float(f)) for f in chol_raw], dtype=float)
    mito_g = np.array([mitochondrial_block(float(cc), p.mito_ic50_nm) for cc in c], dtype=float)
    atp_g = np.array(
        [aten_floor_factor(float(m), p.atp_floor, p.kexpz_atp * float(m)) for m in mito_g],
        dtype=float,
    )
    gsh_g = np.array([redox_state(float(cc), p.redox_ic50_nm)[1] for cc in c], dtype=float)
    imm_g = np.array(
        [
            immune_hazard(float(cc), p.immune_ic50_nm) if p.immune_ic50_nm is not None else 0.0
            for cc in c
        ],
        dtype=float,
    )

    regen_g = (
        np.ones_like(t, dtype=float)
        if proliferation_signal is None
        else regeneration_scale(
            proliferation_signal,
            p.regen_pathway_gain,
            p.regen_pathway_min,
            p.regen_pathway_max,
        )
    )

    def stress_at(sol_t: float, imm: float) -> float:
        chol = float(np.interp(sol_t, t, chol_g))
        atp_frac = float(np.interp(sol_t, t, atp_g))
        gsh = float(np.interp(sol_t, t, gsh_g))
        return combined_stress(chol, atp_frac, gsh, immune=imm, immune_weight=p.immune_weight)

    def regen_at(sol_t: float) -> float:
        return float(np.interp(sol_t, t, regen_g))

    def immune_at(sol_t: float) -> float:
        return float(np.interp(sol_t, t, imm_g))

    def rhs(sol_t: float, y: NDArray) -> NDArray:
        dead = float(y[0])
        imm = float(y[1])
        d_dead = _death_rhs(dead, stress_at(sol_t, imm), p, regen_at(sol_t))
        d_imm = p.immune_recruit_1h * immune_at(sol_t) * (1.0 - imm) - p.immune_decay_1h * imm
        return np.array([d_dead, d_imm], dtype=float)

    y0 = np.array([0.0, 0.0], dtype=float)
    t_eval = np.linspace(t[0], t[-1], n_eval)
    sol = solve_ivp(rhs, (t[0], t[-1]), y0, t_eval=t_eval, method="LSODA", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"liver solve failed: {sol.message}")
    dead = sol.y[0]
    imm = sol.y[1]

    chol: NDArray = np.interp(t_eval, t, chol_g)
    atp: NDArray = np.interp(t_eval, t, atp_g)
    gsh: NDArray = np.interp(t_eval, t, gsh_g)
    stress: NDArray = np.array(
        [stress_at(float(tp), float(np.interp(tp, t_eval, imm))) for tp in t_eval], dtype=float
    )
    free: NDArray = np.interp(t_eval, t, c)
    regen_traj: NDArray = np.interp(t_eval, t, regen_g)

    alt: NDArray = p.alt_uln_u_l * (1.0 + p.alt_release_per_dead * dead)
    ast: NDArray = p.ast_uln_u_l * (1.0 + p.ast_release_per_dead * dead)
    bili_raw: NDArray = p.bilirubin_uln_mg_dl * (
        1.0 + p.bili_rise_per_bsep * chol + p.bili_rise_per_dead * dead
    )
    bili: NDArray = np.minimum(bili_raw, p.bilirubin_uln_mg_dl * p.bile_rise_max_fold)
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
        immune=imm,
        dead_frac=dead,
        regen_scale=regen_traj,
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
    "BileAcidParams",
    "LiverParams",
    "LiverStress",
    "LiverTrajectory",
    "aten_floor_factor",
    "bile_acid_stress",
    "bsep_ki_from_ic50_nm",
    "combined_stress",
    "dili_grade",
    "immune_hazard",
    "inhibition",
    "liver_params_from_panel",
    "mitochondrial_block",
    "redox_state",
    "regeneration_scale",
    "simulate_gcdca_pbk",
    "simulate_liver",
]
