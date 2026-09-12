"""Off-target resolution: production-model heads re-score the safety panel.

For a novel molecule scored by ADMET-AI but carrying no measured affinity,
the safety panel of ``targets.safety_panel`` is class-typical by construction.
This module supplies the *two-path* replacement for the bare binary hERG sieve
(doc/12 §1 path B and path A, decision record D10):

- **path B — head-as-specialist:** a single corpus-calibrated, monotone
  probability -> affinity curve maps each classifier head that measures the
  same biological interaction named by a panel site onto a continuous KD.
  Today that is the hERG head and the Veith CYP2D6/3A4/2C9 inhibition heads;
  the P-glycoprotein head is a *substrate* classifier (substrate probability
  is not inhibition potency) so the P-gp site stays on its class prior, as do
  BSEP/OATP1B1/MRP/MRC, the nuclear receptors and any primary-target site —
  those await the sequence-DTI resolver of path A,
- **path A — general resolver seam:** ``resolve_admet_panel`` is the single
  hook the future DeepDTA-family estimator (doc/12 §1 row 2b, P5) will fill in
  for the remaining sites; the G4 factory rule requires any such replacement
  to ship a corpus-calibration + equivalence case before it is admitted.

The curve itself::

    log10 KD(p) = log10(KD_weak) - (log10(KD_weak) - log10(KD_prior)) * p

with ``p`` the clamped head probability, ``KD_weak`` the weak-interaction floor
(1 mM) and ``KD_prior`` the site's conservative class prior.  It is strictly
decreasing in ``p``, continuous, and *never more potent than the panel prior*:

- ``p = 1`` reproduces the panel prior exactly (flagship conservative anchor),
- ``p = 0`` reproduces the weak floor exactly,
- ``p = 0.5`` (the classifier threshold) lands on the geometric mean of the two
  anchors — the corpus-typical weak potency of the hERG Central corpus
  (doc/08 R-2: median % inhibition at 1 uM is small, so the typical IC50 sits
  around the log-midpoint between a strong blocker and the floor).

All invariants are pinned in ``tests/test_target_resolver.py`` and the L2
validation case ``case_herg_calibration`` (doc/08 R-8).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from drugos.target.targets import Target

KD_WEAK_NM = 1.0e6  # 1 mM weak-interaction floor (doc/05 2.2, formerly _HERG_WEAK_KD_NM)

# Panel sites whose ADMET-AI head is an inhibition classifier of the same
# interaction named by the site.  hERG is resolved by the pipeline sieve (its
# measured ``qt_ic50_nm`` override must stay absolute), so it is not listed
# here; P-gp and every sequence/transporter/nuclear-receptor site are excluded
# on purpose (see module docstring) and stay on their class priors.
_SITE_TO_HEAD: dict[str, str] = {
    "CYP2D6 inhibition": "cyp2d6_inhibitor",
    "CYP3A4 inhibition": "cyp3a4_inhibitor",
    "CYP2C9 inhibition": "cyp2c9_inhibitor",
}


def kd_from_score(
    prob: float,
    prior_kd_nm: float,
    weak_kd_nm: float = KD_WEAK_NM,
) -> float:
    """Map a clamped head probability in [0, 1] onto a continuous KD (nM).

    Monotone (never more potent than the prior), continuous, with the exact
    anchors ``p=1 -> prior_kd_nm``, ``p=0 -> weak_kd_nm`` (module docstring).
    """
    p = min(max(prob, 0.0), 1.0)
    if prior_kd_nm <= 0 or weak_kd_nm <= 0:
        raise ValueError("KD anchors must be positive")
    lo = math.log10(weak_kd_nm)
    hi = math.log10(prior_kd_nm)
    return float(10.0 ** (lo - (lo - hi) * p))


def herg_scored_kd_nm(
    prob: float | None,
    prior_kd_nm: float,
    weak_kd_nm: float = KD_WEAK_NM,
) -> float | None:
    """hERG KD from the ADMET-AI hERG head, or ``None`` for a missing head."""
    return kd_from_score(prob, prior_kd_nm, weak_kd_nm) if prob is not None else None


def _rebind(
    targets: Sequence[Target],
    name: str,
    kd_nm: float,
    reference: str,
    low_confidence: bool = True,
) -> tuple[Target, ...]:
    out: list[Target] = []
    for t in targets:
        if t.name != name:
            out.append(t)
            continue
        out.append(
            Target(
                name=t.name,
                kd_nm=kd_nm,
                kon_nm_h=t.kon_nm_h,
                r0_nm=t.r0_nm,
                rho_h=t.rho_h,
                kint_h=t.kint_h,
                low_confidence=low_confidence,
                reference=reference,
            )
        )
    return tuple(out)


def resolve_admet_panel(panel: Sequence[Target], admet: object) -> tuple[Target, ...]:
    """Re-score the panel sites backed by an ADMET-AI inhibition head.

    A missing ``admet`` or a head that returned ``None`` leaves those sites on
    their class priors; every other site is re-bound to
    ``kd_from_score(head_prob, site.prior)`` under the *never-more-potent*
    ceiling, flagged ``low_confidence=True`` (ML-derived, not measured) with a
    traceable head reference.
    """
    if admet is None:
        return tuple(panel)
    out: list[Target] = []
    for t in panel:
        head = _SITE_TO_HEAD.get(t.name)
        if head is not None:
            prob = getattr(admet, head, None)
            if prob is not None:
                kd = kd_from_score(float(prob), t.kd_nm)
                out.append(
                    Target(
                        name=t.name,
                        kd_nm=kd,
                        kon_nm_h=t.kon_nm_h,
                        r0_nm=t.r0_nm,
                        rho_h=t.rho_h,
                        kint_h=t.kint_h,
                        low_confidence=True,
                        reference=(
                            f"ADMET-AI {head.split('_')[0].upper()} head re-scored "
                            "panel prior (corpus-calibrated monotone P->KD, doc/12 D10)"
                        ),
                    )
                )
                continue
        out.append(t)
    return tuple(out)


def bind_site(
    targets: Sequence[Target],
    name: str,
    kd_nm: float,
    reference: str,
    low_confidence: bool = True,
) -> tuple[Target, ...]:
    """Rebind one named site to ``kd_nm`` everywhere it appears.

    ``low_confidence`` defaults to ``True`` (a surrogate or class estimate);
    pass ``low_confidence=False`` for an absolute measured override so the
    trust record does not disclose a measured anchor as a class prior.
    """
    return _rebind(targets, name, kd_nm, reference, low_confidence)


__all__ = [
    "KD_WEAK_NM",
    "bind_site",
    "herg_scored_kd_nm",
    "kd_from_score",
    "resolve_admet_panel",
]
