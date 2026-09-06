"""Composite toxicity scoring with three-line evidence fusion (Stage 5, D19).

Implements doc/05 5.2: an endpoint risk model combining (i) mechanistic
organ output from Stage 4 (CTCAE-style grades), (ii) an exposure-ratio line

    log10( Cmax_unbound / in-vitro IC50 )

fitted to published QST/toxicity ROC anchors (literature AUC ~0.91-0.96,
doc/05 5.2), and (iii) structural ADMET-AI flags.  Lines are fused in
log-odds space with per-line trust weights and a reported uncertainty
interval, following the Baldwin 2017 logistic/odds-ratio meta-analysis
framing (doc/05 5.2 evidence table).
"""

from __future__ import annotations

import enum
import math
from collections.abc import Mapping, Sequence, Set
from dataclasses import dataclass

from scipy.stats import beta as beta_dist

from drugos.clinical.biomarkers import SEVERITY_LABELS
from drugos.pk.admet import AdmetOutput

_MIN_P = 1e-6


class Endpoint(enum.StrEnum):
    """Toxicity endpoints scored by the composite model."""

    DILI = "dili"
    QT = "qt"
    AKI = "aki"
    CNS = "cns"


class EvidenceKind(enum.StrEnum):
    """The three evidence lines of doc/05 5.2."""

    MECHANISTIC = "mechanistic"
    EXPOSURE = "exposure_ratio"
    STRUCTURAL = "structural"


@dataclass(frozen=True, slots=True)
class EndpointRule:
    """Model parameters for one endpoint's three evidence lines.

    ``grade_probs`` is the mechanistic grade -> risk mapping (grade 0..4);
    ``ratio_center`` is the log10(Cmax_unbound/IC50) at 50% exposure risk and
    ``ratio_slope`` its steepness; ``structural_field`` is the AdmetOutput
    attribute that supplies the structural probability (``None`` when the
    structural line is unavailable for the endpoint).
    """

    label: str
    prior: float
    grade_probs: tuple[float, float, float, float, float]
    ratio_center: float
    ratio_slope: float
    structural_field: str | None


RULES: dict[Endpoint, EndpointRule] = {
    Endpoint.DILI: EndpointRule(
        label="Drug-induced liver injury",
        prior=0.25,
        grade_probs=(0.03, 0.06, 0.18, 0.45, 0.80),
        ratio_center=-0.5,
        ratio_slope=1.2,
        structural_field="DILI",
    ),
    Endpoint.QT: EndpointRule(
        label="QT prolongation / TdP",
        prior=0.20,
        grade_probs=(0.02, 0.06, 0.18, 0.45, 0.85),
        ratio_center=-1.3,
        ratio_slope=1.2,
        structural_field="hERG",
    ),
    Endpoint.AKI: EndpointRule(
        label="Acute kidney injury",
        prior=0.22,
        grade_probs=(0.02, 0.06, 0.18, 0.45, 0.80),
        ratio_center=-0.5,
        ratio_slope=1.2,
        structural_field=None,
    ),
    Endpoint.CNS: EndpointRule(
        label="CNS liability",
        prior=0.20,
        grade_probs=(0.03, 0.06, 0.18, 0.45, 0.75),
        ratio_center=-0.5,
        ratio_slope=1.2,
        structural_field="BBB",
    ),
}

MECHANISTIC_WEIGHT = 1.2
EXPOSURE_WEIGHT = 1.0
STRUCTURAL_WEIGHT = 0.7

RISK_GRADE_RANGE: tuple[float, float, float, float] = (0.35, 0.5, 0.7, 0.9)


@dataclass(frozen=True, slots=True)
class Evidence:
    """One evidence line contributing to an endpoint risk.

    ``probability`` is the line's raw P(endpoint toxicity) and
    ``baseline_probability`` the probability of the *neutral* reference
    condition of that line; ``fuse`` applies the log-odds ratio
    logit(p) - logit(base) scaled by ``weight``, so the prior is only
    moved when a line departs from its neutral reference.
    """

    kind: EvidenceKind
    endpoint: Endpoint
    probability: float
    weight: float
    baseline_probability: float = 0.5
    note: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind.value,
            "endpoint": self.endpoint.value,
            "probability": round(self.probability, 4),
            "weight": self.weight,
            "baseline_probability": round(self.baseline_probability, 4),
            "note": self.note,
        }


@dataclass(frozen=True, slots=True)
class EndpointRisk:
    """Fused risk for one endpoint."""

    endpoint: Endpoint
    label: str
    risk: float
    ci_lo: float
    ci_hi: float
    grade: int
    severity: str
    driver: EvidenceKind
    evidence: tuple[Evidence, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "endpoint": self.endpoint.value,
            "label": self.label,
            "risk": round(self.risk, 4),
            "ci": [round(self.ci_lo, 4), round(self.ci_hi, 4)],
            "grade": self.grade,
            "severity": self.severity,
            "driver": self.driver.value,
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(frozen=True, slots=True)
class ToxicityReport:
    """The full composite toxicity output consumed by the report layer."""

    risks: tuple[EndpointRisk, ...]

    def by_endpoint(self, endpoint: Endpoint) -> EndpointRisk | None:
        for r in self.risks:
            if r.endpoint is endpoint:
                return r
        return None

    def overall_risk(self) -> float:
        return max((r.risk for r in self.risks), default=0.0)

    def overall_driver(self) -> str:
        worst = max(self.risks, key=lambda r: r.risk, default=None)
        if worst is None:
            return "no_evidence"
        return worst.endpoint.value

    def to_dict(self) -> dict[str, object]:
        return {
            "overall_risk": round(self.overall_risk(), 4),
            "overall_driver": self.overall_driver(),
            "endpoints": [r.to_dict() for r in self.risks],
        }


def clamp_prob(p: float) -> float:
    """Clamp a probability to the open (0,1) interval used by log-odds math."""
    return max(_MIN_P, min(1.0 - _MIN_P, p))


def logit(p: float) -> float:
    return math.log(clamp_prob(p) / (1.0 - clamp_prob(p)))


def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def mechanistic_evidence(
    endpoint: Endpoint, grade: int, weight: float = MECHANISTIC_WEIGHT
) -> Evidence:
    """Evidence line from the Stage-4 organ grade (doc/05 5.2 line 1)."""
    rule = RULES[endpoint]
    grade = min(4, max(0, grade))
    prob = rule.grade_probs[grade]
    return Evidence(
        kind=EvidenceKind.MECHANISTIC,
        endpoint=endpoint,
        probability=prob,
        weight=weight,
        baseline_probability=rule.grade_probs[0],
        note=(f"mechanistic organ QST grade {grade} (severity: {SEVERITY_LABELS[grade]})"),
    )


def exposure_evidence(
    endpoint: Endpoint,
    cmax_unbound_nm: float,
    ic50_nm: float,
    weight: float = EXPOSURE_WEIGHT,
) -> Evidence:
    """Exposure-ratio line: log10(Cmax_unbound / in-vitro IC50) ROC anchor."""
    if ic50_nm <= 0:
        raise ValueError("ic50_nm must be positive")
    if cmax_unbound_nm < 0:
        raise ValueError("cmax_unbound_nm must be non-negative")
    rule = RULES[endpoint]
    ratio = math.log10(cmax_unbound_nm / ic50_nm)
    prob = sigmoid(rule.ratio_slope * (ratio - rule.ratio_center))
    return Evidence(
        kind=EvidenceKind.EXPOSURE,
        endpoint=endpoint,
        probability=prob,
        weight=weight,
        note=f"exposure ratio Cmax_unbound/IC50 = 10^{ratio:.2f}",
    )


def structural_evidence(
    endpoint: Endpoint,
    admet: AdmetOutput,
    weight: float = STRUCTURAL_WEIGHT,
) -> Evidence:
    """Structural ADMET-AI line (doc/05 5.2 line 3); absent when no field."""
    rule = RULES[endpoint]
    field = rule.structural_field
    if field is None:
        raise ValueError(f"{endpoint.value} has no structural field")
    value = getattr(admet, field)
    if value is None:
        raise ValueError(f"ADMET prediction field '{field}' is missing")
    return Evidence(
        kind=EvidenceKind.STRUCTURAL,
        endpoint=endpoint,
        probability=clamp_prob(float(value)),
        weight=weight,
        note=f"structural ADMET-AI flag ({field})",
    )


def fuse(prior: float, evidence: Sequence[Evidence]) -> tuple[float, float, float]:
    """Fuse prior with weighted log-odds-ratio updates.

    Returns (risk, ci_lo, ci_hi).  The credible interval comes from a Beta
    posterior whose concentration ``k = 1 + sum(weights)`` shrinks as more
    (trusted) evidence lines participate.
    """
    lgt = logit(prior)
    k = 1.0
    for ev in evidence:
        delta = logit(ev.probability) - logit(ev.baseline_probability)
        lgt += ev.weight * delta
        k += max(ev.weight, 0.0)
    risk = sigmoid(lgt)
    a = 1.0 + risk * k
    b = 1.0 + (1.0 - risk) * k
    lo, hi = beta_dist.ppf(0.025, a, b), beta_dist.ppf(0.975, a, b)
    return risk, max(0.0, float(lo)), min(1.0, float(hi))


def risk_grade(risk: float) -> int:
    """Map a fused risk to the 0..4 severity ladder."""
    return min(4, sum(1 for t in RISK_GRADE_RANGE if risk >= t))


def score_toxicity(
    mechanistic_grades: Mapping[Endpoint, int],
    cmax_unbound_nm: float | None,
    ic50_nm: Mapping[Endpoint, float | None],
    admet: AdmetOutput | None,
    structural_mask: Set[Endpoint] | None = None,
) -> ToxicityReport:
    """Run every endpoint through the three evidence lines and fuse.

    ``structural_mask`` suppresses the structural ADMET-AI line for the listed
    endpoints (doc/05 5.2: the CNS line is only fused when the molecule is
    CNS-anchored, so the pipeline masks it otherwise).
    """
    masked = structural_mask if structural_mask is not None else frozenset()
    risks: list[EndpointRisk] = []
    for endpoint, rule in RULES.items():
        evidence: list[Evidence] = []
        grade = mechanistic_grades.get(endpoint)
        if grade is not None:
            evidence.append(mechanistic_evidence(endpoint, grade))
        ic50 = ic50_nm.get(endpoint)
        if cmax_unbound_nm is not None and ic50 is not None and ic50 > 0:
            evidence.append(exposure_evidence(endpoint, cmax_unbound_nm, ic50))
        if admet is not None and rule.structural_field is not None and endpoint not in masked:
            field_value = getattr(admet, rule.structural_field)
            if field_value is not None:
                evidence.append(structural_evidence(endpoint, admet))
        risk, lo, hi = fuse(rule.prior, evidence)
        driver = max(
            evidence,
            key=lambda e: abs(e.weight * logit(e.probability)),
            default=None,
        )
        driver_kind = EvidenceKind.MECHANISTIC if driver is None else driver.kind
        risks.append(
            EndpointRisk(
                endpoint=endpoint,
                label=rule.label,
                risk=risk,
                ci_lo=lo,
                ci_hi=hi,
                grade=risk_grade(risk),
                severity=SEVERITY_LABELS[risk_grade(risk)],
                driver=driver_kind,
                evidence=tuple(evidence),
            )
        )
    return ToxicityReport(risks=tuple(risks))


__all__ = [
    "Endpoint",
    "EndpointRisk",
    "EndpointRule",
    "Evidence",
    "EvidenceKind",
    "EXPOSURE_WEIGHT",
    "MECHANISTIC_WEIGHT",
    "RULES",
    "RISK_GRADE_RANGE",
    "STRUCTURAL_WEIGHT",
    "ToxicityReport",
    "clamp_prob",
    "exposure_evidence",
    "fuse",
    "logit",
    "mechanistic_evidence",
    "risk_grade",
    "score_toxicity",
    "sigmoid",
    "structural_evidence",
]
