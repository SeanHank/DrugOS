"""Off-target resolver unit tests (drugos.target.resolver + pipeline wiring).

Pins the corpus-calibrated monotone P -> KD curve (doc/12 D10, path B) and the
specialist-head panel re-scoring seam (path A): never more potent than the
panel prior, exact anchors, strict monotonicity, measured overrides absolute,
and every non-inhibition site untouched.
"""

from __future__ import annotations

import math
from dataclasses import replace

import pytest

import drugos.pipeline as pl
from drugos.pk.admet import AdmetOutput
from drugos.target.resolver import (
    KD_WEAK_NM,
    bind_site,
    herg_scored_kd_nm,
    kd_from_score,
    resolve_admet_panel,
)
from drugos.target.targets import Target, safety_panel

_PRIOR = 2.0


def _warfarin() -> pl.RunSpec:
    return pl.spec_from_benchmark_data(pl.benchmark_data("warfarin"))


def test_kd_from_score_anchors() -> None:
    assert kd_from_score(0.0, _PRIOR) == pytest.approx(KD_WEAK_NM)
    assert kd_from_score(1.0, _PRIOR) == pytest.approx(_PRIOR)
    mid = math.sqrt(_PRIOR * KD_WEAK_NM)
    assert kd_from_score(0.5, _PRIOR) == pytest.approx(mid)


def test_kd_from_score_strictly_monotone_and_conservative() -> None:
    values = [kd_from_score(p, _PRIOR) for p in (0.0, 0.1, 0.25, 0.4, 0.5, 0.75, 0.9, 0.99, 1.0)]
    for lo, hi in zip(values, values[1:], strict=False):
        assert lo > hi
        assert lo <= KD_WEAK_NM
    assert all(v >= _PRIOR for v in values)  # never more potent (smaller KD) than the prior
    assert all(v > 0.0 for v in values)


def test_kd_from_score_clamps_probability() -> None:
    assert kd_from_score(-1.0, _PRIOR) == pytest.approx(KD_WEAK_NM)
    assert kd_from_score(1.7, _PRIOR) == pytest.approx(_PRIOR)


def test_kd_from_score_rejects_bad_anchors() -> None:
    with pytest.raises(ValueError):
        kd_from_score(0.5, 0.0)
    with pytest.raises(ValueError):
        kd_from_score(0.5, 2.0, weak_kd_nm=0.0)


def test_herg_scored_kd_nm_passes_none_through() -> None:
    assert herg_scored_kd_nm(None, _PRIOR) is None
    assert herg_scored_kd_nm(0.9, _PRIOR) == pytest.approx(kd_from_score(0.9, _PRIOR))


def test_resolve_admet_panel_without_admet_is_identity() -> None:
    panel = safety_panel()
    assert resolve_admet_panel(panel, None) == panel


def test_resolve_admet_panel_rescores_cyp_inhibition_sites() -> None:
    admet = AdmetOutput(
        smiles="c1ccccc1",
        hERG=0.95,
        cyp2d6_inhibitor=0.8,
        cyp3a4_inhibitor=0.6,
        cyp2c9_inhibitor=None,
        Pgp=0.9,
    )
    out = resolve_admet_panel(safety_panel(), admet)
    for t in out:
        if t.name == "CYP2D6 inhibition":
            assert t.kd_nm == pytest.approx(kd_from_score(0.8, 15.0e3))
            assert t.low_confidence
            assert "CYP2D6 head" in t.reference and "D10" in t.reference
        elif t.name == "CYP3A4 inhibition":
            assert t.kd_nm == pytest.approx(kd_from_score(0.6, 12.0e3))
        elif t.name == "CYP2C9 inhibition":
            # Missing head probability leaves the class prior untouched.
            assert t.kd_nm == pytest.approx(15.0e3)


def test_resolve_admet_panel_never_touches_pgp_and_other_sites() -> None:
    admet = AdmetOutput(smiles="c1ccccc1", hERG=0.99, Pgp=0.99, DILI=0.9)
    before = {t.name: (t.kd_nm, t.reference, t.low_confidence) for t in safety_panel()}
    out = resolve_admet_panel(safety_panel(), admet)
    after = {t.name: (t.kd_nm, t.reference, t.low_confidence) for t in out}
    for name in (
        "P-gp (MDR1)",
        "hERG (Kv11.1)",
        "BSEP (cholestasis)",
        "OATP1B1",
        "Mitochondrial complex I",
        "Estrogen receptor",
    ):
        assert after[name] == before[name]


def test_measured_override_bypasses_curve() -> None:
    spec = replace(_warfarin(), qt_ic50_nm=12.5)
    assert pl._effective_herg_kd_nm(spec) == pytest.approx(12.5)


def test_pipeline_effective_herg_branches() -> None:
    base = replace(_warfarin(), qt_ic50_nm=None)
    # No ADMET and a missing head both fall back to the panel class prior.
    assert pl._effective_herg_kd_nm(base) == pytest.approx(2.0)
    spec = replace(base, admet=AdmetOutput(smiles="CN", hERG=None))
    assert pl._effective_herg_kd_nm(spec) == pytest.approx(2.0)
    # With a head probability the calibrated curve applies.
    spec = replace(base, admet=AdmetOutput(smiles="CN", hERG=0.75))
    assert pl._effective_herg_kd_nm(spec) == pytest.approx(kd_from_score(0.75, 2.0))
    # No hERG site in the panel -> no effective KD, panel left untouched.
    no_herg = replace(base, panel=(Target(name="sodium channel", kd_nm=5.0e4),))
    spec = replace(no_herg, admet=AdmetOutput(smiles="CN", hERG=0.9))
    assert pl._effective_herg_kd_nm(spec) is None
    assert pl._herg_sieved_panel(spec) == spec.panel


def test_pipeline_sieved_panel_no_override_binds_prior() -> None:
    base = replace(_warfarin(), qt_ic50_nm=None)
    sieved = pl._herg_sieved_panel(base)
    assert next(t.kd_nm for t in sieved if "hERG" in t.name) == pytest.approx(2.0)


def test_bind_site() -> None:
    panel = safety_panel()
    out = bind_site(panel, "hERG (Kv11.1)", 9.9, "test ref")
    assert sum(1 for t in out if "hERG" in t.name) == 1
    herg = next(t for t in out if "hERG" in t.name)
    assert herg.kd_nm == pytest.approx(9.9)
    assert herg.reference == "test ref"
    assert herg.low_confidence
    # Unmatched name returns the panel unchanged.
    assert bind_site(panel, "not a site", 3.0, "x") == panel
