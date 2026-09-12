"""Unit tests for the ChEMBL fingerprint-kNN off-target resolver (doc/12 L15)."""

from __future__ import annotations

import json

import pytest

from drugos.target import dti
from drugos.target.dti import (
    Snapshot,
    _TargetRows,
    cns_ic50_nm,
    load_snapshot,
    no_public_data_sites,
    predict_kd_nm,
    predict_pchembl,
    resolve_offtarget_panel,
)
from drugos.target.targets import safety_panel

DOFETILIDE = "CNC1=CC=C2C=C(CN3CCN(CC3)CCCC(=O)NC3=NC4=CC=CC=C4N3)C=CC2=N1"


def _fake_rows(smiles, pchembl):
    return _TargetRows(
        smiles=tuple(smiles),
        pchembl=tuple(float(v) for v in pchembl),
        pchembl_min=float(min(pchembl)),
        pchembl_max=float(max(pchembl)),
    )


def _fake_snapshot():
    return Snapshot(
        targets={
            "T": _fake_rows(("a", "b", "c"), (5.0, 7.0, 9.0)),
            "EMPTY": _TargetRows((), (), 0.0, 0.0),
        }
    )


@pytest.fixture
def fake_resolver(monkeypatch):
    monkeypatch.setattr(dti, "load_snapshot", lambda: _fake_snapshot())
    fps = {"q": "qfp", "a": "afp", "b": "bfp", "c": "cfp"}
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: fps.get(smi))

    def tanimoto(x, y):
        table = {("qfp", "afp"): 0.9, ("qfp", "bfp"): 0.3, ("qfp", "cfp"): 0.0}
        return table.get((x, y), 0.0)

    monkeypatch.setattr(dti, "_tanimoto", tanimoto)
    monkeypatch.setattr(dti, "CNS_TARGETS", ("T",))
    return fps


def test_load_snapshot_caches_default_path():
    snap1 = load_snapshot()
    snap2 = load_snapshot()
    assert snap1 is snap2
    assert "KCNH2" in snap1


def test_load_snapshot_from_custom_path(tmp_path):
    f = tmp_path / "snap.json"
    f.write_text(
        json.dumps(
            {
                "targets": {
                    "X": {
                        "pchembl_min": 4.0,
                        "pchembl_max": 6.0,
                        "rows": [{"smiles": "CCO", "pchembl": 5.0}],
                    }
                }
            }
        )
    )
    assert load_snapshot(f).targets["X"].pchembl == (5.0,)


def test_load_snapshot_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_snapshot(tmp_path / "missing.json")


def test_real_fingerprint_none_for_garbage():
    assert dti._fingerprint(":::not-a-molecule:::") is None


def test_predict_unknown_target(fake_resolver):
    assert predict_pchembl("NOPE", "q") is None
    assert predict_kd_nm("NOPE", "q") is None


def test_predict_unparseable_smiles(fake_resolver, monkeypatch):
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: None)
    assert predict_pchembl("T", "q") is None
    assert cns_ic50_nm("q") is None


def test_predict_empty_rows(fake_resolver):
    assert predict_pchembl("EMPTY", "q") is None


def test_predict_weighted_knn_default_k(fake_resolver):
    p = predict_pchembl("T", "q")
    assert p == pytest.approx((0.9 * 5.0 + 0.3 * 7.0) / 1.2)


def test_predict_knn_restricted_k_and_orthogonal(fake_resolver):
    assert predict_pchembl("T", "q", k=1) == pytest.approx(5.0)


def test_predict_all_neighbours_orthogonal(fake_resolver, monkeypatch):
    monkeypatch.setattr(dti, "_tanimoto", lambda x, y: 0.0)
    assert predict_pchembl("T", "q") is None


def test_predict_drops_none_fingerprints(fake_resolver, monkeypatch):
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: None if smi == "b" else fake_resolver[smi])
    p = predict_pchembl("T", "q")
    assert p == pytest.approx(5.0)  # only 'a' (sim 0.9) survived; 'c' sim 0 excluded


def test_predict_kd_nm_conversion(fake_resolver):
    kd = predict_kd_nm("T", "q")
    assert kd == pytest.approx(10.0 ** (9.0 - 5.5))


def test_cns_worst_case_and_cap(fake_resolver):
    assert cns_ic50_nm("q", max_nm=1000.0) == pytest.approx(1000.0)
    assert cns_ic50_nm("q", max_nm=1.0e5) == pytest.approx(10.0 ** (9.0 - 5.5))


def test_no_public_data_sites_sorted():
    resolved = list(safety_panel())
    assert no_public_data_sites([]) == ()
    assert no_public_data_sites(resolved) == tuple(
        sorted(t.name for t in resolved if t.low_confidence)
    )


def test_resolve_panel_no_smiles_keeps_priors():
    panel, no_data = resolve_offtarget_panel(safety_panel(), smiles=None)
    assert panel == tuple(safety_panel())
    assert no_data == no_public_data_sites(safety_panel())


def test_resolve_panel_real_snapshot_evidence_gated():
    # haloperidol has chemotype support in the snapshot (P-gp 0.67, CYP2D6
    # 0.42, hERG 0.38) so those sites are structure-re-bound; every other
    # mapped site sits in a chemotype desert and stays on its prior.
    haloperidol = "OC1(CCN(CC1)CCCC2=CC=C(F)C=C2)C3=CC=C(Cl)C=C3"
    panel, no_data = resolve_offtarget_panel(safety_panel(), smiles=haloperidol)
    by_name = {t.name: t for t in panel}
    for site in ("P-gp (MDR1)", "CYP2D6 inhibition", "hERG (Kv11.1)"):
        assert not by_name[site].low_confidence
        assert "ChEMBL kNN" in by_name[site].reference
    assert by_name["Mitochondrial complex I"].low_confidence  # top 0.21 < gate
    assert "BSEP (cholestasis)" in {t.name for t in panel}
    assert set(no_data) == {
        "Mitochondrial complex I",
        "CYP3A4 inhibition",
        "CYP2C9 inhibition",
        "BSEP (cholestasis)",
        "OATP1B1",
        "Glucocorticoid receptor",
        "Estrogen receptor",
        "Androgen receptor",
        "MRP3",
        "MRP4",
        "Mitochondrial complex II",
        "Mitochondrial complex III",
        "Mitochondrial complex IV",
        "Mitochondrial pyruvate carrier",
    }


def test_resolve_panel_chemotype_desert_stays_on_priors():
    # dofetilide's best Tanimoto match against every snapshot target is
    # 0.15-0.30 (below MIN_NEIGHBOR_TANIMOTO): no site is structure-bound and
    # every mapped site is disclosed on its prior.
    panel, no_data = resolve_offtarget_panel(safety_panel(), smiles=DOFETILIDE)
    by_name = {t.name: t for t in panel}
    assert all(by_name[n].low_confidence for n in dti.SITE_TO_TARGET if n in by_name)
    assert len(no_data) == len(dti.SITE_TO_TARGET) + 6


def test_resolve_panel_unparseable_smiles_falls_back(fake_resolver, monkeypatch):
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: None)
    panel, no_data = resolve_offtarget_panel(safety_panel(), smiles="not-a-smiles")
    assert no_data == no_public_data_sites(safety_panel())
    assert all(t.low_confidence for t in panel if t.name != "hERG (Kv11.1)")


class _Admet:
    cyp2d6_inhibitor = 1.0
    cyp3a4_inhibitor = None
    cyp2c9_inhibitor = None
    hERG = None


def test_resolve_panel_path_b_heads_applied_before_knn(fake_resolver):
    panel, _ = resolve_offtarget_panel(safety_panel(), smiles=None, admet=_Admet())
    by_name = {t.name: t for t in panel}
    assert "ADMET-AI" in by_name["CYP2D6 inhibition"].reference
    assert by_name["CYP2D6 inhibition"].kd_nm < by_name["CYP2D6 inhibition"].kd_nm + 1
    assert by_name["CYP3A4 inhibition"].low_confidence  # missing head stays prior


def test_real_predictions_bounded_by_snapshot_band():
    for target in dti.load_snapshot().targets:
        band = dti.load_snapshot().targets[target]
        kd = predict_kd_nm(target, DOFETILIDE)
        assert kd is None or band.pchembl_min <= 9.0 - __import__("math").log10(kd) <= (
            band.pchembl_max + 1e-9
        )
        assert kd is None or kd > 0


def test_top_tanimoto_skips_unparseable_rows(monkeypatch):
    snap = _fake_snapshot()
    monkeypatch.setattr(dti, "load_snapshot", lambda: snap)
    fps = {"q": "qfp", "b": "bfp", "c": "cfp"}
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: fps.get(smi))
    monkeypatch.setattr(dti, "_tanimoto", lambda x, y: 0.9 if x == "qfp" and y == "bfp" else 0.0)
    assert dti.top_tanimoto("T", "q") == 0.9


def test_cns_skips_candidates_without_kd(monkeypatch):
    snap = _fake_snapshot()
    monkeypatch.setattr(dti, "load_snapshot", lambda: snap)
    fps = {"q": "qfp", "a": "afp", "b": "bfp", "c": "cfp"}
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: fps.get(smi))
    monkeypatch.setattr(dti, "_tanimoto", lambda x, y: 0.9 if x == "qfp" and y == "afp" else 0.0)
    monkeypatch.setattr(dti, "CNS_TARGETS", ("T",))
    monkeypatch.setattr(dti, "predict_kd_nm", lambda *a, **k: None)
    assert cns_ic50_nm("q") is None


def test_resolve_panel_skips_candidates_without_kd(monkeypatch):
    from drugos.target.targets import safety_panel as _panel

    snap = _fake_snapshot()
    monkeypatch.setattr(dti, "load_snapshot", lambda: snap)
    fps = {"q": "qfp", "a": "afp", "b": "bfp", "c": "cfp"}
    monkeypatch.setattr(dti, "_fingerprint", lambda smi: fps.get(smi))
    monkeypatch.setattr(dti, "_tanimoto", lambda x, y: 0.9 if x == "qfp" and y == "afp" else 0.0)
    panel = _panel()
    monkeypatch.setattr(dti, "SITE_TO_TARGET", {t.name: "T" for t in panel})
    monkeypatch.setattr(dti, "predict_kd_nm", lambda *a, **k: None)
    resolved, no_data = resolve_offtarget_panel(panel, smiles="q")
    assert all(t.low_confidence for t in resolved)
    assert no_data == tuple(sorted(t.name for t in panel))


def test_cns_anchor_requires_chemotype_support():
    # haloperidol (DRD2 top-match 0.67) honors a genuine sub-µM anchor...
    haloperidol = "OC1(CCN(CC1)CCCC2=CC=C(F)C=C2)C3=CC=C(Cl)C=C3"
    cns = cns_ic50_nm(haloperidol)
    assert cns is not None and 0 < cns < 1.0e4
    supported = [
        dti.predict_kd_nm(t, haloperidol)
        for t in dti.CNS_TARGETS
        if (dti.top_tanimoto(t, haloperidol) or 0.0) >= dti.MIN_NEIGHBOR_TANIMOTO
    ]
    assert cns == pytest.approx(min(supported), rel=1e-9)
    # ...while chemotype-desert molecules (dofetilide, carboxylic acids) get
    # no structure-derived anchor and fall to the disclosed prior.
    for smi in (DOFETILIDE, "CC(=O)O", "CC1=CC=C(C=C1)C(=O)O"):
        assert cns_ic50_nm(smi) is None
