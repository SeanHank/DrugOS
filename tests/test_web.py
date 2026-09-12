"""Playground API tests (drugos.web.app) — user-checkpoint #5 completeness."""

from __future__ import annotations

import runpy
from pathlib import Path

import flask
import pytest

from drugos.web import app as web_app
from drugos.web.app import _float, _int, create_app


@pytest.fixture()
def client() -> object:
    return create_app().test_client()


def test_index_page(client: object) -> None:
    r = client.get("/")
    assert r.status_code == 200
    body = r.get_data(as_text=True)
    for needle in (
        "Run pipeline",
        "MULTISCALE DRUG RESPONSE",
        "uncertainty (D21)",
        "population (D22)",
        "sensitivity (D23)",
        "c_brain",
        "CNS · brain free exposure",
        "Uncertainty ensemble (D21)",
        "Virtual population (D22)",
        "Sensitivity (D23)",
        "Target engagement",
        "Pathway signaling",
        "occ_rows",
        "pathway_html",
    ):
        assert needle in body


def test_benchmarks_api(client: object) -> None:
    r = client.get("/api/benchmarks")
    assert r.status_code == 200
    data = r.get_json()
    assert "dofetilide" in data["benchmarks"]


def test_run_default(client: object) -> None:
    r = client.post("/api/run", json={"benchmark": "dofetilide", "dose": 0.5})
    assert r.status_code == 200
    data = r.get_json()
    # Full fidelity auto-engages a TMDD sink at the hERG site, which lowers a
    # 0.5 mg dofetilide run from "High" to "Elevated composite risk" — the
    # cardiac endpoint is still the flagged one (re-baselined under the default
    # mandated realism, doc/06 §6).
    assert data["clinical"]["verdict"].startswith(
        ("High composite risk", "Elevated composite risk")
    )
    assert "cns" in data["organ"]
    assert data["organ"]["cardiac"]["tdpr_band"] != ""
    assert data["occupancy"]["primary"] > 0.0
    assert data["occupancy"]["ranked"][0]["time_at_target_h"] > 0.0
    assert data["pathway"]["available"] is True
    assert data["pathway"]["readout_peak"] > 0.0


def test_run_unknown_benchmark(client: object) -> None:
    r = client.post("/api/run", json={"benchmark": "nope"})
    assert r.status_code == 400


def test_run_reliability_and_measured_payload(client: object) -> None:
    # DISCLAIMER §2: default (benchmark) run reports the validated regime.
    r = client.post("/api/run", json={"benchmark": "acetaminophen"})
    assert r.status_code == 200
    reliability = r.get_json()["trust"]["reliability"]
    assert reliability["regime"] == "validated_in_range_on_label"
    assert reliability["band_cv"] == 0.2
    # a measured-true-parameters override is accepted through the same endpoint
    r2 = client.post(
        "/api/run",
        json={
            "benchmark": "acetaminophen",
            "measured": {"fup": 0.75, "cl_hep_l_h": 22.0, "cl_renal_l_h": 2.0},
            "empirical": {"plasma_cmax_mg_l": 12.0},
        },
    )
    assert r2.status_code == 200
    data2 = r2.get_json()
    assert data2["trust"]["reliability"]["regime"] == "measured_in_range_on_label"
    assert data2["trust"]["empirical_agreement"] is not None
    obs = data2["trust"]["empirical_agreement"]["observations"][0]
    assert obs["endpoint"] == "plasma_cmax_mg_l"


def test_run_robustness_flags(client: object) -> None:
    r = client.post(
        "/api/run",
        json={
            "benchmark": "dofetilide",
            "dose": 0.5,
            "with_uncertainty": True,
            "with_population": True,
            "with_sensitivity": True,
            "n_unc": 3,
            "n_pop": 2,
            "n_sobol": 2,
        },
    )
    assert r.status_code == 200
    data = r.get_json()
    assert data["uncertainty"]["n_runs"] == 3
    assert "band_90" in data["uncertainty"]
    assert "brain_free_nm" in data["uncertainty"]["band_90"]
    assert data["population"]["n_individuals"] == 2
    assert "incidence_grade_ge_1" in data["population"]
    assert data["sensitivity"]["sobol"]["design_n"] == 2


def test_run_smiles_with_robustness(client: object) -> None:
    r = client.post(
        "/api/run",
        json={
            "smiles": "CC(=O)Nc1ccc(O)cc1",
            "dose": 20000,
            "with_uncertainty": True,
            "n_unc": 2,
        },
    )
    assert r.status_code == 200
    data = r.get_json()
    assert data["clinical"]["verdict"].startswith("High composite risk")


def test_run_empty_payload_defaults(client: object) -> None:
    r = client.post("/api/run", data="not json", content_type="application/json")
    assert r.status_code == 200  # no payload -> acetaminophen default
    r = client.post("/api/run", json=None)
    assert r.status_code == 200
    assert "acetaminophen" in r.get_json()["manifest"]["name"]


def test_run_malformed_input_is_explicit_400(client: object) -> None:
    r = client.post("/api/run", json={"benchmark": "warfarin", "dose": "abc"})
    assert r.status_code == 400
    assert "error" in r.get_json()


def test_coercion_helpers() -> None:
    assert _float(None) is None
    assert _float("5.5") == 5.5
    assert _float(4) == 4.0
    with pytest.raises(ValueError):
        _float("abc")
    with pytest.raises(ValueError):
        _float([1])
    with pytest.raises(ValueError):
        _float(True)
    assert _int(None, 7) == 7
    assert _int("9", 7) == 9
    with pytest.raises(ValueError):
        _int("abc", 7)
    with pytest.raises(ValueError):
        _int("0", 7)
    with pytest.raises(ValueError):
        _int(True, 7)
    assert _int("3", 7) == 3


def test_run_no_pathway_flag(client: object) -> None:
    r = client.post("/api/run", json={"benchmark": "warfarin", "no_pathway": True})
    assert r.status_code == 200
    data = r.get_json()
    assert data["pathway"]["available"] is False


def test_module_main_guard_false_branch() -> None:
    # ``drugos.web.app`` is imported as a library here (not __main__), so the
    # ``if __name__ == "__main__"`` guard's False branch is taken at import.
    assert web_app.__name__ == "drugos.web.app"


def test_module_main_guard_true_branch(monkeypatch: pytest.MonkeyPatch) -> None:
    # Executing the module as ``__main__`` takes the True branch; the real
    # server is never started because Flask.run is replaced by a recorder.
    calls: list[dict[str, object]] = []

    def fake_run(self: object, **kwargs: object) -> None:
        calls.append(kwargs)

    monkeypatch.setattr(flask.Flask, "run", fake_run)
    runpy.run_path(str(Path(web_app.__file__)), run_name="__main__")
    assert calls  # the guard body ran: create_app().run() was invoked
