"""CLI tests (drugos.cli, python -m drugos)."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pytest

import drugos.cli as cli
from drugos.cli import (
    _smiles_spec,
    main,
    spec_from_cli,
    write_report_directory,
)
from drugos.pipeline import RunSpec, run_pipeline


def test_version_flag(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "drugos" in capsys.readouterr().out


def test_main_dispatch_run(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = Mock(return_value=3)
    monkeypatch.setattr(cli, "cmd_run", fake)
    assert main(["run", "--benchmark", "warfarin"]) == 3
    assert fake.call_args.args[0].benchmark == "warfarin"


def test_main_dispatch_benchmarks(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = Mock(return_value=7)
    monkeypatch.setattr(cli, "cmd_benchmarks", fake)
    assert main(["benchmarks"]) == 7


def test_main_dispatch_serve(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = Mock(return_value=9)
    monkeypatch.setattr(cli, "cmd_serve", fake)
    assert main(["serve", "--port", "9090"]) == 9
    assert fake.call_args.args[0].port == 9090


def test_main_dispatch_study(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = Mock(return_value=5)
    monkeypatch.setattr(cli, "cmd_study", fake)
    assert main(["study", "--n-unc", "2"]) == 5


def test_main_dispatch_version(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = Mock(return_value=4)
    monkeypatch.setattr(cli, "cmd_version", fake)
    assert main(["version"]) == 4


def test_cmd_version(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.cmd_version() == 0
    assert capsys.readouterr().out == f"drugos {cli.__version__}\n"


def test_main_unknown_dispatch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli, "_dispatch", lambda command: None)
    assert main(["run", "--benchmark", "warfarin"]) == 0


def test_dispatch_lookup() -> None:
    assert cli._dispatch("run") is cli.cmd_run
    assert cli._dispatch("benchmarks") is cli.cmd_benchmarks
    assert cli._dispatch("serve") is cli.cmd_serve
    assert cli._dispatch("study") is cli.cmd_study
    assert cli._dispatch("version") is cli.cmd_version
    assert cli._dispatch("nope") is None


def test_main_module_import_as_library() -> None:
    """Normal import must not run the entry guard (False branch)."""
    import importlib

    importlib.import_module("drugos.__main__")


def test_main_module_entry(monkeypatch: pytest.MonkeyPatch) -> None:
    import runpy
    import sys

    monkeypatch.setattr(sys, "argv", ["drugos", "benchmarks"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module("drugos", run_name="__main__")
    assert exc.value.code == 0


def test_cli_main_entry_guard(monkeypatch: pytest.MonkeyPatch) -> None:
    import runpy
    import sys

    monkeypatch.setattr(sys, "argv", ["drugos", "benchmarks"])
    entry = Path(__file__).resolve().parents[1] / "src" / "drugos" / "cli.py"
    with pytest.raises(SystemExit) as exc:
        runpy.run_path(str(entry), run_name="__main__")
    assert exc.value.code == 0


def test_benchmark_lookup() -> None:
    assert cli._benchmark_lookup("warfarin") is not None
    assert cli._benchmark_lookup("not-a-compound") is None


def test_spec_from_cli_benchmark(fast_warfarin: RunSpec) -> None:
    spec = spec_from_cli(benchmark="warfarin")
    assert spec.name == "warfarin"
    assert spec.profile.weight_kg == 70.0
    spec = spec_from_cli(
        benchmark="warfarin",
        dose=25.0,
        route="iv_bolus",
        sex="female",
        age=55.0,
        height=160.0,
        weight=60.0,
        no_pathway=True,
    )
    assert spec.dose_plan.events[0].dose_mg == 25.0
    assert spec.dose_plan.events[0].route.value == "iv_bolus"
    assert spec.profile.sex.value == "female"
    assert spec.include_pathway is False
    with pytest.raises(ValueError):
        spec_from_cli(benchmark="does-not-exist")


def test_spec_from_cli_default_benchmark() -> None:
    assert spec_from_cli().name == "acetaminophen"


def _fake_admet(**kw: Any) -> SimpleNamespace:
    base = dict(fup_plasma=0.4, cl_int_hep_ml_min_kg=8.0, DILI=0.1, hERG=0.1, BBB=0.1)
    base.update(kw)
    return SimpleNamespace(**base)


def _fake_mol(mw: float = 150.0, name: str = "custom") -> Any:
    return SimpleNamespace(mw=mw, name=name)


def test_smiles_spec_success(monkeypatch: pytest.MonkeyPatch) -> None:
    import importlib

    ps = importlib.import_module("drugos.inputs.parse_structure")
    admet_mod = importlib.import_module("drugos.pk.admet")
    monkeypatch.setattr(ps, "parse_structure", lambda s, name: _fake_mol())
    monkeypatch.setattr(admet_mod, "predict_admet", lambda s: _fake_admet())
    spec = _smiles_spec(
        SimpleNamespace(
            smiles="CC",
            dose=20.0,
            route="oral",
            sex="male",
            age=40.0,
            height=170.0,
            weight=70.0,
            no_pathway=False,
        )
    )
    assert spec.name == "custom"
    assert spec.mw == 150.0
    assert spec.fup == 0.4
    assert spec.dose_plan.events[0].dose_mg == 20.0
    assert spec.admet is not None


def test_smiles_spec_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    import importlib

    ps = importlib.import_module("drugos.inputs.parse_structure")
    admet_mod = importlib.import_module("drugos.pk.admet")
    ns = SimpleNamespace(
        smiles="CC",
        dose=10.0,
        route=None,
        sex="male",
        age=40.0,
        height=170.0,
        weight=70.0,
        no_pathway=False,
    )
    monkeypatch.setattr(ps, "parse_structure", lambda s, name: _fake_mol(0.0))
    monkeypatch.setattr(admet_mod, "predict_admet", lambda s: _fake_admet())
    with pytest.raises(ValueError):
        _smiles_spec(ns)
    monkeypatch.setattr(ps, "parse_structure", lambda s, name: _fake_mol())
    monkeypatch.setattr(admet_mod, "predict_admet", lambda s: _fake_admet(fup_plasma=None))
    with pytest.raises(ValueError):
        _smiles_spec(ns)
    monkeypatch.setattr(admet_mod, "predict_admet", lambda s: _fake_admet(fup_plasma=0.0))
    with pytest.raises(ValueError):
        _smiles_spec(ns)
    monkeypatch.setattr(
        admet_mod, "predict_admet", lambda s: _fake_admet(cl_int_hep_ml_min_kg=None)
    )
    with pytest.raises(ValueError):
        _smiles_spec(ns)
    ns_no_pathway = SimpleNamespace(**{**ns.__dict__, "no_pathway": True})
    monkeypatch.setattr(admet_mod, "predict_admet", lambda s: _fake_admet())
    spec = _smiles_spec(ns_no_pathway)
    assert spec.include_pathway is False


def test_cmd_run_json(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, fast_warfarin: RunSpec
) -> None:
    monkeypatch.setattr(cli, "spec_from_cli", lambda **kw: fast_warfarin)
    rc = cli.cmd_run(
        SimpleNamespace(
            benchmark="warfarin",
            smiles=None,
            dose=None,
            route=None,
            sex="male",
            age=40.0,
            height=170.0,
            weight=70.0,
            no_pathway=False,
            sc_im_ka=None,
            format="json",
            out=None,
        )
    )
    assert rc == 0
    data = json.loads(capsys.readouterr().out)
    assert data["manifest"]["name"] == "warfarin"


def test_cmd_run_error(capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(**kw: Any) -> RunSpec:
        raise ValueError("bad compound")

    monkeypatch.setattr(cli, "spec_from_cli", boom)
    ns = SimpleNamespace(
        benchmark="x",
        smiles=None,
        dose=None,
        route=None,
        sex="male",
        age=40.0,
        height=170.0,
        weight=70.0,
        no_pathway=False,
        sc_im_ka=None,
        format="json",
        out=None,
    )
    assert cli.cmd_run(ns) == 2
    assert "error: bad compound" in capsys.readouterr().err


def test_cmd_run_markdown(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    fast_warfarin: RunSpec,
) -> None:
    monkeypatch.setattr(cli, "spec_from_cli", lambda **kw: fast_warfarin)
    out_dir = str(tmp_path / "reports")
    ns = SimpleNamespace(
        benchmark="warfarin",
        smiles=None,
        dose=None,
        route=None,
        sex="male",
        age=40.0,
        height=170.0,
        weight=70.0,
        no_pathway=False,
        sc_im_ka=None,
        format="markdown",
        out=out_dir,
    )
    assert cli.cmd_run(ns) == 0
    assert (tmp_path / "reports" / "report.md").exists()
    assert "wrote json report" in capsys.readouterr().out


def test_spec_from_cli_sc_im_ka() -> None:
    plain = cli.spec_from_cli(benchmark="warfarin")
    assert plain.sc_im_ka_per_h is None
    sc = cli.spec_from_cli(benchmark="warfarin", sc_im_ka=0.1)
    assert sc.sc_im_ka_per_h == 0.1


def test_cmd_benchmarks(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.cmd_benchmarks() == 0
    names = capsys.readouterr().out.split()
    assert "warfarin" in names and "dofetilide" in names


def test_cmd_serve(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    app = SimpleNamespace(run=Mock())
    fake_module = SimpleNamespace(create_app=Mock(return_value=app))
    monkeypatch.setattr(cli.importlib, "import_module", lambda name: fake_module)
    ns = SimpleNamespace(host="0.0.0.0", port=1234)
    assert cli.cmd_serve(ns) == 0
    assert "http://0.0.0.0:1234" in capsys.readouterr().out
    assert fake_module.create_app.call_args.kwargs == {"host": "0.0.0.0", "port": 1234}
    app.run.assert_called_once()
    assert app.run.call_args.kwargs["port"] == 1234


def test_write_report_directory(tmp_path: Path, fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    write_report_directory(str(tmp_path / "deep" / "dir"), result)
    assert (tmp_path / "deep" / "dir" / "report.json").exists()


def _study_args(**kw: Any) -> SimpleNamespace:
    base = dict(
        benchmark="warfarin",
        dose=None,
        route=None,
        sex="male",
        age=40.0,
        height=170.0,
        weight=70.0,
        no_pathway=False,
        n_unc=1,
        n_pop=1,
        n_sobol=1,
        seed=7,
        format="json",
        out=None,
    )
    base.update(kw)
    return SimpleNamespace(**base)


@pytest.fixture(scope="module")
def study_artifacts(fast_warfarin: RunSpec) -> tuple[Any, Any, Any, Any]:
    return cli._study_objects(fast_warfarin, _study_args())


def _patch_study(monkeypatch: pytest.MonkeyPatch, artifacts: tuple[Any, Any, Any, Any]) -> None:
    monkeypatch.setattr(cli, "_study_objects", lambda spec, args: artifacts)


def test_study_contract(
    fast_warfarin: RunSpec,
    monkeypatch: pytest.MonkeyPatch,
    study_artifacts: tuple[Any, Any, Any, Any],
) -> None:
    _patch_study(monkeypatch, study_artifacts)
    study = cli._study_contract(fast_warfarin, _study_args())
    assert set(study) == {"compound", "dose_mg", "uncertainty", "population", "sensitivity"}
    assert study["compound"] == "warfarin"
    assert study["uncertainty"]["n_runs"] == 1
    assert study["population"]["n_individuals"] == 1
    assert study["sensitivity"]["sobol"]["design_n"] == 1
    assert len(study["sensitivity"]["drivers_overall_risk"]) == 6


def test_study_markdown(
    fast_warfarin: RunSpec,
    monkeypatch: pytest.MonkeyPatch,
    study_artifacts: tuple[Any, Any, Any, Any],
) -> None:
    _patch_study(monkeypatch, study_artifacts)
    md = cli._study_markdown(fast_warfarin, _study_args())
    assert md.startswith("# DrugOS robustness study — warfarin")
    assert "## D21 Uncertainty" in md
    assert "## D22 Virtual population" in md
    assert "## D23 Sensitivity" in md
    assert "Top local drivers" in md


def test_cmd_study_stdout_json(
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    fast_warfarin: RunSpec,
    study_artifacts: tuple[Any, Any, Any, Any],
) -> None:
    _patch_study(monkeypatch, study_artifacts)
    monkeypatch.setattr(cli, "spec_from_cli", lambda **kw: fast_warfarin)
    assert cli.cmd_study(_study_args()) == 0
    data = json.loads(capsys.readouterr().out)
    assert "uncertainty" in data


def test_cmd_study_markdown_stdout(
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    fast_warfarin: RunSpec,
    study_artifacts: tuple[Any, Any, Any, Any],
) -> None:
    _patch_study(monkeypatch, study_artifacts)
    monkeypatch.setattr(cli, "spec_from_cli", lambda **kw: fast_warfarin)
    assert cli.cmd_study(_study_args(format="markdown")) == 0
    out = capsys.readouterr().out
    assert out.startswith("# DrugOS robustness study")


def test_cmd_study_out_dir(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    fast_warfarin: RunSpec,
    study_artifacts: tuple[Any, Any, Any, Any],
) -> None:
    _patch_study(monkeypatch, study_artifacts)
    monkeypatch.setattr(cli, "spec_from_cli", lambda **kw: fast_warfarin)
    out = str(tmp_path / "study")
    assert cli.cmd_study(_study_args(out=out)) == 0
    assert (tmp_path / "study" / "study.json").exists()
    assert (tmp_path / "study" / "study.md").exists()
    assert "wrote json study" in capsys.readouterr().out


def test_cmd_study_error(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    def boom(**kw: Any) -> RunSpec:
        raise ValueError("unknown benchmark")

    monkeypatch.setattr(cli, "spec_from_cli", boom)
    assert cli.cmd_study(_study_args()) == 2
    assert "error: unknown benchmark" in capsys.readouterr().err
