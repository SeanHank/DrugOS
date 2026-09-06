"""Shared fast fixtures: tiny-horizon specs that keep the suite quick."""

from __future__ import annotations

from typing import Any

import pytest

from drugos.pipeline import RunSpec, _BenchmarkLike, spec_from_benchmark_data


def _bench(
    name: str,
    smiles: str,
    fup: float,
    bp: float,
    dose: float,
    log_p: float,
    acids: list[float],
    bases: list[float],
    published: dict[str, tuple[float, float]],
) -> _BenchmarkLike:
    return _BenchmarkLike(
        name=name,
        smiles=smiles,
        fup=fup,
        bp=bp,
        route="oral",
        dose_mg=dose,
        tmax_h=6.0,
        n_eval=61,
        log_p=log_p,
        pka_acids=acids,
        pka_bases=bases,
        published=published,
    )


@pytest.fixture(scope="module")
def fast_warfarin() -> RunSpec:
    bench = _bench(
        "warfarin",
        "CC(=O)c1ccc(cc1)C(C(=O)O)",
        0.01,
        1.15,
        5.0,
        3.1,
        [5.0],
        [],
        {"cl_plasma_l_h": (0.15, 0.2), "urine_fraction": (0.8, 0.9)},
    )
    return spec_from_benchmark_data(bench)


@pytest.fixture(scope="module")
def fast_apap() -> RunSpec:
    bench = _bench(
        "acetaminophen",
        "CC(=O)Nc1ccc(O)cc1",
        0.75,
        1.0,
        1000.0,
        0.5,
        [9.9],
        [],
        {"cl_plasma_l_h": (20.0, 25.0), "urine_fraction": (0.0, 0.05)},
    )
    return spec_from_benchmark_data(bench)


@pytest.fixture(scope="module")
def fast_apap_od() -> RunSpec:
    bench = _bench(
        "acetaminophen",
        "CC(=O)Nc1ccc(O)cc1",
        0.75,
        1.0,
        20000.0,
        0.5,
        [9.9],
        [],
        {"cl_plasma_l_h": (20.0, 25.0), "urine_fraction": (0.0, 0.05)},
    )
    return spec_from_benchmark_data(bench)


@pytest.fixture()
def no_safety_panel(fast_warfarin: RunSpec) -> RunSpec:
    from dataclasses import replace

    return replace(fast_warfarin, qt_ic50_nm=None, panel=())


@pytest.fixture()
def non_herg_panel(fast_warfarin: RunSpec) -> RunSpec:
    from dataclasses import replace

    from drugos.target.targets import Target

    return replace(
        fast_warfarin,
        qt_ic50_nm=None,
        panel=(
            Target(
                name="sodium channel",
                kd_nm=5.0e4,
                r0_nm=1.0,
                reference="synthetic fixture",
            ),
        ),
    )


@pytest.fixture()
def default_panel_no_override(fast_warfarin: RunSpec) -> RunSpec:
    from dataclasses import replace

    return replace(fast_warfarin, qt_ic50_nm=None)


@pytest.fixture()
def cns_anchored_apap(fast_apap: RunSpec) -> RunSpec:
    from dataclasses import replace

    return replace(fast_apap, cns_ic50_nm=2.0e4, dili_ic50_nm=5.0e5, qt_ic50_nm=1.0e6)


@pytest.fixture()
def study_args() -> Any:
    from types import SimpleNamespace

    return SimpleNamespace(
        benchmark="warfarin",
        dose=None,
        route=None,
        sex="male",
        age=40.0,
        height=170.0,
        weight=70.0,
        no_pathway=False,
        n_unc=3,
        n_pop=3,
        n_sobol=2,
        seed=7,
        format="json",
        out=None,
    )
