"""Smoke tests for the DrugOS development environment."""

import tomllib
from pathlib import Path

from drugos import __version__
from drugos.version import PACKAGE_VERSION


def test_project_version_matches_pyproject() -> None:
    project_root = Path(__file__).resolve().parents[1]
    with (project_root / "pyproject.toml").open("rb") as fh:
        pyproject_version = tomllib.load(fh)["project"]["version"]
    assert __version__ == pyproject_version
    assert PACKAGE_VERSION == __version__


def test_version_scheme() -> None:
    year, month, rev = (int(p) for p in __version__.split("."))
    assert year == 2026
    assert 1 <= month <= 12
    assert rev >= 0
