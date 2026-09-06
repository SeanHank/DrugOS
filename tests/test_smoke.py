"""Smoke tests for the DrugOS development environment."""

from drugos import __version__
from drugos.version import PACKAGE_VERSION


def test_project_version() -> None:
    assert __version__ == "2026.9.0"
    assert PACKAGE_VERSION == __version__


def test_version_scheme() -> None:
    year, month, rev = (int(p) for p in __version__.split("."))
    assert year == 2026
    assert 1 <= month <= 12
    assert rev >= 0
