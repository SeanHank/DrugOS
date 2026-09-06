"""Coverage for the release tooling in ``scripts/release.py``.

The fallback audit (G5) is itself a gate: these tests pin the AST detector so
that a regression in the tooling cannot silently re-allow a swallow handler.
The tool is driven on the real package (must match the pinned allowlist) and on
tiny synthetic trees so the detector's judgment is asserted directly.
"""

import ast
import importlib.util
from collections.abc import Sequence
from pathlib import Path

import pytest

_RELEASE = Path("scripts/release.py").resolve()
_SPEC = importlib.util.spec_from_file_location("release_tool", _RELEASE)
assert _SPEC is not None and _SPEC.loader is not None
release = importlib.util.module_from_spec(_SPEC)
assert isinstance(_SPEC.loader, importlib.abc.Loader)
_SPEC.loader.exec_module(release)


def _handler(body_source: str) -> ast.ExceptHandler:
    """Build an ExceptHandler AST node with the given body statements."""
    tree = ast.parse(
        "try:\n    1 / 0\nexcept Exception:\n"
        + "\n".join(f"    {line}" for line in body_source.splitlines())
    )
    return next(n for n in ast.walk(tree) if isinstance(n, ast.ExceptHandler))


@pytest.mark.parametrize(
    ("body", "explicit"),
    [
        ("raise ValueError('boom')", True),
        ("return jsonify({'error': 'boom'}), 400", True),
        ("return 2", True),
        ("pass", False),
        ("return None", False),
        ("result = 5", False),
        ("if bad:\n        raise ValueError('boom')", True),
        ("if not bad:\n        return", False),
    ],
)
def test_surfaces_explicit_error(body: str, explicit: bool) -> None:
    node = _handler(body)
    assert release._surfaces_explicit_error(node) is explicit


def test_except_handlers_real_package_matches_allowlist() -> None:
    violations = release.fallback_audit()
    assert violations == [], "\n".join(violations)


def test_except_handlers_scans_tmp_package(tmp_path: Path) -> None:
    pkg = tmp_path / "sample_pkg"
    pkg.mkdir()
    (pkg / "mod.py").write_text(
        "def f():\n"
        "    try:\n"
        "        return 1\n"
        "    except Exception:\n"
        "        return 2\n"
        "def g():\n"
        "    try:\n"
        "        x = 1\n"
        "    except Exception:\n"
        "        raise RuntimeError('boom')\n",
        encoding="utf-8",
    )
    found = release._except_handlers(pkg)
    assert found == [("mod.py", 4), ("mod.py", 9)]


def test_versions_consistent(monkeypatch: pytest.MonkeyPatch) -> None:
    assert release.read_version() == release.pyproject_version()
    release.ensure_versions_consistent()


def test_versions_mismatch_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(release, "pyproject_version", lambda: "9999.1.1")
    with pytest.raises(SystemExit):
        release.ensure_versions_consistent()


def test_validation_counts_parses_report() -> None:
    passed, total = release.validation_counts()
    assert total >= 1
    assert 0 <= passed <= total


def test_next_version_modes() -> None:
    today = type("D", (), {"year": 2026, "month": 9})()
    assert release.next_version("2026.9.0", "none", today) == "2026.9.0"
    assert release.next_version("2026.9.0", "auto", today) == "2026.9.1"
    assert release.next_version("2026.8.0", "auto", today) == "2026.9.0"
    assert release.next_version("2026.8.3", "revision", today) == "2026.8.4"
    assert release.next_version("2026.8.3", "month", today) == "2026.9.0"
    assert release.next_version("2026.8.3", "year", today) == "2026.0.0"


def test_gate_log_has_header_and_passed_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    log = tmp_path / "quality_gate.log"
    monkeypatch.setattr(release, "GATE_LOG", log)
    monkeypatch.setattr(release, "ensure_versions_consistent", lambda: None)
    monkeypatch.setattr(release, "verify_data_checksums", lambda: None)
    monkeypatch.setattr(release, "fallback_audit", lambda: [])
    commands: list[list[str]] = []

    def fake_run_cmd(log_path: Path, label: str, cmd: Sequence[str], **_kwargs: object) -> int:
        commands.append(list(cmd))
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"== {label} ==\n")
        return 0

    monkeypatch.setattr(release, "_run_cmd", fake_run_cmd)
    release.run_gates(use_xdist=False)
    text = log.read_text(encoding="utf-8")
    assert text.startswith("DrugOS release quality gate")
    assert "started " in text
    assert "== G1 lint (ruff check) ==" in text
    assert "G4 validation suite" in text
    assert "fallback audit: clean" in text
    assert text.endswith("ALL GATES PASSED\n")
    g3 = next(c for c in commands if len(c) >= 3 and c[2] == "pytest")
    assert g3 == [release.python(), "-m", "pytest"]


def test_run_gates_uses_xdist_when_enabled(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    log = tmp_path / "quality_gate.log"
    monkeypatch.setattr(release, "GATE_LOG", log)
    monkeypatch.setattr(release, "ensure_versions_consistent", lambda: None)
    monkeypatch.setattr(release, "verify_data_checksums", lambda: None)
    monkeypatch.setattr(release, "fallback_audit", lambda: [])
    commands: list[list[str]] = []
    monkeypatch.setattr(
        release,
        "_run_cmd",
        lambda log_path, label, cmd, **_kwargs: commands.append(list(cmd)) or 0,
    )
    release.run_gates(use_xdist=True)
    g3 = next(c for c in commands if len(c) >= 3 and c[2] == "pytest")
    assert g3 == [release.python(), "-m", "pytest", "-n", "auto"]
