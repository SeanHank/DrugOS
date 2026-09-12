"""Coverage for the release tooling in ``scripts/release.py``.

The fallback audit (G5) and the marker audit (G6) are themselves gates: these
tests pin the AST/grepper detectors so a regression in the tooling cannot
silently re-allow a swallowed handler or a retained marker. The tool is driven
on the real package (G5's handler inventory reconciles; G6's marker audit is
clean) and on tiny synthetic trees so the detectors' judgment is asserted
directly.
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


def test_release_requires_explicit_version() -> None:
    # No auto-bump: the release command must be told exactly which version.
    with pytest.raises(SystemExit):
        release.main(["release", "--no-gates"])


def test_release_rejects_malformed_version() -> None:
    with pytest.raises(SystemExit):
        release.main(["release", "--version", "x.y.z", "--no-gates"])


def test_release_accepts_wellformed_version() -> None:
    # A well-formed --version passes the format guard; the gate path is what
    # the release command drives next (exercised by the integration gates).
    assert (
        release.cmd_release(
            type(
                "A",
                (),
                {
                    "version": "2026.9.1",
                    "dry_run": True,
                    "no_gates": True,
                    "no_sync": True,
                    "no_xdist": True,
                },
            )()
        )
        == 0
    )


def test_hard_marker_pattern_detects_markers() -> None:
    # Marker literals are assembled from fragments so this detector test does
    # not itself trip the very markers it is pinning.
    f_todo = "T" + "ODO"
    f_nie = "Not" + "Implemented" + "Error"
    f_unim = "unim" + "plemented"
    f_nimp = "not im" + "plemented"
    assert release.HARD_MARKER_PATTERN.search(f"# {f_todo} wire this") is not None
    assert release.HARD_MARKER_PATTERN.search(f"raise {f_nie}") is not None
    assert release.HARD_MARKER_PATTERN.search(f"{f_unim} lane") is not None
    assert release.HARD_MARKER_PATTERN.search(f_nimp) is not None
    assert release.HARD_MARKER_PATTERN.search("an honest measurement") is None


def test_marker_tokens_in_covers_full_set() -> None:
    # Every scanned file, code and docs alike, is checked against the full
    # marker set; the audit keeps no allowance file and no whitelist.
    f_def = "defe" + "rred"
    assert release.marker_tokens_in(f_def + " until the next phase") == [f_def]
    f_todo = "TO" + "DO"
    f_nie = "Not" + "Implemented" + "Error"
    assert release.marker_tokens_in(f"# {f_todo} later") == [f_todo]
    assert release.marker_tokens_in(f_nie) == [f_nie]
    assert release.marker_tokens_in("plain prose") == []


def test_no_deferral_audit_real_scope_clean() -> None:
    violations = release.no_deferral_audit()
    assert violations == [], "\n".join(violations)


def test_marker_tokens_in_catches_inflections() -> None:
    # Marker literals are assembled from fragments so this detector test does
    # not itself trip the very markers it is pinning.
    f_stub = "s" + "tub"
    f_simpl = "simplif" + "ication"
    f_deferral = "defer" + "ral"
    f_catalog = "cata" + "log"
    hits = release.marker_tokens_in(f"those {f_stub}s, their {f_simpl}s, {f_deferral}s")
    assert f_stub in hits and f_simpl in hits and f_deferral in hits
    assert release.marker_tokens_in(f"a {f_catalog}ed americanized entry") == [f_catalog]
    assert release.marker_tokens_in("this work is real") == []


def test_marker_tokens_in_recognizes_g2_typing_path() -> None:
    # ``stubs/`` is the committed G2 typing-declarations directory (a realized
    # artifact); its literal path must not be misread as the marker word.
    f_stub = "s" + "tub"
    f_ph = "place" + "holder"
    for line in (
        f"typed by ``{f_stub}s/libsbml/__init__.pyi`` (doc/06, G2)",
        f"the typing packages under `{f_stub}s/` cover the untyped wheel",
        f"consumed subset typed by ``{f_stub}s/rpy2/robjects.pyi``",
    ):
        assert release.marker_tokens_in(line) == [], line
    assert release.marker_tokens_in(f"a {f_stub}bed-out {f_ph} remains") != []


def test_unrealized_markers_in_detects_status_vocab() -> None:
    # Status-word literals are assembled from fragments so this G7 detector
    # test carries none of the wording it is pinning open in the test file.
    f_planned = "p" + "lanned"
    f_blocked = "bl" + "ocked"
    f_partial = "parti" + "al"
    f_pending = "pen" + "ding"
    f_candidate = "can" + "didate"
    f_future = "fu" + "ture"
    f_roadmap = "road" + "map"
    f_not_downloaded = "not " + "downloaded"
    f_out_of_scope = "out of " + "scope"
    f_phase8 = "Phase-" + "8"
    f_release_track = "Release-" + "Track"
    f_p5 = "P" + "5"
    assert f_planned in release.unrealized_markers_in("the lane is " + f_planned)
    assert f_blocked in release.unrealized_markers_in("the dependency is " + f_blocked)
    assert f_partial in release.unrealized_markers_in("a " + f_partial + "ly wired seam")
    assert f_pending in release.unrealized_markers_in(f_pending + " dataset")
    assert f_candidate in release.unrealized_markers_in("an upstream " + f_candidate + " vendor")
    assert f_future in release.unrealized_markers_in(f_future + " work")
    assert f_roadmap in release.unrealized_markers_in("the " + f_roadmap)
    assert "phrase 'not downloaded'" in release.unrealized_markers_in(
        f_not_downloaded + " (see row 9)"
    )
    assert "phrase 'out of scope'" in release.unrealized_markers_in(f_out_of_scope)
    assert "phrase 'phase 8'" in release.unrealized_markers_in(f_phase8 + " coupling")
    assert "phrase 'release track'" in release.unrealized_markers_in(f_release_track + " register")
    assert "phrase 'P5'" in release.unrealized_markers_in(f_p5 + "-P9 lanes")


def test_unrealized_markers_in_ignores_scientific_prose() -> None:
    # Scientific language that legitimately uses the substrings is untouched:
    # the gate pins the status sense, not the channel-block physiology.
    assert release.unrealized_markers_in("hERG channel block reserves IKr") == []
    assert release.unrealized_markers_in("beta-blockers and non-blockers") == []
    assert release.unrealized_markers_in("the fractional IKr conductance") == []
    assert release.unrealized_markers_in("shipped and wired in 2026.9.1") == []
    assert release.unrealized_markers_in("blockade of hERG") == []
    assert release.unrealized_markers_in("forwards to the next stage") == []
    assert release.unrealized_markers_in("the regulator publishes guidance") == []


def test_docs_truth_audit_exempts_citation_surface(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The citation/license surface is read-only attribution (the G6/G7 red
    # line): a cited work's real title may use any wording without becoming a
    # claim about this product, and the gate must not require it to be
    # rewritten.
    f_planned = "p" + "lanned"
    f_out_of_scope = "out of " + "scope"
    doc = tmp_path / "chip.md"
    doc.write_text(
        "# Chip\n"
        f"The lane is {f_planned}, an {f_out_of_scope} seam.\n"
        "Every row below is wired.\n"
        "\n"
        "## References\n"
        f"Braun, C. (2026). The {f_planned} readout. J. C:\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(release, "ROOT", tmp_path)
    monkeypatch.setattr(release, "_doc_scope_files", lambda: [("doc", doc)])
    monkeypatch.setattr(release, "verify_doc_claims", lambda: [])
    violations = release.docs_truth_audit()
    assert any(":2:" in v for v in violations), violations
    assert all(":6:" not in v for v in violations), violations


def test_docs_truth_audit_real_scope_clean() -> None:
    # G7 across the whole shipped docs set: every document describes only
    # realized, wired behaviour, and every back-ticked claim in doc/12
    # resolves to a real src/drugos symbol, registered validation case, or
    # manifest-pinned file.
    violations = release.docs_truth_audit()
    assert violations == [], "\n".join(violations)


def test_verify_doc_claims_real_register_clean() -> None:
    violations = release.verify_doc_claims()
    assert violations == [], f"{len(violations)} unbacked claims:\n" + "\n".join(violations)


def test_gate_log_has_header_and_passed_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    log = tmp_path / "quality_gate.log"
    monkeypatch.setattr(release, "GATE_LOG", log)
    monkeypatch.setattr(release, "ensure_versions_consistent", lambda: None)
    monkeypatch.setattr(release, "verify_data_checksums", lambda: None)
    monkeypatch.setattr(release, "fallback_audit", lambda: [])
    monkeypatch.setattr(release, "no_deferral_audit", lambda: [])
    monkeypatch.setattr(release, "docs_truth_audit", lambda: [])
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
    assert "marker audit: clean" in text
    assert "docs-truth audit: clean" in text
    assert text.endswith("ALL GATES PASSED\n")
    g3 = next(c for c in commands if len(c) >= 3 and c[2] == "pytest")
    assert g3 == [release.python(), "-m", "pytest"]


def test_run_gates_uses_xdist_when_enabled(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    log = tmp_path / "quality_gate.log"
    monkeypatch.setattr(release, "GATE_LOG", log)
    monkeypatch.setattr(release, "ensure_versions_consistent", lambda: None)
    monkeypatch.setattr(release, "verify_data_checksums", lambda: None)
    monkeypatch.setattr(release, "fallback_audit", lambda: [])
    monkeypatch.setattr(release, "no_deferral_audit", lambda: [])
    monkeypatch.setattr(release, "docs_truth_audit", lambda: [])
    commands: list[list[str]] = []
    monkeypatch.setattr(
        release,
        "_run_cmd",
        lambda log_path, label, cmd, **_kwargs: commands.append(list(cmd)) or 0,
    )
    release.run_gates(use_xdist=True)
    g3 = next(c for c in commands if len(c) >= 3 and c[2] == "pytest")
    assert g3 == [release.python(), "-m", "pytest", "-n", "auto"]
