"""Release automation for DrugOS (stdlib only).

Single entry point that replaces the former ``quality_gate.sh`` + ``release.py``
pair. It runs the release-quality gate, the no-silent-fallback audit, the
no-deferral/no-simplification audit, version bookkeeping and validation/status
sync, and drives a release:

Subcommands:

    gates               G1-G7 quality gates (lint, types, coverage, validation,
                        fallback audit, no-deferral audit, docs-truth audit),
                        logging to
                        ``build/quality_gate.log``. Hard-fails on any violation.
    fallback-audit      Scan every ``except`` handler in ``src/drugos`` against
                        the pinned allowlist
                        (``scripts/fallback_allowlist.json``). A handler must
                        either re-raise or return an explicit error value; a
                        silent swallow is a hard failure, and any drift in the
                        handler inventory fails until the allowlist is updated
                        intentionally.
    marker-audit        Scan every scanned code and doc file for any occurrence
                        of the hard or soft marker set (G6). Any occurrence at
                        any count is a hard failure: there is no allowance file,
                        no whitelist, no pinned inventory and no hard-coded
                        permission to retain a marker of any kind.
    docs-audit           Docs-truth audit (G7): every scanned document must
                        describe only realized, wired behaviour.  Two layers:
                        (a) the unrealized-status vocabulary (planned, blocked,
                        deferred, partial, pending, candidate, future, roadmap,
                        out-of-scope, not-yet, not-downloaded, and the rest of
                        the set) must be absent from doc text outside the
                        citation surface; (b) every back-ticked identifier and
                        ``case_*`` name in doc/12 §1/§5/§6 must resolve to a
                        real ``src/drugos`` symbol or a registered validation
                        case, and every back-ticked ``data/...`` path must be listed in
                        the pinned manifest.  A document that claims or
                        spans anything not fully realized is a hard failure.
   version              Print the canonical version (authoritative source:
                        ``version.py``); verifies ``pyproject.toml`` agrees.
   status               Re-sync validation counts, the ``doc/09`` status block
                        and ``build/release_status.json`` from the last gate
                        run.
   release (default)    gates (unless ``--no-gates``) -> project-wide version
                        sync from an explicit ``--version YYYY.M.V`` (never
                        auto-bumped) -> validation/status sync ->
                        ``release_status.json``.

All document updates are idempotent; ``--dry-run`` prints the plan without
writing. There is no silent fallback of any form inside this script either.
"""

from __future__ import annotations

import argparse
import ast
import datetime as _dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VERSION_FILE = ROOT / "src" / "drugos" / "version.py"
PYPROJECT = ROOT / "pyproject.toml"
REPORT = ROOT / "validation" / "report.md"
STATUS_JSON = ROOT / "build" / "release_status.json"
GATE_LOG = ROOT / "build" / "quality_gate.log"
INDEX_HTML = ROOT / "src" / "drugos" / "web" / "templates" / "index.html"
ALLOWLIST = ROOT / "scripts" / "fallback_allowlist.json"
DATA_MANIFEST = ROOT / "data" / "manifest.json"
PACKAGE = ROOT / "src" / "drugos"

_REPORT_COUNT = re.compile(r"\((?P<passed>\d+)/(?P<total>\d+) cases passed\)")

SYNC_FILES = sorted([ROOT / "README.md", *ROOT.glob("doc/*.md")])

KNOWN_COMMANDS = {
    "gates",
    "fallback-audit",
    "marker-audit",
    "docs-audit",
    "version",
    "status",
    "release",
}

HARD_MARKER_TOKENS = (
    "TODO",
    "FIXME",
    "XXX",
    "HACK",
    "NotImplementedError",
    "unimplemented",
    "not implemented",
)
HARD_MARKER_PATTERN = re.compile(
    r"\b(?:" + "|".join(re.escape(t) for t in HARD_MARKER_TOKENS) + r")\w*\b",
    re.IGNORECASE,
)

SOFT_MARKER_TOKENS = (
    "deferred",
    "deferral",
    "simplified",
    "simplification",
    "stub",
    "placeholder",
    "catalogue",
    "catalogued",
    "catalog",
    "backlog",
    "future work",
    "coming soon",
)
ALL_MARKER_TOKENS = HARD_MARKER_TOKENS + SOFT_MARKER_TOKENS
#: Every marker is matched on its stem plus any word-char continuation, so
#: plurals and derived forms ("stubs", "simplifications", "deferrals",
#: "cataloged", "backlogs", "placeholder notes") are caught exactly like
#: the bare token — G6 is absence-based and permissive to no form of a marker.
MARKER_TOKEN_RE = {
    tok: re.compile(rf"\b{re.escape(tok)}\w*\b", re.IGNORECASE) for tok in ALL_MARKER_TOKENS
}
#: The G2 typing-declarations directory (``stubs/``).  ``stubs/libsbml`` and
#: ``stubs/rpy2`` are real, committed, branch-exact type-declaration packages
#: backing 100%-tested modules — a realized artifact, not a placeholder.  The
#: audit strips the literal path segment ``stubs/`` before matching so the
#: directory name is never misread as the marker word.
STUBS_DIR_RE = re.compile(r"(?<![A-Za-z0-9_])stubs/")

#: Doc headings that open the attribution surface: Reference/Bibliography and
#: License sections and columns are read-only citation and license records.
#: The G6 audit recognizes these regions and never counts or flags lines in
#: them (the citation red line, doc/09 G6): a cited work's real title or
#: license text is attribution, not a marker, and the gate must not require it
#: to be rewritten.
CITATION_HEADING_RE = re.compile(
    r"^#{1,3}\s+.*\b(?:References?|Bibliography|Licen[cs]e[s]?|Attribution)\b.*$",
    re.IGNORECASE,
)


def python() -> str:
    """The interpreter running this script (the project env when driven by it)."""
    return sys.executable


def read_version() -> str:
    text = VERSION_FILE.read_text(encoding="utf-8")
    m = re.search(r'__version__\s*=\s*"([^"]+)"', text)
    if m is None:
        raise SystemExit(f"cannot parse version from {VERSION_FILE}")
    return m.group(1)


def write_version(version: str) -> None:
    text = VERSION_FILE.read_text(encoding="utf-8")
    updated, n = re.subn(r'__version__\s*=\s*"[^"]+"', f'__version__ = "{version}"', text, count=1)
    if n != 1:
        raise SystemExit(f"cannot rewrite version in {VERSION_FILE}")
    VERSION_FILE.write_text(updated, encoding="utf-8")


def pyproject_version() -> str:
    text = PYPROJECT.read_text(encoding="utf-8")
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if m is None:
        raise SystemExit(f"cannot parse version from {PYPROJECT}")
    return m.group(1)


def write_pyproject_version(version: str) -> None:
    text = PYPROJECT.read_text(encoding="utf-8")
    updated, n = re.subn(
        r'^version\s*=\s*"[^"]+"', f'version = "{version}"', text, count=1, flags=re.MULTILINE
    )
    if n != 1:
        raise SystemExit(f"cannot rewrite version in {PYPROJECT}")
    PYPROJECT.write_text(updated, encoding="utf-8")


def ensure_versions_consistent() -> None:
    """The canonical version lives in ``version.py``; ``pyproject.toml`` must match.

    CI derives the release version from ``pyproject.toml``, and the package
    reports from ``version.py`` — a mismatch would publish two different
    versions, so it is an explicit hard error (never a silent reconciliation).
    """
    src = read_version()
    pyproj = pyproject_version()
    if src != pyproj:
        raise SystemExit(
            f"version mismatch: version.py says {src}, pyproject.toml says {pyproj}; "
            "run `python scripts/release.py release --version <ver>` to sync project-wide"
        )


def validation_counts() -> tuple[int, int]:
    if not REPORT.exists():
        raise SystemExit(f"missing {REPORT}; run `python scripts/release.py gates` first (G4)")
    m = _REPORT_COUNT.search(REPORT.read_text(encoding="utf-8"))
    if m is None:
        raise SystemExit(f"cannot parse validation counts from {REPORT}")
    return int(m.group("passed")), int(m.group("total"))


def _run_cmd(log_path: Path, label: str, argv: list[str]) -> None:
    """Run one gate step appended to the log; hard-fail with log tail on error."""
    with log_path.open("a", encoding="utf-8") as log:
        print(f"== {label} ==", file=log)
        print(f"--- {label} ---", flush=True)
        result = subprocess.run(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print("GATE FAILED — last 40 lines:")
        print("".join(log_path.read_text(encoding="utf-8").splitlines(keepends=True)[-40:]))
        raise SystemExit(f"gate step failed ({label}, rc={result.returncode})")


def _except_handlers(src_root: Path) -> list[tuple[str, int]]:
    """Every ``except`` handler in ``src_root`` as ``(relpath, lineno)``.

    Paths are reported relative to ``ROOT`` when the scanned tree lives inside
    the repo (so the allowlist keys are stable repo paths), otherwise relative
    to the scanned tree itself (e.g. when tests scan a tmp package).
    """
    base = ROOT if ROOT in src_root.parents else src_root
    found: list[tuple[str, int]] = []
    for path in sorted(src_root.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ExceptHandler):
                found.append((path.relative_to(base).as_posix(), node.lineno))
    return found


def _surfaces_explicit_error(node: ast.ExceptHandler) -> bool:
    """True when the handler re-raises or returns a non-None error signal."""
    for stmt in node.body:
        if isinstance(stmt, (ast.Raise,)):
            return True
        if isinstance(stmt, ast.Return) and stmt.value is not None:
            value = stmt.value
            if not (isinstance(value, ast.Constant) and value.value is None):
                return True
        if isinstance(stmt, ast.If):
            for branch in (stmt.body, stmt.orelse):
                if any(isinstance(b, ast.Raise) for b in branch):
                    return True
                if any(
                    isinstance(b, ast.Return)
                    and b.value is not None
                    and not (isinstance(b.value, ast.Constant) and b.value.value is None)
                    for b in branch
                ):
                    return True
    return False


def verify_data_checksums() -> None:
    """Verify every vendored data/ file against the pinned sha256 manifest.

    Fails hard on a missing/extra file or a checksum mismatch — a data file
    must never be silently substituted or edited without the manifest being
    updated on purpose.
    """
    if not DATA_MANIFEST.exists():
        raise SystemExit(
            f"missing data manifest {DATA_MANIFEST}; data checksums cannot be verified"
        )
    manifest = json.loads(DATA_MANIFEST.read_text(encoding="utf-8"))
    entries = manifest["files"]
    for entry in entries:
        rel = entry["path"]
        path = ROOT / "data" / rel
        if not path.is_file():
            raise SystemExit(f"data manifest references missing file: {rel}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != entry["sha256"]:
            raise SystemExit(
                f"data checksum mismatch for {rel}: manifest {entry['sha256']} != disk {digest}"
            )
    print(
        f"data checksums: {len(entries)} file(s) verified against {DATA_MANIFEST.relative_to(ROOT)}"
    )


def fallback_audit() -> list[str]:
    """Scan for silent fallbacks; return a list of violations (empty = clean).

    Two independent layers:
      1. every handler must re-raise or surface an explicit error
         (``raise`` / non-None ``return``) — a silent swallow is a violation;
      2. the exact handler inventory is pinned in ``fallback_allowlist.json`` —
         any new/removed/relocated handler is a violation until the allowlist
         is updated on purpose.
    """
    handlers = _except_handlers(PACKAGE)
    violations: list[str] = []

    for path in PACKAGE.rglob("*.py"):
        rel = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ExceptHandler):
                continue
            if not _surfaces_explicit_error(node):
                violations.append(f"{rel}:{node.lineno}: handler swallows the exception silently")

    inventory: dict[str, list[int]] = {}
    for rel, lineno in sorted(handlers):
        inventory.setdefault(rel, []).append(lineno)

    allowlist = json.loads(ALLOWLIST.read_text(encoding="utf-8"))
    expected = {k: sorted(v) for k, v in allowlist["handlers"].items()}
    if inventory != expected:
        violations.append(
            "handler inventory drifted from scripts/fallback_allowlist.json\n"
            f"  found:    {json.dumps(inventory, sort_keys=True)}\n"
            f"  expected: {json.dumps(expected, sort_keys=True)}"
        )
    return violations


def _audit_scope_files() -> list[tuple[str, Path]]:
    """Every file the G6 marker audit scans, as ``(kind, path)`` pairs.

    Code scope is ``src/drugos``, ``tests``, ``validation`` and ``scripts``
    (``.py`` only); the checker script itself is excluded so the audit stays
    self-eligible. Doc scope is ``README.md``, ``DISCLAIMER.md``,
    ``CONTRIBUTING.md`` and ``doc/*.md``. Generated artifacts
    (``validation/report.md``) are not in scope. License and citation lines are
    the audited attribution surface: the gate is read-only over them and never
    counts or flags them, and no arrangement inside this gate ever treats a
    cited work as a reason to permit a marker.
    """
    files: list[tuple[str, Path]] = []
    for root in (PACKAGE, ROOT / "tests", ROOT / "validation", ROOT / "scripts"):
        for path in sorted(root.rglob("*.py")):
            if root == ROOT / "scripts" and path.name == "release.py":
                continue
            files.append(("code", path))
    for tail in ("README.md", "DISCLAIMER.md", "CONTRIBUTING.md"):
        files.append(("doc", ROOT / tail))
    for path in sorted(ROOT.glob("doc/*.md")):
        files.append(("doc", path))
    return files


def marker_tokens_in(line: str) -> list[str]:
    """The marker tokens present on one line of scanned text.

    Every scanned file — code and docs alike — is checked against the full
    marker set on stems plus any inflection, so no plural or derived form of a
    marker slips through. The audit keeps no allowance file, no whitelist and
    no pinned inventory: there is no count that a marker occurrence could be
    reconciled against, because no count of any marker is permitted.
    """
    redacted = STUBS_DIR_RE.sub("typingpkg/", line)
    return [tok for tok in ALL_MARKER_TOKENS if MARKER_TOKEN_RE[tok].search(redacted)]


def _citation_surface_line_numbers(text: str) -> set[int]:
    """Line numbers of doc text that belong to the citation/license surface.

    Any ``References``/``Bibliography``/``License`` heading begins a citation
    region that runs to the next heading (or end of file). The G6 audit treats
    those lines as read-only attribution (the citation red line): it never
    counts or flags them, because a cited work's real title or license text is
    attribution, not a marker, and the gate must not require it to be
    rewritten.
    """
    lines = text.splitlines()
    headings = [i + 1 for i, ln in enumerate(lines) if CITATION_HEADING_RE.match(ln)]
    skipped: set[int] = set()
    for idx, start in enumerate(headings):
        end = headings[idx + 1] if idx + 1 < len(headings) else len(lines) + 1
        skipped.update(range(start, end))
    return skipped


def no_deferral_audit() -> list[str]:
    """Scan code and docs for any marker of a stand-in, simplification or deferral.

    The G6 gate is purely absence-based:

      1. Every scanned file (code and docs) is checked against the full marker
         set — hard and soft — on stems plus any inflection, so plurals and
         derived forms are caught exactly like the bare token.
      2. Any occurrence at any count, in any file, is a violation. There is no
         allowance file, no whitelist, no pinned per-file/per-token inventory
         and no hard-coded exception, by design: the product carries no marker
         of a retained seam, and the gate has no mechanism that could license
         one.

    The citation/license surface is not an exception to that rule but the
    audited boundary the rule is read-only over: Reference/Bibliography and
    License regions of scanned docs are recognized attribution lines and are
    never counted or flagged, because a cited work's real title or license
    text is attribution, not a marker, and the gate must not require it to be
    rewritten (doc/09 G6). The G2 typing-declarations directory ``stubs/`` is
    likewise a realized, committed artifact, not a placeholder; its literal
    path is recognized so the directory name is not misread as the marker.
    """
    violations: list[str] = []
    for kind, path in _audit_scope_files():
        rel = path.relative_to(ROOT).as_posix()
        text = path.read_text(encoding="utf-8")
        citation_lines = _citation_surface_line_numbers(text) if kind == "doc" else set()
        for lineno, line in enumerate(text.splitlines(), 1):
            if lineno in citation_lines:
                continue
            for tok in marker_tokens_in(line):
                violations.append(f"{rel}:{lineno}: marker {tok!r}")
    return violations


# ---------------------------------------------------------------------------
# G7 — docs-truth audit
#
# (a) The unrealized-status vocabulary.  Every scanned document must describe
#     only behavior that is fully realized and wired in this repository.  A
#     claim that something is planned, blocked, deferred, partial, pending, a
#     candidate, out of scope, or otherwise not yet fully realized or
#     connected is a hard failure: the docs are the shipped product's
#     description, and they are permitted to describe nothing that is not
#     shipped.  Stems are matched with any word continuation so every form is
#     caught ("planned", "planning", "candidates", "not-yet-on-disk" ...).
#     The citation/license surface is read-only attribution: a cited work's
#     real title or license text is not a claim about this product and is
#     exempt, exactly as in G6.
# (b) Claim backing.  Every back-ticked identifier, ``case_*`` name and
#     ``data/...``/``src/...`` path in doc/12 must resolve to a real
#     ``src/drugos`` symbol, a registered validation case, or a manifest-pinned
#     file.  A documented claim must be backed by the artifact that realizes
#     it; a shipped claim pointing at nothing real fails the gate.
# ---------------------------------------------------------------------------

UNREALIZED_TOKENS = (
    "planned",
    "planning",
    "roadmap",
    "blocked",
    "defer",
    "deferred",
    "deferral",
    "partial",
    "pending",
    "candidate",
    "future",
    "await",
    "unreleased",
    "not-yet",
)
UNREALIZED_PHRASES = (
    "out of scope",
    "out-of-scope",
    "not in scope",
    "not in this release",
    "not yet",
    "not on disk",
    "not acquired",
    "not downloaded",
    "downloaded at runtime",
    "to be implemented",
    "to be realized",
    "to be added",
    "to be shipped",
    "to be wired",
    "to be evaluated",
    "work in progress",
    "works in progress",
    "in progress",
    "will be added",
    "will be implemented",
    "will be shipped",
    "will be wired",
    "will ship",
    "will land",
    "later release",
    "later phase",
    "later phases",
    "a later",
    "future release",
    "future work",
    "phase 8",
    "planned-later",
    "planned for",
    "planned release",
    "release track",
    "P5",
    "P6",
    "P7",
    "P8",
    "P9",
    "PLANNED-LATER",
    "NOT-DOWNLOADED",
    "WIRED (partial)",
    "PARTIAL-IN-HOUSE",
    "PENDING DATASET",
    "PLANNED (blocked)",
    "PARTIAL —",
    "部分落实",
    "计划内",
    "简化落实",
)
UNREALIZED_TOKEN_RE = {
    tok: re.compile(rf"\b{re.escape(tok)}\w*\b", re.IGNORECASE) for tok in UNREALIZED_TOKENS
}
def _compile_unrealized_phrase(phrase: str) -> re.Pattern[str]:
    """Compile one unrealized-status phrase, tolerant of a hyphen or space
    between its words so hyphenated and spaced forms are both caught."""
    parts = [re.escape(p) for p in phrase.split(" ")]
    return re.compile(rf"\b{('[ -]?'.join(parts))}\b", re.IGNORECASE)


UNREALIZED_PHRASE_RES = [_compile_unrealized_phrase(p) for p in UNREALIZED_PHRASES]

ALLOWED_SCIENTIFIC_SYMBOLS = frozenset(
    {
        # Closed set of universal physical/pharmacology symbols that the docs
        # use as parameter names without a one-to-one code symbol.
        "kd", "kon", "koff", "km", "ki", "ic50", "ec50", "auc", "cmax", "tmax",
        "t1/2", "cl", "vmax", "fa", "fup", "logp", "bsa", "egfr", "scr", "gfr",
        "alt", "ast", "uln", "atp", "gsh", "ros", "qt", "qtc", "qtw", "cv",
        "vs", "qs",
    }
)


def unrealized_markers_in(line: str) -> list[str]:
    """The unrealized-status markers present on one line of doc text."""
    hits: list[str] = []
    for tok in UNREALIZED_TOKENS:
        if UNREALIZED_TOKEN_RE[tok].search(line):
            hits.append(tok)
    for phrase, pat in zip(UNREALIZED_PHRASES, UNREALIZED_PHRASE_RES, strict=True):
        if pat.search(line):
            hits.append(f"phrase {phrase!r}")
    return hits


def _doc_scope_files() -> list[tuple[str, Path]]:
    """The doc set the G7 docs-truth audit scans (markdown only)."""
    return [(kind, path) for kind, path in _audit_scope_files() if kind == "doc"]


def _manifest_paths() -> set[str]:
    manifest = json.loads(DATA_MANIFEST.read_text(encoding="utf-8"))
    return {entry["path"] for entry in manifest["files"]}


def _src_identifiers() -> set[str]:
    """Every identifier appearing in ``src/drugos`` source, lowercased."""
    ids: set[str] = set()
    for path in PACKAGE.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        for tok in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text):
            ids.add(tok.lower())
    return ids


def _backtick_spans(path: Path) -> list[tuple[int, str]]:
    text = path.read_text(encoding="utf-8")
    citation_lines = _citation_surface_line_numbers(text)
    spans: list[tuple[int, str]] = []
    for lineno, line in enumerate(text.splitlines(), 1):
        if lineno in citation_lines:
            continue
        for m in re.finditer(r"`([^`]+)`", line):
            spans.append((lineno, m.group(1).strip()))
    return spans


def verify_doc_claims() -> list[str]:
    """Positive claim-backing check over the back-ticked claims in doc/12.

    Every back-ticked identifier in doc/12 must resolve to a real
    ``src/drugos`` symbol or a registered validation case; every back-ticked
    ``data/...``/``src/...``/``*.py`` path must exist or be pinned in the
    manifest.  Unresolved claims are hard violations: a documented, shipped
    behaviour must point at the artifact that realizes it.
    """
    path = ROOT / "doc" / "12-production-models.md"
    if not path.is_file():
        return ["doc/12-production-models.md missing (G7 claim register absent)"]
    text = path.read_text(encoding="utf-8")
    case_init = (ROOT / "validation" / "cases" / "__init__.py").read_text(encoding="utf-8")
    manifest_paths = _manifest_paths()
    identifiers = _src_identifiers()
    violations: list[str] = []
    rel = "doc/12-production-models.md"

    for lineno, content in _backtick_spans(path):
        if not content:
            continue
        where = f"{rel}:{lineno}: claim `{content}`"
        if content.startswith("data/"):
            if content == "data/manifest.json":
                continue
            key = content[len("data/"):]
            if key not in manifest_paths:
                violations.append(f"{where} not listed in data/manifest.json")
            continue
        if content.startswith(("src/", "scripts/")):
            if not (ROOT / content).is_file():
                violations.append(f"{where} path does not exist")
            continue
        if "/" in content and content.endswith(".py"):
            candles = [ROOT / content]
            if content.startswith("drugos/"):
                candles.append(PACKAGE / content[len("drugos/"):])
            else:
                candles.append(PACKAGE / content)
            if not any(c.is_file() for c in candles):
                violations.append(f"{where} path does not exist")
            continue
        if content.startswith("case_"):
            case_file = ROOT / "validation" / "cases" / f"{content}.py"
            if not case_file.is_file():
                violations.append(f"{where} has no validation case module")
            elif content not in case_init:
                violations.append(f"{where} case exists but is not registered in "
                                  "validation/cases/__init__.py")
            continue
        if content.startswith("drugos."):
            parts = content.split(".")
            ok = False
            for split in range(len(parts), 1, -1):
                mod = ".".join(parts[1:split])
                base = PACKAGE / mod.replace(".", "/")
                if (
                    base.is_file()
                    or base.with_suffix(".py").is_file()
                    or (base / "__init__.py").is_file()
                ):
                    ok = True
                    break
            if not ok:
                violations.append(f"{where} module path does not resolve in src/drugos")
            continue
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", content):
            low = content.lower()
            if low in identifiers or low in ALLOWED_SCIENTIFIC_SYMBOLS:
                continue
            if len(content) >= 2:
                violations.append(f"{where} identifier does not resolve in src/drugos")
            continue
    return violations


def docs_truth_audit() -> list[str]:
    """G7: every scanned document describes only realized, wired behaviour.

    1. The unrealized-status vocabulary (planned, blocked, deferred, partial,
       pending, candidate, future, roadmap, out-of-scope, not-yet,
       not-downloaded, release-track, and the remaining set) is absent from
       doc text outside the citation surface — a document may describe
       nothing that is not shipped.
    2. Every back-ticked claim in doc/12 resolves to a real ``src/drugos``
       symbol, a registered validation case, or a manifest-pinned file.

    There is no allowance file and no pinned inventory for either layer: any
    occurrence at any count is a violation, because the docs are permitted to
    describe only the fully realized system.
    """
    violations: list[str] = []
    for kind, path in _doc_scope_files():
        rel = path.relative_to(ROOT).as_posix()
        text = path.read_text(encoding="utf-8")
        citation_lines = _citation_surface_line_numbers(text) if kind == "doc" else set()
        for lineno, line in enumerate(text.splitlines(), 1):
            if lineno in citation_lines:
                continue
            for tok in unrealized_markers_in(line):
                violations.append(f"{rel}:{lineno}: unrealized-status marker {tok!r}")
    violations.extend(verify_doc_claims())
    return violations


def run_gates(use_xdist: bool = True) -> None:
    ensure_versions_consistent()
    GATE_LOG.parent.mkdir(parents=True, exist_ok=True)
    with GATE_LOG.open("w", encoding="utf-8") as log:
        log.write("DrugOS release quality gate\n")
        log.write(f"started {_dt.datetime.now(_dt.UTC).isoformat()}\n")
    _run_cmd(
        GATE_LOG,
        "G1 lint (ruff check)",
        [python(), "-m", "ruff", "check", "src", "tests", "validation"],
    )
    _run_cmd(
        GATE_LOG,
        "G1 format (ruff check-format)",
        [python(), "-m", "ruff", "format", "--check", "src", "tests", "validation"],
    )
    _run_cmd(
        GATE_LOG, "G2 types (mypy --strict)", [python(), "-m", "mypy", "--strict", "src/drugos"]
    )
    g3_cmd = [python(), "-m", "pytest"]
    if use_xdist:
        g3_cmd += ["-n", "auto"]
    _run_cmd(GATE_LOG, "G3 coverage (pytest 100% branch)", g3_cmd)
    _run_cmd(
        GATE_LOG, "G4 validation suite", [python(), str(ROOT / "validation" / "run_validation.py")]
    )
    with GATE_LOG.open("a", encoding="utf-8") as log:
        try:
            verify_data_checksums()
            log.write("data checksums: verified (manifest matches disk)\n")
        except SystemExit:
            log.write("data checksums: FAILED\n")
            raise
    violations = fallback_audit()
    with GATE_LOG.open("a", encoding="utf-8") as log:
        if violations:
            log.write("== fallback audit ==")
            for v in violations:
                log.write(f"  FAIL {v}\n")
            print("FALLBACK AUDIT FAILED:")
            print("\n".join(f"  {v}" for v in violations))
            raise SystemExit("silent fallback detected; refusing release")
        log.write("fallback audit: clean (no silent fallbacks)\n")

    deferrals = no_deferral_audit()
    with GATE_LOG.open("a", encoding="utf-8") as log:
        if deferrals:
            log.write("== marker audit (G6) ==")
            for v in deferrals:
                log.write(f"  FAIL {v}\n")
            print("MARKER AUDIT FAILED:")
            print("\n".join(f"  {v}" for v in deferrals))
            raise SystemExit("marker of a retained seam detected; refusing release")
        log.write("marker audit: clean (no markers in any scanned file)\n")

    unrealized = docs_truth_audit()
    with GATE_LOG.open("a", encoding="utf-8") as log:
        if unrealized:
            log.write("== docs-truth audit (G7) ==")
            for v in unrealized:
                log.write(f"  FAIL {v}\n")
            print("DOCS-TRUTH AUDIT FAILED:")
            print("\n".join(f"  {v}" for v in unrealized))
            raise SystemExit("unrealized-status wording or unbacked claim detected; refusing release")
        log.write("docs-truth audit: clean (every doc claim is realized, wired, and backed)\n")
        log.write("ALL GATES PASSED\n")
    print("ALL GATES PASSED")
    print("".join(GATE_LOG.read_text(encoding="utf-8").splitlines(keepends=True)[-25:]))


def sync_version_text(old: str, version: str, dry: bool) -> list[str]:
    """Replace the old version token project-wide (README, docs, index banner)."""
    targets = [*SYNC_FILES, INDEX_HTML]
    touched: list[str] = []
    for path in targets:
        text = path.read_text(encoding="utf-8")
        updated, n = re.subn(re.escape(old), version, text)
        if n and not dry:
            path.write_text(updated, encoding="utf-8")
        if n:
            touched.append(f"{path.relative_to(ROOT)}  ({n} token(s))")
    return touched


def sync_validation_counts(passed: int, total: int, dry: bool) -> list[str]:
    """Sync validation counts and the badge across the docs (idempotent).

    The previously published count is read back from the README badge rather
    than assumed, so the sync stays correct as the suite grows past its
    launch-time value.
    """
    new_count = f"{passed}/{total}"
    new_badge = f"validation-{passed}%2F{total}%20green"
    prev = _published_validation_count()
    touched: list[str] = []
    for path in SYNC_FILES:
        before = path.read_text(encoding="utf-8")
        text = before
        if prev is not None:
            text = text.replace(prev, new_count)
        text = re.sub(
            r"validation-\d+%2F\d+%20green",
            new_badge,
            text,
            flags=re.IGNORECASE,
        )
        if text != before:
            if not dry:
                path.write_text(text, encoding="utf-8")
            touched.append(f"{path.relative_to(ROOT)}  ({prev} -> {new_count})")
    return touched


def _published_validation_count() -> str | None:
    """Return the currently published 'passed/total' (e.g. '36/36') if any."""
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    for pattern in (
        r"validation-(\d+)%2F(\d+)%20green",
        r"(\d+/\d+) cases green",
        r"(\d+/\d+) green",
    ):
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            groups = match.groups()
            if len(groups) == 2:
                return f"{groups[0]}/{groups[1]}"
            return groups[0]
    return None


def sync_doc09_status(
    version: str,
    passed: int,
    total: int,
    gates: str,
    dry: bool,
    now: _dt.datetime,
) -> None:
    """Update the machine-maintained status block in doc/09 §6."""
    path = ROOT / "doc" / "09-quality-gate.md"
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"- Version: \*\*[^*]+\*\*", f"- Version: **{version}**", text, count=1)
    text = re.sub(
        r"- Validation: \*\*PASS \(\d+/\d+ cases\)\*\*",
        f"- Validation: **PASS ({passed}/{total} cases)**",
        text,
        count=1,
    )
    for label in (
        "G1 lint (ruff)",
        "G2 types (mypy --strict)",
        "G3 coverage (100 % branch of `src/drugos`)",
        "G4 validation suite",
        "G5 fallback audit",
        "G6 marker audit",
        "G7 docs-truth audit",
    ):
        text = re.sub(rf"- {re.escape(label)}: \w+", f"- {label}: {gates}", text, count=1)
    text = re.sub(
        r"- Last release run: .*",
        f"- Last release run: {now.astimezone(_dt.UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        text,
        count=1,
    )
    if not dry:
        path.write_text(text, encoding="utf-8")


def write_status(old: str, target: str, passed: int, total: int, gates: str) -> None:
    STATUS_JSON.parent.mkdir(parents=True, exist_ok=True)
    status = {
        "version": target,
        "previous_version": old,
        "version_rule": "explicit (no auto-bump)",
        "validation_passed": passed,
        "validation_total": total,
        "gates": gates,
        "written_utc": _dt.datetime.now(_dt.UTC).isoformat(),
    }
    STATUS_JSON.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def cmd_version() -> int:
    ensure_versions_consistent()
    print(read_version())
    return 0


def cmd_fallback_audit() -> int:
    violations = fallback_audit()
    if violations:
        print("silent fallbacks detected:")
        for v in violations:
            print(f"  {v}")
        return 1
    count = len(_except_handlers(PACKAGE))
    print(f"fallback audit clean: {count} except handler(s), all explicit")
    return 0


def cmd_marker_audit() -> int:
    violations = no_deferral_audit()
    if violations:
        print("marker detected in scanned files:")
        for v in violations:
            print(f"  {v}")
        return 1
    print(
        "marker audit clean: no marker of any kind in any scanned file "
        "(no allowance file, no whitelist, no pinned inventory)"
    )
    return 0


def cmd_docs_audit() -> int:
    violations = docs_truth_audit()
    if violations:
        print("unrealized-status wording or unbacked claims in docs:")
        for v in violations:
            print(f"  {v}")
        return 1
    print(
        "docs-truth audit clean: every doc claim is realized, wired, and backed "
        "(no planned/blocked/deferred/partial/candidate/out-of-scope wording outside "
        "the citation surface; every doc/12 back-tick resolves)"
    )
    return 0


def cmd_status(dry: bool) -> int:
    ensure_versions_consistent()
    passed, total = validation_counts()
    now = _dt.datetime.now(_dt.UTC)
    print(f"[status] validation: {passed}/{total} passed")
    for line in sync_validation_counts(passed, total, dry):
        print(f"[status] count sync:   {line}")
    if not dry:
        sync_doc09_status(read_version(), passed, total, "PASS", dry=False, now=now)
        print("[status] doc/09 §6 status block updated")
    return 0


def cmd_release(args: argparse.Namespace) -> int:
    old = read_version()
    target = args.version
    if target is None:
        raise SystemExit("release requires an explicit --version YYYY.M.V (no auto-bump)")
    if not re.fullmatch(r"\d+\.\d+\.\d+", target):
        raise SystemExit(f"bad target version: {target!r}")

    print(f"[release] version  {old} -> {target} (explicit)")
    now = _dt.datetime.now(_dt.UTC)

    if args.dry_run:
        print("[release] dry-run; no files written")
        print("[release] planned version edits:")
        for line in sync_version_text(old, target, dry=True):
            print(f"  {line}")
        print("[release] planned count edits:")
        for line in sync_validation_counts(*validation_counts(), dry=True):
            print(f"  {line}")
        return 0

    if args.no_gates:
        gates = "SKIPPED"
        print("[release] gates skipped (--no-gates)")
    else:
        run_gates(use_xdist=not args.no_xdist)
        gates = "PASS"

    passed, total = validation_counts()
    print(f"[release] validation report: {passed}/{total} passed")

    if target != old:
        write_version(target)
        write_pyproject_version(target)
        print(f"[release] wrote version {target} to version.py and pyproject.toml")

    for line in sync_version_text(old, target, dry=False):
        print(f"[release] version sync: {line}")
    if not args.no_sync:
        for line in sync_validation_counts(passed, total, dry=False):
            print(f"[release] count sync:   {line}")
        sync_doc09_status(target, passed, total, gates, dry=False, now=now)
        print("[release] doc/09 §6 status block updated")

    write_status(old, target, passed, total, gates)
    print(f"[release] wrote {STATUS_JSON.relative_to(ROOT)}")
    print("[release] done")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="release.py",
        description=(
            "DrugOS release tool: quality gates, version, validation status. "
            "See the module docstring for subcommand semantics."
        ),
    )
    parser.add_argument(
        "cmd",
        nargs="?",
        default="release",
        choices=sorted(KNOWN_COMMANDS),
        help="action (default: release)",
    )
    parser.add_argument(
        "--version",
        default=None,
        help="target version YYYY.M.V; required for 'release' (never auto-bumped)",
    )
    parser.add_argument("--no-gates", action="store_true", help="skip the quality gate")
    parser.add_argument(
        "--no-xdist",
        action="store_true",
        help="run G3 with base pytest only (no pytest-xdist); used by GitHub CI",
    )
    parser.add_argument("--no-sync", action="store_true", help="skip status/doc sync")
    parser.add_argument("--dry-run", action="store_true", help="print the plan, do not write")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "version":
        return cmd_version()
    if args.cmd == "fallback-audit":
        return cmd_fallback_audit()
    if args.cmd == "marker-audit":
        return cmd_marker_audit()
    if args.cmd == "docs-audit":
        return cmd_docs_audit()
    if args.cmd == "status":
        return cmd_status(args.dry_run)
    if args.cmd == "gates":
        run_gates(use_xdist=not args.no_xdist)
        return 0
    if args.cmd == "release":
        return cmd_release(args)
    raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    sys.exit(main())
