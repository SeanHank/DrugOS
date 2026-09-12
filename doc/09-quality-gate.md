# Quality Gate (Release-Enforced)

The DrugOS repository is gated by six mandatory checks (G1-G6). **No release,
commit to the protected branches, or "all green" claim is valid while any gate
fails.** The gates are deliberately strict: a mechanistic modeling platform that
feeds safety reasoning must fail loudly rather than excuse itself.

The current gate suite is defined by four quality gates (G1-G4) plus two
audits: the no-silent-fallback audit (G5) and the marker audit (G6), run
locally by `scripts/release.py gates` (the merged successor of the former
`quality_gate.sh`) and in CI by the `ci.yml` workflow.

## 1. The Gates

| Gate | Tool | Rule |
|---|---|---|
| G1 | `ruff check` + `ruff format --check` | zero findings, zero violations of configured rules, zero `format` diffs |
| G2 | `mypy --strict` | zero type errors across `src/drugos` on the pinned `py312` target |
| G3 | `pytest` + `coverage` | 100 % branch coverage of `src/drugos`; `fail_under = 100`. GitHub CI runs **base pytest only** (`python scripts/release.py gates --no-xdist`); local runs may add `-n auto` (pytest-xdist) for speed |
| G4 | validation suite | `validation/` runs end-to-end, every case must **pass**, and `validation/report.md` must be regenerated cleanly |
| G5 | fallback audit | every `except` handler in `src/drugos` must surface an explicit error; the handler inventory is pinned in `scripts/fallback_allowlist.json` and any drift or silent swallow fails the release |
| G6 | marker audit | code and docs carry no marker of a stand-in yet to be realized: any occurrence of the audited marker set at any count in any audited file is a release-blocking violation, with no allowance file, no whitelist and no pinned inventory |

### G1 — Lint and Format

- Configured in `pyproject.toml` (`[tool.ruff]`): `select = ["E","F","W","I","UP","B"]`,
  `ignore = ["B008"]`, `line-length = 100`.
- Run: `ruff check src tests validation`, `ruff format --check src tests validation`.
- **Exclusion rule.** Blanket suppressions are forbidden: no `# noqa` without a
  targeted rule code and a one-line justification, no module-level
  `# flake8: noqa`, no config-level `extend-ignore` of categories. Adding an
  ignore requires a design-doc amendment, not an inline patch.

### G2 — Types

- `mypy --strict src/drugos` with `python_version = "3.12"`.
- **Exclusion rule.** No `# type: ignore`, no `# mypy: ignore-errors`, no
  `# type: ignore[...]` with an empty rationale, no `Any`-escaping shortcuts.
  Where a third-party package ships no typing declarations, the gap is
  recorded in
  `doc/06-technology-stack.md` and the typing packages under `stubs/` cover
  it (G2 stays branch-exact for both mypy and coverage).

### G3 — Coverage

- `pytest` runs the full `tests/` tree; in GitHub CI this is **base pytest**
  (`release.py gates --no-xdist`), while local runs may append `-n auto`
  (pytest-xdist) to parallelize. `coverage` measures `src/drugos` with
  `fail_under = 100` (statement, function and branch totals recorded in the
  gate log).
- **Exclusion rule (strictest in the suite).** Coverage exclusions of any form
  are disallowed:
  - no `# pragma: no cover`,
  - no `[tool.coverage.run] omit` / `exclude_lines` derived from the source,
  - no `[tool.coverage.report] exclude` in the *project* config,
  - no diff-limited or "only-changed-lines" reporting for the release gate,
  - no separate "debugging-only" files that are exempt.
  If the 100 % bar cannot be met for a module, the module is redesigned to be
  testable or removed from the release; the bar is not lowered.
  (The one structural exception is the *reported* HTML/terminal annotation —
  measurement still includes everything under `src/drugos`.)

### G4 — Scientific Validation

- `validation/` is not an approximation of the pipeline: it is the check that
  the pipeline agrees with reality. Each validation case encodes a published
  observation (PK, occupancy, pathway, organ, clinical) and a quantitative
  pass criterion with its source citation.
- Rules:
  - the suite must **compile and run** with `python scripts/release.py gates`;
  - **all** cases must pass (no `xfail`/skip-without-documentation, no "known
    failures" list);
  - every run regenerates `validation/report.md` from the same code path that
    produced it (single source of truth);
  - new model features must add a validation case before merge, not after.

### G5 — No Silent Fallback

- **Policy.** There is no silent fallback of any form in `src/drugos`: an
  `except` handler must either re-raise, raise a new error, or return an
  explicit error value; a handler that swallows the exception and continues
  with a default is a release-blocking violation.
- `scripts/release.py gates` runs an AST scan of every handler in the package.
  Each handler must satisfy the explicit-error rule *and* the exact handler
  inventory must match `scripts/fallback_allowlist.json` — any new, removed or
  relocated handler fails the gate until the allowlist is updated on purpose
  (the allowlist documents why each remaining handler is permitted).
- The same no-silent-fallback rule applies to user-facing coercion: malformed
  CLI/HTTP inputs raise explicit errors (a CLI error + nonzero exit, or an HTTP
  400) instead of silently substituting a default (see `web/app.py` and
  `cli.py`), and missing ADMET data columns raise rather than degrade silently.

### G6 — Marker Audit

- **Policy.** Every audited code file and every audited doc file must be free
  of any marker of a stand-in that is not realized and not final.
- The audited marker set is defined by the hard and soft marker constants in
  `scripts/release.py` and is reviewable there. No scanned document lists the
  set, because a scanned file may not contain a marker at all: the definition
  lives only in the self-excluded checker, and it admits no edits at release
  time.
- **Stem + inflection.** Markers are matched on their stem plus any word-char
  continuation, so no plural or derived form of a marker (a repeated or
  inflected occurrence of the prohibited vocabulary) bypasses the gate.
- **Zero tolerance, no mechanism of allowance.** Any occurrence of the marker
  set at any count in any audited file is a release-blocking violation. The
  audit keeps no allowance file, no whitelist and no pinned inventory of
  permitted items, and it accepts no exclusions, annotations, pragmas or
  ledger entries that would admit one.
- The audited scope is `src/drugos`, `tests`, `validation`, `scripts`
  (`*.py`, excluding the checker itself), plus `README.md`, `DISCLAIMER.md`,
  `CONTRIBUTING.md` and `doc/*.md`. Generated `validation/report.md` is not in
  scope (it is the G4 output, not a source document).
- The audit is **read-only** and never counts or flags a license or citation
  line: the third-party attribution surface (`doc/10`, `doc/11`, `doc/12`
  License columns, the bibliography) is the citation red line and is never
  touched by the gate, and no marker rule can be waived by citing it.
- **Recognized artifact paths.** The G2 typing-declarations directory
  `stubs/` (`stubs/libsbml`, `stubs/rpy2`) is a committed, branch-exact set
  of type-declaration packages backing 100%-tested modules — a realized
  artifact, not an unbuilt stand-in. The literal path segment is recognized
  by the audit so the directory name is not misread as the prohibited word;
  any use of the word outside that path is audited like any other occurrence.
- `scripts/release.py gates` runs the audit; standalone:
  `python scripts/release.py marker-audit`.

## 2. Local Commands

```bash
/opt/anaconda3/envs/drug_os/bin/python -m ruff check src tests validation
/opt/anaconda3/envs/drug_os/bin/python -m ruff format --check src tests validation
/opt/anaconda3/envs/drug_os/bin/python -m mypy --strict src/drugos
/opt/anaconda3/envs/drug_os/bin/python -m pytest -n auto               # full suite, parallel
/opt/anaconda3/envs/drug_os/bin/python validation/run_validation.py   # regenerates validation/report.md
/opt/anaconda3/envs/drug_os/bin/python scripts/release.py gates       # one-shot: all of the above
# GitHub CI runs the same gate with base pytest only:  gates --no-xdist
python scripts/release.py fallback-audit                                # G5 standalone
python scripts/release.py marker-audit                                  # G6 standalone
```

## 3. Configuration Residency

- Coverage measurement parameters live in `pyproject.toml`
  (`[tool.coverage.run]`, `[tool.coverage.report]`, `[tool.pytest.ini_options]`).
  Reviewers treat any edit that shrinks measured scope as a gate-relaxation and
  require an accompanying rationale in the commit message.
- Mypy settings live in `pyproject.toml` (`[tool.mypy]`).
- `scripts/release.py` (gates / fallback-audit / marker-audit /
  version / release / status) is the single entry point for CI and pre-merge
  checks; local phase passes that bypass it do not count.

## 4. CI Enforcement

- `.github/workflows/ci.yml` runs `python scripts/release.py gates` (G1-G6) on
  every push to the protected branches, on pull requests, and on tagged pushes.
- The `quality` job installs the package plus `.[dev]`, reads the project
  version from `pyproject.toml`, and uploads the gate log and regenerated
  `validation/report.md` as artifacts on any outcome.
- When the gate is green and artifacts build cleanly, CI ships a **GitHub
  Release**:
  - **`build` job** (depends on `quality`): reads the release version from
    `pyproject.toml` (never hard-coded, never inferred from the tag),
    `python -m build` produces the wheel (`dist/*.whl`) and sdist
    (`dist/*.tar.gz`), `twine check` validates the metadata, and the wheel is
    smoke-tested in a clean venv (imports, `py.typed` present, `create_app()`
    reads the packaged template, and the installed `drugos` CLI runs). On
    tagged pushes it verifies the tag equals the `pyproject.toml` version and
    creates a GitHub Release with the artifacts attached (`gh release create`,
    auto-generated notes).
  - CI publishes exclusively to GitHub Releases; no PyPI upload is performed.

## 5. Definition of "Release"

A release (`2026.*`) requires:

1. G1-G6 green in CI **and** on the pinned conda interpreter locally;
2. `validation/report.md` checked in, current, and reproducible;
3. `doc/*.md` describing exactly the shipped behavior (docs are part of the
   gate: a doc that contradicts the code fails review);
4. the parameter manifest / data checksums documented and pinned;
5. no silent fallbacks anywhere in `src/drugos` (the fallback audit is clean and
   the allowlist matches the code);
6. no marker of a stand-in yet to be realized in any audited code or doc file
   (the marker audit is clean; there is no allowance file, no whitelist and
   nothing is pinned).

## 6. Current Status

> This block is machine-maintained by `scripts/release.py` (and verified by
> the gate) — do not edit by hand.

- Version: **2026.9.1**
- Validation: **PASS (38/38 cases)**
- G1 lint (ruff): PASS
- G2 types (mypy --strict): PASS
- G3 coverage (100 % branch of `src/drugos`): PASS
- G4 validation suite: PASS
- G5 fallback audit: PASS
- G6 marker audit: PASS
- Last release run: 2026-09-12 04:49 UTC