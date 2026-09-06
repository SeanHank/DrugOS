# Contributing

Thanks for considering a contribution to DrugOS. This project is a
release-grade research platform: contributions are welcome, and they are judged
by the same four gates that protect every release.

## Code of conduct

Be constructive, be specific, and assume good faith. Mechanistic-modeling
debates are scientific debates; argue from equations, data, and validation
cases, not from authority.

## Ground rules

- **A change is not merged until all four gates are green**:
  G1 ruff, G2 `mypy --strict`, G3 100 % branch coverage of `src/drugos`,
  G4 the full validation suite. `python scripts/release.py gates` runs them all
  (plus the no-silent-fallback audit, G5).
- **Docs are part of the gate.** If code behavior changes, the corresponding
  section of `doc/*.md` changes in the same PR. A doc that contradicts the code
  fails review.
- **Blanket suppressions are forbidden.** No `# type: ignore`, no bare
  `# noqa`, no `# pragma: no cover`, no coverage `omit`/`exclude` for source.
  If a rule genuinely cannot be satisfied, amend the design doc — don't patch
  around the gate.
- **New model features add a validation case before merge, not after.** The
  feature ships with the empirical/analytic check that pins its behavior.

## How to propose a change

1. **Open an issue first** for anything that moves the model or the gates.
   Small bug fixes and tests can skip straight to a pull request.
2. Branch from the default branch; keep the change small and focused.
3. Run the gates before pushing:
   ```bash
   python scripts/release.py gates
   ```
4. In the PR description, state what changed, why, and which validation case
   (existing or new) demonstrates it.

## Running the gates

```bash
python -m ruff check src tests validation
python -m ruff format --check src tests validation
python -m mypy --strict src/drugos
python -m pytest -n auto          # full suite, parallel
python validation/run_validation.py # regenerates validation/report.md
python scripts/release.py gates    # one-shot: all of the above
python scripts/release.py fallback-audit                            # every except handler must surface an error
```

## Adding a validation case

- Place the case in `validation/cases/`. Each case encodes a published
  observation (PK, occupancy, pathway, organ, clinical) and a quantitative pass
  criterion with its source citation.
- The suite must run end-to-end with `validation/run_validation.py` and the
  regenerated `validation/report.md` must be committed.
- **All cases must pass.** There are no `xfail` entries and no "known
  failures" list.

## Releasing

Releases are `YYYY.M.V` (e.g. `2026.9.0`). `scripts/release.py` runs the four
gates + fallback audit, bumps the version project-wide, syncs the status
numbers in `README.md` and `doc/*.md`, and writes `build/release_status.json`:

```bash
python scripts/release.py --bump auto
```

See `doc/09-quality-gate.md` for the full release contract.