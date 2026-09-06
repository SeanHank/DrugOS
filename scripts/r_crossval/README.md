# R cross-validation harness (`scripts/r_crossval`)

Optional R-language integration point documented in `doc/06` (rpy2 was listed
there as "optional R support, not required for the baseline release"). This
directory turns that into an executable R cross-validation scaffold.

## Why an R bridge at all

DrugOS computes the PK/organ curves in Python (numpy/scipy). The
pharmacometrics community standard for model-based cross-validation is R
ecosystem tooling (`nlmixr2`, `nonmem`-style VPC, `ospsuite`). Instead of only
claiming "the published bands match", this harness re-fits the *simulated*
individual curve in R with an independent estimator and checks it against the
same published pass-bands DrugOS validates against (`data/benchmarks/
published_pk.json`), producing an independent, reproducible verdict.
Targeting "100 % realism" = both code bases must independently reach the same
model-measurement agreement.

## Layout

```
scripts/r_crossval/
  run.py     # python loader: dumps a DrugOS run to TSV, then (optionally) calls R
  run.r      # base-R script: independent nls fit + AUC/band verdict -> markdown report
  README.md  # this file
```

## Requirements

- Python side: the installed `drugos` package (any route).
- R side: base R **only** (`stats`, `utils`), no CRAN installs required. For a
  deeper `nlmixr2`/VPC pass, install `nlmixr2` separately — `run.r` auto-detects
  it and adds a Newton–Raphson 2-compartment re-fit when present.

## Usage

```bash
# 1) generate the simulated curve TSVs (no R needed):
/opt/anaconda3/envs/drug_os/bin/python scripts/r_crossval/run.py --benchmark warfarin --dose 5

# 2) run the R cross-validation verdict:
Rscript scripts/r_crossval/run.r \
  --sim /tmp/r_crossval/warfarin/sim.tsv \
  --bands /tmp/r_crossval/warfarin/bands.tsv \
  --out /tmp/r_crossval/warfarin/report.md

# or both at once (run.py calls Rscript if it is on PATH):
/opt/anaconda3/envs/drug_os/bin/python scripts/r_crossval/run.py --benchmark dofetilide --dose 0.5 --with-r
```

## What the R script declares

1. Trapezoidal `AUC_t` on the simulated curve (base `integrate` over an
   `approxfun`, relaxed tolerance — independent of the Python NCA).
2. Terminal `t_1/2` from a robust log-linear fit on the tail (last 25 % of
   time > 5 % Cmax), plus an `nls()` mono-exponential alternative when
   `nlmixr2` is not installed.
3. $\lambda_z$ slope estimate with standard error; a simple 90 % band on the
   tail prediction.
4. Band verdict: `PASS / WARN / FAIL` when the R-estimated metrics leave the
   published band from `bands.tsv` (the same values DrugOS validates on).
5. A markdown report with the fitted equations, residuals, and the verdict —
   a machine-checkable cross-check that the Python result is not a numerical
   artifact.

## Licensing of the band values

The TSV bands come from `data/benchmarks/published_pk.json` (published numeric
ranges; academic citation use per `data/README.md`). This harness does not
redistribute third-party text; it only reads the same vendored reference the
Python gates already verify.

## Status

- Implemented: R script + Python loader (base-R only).
- Not yet executed in this environment: **R / Rscript are not installed**
  (`which R` → not found). CI does not depend on this harness; it is an
  optional, off-gate cross-validation lane. Run it in any R-enabled
  environment to produce `report.md`.