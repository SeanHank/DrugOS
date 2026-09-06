# Technology Stack & Dependencies

## 1. Language & Runtime

- **Canonical Python interpreter (project-wide):** `/opt/anaconda3/envs/drug_os/bin/python`
  - All development, execution, and CI-facing commands use this interpreter (or the `drug_os` Conda environment; `conda run -n drug_os ...`).
  - The interpreter is pinned to Python 3.11+; **do not** call a system/global `python`.
- **Required R** integration via `rpy2` (>= 3.6, hard dependency in
  `pyproject.toml`) — every pipeline run verifies the clearance/AUC estimator
  against the **original literature equations in R**
  (`src/drugos/rbridge/literature_pk.R`: Wagner 1976, Gibaldi & Perrier 1982,
  Greenblatt & Koch-Weser 1975, Rowland & Tozer 2010). The bridge lives in
  `src/drugos/rbridge/` (100%-branch tested, `stubs/rpy2/robjects.pyi` for
  mypy) and is exercised by validation case R-1. A missing R runtime is a hard
  error — never a silent fallback. For local R: R >= 4 are supported
  (`/usr/local/bin/R` macOS Homebrew; conda `r-base`). The standalone,
  off-gate R cross-validation harness `scripts/r_crossval/` (base-R) remains
  as an independent second estimator lane.
- Runtime requirement: `R` / `Rscript` must be on `PATH` at import time of
  `drugos.rbridge`; verdicts stream to the `r_verify` contract block and the
  web report.

## 2. Core Stack by Concern

### 2.1 Chemistry

| Library | Purpose |
|---|---|
| RDKit | SMILES/InChI/SDF parsing, canonicalization, descriptor computation (LogP, LogD, pKa, TPSA, HBD/HBA, ESOL) |
| Open Babel (optional) | format conversion fallback |

### 2.2 ADME / ML Prediction

| Library | Purpose |
|---|---|
| admet-ai (open source) | GNN-based ADME/T prediction from structure, DrugBank-reference percentiles |
| chemprop (underlying) | training/fine-tuning structure->property GNNs |
| torch / torch-geometric | ML model inference and custom DTI/DTA models (AttentionDTA-style) |
| scikit-learn / xgboost | descriptor-based QSAR classifiers for toxicity risk fusion |

### 2.3 ODE / Simulation Core

| Library | Purpose |
|---|---|
| SciPy (`solve_ivp`, LSODA/Radau) | ODE solves for compartmental PBPK, occupancy, pathway, and organ models (PySB/pysb-pkpd dropped in 2026.9.x — rule-based sugar replaced by direct ODE implementation) |
| OpenPKPD | population PK/PD estimation/simulation alternative; SBML import |
| myokit (optional) | cardiac cell models (QT/hERG axis) |
| numba / numba jit | optional speedup of inner loops |

### 2.4 Data Engineering & Storage

| Library | Purpose |
|---|---|
| pandas / polars | tabular data (physiology tables, drug-target data) |
| numpy | numeric substrate |
| sqlite3 / DuckDB | local caches, small-DB handles |
| pydantic | data-contract validation (unit + schema checks across stage boundaries) |
| Parquet | serialized derived datasets |

### 2.5 Reporting & Visualization

| Library | Purpose |
|---|---|
| matplotlib / plotly | concentration-time, occupancy, pathway, biomarker plots; interactive HTML reports |
| jinja2 | markdown/HTML report templates |
| JSON schema | structured output contract (`report/schema.json`) |

### 2.6 Optional External Integrations

| Tool | Use |
|---|---|
| Open Systems Pharmacology Suite (PK-Sim/MoBi, Windows or d4k) | cross-validation of Stage 1 against the validated library models; physiology-table source |
| `ospsuite` R package | scripted OSPOps interaction where desired |
| Docking (e.g., AutoDock Vina / smina) | pose generation for structure-based affinity in Stage 2 |
| AlphaFold DB | target structures when PDB absent |

## 3. Recommended Environment

Prerequisites: Conda environment `drug_os` created on
`/opt/anaconda3/envs/drug_os`. Example bootstrap:

```
/opt/anaconda3/envs/drug_os/bin/python -m pip install -r requirements.lock
```

```
# requirements (core)
python>=3.11
rdkit
numpy, scipy, pandas
matplotlib, plotly, jinja2
pydantic
admet-ai
torch (cpu ok for inference)
```

## 3.1 Project Versioning

- The whole project follows the `YYYY.M.V` scheme (year, month, intra-month revision starting at 0) — see `01-project-overview.md` section 0.1.
- The version is emitted by `src/drugos/version.py` and mirrored in `pyproject.toml`; package and report versions must match. Baseline: **2026.9.0**.

## 4. Data Files Layout

```
DrugOS/
  doc/                     # design documents (this tree)
  LICENSE                  # AGPLv3
  data/                    # caches, manifests
    manifest.yaml
    physiology/            # OSP/Willmann-derived tables
    admet/                 # admet-ai + TDC model artifacts
    targets/               # drugbank/chEMBL/bindingDb panels
    pathways/              # KEGG/Reactome-derived graphs
    toxicity/              # DILI/hERG/in-vitro IC50 panels
  src/drugos/
    version.py             # YYYY.M.V project version
    inputs/ pk/ target/ pathway/
    organ/ clinical/ report/
  tests/
  notebooks/
  README.md
```

## 5. Version Pinning & Reproducibility

- `requirements.lock` via pip freeze; ML weights pinned by checksum.
- Every simulation run emits a run-manifest (drug params, human params, model versions, solver settings, seeds) for full reproducibility.
- Deterministic seeds for any stochastic/ML components.
- The project version and interpreter in use (`/opt/anaconda3/envs/drug_os/bin/python --version`) are recorded in every run-manifest and report footer.
- License: AGPLv3 (author: Sean Hank). All contributions are subject to AGPLv3.