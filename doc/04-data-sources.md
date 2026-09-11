# Data Sources & Data Engineering

This document catalogs the external data sources used across the pipeline and how they are ingested, cached, and versioned.

## 1. Catalog

### 1.1 Chemical & Physicochemical Properties

| Source | Content | Stage use |
|---|---|---|
| RDKit (in-house computation) | ln descriptors: MW, LogP, LogD(pH), pKa, TPSA, HBD/HBA, rotatable bonds | Stage 1 physchem |
| PubChem | SMILES/InChI canonicalization, aliases | Input normalization |
| ChEMBL | measured bioactivity (IC50/EC50/Ki/Kd) per target | Stage 2 Kd priors, calibration |

### 1.2 ADME / PK Data

| Source | Content | Stage use |
|---|---|---|
| TDC (Therapeutics Data Commons) | 41 ADMET datasets | Training/benchmarking ADMET-AI-style models |
| ADMET-AI (open-source package, Zenodo-archived models) | pre-trained GNN weights + DrugBank reference percentiles | Structure -> ADME prediction |
| OSP PBPK Model Library | validated whole-body PBPK models + evaluation reports for ~40+ compounds (antibiotics, antivirals, sedatives, cardiovascular agents, etc.) | Stage 1 parameter priors, validation benchmarks |
| Physiome Model Repository | organ models + physiological parameter sets | Stage 4 organ parameterization |

### 1.3 Physiology / Virtual Human Data

| Source | Content | Stage use |
|---|---|---|
| OSP physiology database (from PK-Sim) | organ volumes, blood-flow rates, vascular/interstitial/cellular space volumes for human + lab animals | Stage 1 PBPK compartment parameters |
| Barter et al. (2013) CYP immunoquantification | per-isoform hepatic microsomal abundances (CYP1A2/2A6/2B6/2C8/2C9/2C19/2D6/2E1/3A4 = `CYP_ABUNDANCE_PMOL_MG`) | Segment 1 enzyme-kinetics stage (abundance-scaled per-CYP Vmax via `cyp_vmax_mg_h`, consumed by `PBPKModel.cyp_terms`) |
| Willmann et al. (2007) population model equations | covariate->parameter equations (age, sex, height, weight, BMI) | Human parameter resolver |
| ICBP/Physiome data | organ-level function baselines | Stage 4-5 reference ranges |

### 1.4 Target & Pharmacology Data

| Source | Content | Stage use |
|---|---|---|
| DrugBank | approved drug targets (UniProt IDs, gene symbols, mechanism), pharmacology, carriers/transporters | Stage 2 target resolver + safety targets |
| UniProt | target protein sequences (for sequence-based DTI ML) | Stage 2 |
| PDB / AlphaFold DB | experimental / predicted structures | Stage 2 structure-based affinity |
| BindingDB | measured binding affinities and kinetics | Stage 2 Kd/kon/koff parameterization |

### 1.5 Pathway Data

| Source | Content | Stage use |
|---|---|---|
| KEGG | organism-specific pathway graphs (hsa:) with gene/protein nodes | Stage 3 scaffold |
| Reactome | curated signaling/regulatory reactions with causal relations | Stage 3 scaffold |
| PANTHER | pathway terms + protein membership | Stage 3 annotation |
| WikiPathways | community-edited pathway wiring | Stage 3 supplement |

### 1.6 Toxicity Data

| Source | Content | Stage use |
|---|---|---|
| DILI annotations (LiverTox, TDC DILI datasets) | clinical DILI risk labels per drug | Stage 5 model training/validation |
| ChEMBL/HERGdb, OMIT | hERG blockade IC50 | Stage 2 safety off-target + Stage 5 |
| DILIsym literature (in-vitro assay data for BSEP inhibition, ETC inhibition, oxidative stress) | mechanistic toxicity assay parameters | Stage 4 liver model parameterization |
| Cmax/IC50 toxicity-correlation literature (241-drug DILI analyses, ROC AUC 0.91-0.96) | validated exposure-based toxicity formula | Stage 5 composite scoring |

### 1.7 Clinical Endpoint / Biomarker Reference Data

| Source | Content | Stage use |
|---|---|---|
| CTCAE (NCI) | severity grading of adverse events (Grade 1-4) | Stage 5 phenotype grading |
| Clinical lab reference ranges | ALT/AST/bilirubin/creatinine/QTc normal ranges | Stage 5 grading |

## 2. Ingestion & Storage Strategy

### 2.1 Local Caches

- All large external datasets are cached under `data/` with a manifest
  (`data/manifest.json`) recording source, version, retrieval date, sha256
  checksum for every vendored file; `python scripts/release.py gates` verifies
  the checksums and blocks a release on any mismatch.
- The original experimental PK reference corpus for the validation benchmarks
  is vendored at `data/benchmarks/published_pk.json` and is the single source
  of truth for the benchmark pass-bands — `validation/benchmarks/` loads from
  it and fails loudly (never silently defaults) if the file is missing or
  malformed. See `data/README.md`.
- High-value derived assets (compiled physiology tables, target panels,
  pathway graphs) are exported to compact formats (Parquet/JSON) for fast
  pipeline loading.

### 2.2 Licensing

- Open/CC/Academic sources are used where possible (TDC, ChEMBL, UniProt, PDB, Reactome, Physiome, WikiPathways, OSP community models).
- Commercial or restricted sources (full DrugBank distribution, LiverTox) are used under academic license; APIs without bulk download are integrated with rate-limit-aware fetch and cache.

### 2.3 Versioning

- A `data.version` field accompanies every model artifact and dataset; pipeline runs record these in the run manifest so results are reproducible.

## 3. Priority Ingestion Order (baseline 2026.9.1)

1. **Physiology database** (OSP/Willmann tables) — unlocks Stage 1 for any human profile.
2. **ChEMBL + DrugBank target panel** — unlocks Stage 2 for approved/reference drugs.
3. **ADMET-AI package + TDC-derived models** — unlocks structure-only parameterization.
4. **KEGG/Reactome pathway graphs** for the top-N safety-critical pathways (MAPK/ERK, PI3K/AKT/mTOR, apoptosis, bile-acid transport, mitochondrial respiration, oxidative stress, cardiac ion channels).
5. **OSP PBPK Model Library baseline compounds** — the validation corpus for Stage 1.

## 4. Required Quality Gates

| Gate | Check |
|---|---|
| Canonicalization | Every SMILES canonicalized to a single InChIKey; duplicate inputs merged |
| Parameter coverage | Fail loudly if PBPK parameterization is missing critical enzymes/transporters instead of silently defaulting |
| Unit consistency | Concentration in mol/L; volumes in L; affinities in molar units; a unit-check layer on the data contract |
| Version pins | ML model weights and dataset snapshots pinned by checksum |