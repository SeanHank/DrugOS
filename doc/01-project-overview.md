# DrugOS Project Overview

## 0. Project Metadata

| Field | Value |
|---|---|
| Project name | DrugOS |
| Author | Sean Hank |
| License | AGPLv3 (GNU Affero General Public License v3.0; the repository ships a `LICENSE` file) |
| Current version | 2026.9.0 |
| Primary language | Python (see `06-technology-stack.md`) |

### 0.1 Versioning Convention (project-wide)

All artifacts of the project (packages, releases, model snapshots, documents) use the format **`YYYY.M.V`**:

- `YYYY` — the year (e.g., 2026)
- `M` — the month, 1-12
- `V` — the revision number *within* that month, starting from **0**

Rules:
- Every new release within the same month increments `V` (0, 1, 2, ...).
- When the month (or year) rolls over, `V` resets to **0** (e.g., 2026.9.0, 2026.9.1, ..., 2026.10.0).
- The baseline release documented here is **2026.9.0**; all earlier "v-major.minor" conventions are replaced by this scheme.

## 1. Vision

DrugOS is an interdisciplinary project combining **computational science** and **pharmaceutical science** to build a multiscale, mechanism-based model of the human body. Given (a) the chemical structure of a drug, (b) the route of administration, (c) the dose, and (d) a complete set of human physiological parameters, the system predicts:

1. **Drug concentration in the body over time** (systemic and per-organ pharmacokinetics)
2. **Target engagement** (binding kinetics and occupancy at molecular targets)
3. **Physiological indicator changes** (biomarkers, organ function metrics)
4. **Potential toxic side effects** (mechanistic and probabilistic toxicity flags)

## 2. Core Causal Chain

The system operationalizes the following translational chain, from input to clinical phenotype:

```
drug ──▶ in-vivo concentration ──▶ target binding ──▶ signaling pathway
    ──▶ organ function ──▶ clinical phenotype
```

This chain is the organizing principle of the whole pipeline. Each arrow maps to a well-established modeling discipline:

| Chain step | Modeling discipline | Primary literature basis |
|---|---|---|
| drug -> in-vivo concentration | Physiologically Based Pharmacokinetics (PBPK), compartmental PK | PK-Sim / Open Systems Pharmacology; GastroPlus ACAT; ADMET-AI |
| concentration -> target binding | Drug-target interaction (DTI) / binding kinetics, receptor occupancy | AI binding-affinity models; DrugBank targets; target-mediated drug disposition (TMDD) |
| target binding -> signaling pathway | Quantitative Systems Pharmacology (QSP), signaling network ODEs | MET/EGFR QSP models; Physiome Project; pathway databases (KEGG, Reactome) |
| signaling pathway -> organ function | Organ-level physiological / quantitative systems toxicology (QST) models | DILIsym (liver); cardiovascular circulation models; neural nephron (kidney) |
| organ function -> clinical phenotype | Biomarker-to-outcome translation, adverse outcome pathways (AOP), clinical endpoints | QSP-to-biomarker translational reviews; AOP framework |

## 3. Inputs and Outputs

### 3.1 Inputs

> **Requirement: all input items are arbitrarily configurable.** Every one of the four input dimensions — chemical structure, route of administration, dose, and complete human parameters — accepts any valid value the user chooses. No input is fixed, hardcoded, or restricted to a preset list. Sensible defaults are provided where helpful, but every individual parameter may be overridden or fully specified by the user.

| Input | Format | Scope of configurability |
|---|---|---|
| Chemical structure | SMILES / InChI / SDF | Any valid molecule: approved drugs, novel compounds, generic scaffolds |
| Route of administration | e.g., IV bolus, IV infusion, oral, SC, IM, transdermal | Any route; user may also provide absorption parameters directly (e.g., SC/IM depot rate) |
| Dose | amount + unit(s); regimen, frequency, duration; optional profile | Any amount, any schedule (single/repeated/continuous/infusion rate) |
| Human parameters | age, sex, height, weight, body composition, organ blood flows, plasma proteins, enzyme abundances, target expression, organ function status, disease state | Each parameter configurable individually or as a full profile; arbitrary combination of values allowed |

All four dimensions combine arbitrarily (any structure x any route x any dose x any human profile). Configurability is enforced by the input layer (see `03-system-architecture.md`, `06-technology-stack.md`): inputs are validated for plausibility and physical consistency, but never rejected because they fall outside a predefined preset.

### 3.2 Outputs

| Output type | Description | Granularity |
|---|---|---|
| Concentration-time profiles | Plasma + per-organ tissue exposure curves (Cmax, AUC, t1/2, etc.) | Continuous time series |
| Target occupancy | Fractional receptor occupancy / bound complex dynamics over time | Continuous time series |
| Signaling activity | Pathway node activity (phosphorylation, transcription-factor activity, pathway flux) relative to baseline | Continuous time series |
| Physiological indicators | Biomarkers (e.g., ALT/AST, bilirubin, heart rate, QTc, blood glucose) and organ function metrics | Continuous time series or clinical-grade flags |
| Toxicity assessment | Mechanism labels (e.g., DILI, hERG blockade, mitochondrial dysfunction, bile-acid transport inhibition) + probabilistic risk scores | Structured report |

## 4. Guiding Principles

1. **Mechanistic-first, data-driven where needed.** Use mechanistic ODE models for PK, target binding, and signaling; use machine learning only where the mechanism is unknown or too expensive to simulate (structure-to-ADME property prediction, binding-affinity estimation).
2. **Multiscale integration via a common mathematical substrate.** All stages are expressed as coupled ODE systems so they can be composed into a single simulation or run as loosely coupled sub-models through a defined data contract (concentration time series as wallet-passing interface).
3. **Transparency and uncertainty.** Every predicted quantity carries a confidence level; parameter uncertainty propagates through the pipeline (population simulation / virtual patient ensembles).
4. **Open-source driven.** Prefer validated open-source stacks (Open Systems Pharmacology suite, RDKit, SciPy ecosystem) over closed commercial tools. The project is licensed under AGPLv3.
5. **Fully configurable inputs.** All inputs (structure, route, dose, human parameters) are arbitrarily configurable; defaults exist only as conveniences and never constrain a user-defined configuration.

## 5. Scope

**In scope (baseline 2026.9.0):**
- Small-molecule drugs with known or inferable targets
- Intravenous, oral, subcutaneous, intramuscular and transdermal (depot) administration routes
- Liver-centric toxicity (DILI) with secondary coverage of cardiovascular (QT) and kidney endpoints
- Adult healthy and common-disease virtual populations

**Explicitly out of scope (baseline 2026.9.0, deferred):**
- Biologics/antibodies (require FcRn, immunogenicity sub-models)
- Multi-drug interaction networks beyond a single co-administered pair
- 3D spatial (finite-element) organ models; the system uses lumped-compartment organ models (the arbitrary input configurability in section 3.1 applies to model *inputs*; the modeled physiology depth may be expanded in later releases)

## 6. Related Work Recognized in This Design

Document 02 (`02-literature-review.md`) anchors every architectural decision in peer-reviewed and community-standard literature, including but not limited to:

- Open Systems Pharmacology Suite (PK-Sim, MoBi) and its validated whole-body PBPK models
- ADMET-AI (Chemprop-RDKit GNN trained on 41 TDC ADMET datasets)
- AI protein-ligand binding affinity models and DTI prediction reviews
- DILIsym quantitative systems toxicology model and DILI mechanistic sub-models
- Physiome Project / Virtual Physiological Human organ models
- PK/PD Emax and target-occupancy (TO) modeling theory