# System Architecture

## 1. Architectural Pattern

DrugOS is a **six-stage pipeline** of coupled models. The pipeline follows the translational chain and is implemented as a **loosely coupled orchestration** of stage modules with a single shared data contract, so each stage can be upgraded or swapped independently.

```
                        ┌───────────────────────────────────────────────┐
                        │                  INPUT LAYER                  │
                        │  SMILES | Route | Dose | Human Parameters     │
                        └───────────────┬───────────────────────────────┘
                                        ▼
   ┌────────────────────────────────────────────────────────────────────┐
   │  STAGE 1  PK / PD Parameterization + PBPK Simulation              │
   │  [structure -> ADME/physchem params] -> [whole-body PBPK ODEs]    │
   │  Output: plasma + organ concentration-time profiles               │
   └───────────────────────────────┬────────────────────────────────────┘
                                    ▼
   ┌────────────────────────────────────────────────────────────────────┐
   │  STAGE 2  Target Identification + Binding Kinetics                │
   │  [targets (DrugBank/sequence AI)] + [micromolar -> nM affinities] │
   │  Output: target occupancy / bound-complex dynamics per tissue     │
   └───────────────────────────────┬────────────────────────────────────┘
                                    ▼
   ┌────────────────────────────────────────────────────────────────────┐
   │  STAGE 3  Signaling Pathway Transduction (QSP)                    │
   │  [pathway graph (KEGG/Reactome)] -> [mass-action/Hill ODE module] │
   │  Output: pathway node activity vs baseline                         │
   └───────────────────────────────┬────────────────────────────────────┘
                                    ▼
   ┌────────────────────────────────────────────────────────────────────┐
   │  STAGE 4  Organ Function Models (QST / physiology)               │
   │  [liver, cardiovascular, kidney organ models]                    │
   │  Output: organ-level state variables (ALT, cardiac output, GFR...)│
   └───────────────────────────────┬────────────────────────────────────┘
                                    ▼
┌────────────────────────────────────────────────────────────────────┐
    │  STAGE 5  Clinical Phenotype + Toxicity Assessment                │
    │  [biomarker translation] + [toxicity scoring (mechanistic+AOP)]   │
    │  Output: physiological indicators + adverse-event risk report     │
    └───────────────────────────────┬────────────────────────────────────┘
                                     ▼
   ┌─────────────────────────────────────────────────────────────────────┐
   │  PHASE 6  Uncertainty & Decision Layers (D21-D24)                   │
   │  [parameter ensemble (D21)] [virtual cohort (D22)]                 │
   │  [global/local sensitivity + drivers (D23)] [prospective rerun (D24)]│
   └───────────────────────────────┬─────────────────────────────────────┘
                                     ▼
   ┌─────────────────────────────────────────────────────────────────────┐
   │  OUTPUT LAYER (report + web playground)                             │
   │  [JSON/markdown/HTML report] [Flask /api/run playground + CNS panel]│
   └─────────────────────────────────────────────────────────────────────┘
```

## 2. Component Breakdown

### 2.1 Input Layer

- **Chemical structure parser** — RDKit converts SMILES/InChI/SDF into a canonical molecule; computes physicochemical descriptors (MW, LogP, LogD, pKa, TPSA, HBD/HBA, rotatable bonds).
- **Human parameter resolver** — resolves a sparse user profile (age, sex, height, weight, optionally organ function, disease state) against the physiology database to fill the full set of PBPK/PBPK-population parameters (organ volumes, blood flows, plasma protein concentrations, enzyme abundances, target expression levels). Uses the OSP physiology basis and Willmann-style allometric/covariate equations.
- **Dose/regimen parser** — single dose or repeated dosing schedule, route (IV bolus, IV infusion with rate, oral with food-state flag).

### 2.2 Stage 1 — PK/PD Parameterization + PBPK

- **ADME predictor** — ADMET-AI GNN (Chemprop-RDKit) + RDKit descriptors. Outputs: solubility, permeability class, plasma-protein binding (fraction unbound), intrinsic clearance (metabolic stability), half-life, bioavailability flags. Predicted signals are consumed by the PBPK layer directly: `fup_plasma` (binding), `cl_int_hep_ml_min_kg` (hepatic clearance), HIA or `bioavailable_Ma` (the oral fraction-absorbed gate, `_admet_fa`), and `logS` (a solubility-limited dissolution cap on the intestinal lumen pool, `AbsorptionParams.solubility_mg_ml`). Per-isoform hepatic abundances (`CYP_ABUNDANCE_PMOL_MG`, Barter et al. 2013) are carried on the physiology and consumed by the per-CYP enzyme-kinetics stage (`PBPKModel.cyp_terms`, abundance-scaled Vmax via `cyp_vmax_mg_h`).
- **PBPK model builder** — constructs the whole-body compartment graph (blood/plasma, liver, kidney, gut, lung, heart, brain, muscle, adipose, skin, remainder) with blood-flow distribution, tissue-partitioning (Rodgers/Rowland or Poulin-Theil method), and first-order absorption for oral or bolus/infusion input for IV.
- **Simulation engine** — ODE solver (SciPy LSODA) over the coupled compartment system; optionally delegate to the OSP engine (PK-Sim on Windows/WINE) or re-implement the same equations in Python.

### 2.3 Stage 2 — Target Identification + Binding

- **Target resolver** — DrugBank (approved drugs) or sequence/structure-based DTI ML (e.g., AttentionDTA, MINDG-style) for novel molecules; plus off-target panels (hERG, key CYP enzymes, BSEP, mitochondrial proteins) as curated safety targets.
- **Binding kinetics estimator** — literature kinetics (kon/koff/Kd) where available; otherwise structure-based ML affinity (CORDIAL/IPBind-class) mapped to Kd; receptor density from tissue-specific expression data.
- **Target occupancy module** — reversible binding ODEs with target turnover (ρ); TMDD treatment when target abundance is comparable to dose.

### 2.4 Stage 3 — Signaling Pathway Transduction

- **Pathway assembler** — retrieves the perturbation-relevant sub-network per target (KEGG/Reactome/PANTHER), prunes to a tractable graph (~10-100 nodes), and emits an ODE model (direct SciPy implementation) using mass-action and Hill kinetics.
- **Baseline state** — steady-state initialization of the pathway from physiological/expression data, modeling the chosen human profile (healthy vs. disease state).

### 2.5 Stage 4 — Organ Function Models

- **Liver model (QST)** — bile-acid transport inhibition, mitochondrial/ETC dysfunction, oxidative stress, hepatocyte death; outputs ALT/AST/bilirubin trajectories and DILI risk (DILIsym-style sub-model).
- **Cardiovascular model** — lumped circulation (Physiome-style circulation + ventricular interaction) to translate chronotropy/inotropy/QT effects into HR, BP, and QTc.
- **Kidney model** — nephron-level model for GFR and serum creatinine to reflect excretion and nephrotoxicity.
- Organ models exchange state with the PBPK layer (organ exposure drives toxicity; organ dysfunction feeds back into clearance).

### 2.6 Stage 5 — Clinical Phenotype + Toxicity

- **Biomarker translator** — maps organ-model outputs to clinical-grade biomarkers against reference ranges with CTCAE-like grading (Grade 1-4).
- **Toxicity scorer** — combines: (i) mechanism-based QST outputs, (ii) exposure-ratio scores (Cmax / in-vitro IC50 ratio with ROC AUC ~0.9+), (iii) ADMET-AI structural flags (hERG, AMES, hepatotox). Emits per-endpoint risk probabilities with confidence.
- **CNS brain free exposure** — passive blood-brain partitioning; when the molecule is anchored with an explicit CNS IC50 (`cns_ic50_nm`), a brain-free-exposure/IC50 driven CNS risk line is fused; unanchored compounds keep the 0.20 class prior.
- **Report generator** — structured JSON + human-readable markdown/HTML report with plots.

### 2.7 Phase 6 — Uncertainty & Decision Layers (D21-D24)

- **D21 Parameter ensemble** — multiplicative log-normal perturbation of PK/potency scalars (fixed seed `EnsembleConfig`) → 90% percentile bands for plasma/organ curves and endpoint-risk quantiles with per-member verdict counts.
- **D22 Virtual cohort** — anthropometric sampling (sex/age/BMI, deterministic seed) runs the full pipeline per individual at the same dose → population risk quantiles and grade ≥ 1 / ≥ 2 incidence.
- **D23 Sensitivity & drivers** — one-at-a-time log-sensitivities (`local_sensitivity`) plus Saltelli first/total Sobol indices (`run_sobol_sensitivity`) over the D21 key set, ranking the drivers of each endpoint.
- **D24 Prospective rerun** — reproducibility runbook (fixed-seed determinism) and held-out profile/dose re-validation, registered as validation cases (doc/08).

## 3. Data Contract (Interface Between Stages)

The pipeline's stages are coupled by a **time-indexed data contract** so that any stage can be replaced.  In the shipping implementation the stages exchange a canonical `RunSpec -> RunResult` contract (`drugos/pipeline.py`) whose playable projection is the `/api/run` payload:

```json
{
  "time": [0.0, 0.1, ...],                 // hours
  "plasma_total": [...], "plasma_free": [...],
  "tissues": { "liver": {...}, "kidney": {...}, "heart": {...}, "brain": {...} },
  "occupancy": { "TARGET1": {...}, "TARGET2": {...} },
  "pathways": { "MAPK": {...}, "PI3K_AKT": {...} },
  "organ": { "liver_alt": [...], "gfr": [...], "qtc": [...] , "brain_free_nm": [...] , ...},
  "clinical": { "phenotype": {...}, "toxicity": {...} },
  "robustness": {
    "uncertainty": { "points_q5_q50_q95": {...}, "band_90": {...}, "verdict_counts": {...} },
    "population": { "risk_q5_q50_q95": {...}, "incidence_grade_ge_1": {...}, ... },
    "sensitivity": { "local": {...}, "sobol": { "first": {...}, "total": {...} } }
  }
}
```

All stage modules implement `simulate(input_contract, human_params, drug_params) -> output_contract`.

## 4. Two Execution Modes

1. **Monolithic coupled mode** — a single stiff ODE system integrating PBPK + occupancy + pathway + organ stages (feasible when pathway graphs are small, ~10-30 nodes). All dynamics solved simultaneously; strongest mechanistic coupling.
2. **Sequential mode** — run Stage 1 first, then feed concentration time series into staged simulations. Faster, more modular, allows sub-model swapping; less faithful coupling (e.g., feedback of organ dysfunction on clearance is approximated on a slower timescale).

## 5. Configuration & Reproducibility

- Every run records a **parameter manifest**: drug parameters, human parameters, model versions, ML model versions/hashes, solver settings.
- Configurations are version-controlled and reproducible (deterministic seeds for ML, fixed solver tolerances).

## 6. Module Interface Definitions

```
drugos/
  inputs/
    parse_structure.py         # SMILES/InChI/SDF -> RDKit mol + descriptors
    resolve_human.py           # profile -> full physiology parameter set
    parse_dosing.py            # route/amount/regimen -> input events
  pk/
    admet.py                   # ADMET-AI + descriptor wrappers
    pbpk_build.py              # compartment graph + partition coefficients
    simulate.py                # ODE solve -> concentration-time
  target/
    targets.py                 # DrugBank / DTI-ML target sets
    affinity.py                # literature or ML Kd / kinetics
    occupancy.py               # TMDD / turnover occupancy ODEs
  pathway/
    assemble.py                # pathway DB -> ODE model (SciPy)
    simulate.py                # pathway ODE solve
  organ/
    liver.py                   # QST liver sub-model
    cardiovascular.py          # lumped circulation + QTc
    kidney.py                  # GFR / creatinine model
  clinical/
    biomarkers.py              # organ output -> clinical biomarker + grading
    toxicity.py                # composite toxicity scoring
  organ/
    cns.py                     # brain free exposure (passive barrier)
  robustness/
    uncertainty.py             # D21 parameter ensemble + 90% bands
    population.py              # D22 virtual cohort + incidence
    sensitivity.py             # D23 OAT/Sobol indices + drivers
  pipeline.py                  # RunSpec/RunResult contract + orchestration
  report/
    render.py                  # JSON + markdown/HTML report
  cli.py                       # run / benchmarks / study / serve subcommands
  web/                         # Flask playground, shipped inside the wheel
    app.py                     # /api/run (robustness flags), index route
    templates/index.html       # deep-purple dark theme Playground
  py.typed
  __init__.py
```