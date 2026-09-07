# Implementation Roadmap & Milestones

The roadmap is phased so that each phase is independently valuable and validates against published data before extending the coupling.

## Phase 0 — Foundations (Weeks 1-2)

**Goal:** reproducible skeleton, data ingestion, tests.

- D1: Repository structure (`src/drugos`, `data/`, `tests/`, `notebooks/`), `pyproject.toml` (version `2026.9.0`, AGPLv3, author Sean Hank), `LICENSE` (AGPLv3), `src/drugos/version.py` (the `YYYY.M.V` scheme), CI lint+test. All commands run under `/opt/anaconda3/envs/drug_os/bin/python`.
- D2: Input parsers (SMILES canonicalization, dosing/regimen parser, human-profile parser) — all inputs arbitrarily configurable with defaults as overridable conveniences only.
- D3: Physiology tables ingestion (OSP/Willmann-derived data) with manifest + checksum.
- **Exit criteria:** a canonicalized molecule + a resolved human profile flow end-to-end through a stub pipeline with a validated data contract (pydantic schemas); project builds and tests pass on `/opt/anaconda3/envs/drug_os/bin/python`.

## Phase 1 — Stage 1 ONLY: Concentration-Time Prediction (Weeks 3-6)

**Goal:** plausible plasma + organ concentration-time profiles from structure + human parameters.

- D4: ADME parameterization module (RDKit descriptors + admet-ai).
- D5: Whole-body PBPK builder (compartment graph, partition coefficients, absorption, metabolism, renal).
- D6: ODE simulation engine + PK metric computation.
- D7: Validation against OSP PBPK Model Library compounds + literature plasma profiles (geometric mean fold error < 2 target).
- **Exit criteria:** benchmark compounds (e.g., midazolam, warfarin, ciprofloxacin if in library set) with predicted-vs-observed profiles within ~2-fold across dose range; Cmax, AUC, t1/2 reported.

## Phase 2 — Stage 2: Target Binding (Weeks 7-9)

**Goal:** occupancy and bound-complex kinetics.

- D8: Target resolver (DrugBank/ChEMBL panel; DTI-ML fallback for novel molecules).
- D9: Kd/kinetics parameterization (literature -> structure-ML -> sequence-ML precedence) with confidence levels.
- D10: Occupancy/TMDD ODE module coupled to Stage 1 tissue exposure.
- **Exit criteria:** for reference drugs with published receptor-occupancy curves, predicted occupancy-time matches published data (e.g., time-at-target within published approximate range).

## Phase 3 — Stage 3: Signaling Pathway (Weeks 10-13)

**Goal:** pathway node activity perturbations.

- D11: Pathway assembler (KEGG/Reactome scaffolding, pruning, ODE generation).
- D12: Baseline steady-state initialization from expression/physiology data.
- D13: Perturbation simulation + fold-change outputs; verify emergent Hill/Emax dose-response vs literature (EC50, Hill slope).
- **Exit criteria:** 2-3 pathways (e.g., RTK/MAPK, PI3K/AKT, apoptosis) built and behaving qualitatively correctly on benchmark perturbations; dose-response curve shapes validated.

## Phase 4 — Stage 4: Organ Function (Weeks 14-18)

**Goal:** liver/cardiovascular/kidney endpoints.

- D14: Liver QST sub-model (bile-acid, mitochondrial, oxidative-stress, cell-death) fed by Stage 1 liver exposure + Stage 2 in-vitro parameters.
- D15: Lumped cardiovascular model + hERG/QT axis.
- D16: Kidney GFR/creatinine sub-model + renal clearance coupling.
- D17: Validation on literature case studies (e.g., acetaminophen overdose DILI trajectory; a hERG-active compound QT signal; a nephrotoxic compound GFR decline).
- **Exit criteria:** ALT/AST/bilirubin trajectories reproduce published simulations qualitatively and within tolerance; QT risk ordering across a hERG-IC50-ranked compound set matches literature.

## Phase 5 — Stage 5: Clinical Phenotype + Toxicity (Weeks 19-21)

**Goal:** consolidated report with biomarkers, grading, toxicity risk.

- D18: Biomarker translator + CTCAE-style grading.
- D19: Composite toxicity scorer (mechanistic + exposure-ratio + ADMET flags) with uncertainty.
- D20: Report generator (JSON contract + HTML/markdown with plots).
- **Exit criteria:** full end-to-end pipeline runs on benchmark compounds; toxicity risk ordering agrees with known clinical safety profiles for a held-out compound set (ROC AUC target ~0.85+ on a benchmark like the TDC DILI set when feasible).

## Phase 6 — Robustness & Validation (Weeks 22-26)

**Goal:** uncertainty, population, sensitivity, release. ✅ **Complete** (`drugos/robustness/`).

- D21: Uncertainty propagation (parameter ensemble) + percentile bands. ✅ `uncertainty.py` — log-normal multiplicative ensemble, 90% bands, verdict counts.
- D22: Virtual patient population simulation + population-level incidence. ✅ `population.py` — anthropometric sampling, risk quantiles, grade incidence.
- D23: Sensitivity analysis (local + Sobol) and driver attribution. ✅ `sensitivity.py` — OAT log-sensitivities + Saltelli first/total indices.
- D24: Prospective-style validation study on held-out compounds; documentation; benchmark suite runbook. ✅ registered as validation cases (clinical-grading L2, risk-ordering L3, robustness D21-D24 L1, prospective-fidelity L1) — `validation/cases/`, report regenerated by `validation/run_validation.py` (G4).
- **Exit criteria:** ✅ reproducible benchmark suite passes (25/25 cases green, GMFE table in `validation/report.md`; includes R-1 R-PK, R-2 corpus calibration, R-3 ORd cardiac cross-check, R-4 BBB/CNS partition, R-5 SBML MAPK lane, R-6 CKD-EPI GFR, R-7 bile-acid cholestasis PBK).

## Phase 7 — Web Service & CNS Panel (shipped)

- ✅ CNS organ panel — passive brain free-exposure (`organ/cns.py`); anchored CNS risk line only when `cns_ic50_nm` is explicit, else 0.20 class prior.
- ✅ Web/API service layer — Flask playground (`src/drugos/web/app.py` + `web/templates/index.html`, deep-purple dark theme), `/api/run` with `with_uncertainty` / `with_population` / `with_sensitivity` toggles, CNS brain-exposure plot, robustness panels (90% bands, cohort incidence, local/Sobol driver tables).

## Phase 7b — Production-Validated Model Upgrades (shipped)

Design requirement "upgrade all models to production-validated, downloadable,
open-source models" — full pipeline mapping and decision records in
`doc/12-production-models.md`. Shipped here:

- **Required-R hard dependency** — literature PK estimator executed by R ≥ 4 on
  every run, agreement gated ≤ 2 % (validation R-1, 25/25 cases).
- **ORd cardiac anchor** — O'Hara-Rudy 2011 human ventricular AP model vendored
  (`data/models/ohara-2011.mmt`, BSD-3 Myokit) and run via
  `drugos.organ.cardiac_ap` (validation R-3: baseline APD90 266 ms, dofetilide
  at measured IC50 ΔAPD90 ≈ 115 ms, monotone block ramp, warfarin control 0).
  CiPA-v1 2017 retune vendored for the multi-ionic-block upgrade (P9).
- **CKD-EPI kidney baseline** — the CKD-EPI 2021 race-free creatinine equation
  (published, no code license) is applied as the baseline GFR from a measured
  serum creatinine whenever one is carried on the profile
  (`drugos.organ.kidney.ckdepi_2021_egfr`, Mosteller BSA scaling; R-6).
- **Bile-acid cholestasis PBK** — the de Bruijn & Rietjens 2024 GCDCA
  bile-acid PBK (paper CC BY 4.0; the companion GitHub repo is CC-BY-NC-ND
  and its R code is intentionally not ported) anchors the liver cholestasis
  axis (`drugos.organ.liver.simulate_gcdca_pbk`, BSEP Ki = IC50/2; R-7).
- **Corpus calibration** — ChEMBL measured hERG + hERG Central corpus vendored
  (R-2 calibration, doc/11 rows 6/8).
- P5/P6/P7/P8/P9 upgrades (DTI ML, DILIsym proprietary / DILI mito-redox SBML,
  Physiome nephron/liver SBML, PK-Sim physiology import, ORd main-path QTw)
  tracked in `doc/12` §1/§4.

## Phase 8 — Extensions (Backlog)

- Biologics/antibody PBPK (FcRn, immunogenicity).
- Drug-drug interaction (DDI) networks.
- Additional organ panels (pulmonary, reproductive/endocrine).
- Monolithic coupled-mode stabilization for small-pathway models.

## Milestone Summary

| Milestone | Phase | Demonstrable artefact |
|---|---|---|
| M1 (wk6) | 1 | Plasma + tissue concentration-time for benchmark drugs |
| M2 (wk9) | 2 | Target-occupancy time profiles |
| M3 (wk13) | 3 | Pathway perturbation profiles |
| M4 (wk18) | 4 | Organ-function trajectories + mechanistic toxicity signals |
| M5 (wk21) | 5 | Full 6-stage end-to-end report for benchmark drug |
| M6 (wk26) | 6 | Validated, reproducible, documented production pipeline |
| M7 (shipped) | 7 | Web playground + CNS brain-exposure panel (robustness toggles) |