# Validation Strategy & Risk Assessment

## 1. Validation Strategy

The literature (Section 6 of `02-literature-review.md`) establishes a three-tier validation ladder that the pipeline will follow.

### 1.1 Tier 1 — Benchmark Compounds

Curate a corpus of drugs with rich, published PK, target, pathway, and safety data. Sources: OSP PBPK Model Library, literature PK/PD profiles, DILI annotations (LiverTox/TDC).

**Recommended benchmark set (balanced across mechanisms):**

| Compound | Primary purpose |
|---|---|
| Midazolam | Hepatic CYP3A4 metabolism; PBPK accuracy |
| Warfarin | High protein binding; narrow-therapeutic index |
| Ciprofloxacin | Renal elimination + CYP inhibition |
| Acetaminophen | DILI (mitochondria, oxidative stress) trajectory |
| A hERG-active drug (e.g., dofetilide) | QT/TdP risk axis |
| Imatinib (or similar TKI) | Target-mediated disposition + CYP3A4 clearance |
| A BSEP inhibitor (e.g., a known cholestatic drug) | Bile-acid transport sub-model |

### 1.2 Tier 2 — Stage-Level Validation

Each stage is validated independently against its own literature before coupling:

| Stage | Validation metric | Target |
|---|---|---|
| Stage 1 PK | Predicted-vs-observed plasma profiles; geometric-mean fold error on Cmax/AUC | GMFE <= 2 across dose range |
| Stage 2 occupancy | Published receptor-occupancy / TMDD profiles | qualitative + point-wise match |
| Stage 3 pathway | Dose-response EC50 / Hill slope vs literature readouts | within ~2-3x |
| Stage 4 organ | Published organ-function trajectories (e.g., DILI case studies) | qualitative + alert-level match — implemented: liver DILI dose-response (therapeutic safe vs overdose Hy's Law, cholestasis anchored to the de Bruijn bile-acid PBK R-7), dofetilide Delta-QTc vs published band (ORd-anchored R-3), CKD-EPI-anchored KDIGO GFR/AKI escalation (R-6) |
| Stage 5 toxicity | Risk ordering vs known clinical safety profiles on held-out set | implemented — dofetilide > warfarin QTc, APAP 20g > 1g DILI, CNS class-prior fallback, driver attribution (L3 ordering case) |
| Stage 5 clinical | Analytic point-matches of the fusion equations (exposure line, prior-only posterior, grade ladder, crossing windows) | implemented (L2 grading case), exact |
| Robustness D21-D24 | Determinism + CI bookkeeping of the three engines + prospective rerun | implemented (L1 D21-D24 + D24 prospective case), exact |
| R-1 R literature-PK cross-check | Required-R bridge (`src/drugos/rbridge/literature_pk.R`, Wagner 1976 / Gibaldi & Perrier 1982 / Greenblatt & Koch-Weser 1975 / Rowland & Tozer 2010) re-derives CL/AUC/t½/tipping-fit from the same simulated curve; assert |CL_R − CL_py|/CL_py ≤ 2 % and `r:agree` on all 5 benchmarks | implemented (L3, `case_r_bridge`) |
| R-2 corpus calibration | Vendored measured corpora cross-check: dofetilide ChEMBL hERG IC50 (geomean ≈ 26 nM, outlier ≥ 10 µM excluded) is high-affinity (< 100 nM); model class prior (2 nM) within 20× of measured (conservative); hERG-Central % inhibition at 1 µM is a long tail (median ≈ 8 %, P99.9 ≈ 32 %) → per-compound hERG override is the honest choice | implemented (L3, `case_corpus_calibration`) |
| R-3 ORd cardiac AP cross-check | Production-validated O'Hara-Rudy 2011 human ventricular model (BSD-3 Myokit `data/models/ohara-2011.mmt`) run under fractional IKr block: baseline APD90 physiological (200–350 ms; measured 266 ms), dofetilide at measured-IC50 concentration prolongs strongly (ΔAPD90 ≥ 30 ms; measured ≈ 115 ms), prolongation monotone in block, warfarin control zero — confirming the calibrated encoder's ordering/magnitude class | implemented (L3, `case_cardiac_ap_ord`) |
| R-4 ADMET-AI BBB→CNS partition | When an ADMET-AI call is present, its BBB_Martins classifier head decides the residual brain:plasma unbound class (P ≥ 0.5 → kpu 1.0, P < 0.5 → 0.2), replacing the fixed default class; benchmark runs without ADMET-AI keep kpu 1.0, brain exposure ratio 0.2 for non-penetrants | implemented (L2, `case_admet_bbb_cns`) |
| R-5 Huang/Levchenko SBML MAPK cascade | Production Stage-3 anchor: vendored BioModels BIOMD0000000009 (CC0) parses to 20 reactions/22 species via python-libsbml; drug-free system relaxes to its fully-activated PP_K steady state; signal (occupancy) inhibits upstream E1 monotonically — full occupancy collapses the readout, response is monotone in occupancy | implemented (L2, `case_sbml_mapk_validation`) |
| R-6 CKD-EPI 2021 race-free GFR baseline | Production kidney baseline (doc/12 row 4b): male 60 y Scr 1.0 → 86.2 mL/min/1.73 m² (CKD-2 band); identical Scr yields *lower* female eGFR (64.5, sex correction); BSA-scaled absolute GFR > indexed (98.2 mL/min at 1.87 m²); a Scr-carrying profile drives `gfr_ml_min` through the equation; profiles without Scr keep the untouched default (benchmark invariance) | implemented (L2, `case_ckdepi_2021`) |
| R-7 de Bruijn & Rietjens bile-acid PBK cholestasis | Production liver cholestasis anchor (doc/12 row 4c, paper CC BY 4.0): at 1 µM free-hepatic exposure a ritonavir-class BSEP inhibitor (IC50 0.2 µM → Ki 0.1 µM) accumulates intrahepatic GCDCA ~10.3× (stress 1.00, cholestatic) while an itraconazole-class inhibitor (IC50 10 mM) stays at 1.0× (benign); fold rises strictly with falling Ki (rank-order of clinical cholestasis reproduced); Ki=IC50/2 pinned; `simulate_liver` cholestasis matches the standalone PBK | implemented (L2, `case_liver_cholestasis_pbk`) |

### 1.3 Tier 3 — Prospective-Style Evaluation

- **Held-out compound protocol:** after model freeze, run the full pipeline on compounds excluded from tuning; rank by predicted risk and compare to known outcomes. D24 registers a reproducibility precondition + a held-out profile/dose stability case.
- **Cross-stage consistency checks:** internal checks such as — an Emax dose-response must emerge from Stage 2/3 without force-fitting; increased liver exposure must translate to increased liver-grade signal (monotonicity where mechanism-consistent); DILI risk must fall when `dili_ic50` rises.

### 1.4 CI-Embeddable Checks

- Unit tests per stage on analytic/limiting-case solutions (e.g., one-compartment bolus analytic vs numeric; zero-dose baseline recovery to steady state).
- Golden-file regression tests on benchmark compounds.
- Data-contract schema tests at every stage boundary.
- **G4 rule:** every new model feature must add a validation case before merge; `python validation/run_validation.py` (28/28 cases green, see `validation/report.md`) regenerates the report and fails the gate on any red case.
- **R is a hard runtime dependency:** a pipeline run without R raises (no silent solver-substitution); each `run_pipeline` streams an `r_verify` block into the JSON contract and the report.

## 2. Risk Assessment

| # | Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|---|
| 1 | **ADME ML errors propagate** — admet-ai estimate of clearance/permeability wrong for a novel scaffold | Med | High | Use literature/in-vitro params when available; uncertainty propagation; report parameter-confidence; benchmark GMFE |
| 2 | **Target set incomplete for novel molecules** — missing off-targets causing toxicity underestimation | Med | High | Enforce a mandatory safety off-target panel; flag coverage gaps in report; DTI-ML as additional net |
| 3 | **Pathway model overreach** — too many nodes/species, unidentifiable parameters | Med | Med | Harden pruning rules (<= ~200 reactions in the baseline); parameterize from published QSP studies; avoid over-claiming identity of unknowns |
| 4 | **Organ-model fidelity gaps** — lumped QST/hemodynamic models miss spatial phenomena (e.g., lobular zonation) | Med | Med | Explicitly document lumped approximations; validate against published DILIsym/Physiome trajectories; keep AOP-based linkage |
| 5 | **Coupling feedback** (organ dysfunction -> clearance) missing in sequential mode | Med | Med | Implement slow-timescale outer loop; offer monolithic mode for small pathways; document when feedback materially changes predictions |
| 6 | **Data licensing/API limits** (DrugBank, LiverTox) | Low | Med | Academic licenses; rate-limit-aware caching; free/open alternates |
| 7 | **Interpretability** — ML components opaque, hard to audit for safety claims | Med | High | Prefer mechanistic over ML wherever mechanism is known; ML only for unknown mappings; every ML output labeled with confidence + driving features; DTA/affinity models with attention-level interpretability |
| 8 | **Validation corpus bias** — benchmarks skewed toward liver/hepatotoxics | Med | Med | Balance benchmark set; prioritise cardiovascular + kidney cases in Phase 4 |
| 9 | **Computational cost of ensemble/population simulation** | Low | Med | Numba jit on ODE RHS; aggressive cruft-free models; parallel sampling via multiprocessing |
| 10 | **Scope creep toward "patient diagnosis"** — overclaiming clinical utility | Low | High | This is a research/modeling tool, not regulated software; clear disclaimers; outputs framed as hypothesis-generating |

## 3. Regulatory & Ethics Posture (baseline 2026.9.1)

- DrugOS is a **mechanistic research and education platform**, not a replacement for clinical judgment or a regulated medical device.
- All outputs are accompanied by uncertainty bands; toxicity calls are advisory and hypothesis-generating.
- No patient-identifiable data used; human profiles are parametric and synthetic (all inputs remain arbitrarily configurable by the user, per `01-project-overview.md` section 3.1).

## 4. Acceptance Criteria for Release (baseline 2026.9.1)

1. Benchmark suite (Tier 1 set) passes stage-level and end-to-end checks with documented fold-errors.
2. ROC ordering on a held-out compound set meets the agreed target or the deviation is documented.
3. Every run is reproducible from a pinned manifest (deps, data checksums, seeds).
4. Documentation (all `doc/*.md`) matches the shipped behavior; a runbook exists for the benchmark suite.