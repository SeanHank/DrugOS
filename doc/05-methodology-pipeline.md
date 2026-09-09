# Methodology: End-to-End Prediction Pipeline

This is the core design document. It specifies the mathematical formulation, sub-model choices, and implementation plan for each pipeline stage. Every stage is grounded in the literature reviewed in `02-literature-review.md`.

---

## Stage 1 — Drug -> In-Vivo Concentration (PK)

### 1.1 Goal

Given a structure, route, dose, and human profile, produce `plasma` and per-organ concentration-time profiles.

### 1.2 Step 1A: Structure -> ADME/PhysChem Parameterization

Inputs: canonical SMILES. Outputs: the parameter set required by PBPK.

| Parameter | Derivation |
|---|---|
| Molecular weight, LogP, LogD(pH), pKa, TPSA, HBD/HBA | RDKit descriptors |
| Aqueous solubility (logS) | ADMET-AI / RDKit-ESOL |
| Fraction unbound in plasma (fu) | ADMET-AI (plasma protein binding) |
| Intrinsic clearance (Clint) | ADMET-AI metabolic stability; refined via CYP phenotyping (CYP3A4/2D6/2C9/2C19/1A2) from ML or literature |
| Permeability class / fraction absorbed (Fa) | ADMET-AI bioavailability + PAMPA/Caco-2 models; ACAT absorption for oral |
| Efflux/uptake transporter flags (P-gp, OATP1B1/1B3, BCRP, BSEP) | ChEMBL/ADMET-AI transporter inhibition + literature |
| T1/2, Vss priors | ADMET-AI half-life / volume; overridden by PBPK-computed values |

**Accuracy note.** ADMET-AI ranks first on the TDC leaderboard and is suitable for first-pass parameterization; when measured in-vitro ADME data or validated PBPK models exist (e.g., OSP Model Library), they take precedence.

### 1.3 Step 1B: Human Parameter Resolution

The user provides a sparse profile: age, sex, height, weight, (optional) organ status, disease state. The resolver computes:

- Organ/tissue volumes (vascular, interstitial), blood-flow fractions, cardiac output -> allometry/covariate equations (Willmann et al. 2007; ICBP reference tables).
- Plasma protein concentrations (albumin, AGP) — age/sex adjusted.
- Enzyme abundances (CYP450 per gram liver) — population references with covariance.
- Target expression levels per tissue (from proteomics/expression databases where available).

Output: a full PBPK parameter vector for this individual (used both for point simulation and as the center of a virtual-population distribution for uncertainty quantification).

### 1.4 Step 1C: Whole-Body PBPK Model

**Structure.** A system of compartments including at minimum: lung, heart, liver (with gut), kidney, brain, muscle, adipose, skin, bone, remainder (spleen, pancreas, fat-free/lean tissue lumping). Each compartment consists of sub-spaces where appropriate: vascular, interstitial, and cellular (intracellular) for small molecules, with passive permeability/partition exchange plus active processes.

**Governing equation (per compartment, unbound mass balance).**

```
dA_comp/dt = Q_comp * (C_in - C_out) + net passive flux + active flux - metabolism - excretion
```

where `C_in` is the inflow (arterial) concentration, `C_out = C_free_tissue / Kpu` (partition-corrected), `Q_comp` is organ blood flow obtained from the resolved physiology. Tissue partition coefficients use the Rodgers & Rowland (or Poulin-Theil) methods parameterized by LogP/pKa/LogD, tissue composition, and fu.

**Absorption (oral).** ACAT-style first-order transit: stomach -> small-intestine segments -> colon, with pH-dependent solubility/dissolution, permeability-driven absorption into the portal vein, first-pass hepatic extraction.

**Absorption (SC/IM/transdermal).** Subcutaneous, intramuscular and
transdermal doses enter a first-order `depot` compartment that feeds venous
blood with rate constant `k_depot_absorption` (default ~0.15/h) and routed
bioavailability `F`; the unabsorbed fraction (`1 - F`) is routed to feces, so
mass is conserved. Users may override the absorption rate per run
(e.g. `--sc-im-ka`), giving these routes a blunted, delayed Cmax relative to
the equivalent IV bolus (see §1.6 for the transdermal baseline note).

**Clearance.** Hepatic metabolism via MM/Hill kinetics on unbound liver concentration for each CYP (with abundance-scaled Vmax and literature/ML Km); renal excretion via glomerular filtration of unbound drug (fu * GFR * fraction unbound-driven), plus tubular secretion terms if transporter data exist.

**Tubular transport (renal).** Beyond glomerular filtration of the free
fraction (`fup*GFR`), an active tubular-secretion term (`cl_sec`, OAT/OCT
transporters) adds a first-order flux on the unbound kidney drug; it is off
by default (pure filtration, `RunSpec.cl_sec_l_h = 0`).

**Biliary excretion & enterohepatic recirculation.** A bile pool
(`PBPKModel` ``bile`` compartment) receives free parent from the liver at a
biliary-secretion rate `cl_bil` and empties by `k_bile_emptying` into the
small-intestine lumen; the emptied mass is then subject to the ordinary SI
reabsorption/transit kinetics, so a fraction recirculates via the portal vein
and the rest reaches the feces sink. Off by default (`cl_bil = 0`).

**First-pass gut-wall extraction.** Enterocyte (intestinal) first-pass
metabolism removes a fraction `gut_extraction_eg` of the SI-absorbed flux
before the portal blood, so oral `F = Fa·(1−Eh)·(1−Eg)` (§1.6). Off by
default (`Eg = 0`).

**Saturable hepatic clearance.** Hepatic metabolism can be switched from the
linear `cl_hep` term to a Michaelis-Menten `Vmax/Km` form on the unbound
liver concentration (`RunSpec.hepatic_vmax_mg_h`/`hepatic_km_mg_l`, set
together); below Km the slope equals `Vmax/Km`, so the linear default is the
low-dose limit. Off by default (linear clearance preserved).

**TMDD coupling (forward reference).** If target engagement consumes significant drug (high-affinity, high-abundance target), Stage 2 occupancy fluxes (drug + receptor -> complex -> internalization/clearance) feed back into the tissue mass balances. In monolithic mode this is native; in sequential mode it is approximated iteratively.

**PK metrics.** Cmax, tmax, AUCinf, t1/2, Vss, Cl are computed from the simulated profiles and reported. **Systemic bioavailability `F` is also reported per route** (doc/12 row PK-f): IV = 1; depot routes (SC/IM/transdermal) = the depot bioavailability; oral = `Fa·(1 − Eh)·(1 − Eg)` where `Fa = 1 − feces_fraction` (the colon-transit feces sink), `Eh = CL_h/(Q_h + CL_h)` is the well-stirred first-pass hepatic extraction and `Eg` is the first-pass intestinal extraction (0 by default) — the model already routes absorbed oral drug into the liver compartment, so this is the first-pass-corrected value, `PkMetrics.f_abs`. `Fa` itself follows the permeability-gated absorption constant tuned to `RunSpec.fa` (§1.6), so reported `F` reflects both the drug's absorbability and its extraction.

### 1.5 Implementation

- Formulate as sparse coupled-ODE system; solve with `scipy.integrate.solve_ivp` (LSODA) with tight tolerances for stiffness.
- Optional integration: run the equivalent model in the OSP suite for model-library compounds as a cross-validation; keep equations transparent in Python.
- Population mode: sample virtual individuals (n>=100) from physiological covariance, producing 5-95th percentile bands.

### 1.6 Baseline Implementation Scope (2026.9.0)

The shipped engine implements a faithful *baseline* of the mechanics above.
Baseline simplifications, each recorded in doc/01 §5 (out of scope), are:

- **Clearance** defaults to a lumped single hepatic-metabolic `cl_hep` (from
  the resolved physiology and structure-derived partition) plus renal GFR of
  unbound drug.  Extra realism terms are implemented but **off by default** —
  so every shipped validation run keeps the linear validated baseline:
  active tubular secretion (`RunSpec.cl_sec_l_h`), a saturable Michaelis-
  Menten hepatic term (`hepatic_vmax_mg_h`/`hepatic_km_mg_l`, linear low-dose
  limit), biliary drug excretion with enterohepatic recirculation
  (`cl_bil_l_h`/`bile_emptying_1h`, §1.4) and first-pass gut-wall extraction
  (`gut_extraction_eg`).  Per-CYP MM/Hill richness (consuming the
  `CYP_ABUNDANCE_PMOL_MG` table) is still deferred.
- **Oral absorption** is ACAT-lite: first-order stomach -> small-intestine ->
  colon transit with pH-dependent dissolution; no multi-segment
  small-intestinal dissolution resolution. **The small-intestine absorption
  rate is permeability/Fa-gated**: `_absorption_params` tunes
  `k_si_absorption` to the compound's fraction-absorbed target
  (`absorption_rate_from_fa`, `RunSpec.fa`) — from the benchmark published
  `f_abs` band (or `fa_override`), from the predicted HIA head for ADMET-AI
  SMILES runs (`_admet_fa`), and from the default endogenous rate otherwise —
  so a low-permeability molecule peaks lower/later and loses more to the
  colon-feces sink while a fast one absorbs in the proximal gut.  **Solubility-
  limited dissolution** is gated on the ADMET-AI `logS` for novel molecules:
  the dissolved small-intestine pool is capped at `S * V_gi`
  (`AbsorptionParams.solubility_mg_ml`), so an excess low-solubility dose
  spills undissolved mass forward to the colon/feces sink.  Benchmark runs
  (no ADMET call) are unaffected.
- **SC/IM/transdermal absorption** is a first-order depot (mass-conserving,
  user-tunable `k_depot_absorption` / `F`); no microcirculation/lymphatic
  resolution in the subcutaneous tissue, and transdermal uses the generic
  depot defaults (a skin-permeation multi-layer model is deferred).
- **ADMET -> PK auto-wiring** is first-class: `pipeline.spec_from_admet` builds
  a full `RunSpec` for a novel molecule from `fup_plasma`,
  `cl_int_hep_ml_min_kg` (scaled by body weight), GFR-derived renal clearance
  and the HIA absorption gate — shared verbatim by the CLI
  (`cli._smiles_spec`) and direct `run_pipeline` consumers, so the chain
  "new chemistry -> in-vivo concentration" has exactly one implementation.
- **TMDD** is approximated iteratively via the opt-in `feedback_loop` driver
  (Stage 4 -> Stage 1 PK rescaling, §4.5), not through Stage-2 mass-balance
  flux coupling.

---

## Stage 2 — Concentration -> Target Binding

### 2.1 Goal

Given free drug concentrations at target sites, predict fractional target occupancy and bound-complex kinetics over time, plus off-target engagement.

### 2.2 Step 2A: Target Identification

- **Approved/reference drugs:** DrugBank targets for the given drug (if present in DB) — primary targets + known off-targets (hERG, CYP enzymes, BSEP, transporters, nuclear receptors).
- **Novel molecules:** DTI/DTA machine learning (sequence-based AttentionDTA/MINDG-style or structure-based when a structure is available) ranks candidate targets from the safety-critical panel and mechanism hypotheses.
- **Recommended default off-target safety panel** (baseline 2026.9.0): hERG (Kv11.1), CYP3A4/2D6/2C9 (inhibition), BSEP, MRP3/MRP4, OATP1B1, P-glycoprotein, mitochondrial complex I/II/III/IV/MCT, and the glucocorticoid/sex-hormone receptors (endocrine effects).
- **Baseline scope (2026.9.0).** The shipped engine resolves targets for built-in
  benchmark compounds only and defaults to the off-target safety panel with
  class-typical IC50 priors; the DrugBank resolution kernel and the DTI/DTA-ML
  target ranker for novel molecules are deferred (doc/01 §5).

### 2.3 Step 2B: Binding Kinetics Parameterization

Priority order for Kd / kon / koff of drug-target pairs:

1. **Measured kinetics/affinity** from ChEMBL, BindingDB, or literature (kon, koff, Kd) — highest priority.
2. **Structure-based AI affinity** (CORDIAL/IPBind-class) mapped from predicted binding energy to Kd when a 3D complex is obtainable (docking or AlphaFold target + RDKit pose via docking).
3. **Sequence/descriptor-based DTA** (AttentionDTA-class) otherwise.
4. **Fallback default** for safety targets: use class-typical IC50/EC50 medians for pharmacology-informed priors and flag low confidence.

> **ADMET-AI hERG sieve (2026.9.0).** For novel molecules scored by ADMET-AI
> without a measured hERG potency, the predicted hERG head acts as a *sieve*:
> a predicted non-blocker (`P < 0.5`) demotes the dofetilide-class panel prior
> (2 nM) to a weak-kD floor (1 mM), a predicted blocker keeps the panel prior.
> The effective KD (sieve result, or a per-compound `qt_ic50_nm` override)
> replaces the hERG site in the safety panel once, so Stage-2 occupancy, the
> QT drive and the exposure anchor all stay mutually consistent
> (`_effective_herg_kd_nm` / `_herg_sieved_panel`). Benchmarks (no ADMET-AI
> call) fall back to the panel prior unchanged.

Target abundance (Rtot per tissue) comes from tissue expression/quantitative proteomics tables; store per tissue.

### 2.4 Step 2C: Occupancy Model

Reversible binding with target turnover (Daryaee & Tonge 2019; classic TO model):

```
dR/dt   = ksyn - ρ * R - kon * D_free * R + koff * DR
dD_free/dt = (delivery from tissue PK) - (clearance) - kon*D_free*R + koff*DR
dDR/dt  = kon * D_free * R - koff * DR - kint * DR
```

- `ρ` = target turnover rate; `ksyn = ρ * R0` sets baseline; `kint` = complex internalization rate.
- Fractional occupancy = DR / Rtot.
- **TMDD**: when [drug]/[Rtot] is not large, apply the full TMDD equations or the quasi-steady-state approximation; this couples back into tissue drug mass balances (Stage 1).

### 2.5 Outputs

- Per-target, per-tissue occupancy-time profiles.
- Time-at-target (occupancy-integral) used by Stage 3 as the "exposure signal".

---

## Stage 3 — Target Binding -> Signaling Pathway (QSP)

### 3.1 Goal

Translate occupancy into downstream pathway activity and relative-to-baseline physiological signals.

### 3.2 Step 3A: Pathway Graph Assembly

For each engaged target, assemble the relevant sub-network from KEGG / Reactome / PANTHER (e.g., RTK -> RAS -> RAF -> MEK -> ERK; PI3K -> AKT -> mTOR; apoptosis/caspase; bile-acid/enterohepatic signaling; calcium/ion-channel signaling; mitochondrial respiration coupling). Prune to a **tractable graph** (20-200 reactions) that covers:
- The intended therapeutic route (efficacy pathway)
- The principal toxicity-relevant routes (see Stage 4 organ panels)

> **Baseline scope (2026.9.0).** The shipped kernel is the Huang/Levchenko
> MAPK cascade (R-5): 20 irreversible mass-action reactions over 22 species
> (RTK -> RAS -> RAF -> MEK -> ERK, dual-phosphorylation readout), compiled
> from the in-repo SBML and executed by ``simulate_sbml_pathway``.  A
> ~12-reaction / ~10-species in-repo DSL compiler is retained as the fast
> path, but the pipeline default is the SBML cascade (the 20-200-reaction
> target).  KEGG/Reactome/PANTHER ingestion and larger toxicity-route graphs
> are deferred (doc/01 §5); the compiler contract and ODE export support
> them unchanged.

### 3.3 Step 3B: ODE Generation (direct SciPy implementation)

Encode reactions with:
- **Mass-action kinetics** for receptor-ligand binding, phosphorylation, PP2A/p65 dephosphorylation, degradation, complex formation.
- **Hill-type laws** for transcription-factor activation, gene-regulatory and feed-back terms, and for enzyme-catalyzed steps (Michaelis-Menten).
- **Conservation/steady-state** initialization from physiological baselines (basal receptor density, kinase levels, substrate pools).

```
dx_i/dt = Σ reactions ± (Hill activation/inhibition) - degradation
```

Model geometry mirrors the MET-QSP precedent (~100 species / ~70 ODEs) and is exported to the shared ODE contract.

### 3.4 Step 3C: Perturbation Simulation

- Input: occupancy-time profile from Stage 2 (drug as perturbation term in receptor mass balance).
- Outputs: activity-time traces of key nodes (e.g., pERK fraction, pAKT, p53, RIPK1/RIPK3, mitochondrial membrane potential proxy, bile-acid pool interference index).
- Baseline normalization: report fold-change relative to drug-free steady state.

### 3.5 Validation Points

- Dose-response shapes (Hill/Emax) must emerge naturally: verify EC50 and Hill slope against literature where available (amplification yields EC50 < Kd).
- Known pathway readouts (e.g., phospho-protein biomarker data) compared if published.

---

## Stage 4 — Signaling Pathway -> Organ Function (QST / Physiology)

### 4.1 Goal

Map pathway/toxicity signals to organ-level functional endpoints. The baseline (2026.9.0) covers liver, cardiovascular, and kidney.

> **CYP450 abundance table (Stage 4 upstream data):** `HumanPhysiology` exposes per-isoform hepatic microsomal abundances (`CYP_ABUNDANCE_PMOL_MG`, Barter et al. 2013 central immunoquantified values for CYP1A2/2A6/2B6/2C8/2C9/2C19/2D6/2E1/3A4) and the derived per-isoform liver content (`hepatic_cyp_content_nmol`). This is a *data* table that a future enzyme-kinetics stage consumes; today's lumped `cl_hep` clearances do not read it, so predictions are unchanged.

### 4.2 Liver (DILI) — QST Sub-Model (DILIsym-style)

Mechanistic sub-models:

1. **Bile-acid transport inhibition** — BSEP/NTCP/MRP3/MRP4 inhibition from Stage 2 safety-target data (IC50); enterocyte/hepatocyte bile-acid pools; basolateral and canalicular efflux; predicts cholestasis and biliary biomarkers. **The cholestasis axis is anchored to the production-validated GCDCA bile-acid PBK of de Bruijn & Rietjens 2024** (paper CC BY 4.0; validated vs clinical cholestasis incidence of ~18 marketed drugs, R-7, doc/12 row 4c): free-hepatic drug competitively inhibits BSEP efflux (`Km_app = Km_BSEP·(1 + C_free/Ki)`, `Ki = IC50/2`) and the intrahepatic pool fold-ratio above the 1.5× risk threshold drives cholestasis.
2. **Mitochondrial dysfunction** — ETC complex inhibition (in-vitro IC50), ATP production shortfall, decrease in mitochondrial membrane potential; with adaptive mitogenesis term.
3. **Oxidative stress** — reactive-oxygen-species generation vs. glutathione buffer; GSH depletion; protein/lipid damage.
4. **Hepatocyte death** — apoptosis (caspase-driven via TNF/ligand pathway) + necrosis (ATP/oxidative-threshold-driven); includes regeneration dynamics.
5. **Immune-mediated component** (stub in baseline; deferred).

**Inputs.** PBPK liver tissue exposure (`C_liver(t)`, free), plus in-vitro IC50/assay parameters (BSEP inhibition, ETC inhibition, oxidative stress) from literature/ChEMBL/in-vitro consortia data.
**Outputs.** ALT/AST/total-bilirubin serum trajectories; dead-cell fraction; DILI-grade classification (ALT > 3x ULN; Hy's Law criteria). Follows the fezolinetant DILIsym precedent where PBPK exposure + in-vitro toxicity parameters predict liver signal.

> **Implemented (2026.9.0):** `src/drugos/organ/liver.py` — cholestasis (de Bruijn & Rietjens bile-acid PBK, R-7), mitochondrial (ETC, adaptive mitogenesis), redox/GSH, hepatocyte death (sigmoid kill + regeneration) ODEs; `LiverParams` carries compound-specific IC50s, `liver_params_from_panel()` wires the Stage-2 safety panel. **Pathway→organ closure:** the Stage-3 ERK/MAPK readout fold-change is the hepatocyte `proliferation_signal` that scales `regeneration_1h` via `regeneration_scale()` (damped `regen_pathway_gain=0.5`, clamped to `[regen_pathway_min=0.5, regen_pathway_max=1.5]`, so an off-target pathway suppression can attenuate but never ablate regeneration); the reported total-bilirubin rise is capped at `bile_rise_max_fold`×ULN while the un-clamped excess stays in the mechanical `stress`. **Immune-mediated axis seam:** `combined_stress()` carries an explicit fourth (immune) axis weighted by `LiverParams.immune_weight` (default 0, so every shipped default run keeps the validated cholestasis/ATP/GSH composition exactly); it is the stub endpoint for a future immune-DILI driver. Validated: acetaminophen overdose reproduces ALT > 3x ULN with Hy's Law while therapeutic dosing stays grade 0/1, and the regeneration-coupling + bilirubin ceiling are pinned by L2 cases.

### 4.3 Cardiovascular — Lumped Circulation + Electrical Axis

- **Lumped circulation model** (Physiome-style minimal hemodynamic system with ventricular interaction and valve dynamics): cardiac output, stroke volume, mean arterial pressure, central venous pressure as state variables.
- **hERG/QT axis**: from hERG blockade fraction (Stage 2 off-target), estimate IKr reduction and a QTc-prolongation model (multiplicative / Emax on hERG channel current; literature QTc-hERG relationships). Maps to torsades-de-pointes risk band.
- **Inotropy/chronotropy modulators**: beta/catecholamine pathway effects feed into contractility and rate.

> **Implemented (2026.9.0):** `src/drugos/organ/cardiac.py` — `predict_qtc()` Emax hERG->IKr->QTc axis with TdP banding (450/480/500 ms) and a two-compartment Windkessel (`simulate_hemodynamics()`) for MAP/CVP/CO/SV. The hERG blockade fraction is driven by the **PBPK cardiac (heart) free exposure** (`unbound_tissues["heart"]`), not the hepatic one. **Inotropy/chronotropy coupling (pathway→organ):** the Stage-3 ERK/MAPK amplification ratio gates `CardiacParams.inotropy = chronotropy` via `_cardiac_tone_scale()` (damped `gain=0.2`, capped at 1.3×); the ERK proxy encodes only the *stimulatory* branch (ERK1/2 is downstream of beta-adrenergic E-C coupling), so a baseline or suppressed readout keeps tone exactly 1.0 and the loop `co_fraction` stays neutral for benign drugs, while a genuinely amplified readout raises CO/HR. Sympathetic suppression is explicitly deferred (the ERK readout does not encode it). Validated: 0.5 mg dofetilide peak Delta-QTc 20 ms lies in the published prolongation band vs a negligible-hERG control (L3 case).

### 4.4 Kidney — Nephron + GFR Model

- **Nephron-level model** (Physiome neural-nephron lineage, catalogued P7): glomerular filtration, tubular reabsorption/secretion; clearance coupling with Stage 1 renal elimination.
- **Nephrotoxicity endpoints**: acute kidney injury proxies (GFR decline, tubular injury biomarker (KIM-1 heteromer in practice)) — the baseline keeps serum creatinine + GFR from renal function sub-model.

> **Implemented (2026.9.0):** `src/drugos/organ/kidney.py` — nephron injury sigmoid drives a floored GFR; serum creatinine from the closed-form balance Scr = P/GFR; KDIGO AKI stage from Scr ratio/GFR drop. **The baseline GFR is anchored to the CKD-EPI 2021 race-free creatinine equation** (R-6, doc/12 row 4b) whenever a measured serum creatinine is carried on the profile (`ckdepi_2021_egfr`, BSA-scaled via Mosteller). **KIM-1 tubular biomarker:** the proximal-tubule injury signal is translated to a graded urinary KIM-1 xUNL row (`summarize_clinical_kidney`: `1 + 8·injury`) with CTCAE-style thresholds (1.5/3/5/10 xUNL), so nephrotoxicity is also surfaced by a tubule-specific marker, not only by GFR/creatinine. Validated: exact Scr=P/GFR at zero exposure, monotonic KDIGO escalation, and CKD-EPI reference points / pipeline wiring (L2 cases).

### 4.5 Feedback to PK

Organ dysfunction feeds back to PK: reduced hepatic clearance (liver), reduced GFR (kidney), reduced cardiac output (perfusion-limited distribution). Implement in monolithic mode directly; in sequential mode via a slow-timescale outer loop (recompute PBPK with updated organ parameters once or twice).

> **Implemented (2026.9.0):** `src/drugos/organ/feedback.py` — `feedback_from_results()` maps dead fraction / GFR/CO factors to `OrganFeedback`; `apply_pk_scaling()` returns a PBPK model copy with scaled hepatic/renal clearances and perfusion (hepatic/renal/CO), leaving the original untouched. The cardiac-output factor is **derived from the simulated hemodynamics**: `co_fraction = CO_simulated / CO_reference` (from `CardiacResult.co_l_min`), so a reduced-cardiac-output physiology scales perfusion-limited distribution on the re-run (neutral 1.0 with intact cardiac panels).

---

## Stage 5 — Organ Function -> Clinical Phenotype + Toxicity

### 5.1 Biomarker Translation & Grading

- Map organ outputs to clinical-grade biomarkers with units and reference ranges (ALT, AST, bilirubin, GFR, creatinine, **urinary KIM-1**, HR, BP, QTc).
- Grade each (CTCAE-style Grade 0-4) with severity thresholds; produce time-to-onset and duration.

> **Kidney row set.** `summarize_clinical_kidney` grades GFR (mL/min), serum
> creatinine (xULN) and urinary KIM-1 (xUNL, derived from the nephron injury
> fraction) so AKI is surfaced both by the functional markers and by a
> proximal-tubule injury biomarker (`BIOMARKERS["KIM_1"]`).

> **Spec-pinned grading.** Biomarker grading threads the *spec-driven*
> constants (e.g. `cns_ic50_nm`, the effective hERG KD) into the organ
> summary instead of re-instantiating model defaults (`_grade_organ` carries
> the `CnsParams` used by the simulation), so the reported margins always
> match the simulated exposure-ratio anchors.

### 5.2 Composite Toxicity Scoring

Three independent evidence lines combined with uncertainty:

1. **Mechanistic QST outputs** (Stage 4) — mechanistic DILI/QT/AKI signals, driven by PBPK exposure + in-vitro parameters.
2. **Exposure-ratio score** — validated formula from literature: `Cmax_unbound / in-vitro toxicity_IC50` (e.g., lowest lethal/toxicity concentration across assays); flag if ratio exceeds established thresholds (ROC AUC ~0.91-0.96 in 241-drug DILI analyses).
3. **Structural ADMET flags** — ADMET-AI predictions (hERG blockade probability, AMES, hepatotoxicity probability, DILI probability) as independent structural-molecule evidence.

**Fusion.** Shipped implementation (`drugos/clinical/toxicity.py`): each endpoint derives a probability from its evidence lines —
`P_line = sigmoid(slope * (log10(Cmax_free/IC50) - ratio_center))` for the exposure line, `RULES[].grade_probs[grade]` for the mechanistic line — and the lines are combined in log-odds against a per-endpoint class prior,
`logit(P) = logit(prior) + Σ w_i (logit(P_i) - logit(base_i))`, with a Beta CI shrunk by the sum of evidence weights. Per-endpoint priors DILI 0.25 / QT 0.20 / AKI 0.22 / CNS 0.20; weights MECH 1.2 / EXP 1.0 / STRUCT 0.7.

**Four endpoints (Stage 4 + CNS).**

| Endpoint | Evidence enabled by | Note |
|---|---|---|
| DILI | liver organ + `dili_ic50_nm` (in-vitro override or class prior via BSEP) | acetaminophen class prior gives the DILI anchor |
| QT / TdP | cardiac organ + `qt_ic50_nm` (hERG; dofetilide override vs negligible-hERG defaults) | warfarin/hERG-default compounds stay flat |
| AKI | kidney organ + pharmacophore flag | KDIGO-style GFR/creatinine escalation |
| CNS | `brain_free_exposure` from `organ/cns.py` + `cns_ic50_nm` | **only fused when the molecule is CNS-anchored**; otherwise the 0.20 class prior holds (default 100 µM would wrongly flag APAP) |

Grading: CTCAE-like Grade 0-4 ladders per biomarker (`grade_value`), time-to-onset + duration from threshold crossing (`crossing_interval`). The overall composite risk is resolved to a **low / elevated / High verdict** with the dominant driver named via `overall_driver()`.

### 5.3 Clinical Phenotype Report

Structured output:
- Predicted concentration-time and PK metrics
- Target occupancy & pathway activity profiles
- Predicted physiological indicator trajectories (each graded + CNS brain-free exposure curve)
- Adverse-event risk table per endpoint (DILI, QT prolongation/TdP, nephrotoxicity, CNS effects) with mechanism attribution and confidence
- Sensitivity/driver analysis (which parameters drive the risk margins)

---

## Phase 6 — Uncertainty & Decision Layers (D21-D24)

Shipped implementations in `drugos/robustness/` (all deterministic via fixed seeds):

1. **D21 Parameter ensemble** (`uncertainty.py`) — multiplicative log-normal perturbation (CV default 0.30) of the PK/potency scalars (`cl_hep`, `cl_renal`, `fup`, `qt/dili/cns_ic50`); each member is a full pipeline run; 90% percentile bands (`band()`, q5/q50/q95) for plasma/organ curves, endpoint-risk quantiles, and per-member verdict counts.
2. **D22 Virtual cohort** (`population.py`) — anthropometric sampling (sex/age/height via BMI, deterministic seed); one full pipeline run per individual at the same dose; `risk_quantiles()`, grade ≥ 1 / ≥ 2 `incidence()`, verdict counts.
3. **D23 Sensitivity & drivers** (`sensitivity.py`) — `local_sensitivity` one-at-a-time ±10% normalized log-sensitivities d ln y / d ln p with a ranked `drivers` list; `run_sobol_sensitivity` Saltelli first/total indices (scrambled Sobol design, independent evaluator injectable for CI speed).
4. **D24 Prospective rerun** — reproducibility runbook (identical fixed-seed runs bit-equal) + held-out profile/dose stability; enforced by registration as validation cases (doc/08; see the D21-D24 case in `validation/cases/`).

The Playground exposes D21-D23 as UI toggles (`with_uncertainty` / `with_population`
/ `with_sensitivity`); the CNS brain-free-exposure curve is always plotted in the
organ-panel view and drives the CNS composite risk line only when the molecule is
CNS-anchored (un-gated structural bands are excluded, §5.2).

---

## Monolithic vs Sequential: Recommended Default

- **Default execution = sequential** for the baseline's rigor: run Stage 1 fully, ingest profiles into Stage 2, then Stage 3 and 4, then Stage 5. Simpler to debug and validate per stage.
- **Monolithic execution** optional for small pathway models (single ODE solve); enables full feedback. Enable when the pathway graph <= ~30 nodes.

This staged structure also cleanly maps to test/validation plans in `08-validation-and-risk.md`.