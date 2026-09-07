# Production-Validated Model Integration Matrix

> Design upgrade: every stage of the five-stage pipeline
> (`drug ──▶ in-vivo concentration ──▶ target binding ──▶ signaling pathway
> ──▶ organ function ──▶ clinical phenotype`, `doc/05`) is mapped to the
> **production-validated, open-source, downloadable** model that anchors it —
> with license, source, and integration status.  This is the "upgrade models
> to already production-validated open-source models" design requirement
> made concrete and auditable.  Companion: `doc/10` (dataset comparability),
> `doc/11` (download status), `doc/08` (validation cases incl. R-1/R-2/R-3).

## 0. Rulebook for "production-validated, open-source, downloadable"

A provider must, to qualify for this matrix:

1. be **real**: published, peer-reviewed (or FDA/industry-adopted), with code
   or data obtainable by anyone;
2. be **open-source or openly downloadable**: OSI-approved license (or
   CC-BY/CC0 data) compatible inbound with this project's AGPL-3.0;
3. be **downloadable**: pip/conda-installable, or vendored with a pinned
   sha256 manifest entry under `data/` (`data/manifest.json`);
4. be **verifiable here**: actually run in this repo's validation suite, or —
   if only catalogued for a future phase — blocked by a *named* dependency
   (mapped to a P1–P9 roadmap step in `doc/07`), never silently dropped.

Key constraint carried through every row **G5 (doc/09)**: no silent fallbacks —
if the production model is the source of truth for a stage, the fast analogue
is either *not* present, or is documented as the explicit fast-path encoder
with the production model wired as the validation anchor.

## 1. The matrix

| Stage | Quantity produced | Production-validated model / source | License | Integration status | Runs in this repo |
|---|---|---|---|---|---|
| **1. Drug → in-vivo concentration** | Physicochemical + ADME/T priors | **ADMET-AI** (Swanson et al. 2024, *Nat Mach Intell*; Zenodo weights via PyPI) | BSD-3 (code), model weights CC-BY-4.0 | **WIRED** — Stage-1 ADME primitives (`pk.admet.ADMETPredictor`) | yes (prediction tests + benchmark suite) |
| **1b. PK estimator** | CL, AUC, t½, Vss | **R literature PK** (Wagner 1976; Gibaldi & Perrier 1982; Greenblatt & Koch-Weser 1975; Rowland & Tozer 2010) in `rbridge/literature_pk.R`, executed on **R ≥ 4** (rpy2, required dependency) | BSD-3 (R code, ours) | **WIRED + REQUIRED** — every pipeline run is cross-checked by R (validation R-1, agreement ≤ 2 %) | yes (R-1, gate) |
| **1c. Physiology tables** | Volumes, flows, tissue composition | **Open Systems Pharmacology / PK-Sim** physiology DB (Open-Systems-Pharmacology) | Apache-2.0 | **CATALOGUED** (doc/11 row 4; intended `ospsuite` import; blocked on macOS — OSP engine is .NET/Windows-Linux, no PyPI wheel; P8) | no (see decision record §4) |
| **2. Concentration → target binding** | Off-target occupancy incl. hERG | **ChEMBL** measured hERG IC50 (row-level) + **hERG Central** dose-response corpus (Du et al. 2022; FDA-provenance TDC `Herg` set) | CC0 / CC BY (ChEMBL CC BY-SA 4.0 summary terms; dataverse CC0) | **WIRED** — vendored `data/benchmarks/herg_measured_nm.json` + `data/corpora/herg_central.tsv.gz`; drives R-2 calibration + per-compound hERG override | yes (R-2, gate) |
| **2b. Novel-molecule DTI** | Off-target affinity for un-catalogued chemistry | *(candidate)* **DeepDTA / drug-target ML or ADMET-AI hERG head** | varies | **CATALOGUED** (doc/07 P5) — not wired: ChEMBL panel is source-of-truth today | no |
| **3. Signaling pathway** | Node activity / dose-response | **Huang & Ferrell 1996 ultrasensitive MAPK cascade** (BioModels BIOMD0000000009, CC0) — the canonical ERK cascade, now the pipeline default; **Physiome + Reactome** additional CC BY 4.0 SBML scaffolds catalogued for P7 | CC0 (vendored SBML); CC BY-SA / CC BY 4.0 (future scaffolds) | **WIRED** — `data/models/huang1996-mapk-cascade.xml` + `pathway.sbml_pathway` (python-libsbml); pipeline default `simulate_sbml_pathway`; R-5 pins parse/steady-state/monotone inhibition | yes (R-5, gate) |
| **4. Organ function** | Cardiac repolarisation / QT | **O'Hara-Rudy 2011 (ORd) human ventricular AP model** — validated vs >100 undiseased human hearts (PLoS CB e1002061); BSD-3 Myokit encoding | BSD-3 (Myokit); CC-BY publication | **WIRED** — `data/models/ohara-2011.mmt` + `organ.cardiac_ap` lane; R-3 cross-check on the hERG/QT axis; **CiPA-v1 2017 retune vendored** for the multi-ionic-block upgrade | yes (R-3, gate) |
| **4b. Nephron / kidney** | GFR, AKI grade | **CKD-EPI 2021 race-free equation** (Levey et al., *N Engl J Med* 2021;385:1737) — the standard published clinical eGFR baseline, BSA-scaled per subject; applied whenever a measured serum creatinine is carried on the profile. Deeper *candidate*: **CMR Physiome nephron models** (Layton)/RBF models (CC BY-SA) | CC BY 4.0 (published equation; no code license) | **WIRED** — `organ.kidney.ckdepi_2021_egfr` drives `physiology.gfr_ml_min`; R-6 pins reference eGFR points, BSA scaling, and pipeline wiring; CMR nephron SBML still **CATALOGUED** for P7 | yes (R-6, gate) |
| **4c. Liver / DILI** | Hepatotoxicity grade (cholestasis axis) | **GCDCA bile-acid PBK** (de Bruijn & Rietjens 2024, *Arch. Toxicol.* 98:3077; paper **CC BY 4.0**) — validated against the clinical cholestasis incidence of ~18 marketed drugs; competitive BSEP efflux inhibition by free-hepatic drug drives intrahepatic bile-acid accumulation past the 1.5× risk threshold. Deeper candidates: **DILIsym QSP / Breitwieser et al. 2022 DILI-QSP** (proprietary / catalogued) | CC BY 4.0 (paper equations; repo is CC-BY-NC-ND — code intentionally *not* ported) | **WIRED** — `organ.liver.simulate_gcdca_pbk` anchors `simulate_liver` cholestasis stress; R-7 pins the cholestatic-vs-benign ranking, Ki=IC50/2, and organ→submodel wiring | yes (R-7, gate) |
| **5. Clinical phenotype** | QTc/TdP, DILI, AKI, CNS grading | **ICH E14 / CTCAE ladders + Redfern TdP bands** + QRd-verified cardiac anchor | published thresholds | **WIRED** — `clinical/` ladders; R-3 ORd anchor provides independent ionic confirmation; benchmark L3 cases certify (doc/08) | yes (gate) |

## 2. Status legend

- **WIRED + REQUIRED** — hard dependency; absence = hard error (no fallback,
  G5).  The R bridge is the exemplar.
- **WIRED** — vendored/installed, actively computed in the validation suite and
  gated (sha256 + case).
- **CATALOGUED** — named, licensed, URL-located, mapped to a roadmap step; not
  yet on disk (blocked reason recorded in doc/11 §NOT-DOWNLOADED).

## 3. Experience/intent mapping to pipeline stages

| Pipeline stage | Today (default) | With the production anchor live |
|---|---|---|
| in-vivo concentration | ADMET-AI + PBPK (OSP-style equations) + required-R cross-check | + PK-Sim physiology imported (P8) |
| target binding | ChEMBL measured hERG + class priors | + ChEMBL/DeepDTA DTI for novel molecules (P5) |
| signaling pathway | in-house ODEs (DSL fast path) | + Huang/Levchenko SBML MAPK as the pipeline default (DONE, R-5); Reactome/Physiome scaffolds (P7) |
| organ function | calibrated encoders (QTc Emax, ALT Emax) + CKD-EPI GFR + bile-acid PBK cholestasis | + ORd/IKr ionic cross-check (DONE, R-3); nephron SBML (P7); DILI mito/redox SBML (P6) |
| clinical phenotype | ICH/CTCAE ladders | unchanged (ladders are the published phenotype ground truth) |

## 4. Decision record

- **D1 — required-R instead of optional R (this work).** The literature PK
  estimator is executed by R on every run and agreement is gated ≤ 2 %
  (validation R-1).  Rationale: the "upgrade to production-validated" rule
  says the *primary* estimator should be the literature method; a parallel
  numpy implementation that could silently diverge in the same process is the
  exact failure G5 forbids, so both must always run and agree.
- **D2 — ORd as the cardiac anchor, encoder kept as the fast path.** The
  web/robustness compute path keeps the RdRedfern-calibrated Emax QTc encoder
  (validated against wolf-qtc clinical bands in the cardiac L3 case);
  ORd-2011 (and the vendored CiPA-v1 2017 retune) run as the independent ionic
  anchor (R-3).  This is *not* a fallback: the anchor is mandatory in the
  validation suite, and the encoder is not silently substituted (R-3 asserts
  ordering direction and magnitude class).  A full main-path ORd/QTw
  computation (multi-channel block incl. ICaL/INaL/IKs/IK1 via the CiPA
  retune) is P9.
- **D3 — why FDA/CiPA MATLAB/C was not in-process.** The official FDA/CiPA
  uncertainty code (`github.com/FDA/CiPA`, GPL-3.0) is license-compatible in
  principle (AGPL/GPL) but runs on R/Matlab and does not fit the packaged
  Python wheel; the same equations are available BSD-3 from `myokit/models`
  (`ohara-cipa-v1-2017.mmt`) — that is the vendored artifact.
- **D4 — OSPSuite/PK-Sim blocked on this host.** `ospsuite` ships no PyPI
  wheel; the OSP runtime is .NET and Linux/Windows-only (macOS via Docker).
  Documented not-viable here; physiology equations already mirror the
  Willmann/OSP power-laws (doc/11 row 5) and P8 tracks the Docker route.
- **D5 — ADMET-AI is the production-validated ADME/T frontend.** BSD-3 code +
  Zenodo weights, used widely in pharma AI pipelines; it *is* the primary
  "production-validated, downloadable, open-source" concentration-stage model
  (row 1).  Its hERG/TdP heads feed the Stage-2 off-target prior for novel
  chemistry pending DTI (D2/P5).
- **D6 — Huang/Levchenko SBML as the Stage-3 (pathway) production anchor.**
  The pipeline default is now the vendored BioModels BIOMD0000000009 (CC0)
  parsed through python-libsbml; the DSL cascade remains the auditable fast
  path / equivalence test, not a silent fallback (R-5 pins the production lane;
  the DSL is exercised by the analytic amplification case in parallel, and the
  two are never mixed in one run).  Coupling: target occupancy s(t) scales the
  upstream MAPKKK activator E1 down (`E1_eff = E1·(1−s)`), readout is
  doubly-phosphorylated ERK (PP_K); drug-free system sits at its fully
  activated steady state, occupancy suppresses the readout monotonically.
- **D7 — de Bruijn & Rietjens 2024 GCDCA bile-acid PBK as the liver
  cholestasis production anchor.** DILIsym is proprietary (no download) and
  the Breitwieser DILI-QSP catalog entry remains the deeper P6 target; the
  published, openly licensed (CC BY 4.0) bile-acid PBK of de Bruijn &
  Rietjens reproduces the *validated* human behavior — clinical cholestasis
  incidence of ~18 marketed drugs — via competitive BSEP inhibition.  Its
  companion GitHub repository (Veronique-de-Bruijn/PBK-model-cholestasis) is
  **CC-BY-NC-ND**: we do *not* vendor or re-derive its R code; the equations
  ported into `organ.liver` are the paper's open-access ones, and only the
  published numeric BSEP IC50 anchors are vendored
  (`data/models/bsep_shh_ic50_reference.json`, CC BY 4.0 attribution in the
  manifest).  R-7 pins the cholestatic-vs-benign *ranking* — the model's key
  falsifiable output — plus Ki=IC50/2 and the organ→submodel wiring.
- **D8 — CKD-EPI 2021 race-free as the kidney GFR baseline.** The prior
  age/sex default (125 mL/min) is replaced, whenever a measured serum
  creatinine is present on the profile, by the CKD-EPI 2021 creatinine
  equation (Levey et al., NEJM 2021) — the accepted clinical standard, openly
  published, BSA-scaled per subject.  The depth-3 CMR Physiome nephron SBML
  remains catalogued for P7; GFR baseline and AKI grading do not silently
  fall back (R-6 asserts the benchmark invariance when Scr is absent).

## 5. Gate-keeping notes

- Every `data/` file above is sha256-pinned (`data/manifest.json`); adding any
  new file without a manifest entry fails `gates` (checksum gate).
- New `except` handlers in provider wrappers must register with
  `scripts/fallback_allowlist.json` (G5).
- The ORd lane needs `myokit>=1.39` (pip) **and** SUNDIALS headers
  (`conda install -c conda-forge sundials`, `apt install libsundials-dev`, or
  `brew install sundials`) for myokit's C codegen.  This is a documented
  runtime requirement for the validation suite (doc/06).
- Validation cases introduced by this matrix: **R-1** (R PK agreement),
  **R-2** (corpus/measured-hERG calibration), **R-3** (ORd APD90 cardiac
  cross-check), **R-4** (ADMET-AI BBB→CNS partition), **R-5** (Huang/Levchenko
  SBML MAPK cascade integration), **R-6** (CKD-EPI 2021 race-free GFR
  baseline), **R-7** (de Bruijn & Rietjens bile-acid PBK cholestasis ranking).
  All seven run in the G4 gate.

## 6. References (additions beyond doc/02)

- O'Hara T, Virág L, Varró A, Rudy Y. *Simulation of the undiseased human
  cardiac ventricular action potential...* PLoS Comput Biol 2011;7(5):e1002061.
- Dutta S, Chang KC, et al. (Li Z). *Optimization of an in-silico cardiac cell
  model for proarrhythmia risk assessment.* Front Physiol 2017;8:616.
- Swanson K, Walther P, et al. *ADMET-AI: a machine learning ADMET platform.*
  Nat Mach Intell 2024.
- Du F, et al. *hERG Central.* J Chem Inf Model 2022 (Harvard Dataverse
  doi:10.7910/DVN/7BVDG8).
- Wilhelms M, et al. / CiPA in-silico guidance — multi-channel ORd use (P9).
- Mirams G, et al. Myokit: a framework for cardiac cell modelling
  (myokit.org), BSD-3.
- Huang CY, Ferrell JE Jr. *Ultrasensitivity in the mitogen-activated protein
  kinase cascade.* Proc Natl Acad Sci USA 1996;93(19):10078-83.  SBML encoding
  curated on BioModels as BIOMD0000000009 (CC0).
- Levchenko A, Bruck J, Sternberg PW. *Scaffold proteins may biphasically
  affect the levels of mitogen-activated protein kinase signaling and reduce
  its threshold properties.* Proc Natl Acad Sci USA 2000;97(11):5818-23.
- S. Hoops et al. *COPASI* / libSBML — python-libsbml BSD licensed; SBML L2V4
  parsing for the vendored pathway model.
- Levey AS, Stevens LA, Schmid CH, et al. *A new equation to estimate
  glomerular filtration rate.* Ann Intern Med 2009;150:604-12 (CKD-EPI 2009).
- Inker LA, Eneanya ND, Coresh J, et al. *New creatinine- and
  cystatin-C–based equations to estimate GFR without race.*
  N Engl J Med 2021;385:1737-49 (CKD-EPI 2021, race-free — R-6).
- Mosteller RD. *Simplified calculation of body-surface area.*
  N Engl J Med 1987;317:1098 (BSA scaling for absolute eGFR).
- de Bruijn V, Rietjens IMCM. *Mechanistic insights on the role of the
  Bile Salt Export Pump in drug-induced cholestasis...* Arch Toxicol 2024;98:
  3077-3095. doi:10.1007/s00204-024-03775-6 (CC BY 4.0 — R-7).