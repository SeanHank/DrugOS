# Dataset Requirements, Availability & Scientific Comparability

> Companion to `doc/04-data-sources.md` (catalog), `doc/05-methodology-pipeline.md`
> (what each stage consumes), `doc/11-download-status.md` (what is actually on
> disk today), and `data/README.md` (vendoring + checksum policy).
> Dates are retrieve/cutoff aware; every commercial dataset below requires an
> academic license and rate-limit-aware fetching.

## 0. Scope

This report answers three questions for every dataset the pipeline **needs**
(gap-closing datasets) and every dataset it **already maps to** (current
anchors):

1. **Needed**: which model quantity would it parameterize/validate?
2. **Availability**: official source, license, format, approximate size, access
   route, and whether it is already on disk (see `doc/11`).
3. **Scientific comparability**: does the dataset measure the same construct we
   compute, in the same units, on the same human-relevant scale — and how
   strong is the evidence line it would carry (anchored L3 vs. class-prior L2)?

Stages refer to the five-stage pipeline in `doc/03`/`doc/05`.

## 1. Current anchors (already in the loop) and how comparable they are

### 1.1 Published PK reference corpus — `data/benchmarks/published_pk.json` (VENDORED)

| Quantity | Dataset | Comparability |
|---|---|---|
| CL, t½, Fa, fe, Vss for midazolam / acetaminophen / warfarin / ciprofloxacin / dofetilide | Goodman & Gilman (13th ed.), host USPIs, Bergan 1986, Greenblatt 1984, Smith 1992 | Direct: same measurable PK metrics (units L/h, h, fraction, L). Used as the L2/L3 pass-bands in validation. Highest-comparability anchor in the repo. |

**Science note.** Literature PK bands are inter-study medians with wide ranges
(compound-specific published variability); they validate *pass-bands*, not
point accuracy. This is the correct bar for a mechanism-based simulator.

### 1.2 ADMET-AI GNN package (pip, OSS) — used only on the custom-SMILES path

| Prediction | Package column | Comparability vs. real assay |
|---|---|---|
| fup | `PPBR_AZ` (0-100%) | Plasma protein binding (equilibrium dialysis / ultracentrifugation corpus). Converted to fup; comparable quantity, model-derived, not experimental. |
| Cl_int hep | `Clearance_Hepatocyte_AZ` | Human hepatocyte intrinsic clearance (mL/min/kg); same construct as our `cl_hep` scaling. |
| t½, Vss | `Half_Life_Obach`, `VDss_Lombardo` | Human PK descriptors estimated from structure/corpus; same units as our NCA outputs. |
| hERG, DILI, AMES, BBB, P-gp, HIA, CYP panel | categorical/continuous | Structural flags used as weak independent lines in Stage-5 fusion; comparable to *in-silico* flags, NOT to clinical labels. |

**Science note.** "Comparability" here means the mapping is onto the same
biophysical construct, but the dataset is an ML estimate with forecast
uncertainty. The pipeline treats these as low-weight structural evidence and
never as hard truth; this is the defensible reading.

### 1.3 Hardcoded class-prior physiology / target / organ tables

Rodgers–Rowland tissue composition (11 tissues), Ye 2016 organ volumes/flows,
Willmann-style scaling, the 16-site safety panel Kd/IC50 medians, organ
IC50 priors (BSEP 300 µM IC50 / Ki 150 µM benign default, mito 300 µM,
redox 600 µM, hERG 2 nM, neurotox 100 µM),
CTCAE-type grade ladders, toxicity fusion priors/weights. All are hardcoded
literals (see `doc/11`, rows "not downloaded"; the BSEP anchors are the
vendored de Bruijn & Rietjens reference values, row 19b).

**Science note.** These are *literature-class medians*: comparable in construct
and unit to the corresponding assays/registries, but they carry class-typical
(not compound-specific) information. Their correct scientific role is as
calibrated priors, which is exactly how the code uses them (they are always
overridable on `RunSpec`).

## 2. Gap-closing datasets: needs, availability, comparability

Priority ranking follows `doc/04 §3` (ingestion order). For each row:
**Need → Source / access / license / size → Scientific comparability analysis.**

### P1. OSP / PK-Sim physiology database (organ volumes, flows, composition)

- **Need.** Replaces `REF_VOLUME_L`/`REF_FLOW_L_MIN`/`TISSUE_COMPOSITION` hard
  literals; enables age/sex/BMI-realistic parameter sets (current model has
  only power-law volume scaling + blunt protein/GFR cutoffs — largest realism
  gap per code review).
- **Availability.** Open Systems Pharmacology Suite (open source, Apache-2.0;
  `github.com/Open-Systems-Pharmacology`), physiology database ships with
  PK-Sim; example physiology imports via `ospsuite` (R/Python). Size ~MBs.
- **Comparability.** Same constructs we compute (organ volumes in L, flows in
  L/min, tissue water/lipid fractions) at cohort level. This dataset is the
  *ground source* for two of our cited equations (Willmann/Ye), so ingestion is
  a strict upgrade with zero construct mismatch. Sci-reliability: high.

### P2. ChEMBL bioactivity + DrugBank target pharmacology

- **Need.** Replace the 16-site hardcoded panel Kd/IC50 medians with measured
  binding/functional potency (Stage 2 Kd priors) and approved-target
  mechanism.
- **Availability.** ChEMBL (EMBL-EBI, CC BY-SA 4.0, bulk FTP ~GBs, query API).
  DrugBank (full distribution: academic licence application, ~hundreds MB XML).
- **Comparability.** ChEMBL pChEMBL/pIC50/Ki are directly the quantity our
  `Target.kd_nm` consumes (nM). **Convert via Cheng–Prusoff, never liter-wire
  IC50→Kd** (the current `kd_from_ic50_um` is the acknowledged simplification).
  Sci-reliability: measured data, assay-units caveat → fit with an assay-variance
  model. DrugBank targets define *which* panel sites matter per compound
  (construct match: mechanism of action lists vs. our fixed panel).

### P3. hERG / TdP assay corpus (ChEMBL-assayed hERG + EC50-summary corpora)

- **Need.** Calibrate the hERG→IKr→QTc Emax (currently b50=0.5, ΔQTc max 40 ms
  are arbitrary; the pipeline uses one fixed sigmoid for all compounds).
- **Availability.** ChEMBL (hERG binding/blockade measurements, free);
  plus literature compilations (e.g., published TdP/hERG logistic corpora with
  Cmax/IC50 labels). 
- **Comparability.** hERG patch-clamp IC50 (nM) = the exact anchor the model
  wants (`qt_ic50_nm`). The *scientific limit* is that hERG IC50 correlates with
  (does not equal) TdP risk; a compound-level logistic on
  Cmax_unbound/hERG_IC50 is the clinically validated correlate (CNS/DBQT,
  Redfern 2003 lineage). This is the single highest-value toxicity-calibration
  enabler.

### P4. DILI corpus (LiverTox + TDC DILI labels + fezolinetant/DILIsym-like in-vitro assays)

- **Need.** Replace the hardcoded DILI exposure sigmoid (slope 1.2, center
  −0.5) and the fixed grade-prob ladder with a dataset-fitted relationship and
  an *in-vitro* parameter table (BSEP inhibition, ETC, oxidative stress) that
  the liver model's default IC50s (90/300/600 µM) should derive from.
- **Availability.** LiverTox (NIH/NCATS, free to query, bulk under academic
  terms); TDC DILI (open); OATP/BSEP inhibition literature tables.
- **Comparability.** LiverTox labels are *clinical* DILI causality (the
  endpoints we grade); TDC DILI is a drug-level binary — construct match at the
  compound level, not dose-response. The strongest comparable model is the
  241-drug Cmax/IC50 ROC (AUC 0.91–0.96) our docstring cites but does not
  embed; fitting the same exposure-ratio discriminant to a snapshot is an
  exact, upgradeable comparability.

### P5. KEGG / Reactome / PANTHER pathway wiring

- **Need.** Replace the single hand-coded 6-species `mapk_cascade()` with
  curated networks (Stage 3). Currently arbitrary kinetic constants; no KEGG/
  Reactome/OmniPath loader exists in `src/`.
- **Availability.** KEGG (license required for academic use, `fTP` mirror);
  Reactome (CC BY 4.0, graph + SBML, free); PANTHER (free). Sizes: MBs-tens of
  MBs.
- **Comparability.** Reactome/KEGG store *network topology & causal relations*
  (construct: signal flow) but **almost no kinetic rate constants** — we would
  still need kinetic parameterization (Kholodenko-style literature rates).
  Comparability caveat: topology ≠ dynamics; a loaded network with ARBITRARY
  rates is not more realistic than the current cascade unless rate constants
  come from quantitative biochemistry (Kd, Vmax, kcat from BRENDA/SABIO-RK).

### P6. Physiome organ models (heart AP, nephron, liver)

- **Need.** Replace static/algebraic organ ports: cardiac ion-channel Markov
  APD model (no state dynamics now), nephron tubule transport + creatinine
  kinetic lag, liver DILI lifecycle beyond the current bile-acid + death axes.
- **Partial today.** The liver cholestasis axis is already the **de Bruijn &
  Rietjens 2024 bile-acid PBK** (R-7) and the kidney GFR baseline is the
  **CKD-EPI 2021 race-free equation** (R-6) — both production-validated open
  anchors (doc/12 rows 4b/4c); a full Physiome nephron/liver SBML integration
  remains this P-item.
- **Availability.** Physiome Model Repository (CC BY-SA, SBML, free). Size:
  individual models ~100 KB–10 MB SBML.
- **Comparability.** The benchmark QTc standards here are *human ECG*
  (Thüringer's QTc project, or the E14/ICH derivation), and nephron kinetics
  validate against Scr time-course clinical data. Construct match is good, but
  these are complex systems whose baseline/calibration must be reproduced in
  our units — a multi-week validation exercise, not a drop-in file.

### P7. Population covariates (NHANES / OSP virtual populations)

- **Need.** Replace `population.py`'s parametric guess (age uniform 18–80,
  Height N(177.5/163.5, 7), BMI N(26, 4)) with correlated cross-sectional
  anthropometry + renal/hepatic function covariates.
- **Availability.** NHANES (US CDC, public-domain, bulk SAS/XPT); OSP virtual
  populations (open). Size: tens to hundreds MB.
- **Comparability.** Direct: age, sex, height, weight, BMI, Scr/GFR covariates
  are measured on the very quantities `resolve_human`/`population.py` sample,
  so correlation structure (e.g., height×weight, age×GFR decline) becomes data
  rather than assumption. Highest-comparability population-data fix.

### P8. QTc / ECG clinical baselines

- **Need.** An empirical human QTc/ΔQTc distribution (not a point 415 ms) and
  the ΔQTc↔risk calibration.
- **Availability.** Published cardio-safety databases and the literature
  ΔQTc↔TdP models (ICH E14 summaries). 
- **Comparability.** Our `delta_qtc = 40ms*(1−IKr_frac)` is a single-Emax
  surrogate; a clinical distribution makes the QT risk line an evidence-based
  distribution shift. Caveat: E14 data are pharma-confidential; only
  published aggregates are usable under citation (never redistribute raw XML).

### P9. DILI-on-hERG-free toxicity precedence datasets (SIDER / Offsides / FAERS)

- **Need.** FDA labeling/FEARS-derived endpoint co-occurrence for evidence
  dependencies *between* endpoints (currently assumed independent: weighted
  log-odds sum).
- **Availability.** SIDER (CC BY-NC-SA, free), Offsides (free), FAERS (public
  FDA, requires pruning for duplicate/confounded reports). Size: hundreds MB.
- **Comparability.** These label *drug→AE* associations (epidemiological) — not
  the mechanism-specific quantities we model. Use only as a cross-check of the
  *direction* of endpoint co-occurrence, never as a quantitative parameter
  source (weak comparability).

## 3. Comparability principles enforced by this report

1. **Same construct**: dataset quantity must be the same physical/clinical
   quantity the pipeline computes (µM IC50 → nM Kd only via Cheng–Prusoff;
   mL/min/kg hepatic intrinsic clearance → L/h, never a bare unit flip).
2. **Same time/reference frame**: clinical band datasets validate pass-bands
   (L2/L3), single-study point datasets validate point accuracy; a model can
   only claim the stronger claim if the reference supports it.
3. **Uncertainty carries over**: ML/structural flags enter as low-weight
   evidence lines; measured assay data enter as calibrated priors; only
   clinical outcome corpora can calibrate the Stage-5 risk fusion.
4. **License gates admission**: commercial/academic-gated bulk corpora
   (DrugBank full, LiverTox bulk, KEGG) are documented but never shipped;
   checksum-pinned snapshots are recorded in `data/manifest.json`.

## 4. Recommendation summary (cheapest → biggest comparability gain)

1. **P1 OSP physiology + P7 NHANES** — pure construct-match, replaces our least
   real tables (cheapest, highest trust).
2. **P3 hERG + P4 DILI** exposure-ratio calibration — converts two toxicity
   lines from calibrated priors to data-fit discriminants.
3. **P2 ChEMBL/DrugBank** — kills the "class-median panel" simplification.
4. **P5/P6 pathway/organ** — highest effort; topology without kinetic rates
   does not raise realism (must pair with BRENDA/SABIO-RK rate sets).
5. **R bridge (`scripts/r_crossval`)**: operate the P1/P3 fit lines in R's
   pharmacometric ecosystem (nlmixr2/stats) for independent cross-validation
   and VPC-style checks (see `scripts/r_crossval/README.md`).