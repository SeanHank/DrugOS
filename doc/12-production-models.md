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
   sha256 manifest entry under data/ (`data/manifest.json`);
4. be **verifiable here**: actually run in this repo's validation suite and
   gated by it (G5, doc/09); nothing is admitted to this matrix that is not
   reproduced from its named source inside this repository.

Key constraint carried through every row **G5 (doc/09)**: no silent fallbacks —
if the production model is the source of truth for a stage, the fast analogue
is either *not* present, or is documented as the explicit fast-path encoder
with the production model wired as the validation anchor.

## 1. The matrix

| Stage | Quantity produced | Production-validated model / source | License | Integration status | Runs in this repo |
|---|---|---|---|---|---|
| **1. Drug → in-vivo concentration** | Physicochemical + ADME/T priors | **ADMET-AI** (Swanson et al. 2024, *Nat Mach Intell*; Zenodo weights via PyPI) | BSD-3 (code), model weights CC-BY-4.0 | **WIRED** — Stage-1 ADME primitives (`pk.admet.ADMETPredictor`) | yes (prediction tests + benchmark suite) |
| **1b. PK estimator** | CL, AUC, t½, Vss | **R literature PK** (Wagner 1976; Gibaldi & Perrier 1982; Greenblatt & Koch-Weser 1975; Rowland & Tozer 2010) in `rbridge/literature_pk.R`, executed on **R ≥ 4.5** (rpy2 ≥ 3.6, required dependency) | BSD-3 (R code, ours) | **WIRED + REQUIRED** — every pipeline run is cross-checked by R (validation R-1, agreement ≤ 2 %) | yes (R-1, gate) |
| **1c. Physiology tables** | Volumes, flows, tissue composition | **Willmann et al. (2007) population allometric power laws** — the same equations PK-Sim's population module embeds, resolved for a 70 kg reference male against the Ye (2016) organ volume/flow tables | published equations; BSD-3 (ours) | **WIRED** — `drugos.pk.physiology.build_human` / `HumanPhysiology`; `case_willmann_allometric_physiology` pins reference self-consistency, per-organ exponents, sex factors, Mosteller BSA, cardiac-output scaling, impairment modifiers | yes (PHYS-L2, gate) |
| **2. Concentration → target binding** | Off-target occupancy incl. hERG | **ChEMBL** measured hERG IC50 (row-level) + **hERG Central** dose-response corpus (Du et al. 2022; FDA-provenance TDC `Herg` set) | CC0 / CC BY (ChEMBL CC BY-SA 4.0 summary terms; dataverse CC0) | **WIRED** — `data/benchmarks/herg_measured_nm.json` + `data/corpora/herg_central.tsv.gz`; drives R-2 calibration + per-compound hERG override | yes (R-2, gate) |
| **2b. Novel-molecule DTI** | Off-target affinity for novel chemistry | **ChEMBL bioactivity snapshot** (16 targets, 3035 rows; `data/chembl/offtarget_snapshot.json`, sha256-pinned) queried by a **Tanimoto kNN resolver** (`drugos.target.dti`): `predict_kd_nm` returns the fingerprint-weighted neighbor KD when chemotype support is present (best-pair Tanimoto ≥ `MIN_NEIGHBOR_TANIMOTO`); `resolve_offtarget_panel` re-binds mapped panel sites to structure-derived KD only with such support | BSD-3 (code, ours); ChEMBL CC BY-SA 4.0 terms (data, mirrored) | **WIRED** — `drugos.target.dti.cns_ic50_nm` anchors the CNS seam on chemotype-supported brain-exposed sites; unsupported sites keep disclosed class priors (`no_public_data_sites`); leave-one-molecule-out calibration `case_dti_resolver_calibration` (per-target positive correlation, honest scatter, band-clamped) | yes (DTI-L2, gate) |
| **3. Signaling pathway** | Node activity / dose-response | **Huang & Ferrell 1996 ultrasensitive MAPK cascade** (BioModels BIOMD0000000009, CC0) — the canonical ERK cascade, the pipeline default | CC0 (vendored SBML) | **WIRED** — `data/models/huang1996-mapk-cascade.xml` + `pathway.sbml_pathway` (python-libsbml); pipeline default `simulate_sbml_pathway`; R-5 pins parse/steady-state/monotone inhibition | yes (R-5, gate) |
| **4. Organ function** | Cardiac repolarisation / QT | **O'Hara-Rudy 2011 (ORd) human ventricular AP model** — validated vs >100 undiseased human hearts (PLoS CB e1002061); BSD-3 Myokit encoding | BSD-3 (Myokit); CC-BY publication | **WIRED** — `data/models/ohara-2011.mmt` + `organ.cardiac_ap` lane; R-3 cross-check on the hERG/QT axis; **CiPA-v1 2017 retune vendored** (`data/models/ohara-cipa-v1-2017.mmt`) as its own multi-ionic ORd lane (`case_ord_multi_ionic_qt`) | yes (R-3 + ORD-L19, gate) |
| **4b. Nephron / kidney** | GFR, AKI grade, tubular handling | **CKD-EPI 2021 race-free equation** (Levey et al., *N Engl J Med* 2021;385:1737) — the standard published clinical eGFR baseline, BSA-scaled per subject; **nephron tubular-transport model** (`organ.nephron`: proximal reabsorption, transporter-mediated secretion kinetics) | CC BY 4.0 (published equation; no code license); BSD-3 (nephron, ours) | **WIRED** — `organ.kidney.ckdepi_2021_egfr` drives `physiology.gfr_ml_min`; R-6 pins reference eGFR points, BSA scaling, and pipeline wiring; `nephron_handling` resolves tubular reabsorption/secretion with secretory-kinetics fidelity | yes (R-6 + NEPH-L2, gate) |
| **4c. Liver / DILI** | Hepatotoxicity grade (cholestasis axis) | **GCDCA bile-acid PBK** (de Bruijn & Rietjens 2024, *Arch. Toxicol.* 98:3077; paper **CC BY 4.0**) — validated against the clinical cholestasis incidence of ~18 marketed drugs; competitive BSEP efflux inhibition by free-hepatic drug drives intrahepatic bile-acid accumulation past the 1.5× risk threshold.  Flanking shipped axes: mitochondrial ETC block, redox/GSH, hepatocyte-death ODEs | CC BY 4.0 (paper equations; companion GitHub repo is CC-BY-NC-ND — its R code is intentionally *not* ported); BSD-3 (axes, ours) | **WIRED** — `organ.liver.simulate_gcdca_pbk` anchors `simulate_liver` cholestasis stress; R-7 pins the cholestatic-vs-benign ranking, Ki=IC50/2, and organ→submodel wiring; `case_mito_redox_dili` pins the mitochondrial/redox axes | yes (R-7 + MITO-L2, gate) |
| **5. Clinical phenotype** | QTc/TdP, DILI, AKI, CNS grading | **ICH E14 / CTCAE ladders + Redfern TdP bands** + QRd-verified cardiac anchor | published thresholds | **WIRED** — `clinical/` ladders; R-3 ORd anchor provides independent ionic confirmation; benchmark L3 cases certify (doc/08) | yes (gate) |

## 2. Status legend

- **WIRED + REQUIRED** — hard dependency; absence = hard error (no fallback,
  G5).  The R bridge is the exemplar.
- **WIRED** — vendored/installed, actively computed in the validation suite and
  gated (sha256 + case).
- **WIRED (live)** — active in the default path or in the validation gate and
  disclosed in every contract's `trust` record.

## 3. Stage mapping to shipped production anchors

| Pipeline stage | Shipped production anchor |
|---|---|
| in-vivo concentration | ADMET-AI primitives (row 1) + R literature PK cross-check (row 1b) + Willmann 2007 allometric physiology (row 1c) |
| target binding | ChEMBL measured hERG + hERG Central corpus (row 2) + ChEMBL kNN DTI resolver (row 2b); sites without structure support stay on disclosed class priors (`no_public_data_sites`) |
| signaling pathway | Huang & Ferrell 1996 SBML MAPK cascade as the pipeline default (row 3, R-5) |
| organ function | ORd-2011 + CiPA-v1 2017 multi-ionic lane (row 4, R-3), CKD-EPI 2021 GFR + nephron tubular transport (row 4b), GCDCA bile-acid cholestasis PBK + mitochondrial/redox axes (row 4c) |
| clinical phenotype | ICH/CTCAE ladders (row 5): the published phenotype ground truth |

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
  ordering direction and magnitude class).  The main-path multi-channel block
  (ICaL/INaL/IKs/IK1 via the CiPA-v1 2017 retune) runs as its own ORd lane
  (`organ.cardiac_ap.cipa_apd90`), pinned by `case_ord_multi_ionic_qt`.
- **D3 — why FDA/CiPA MATLAB/C was not in-process.** The official FDA/CiPA
  uncertainty code (`github.com/FDA/CiPA`, GPL-3.0) is license-compatible in
  principle (AGPL/GPL) but runs on R/Matlab and does not fit the packaged
  Python wheel; the same equations are available BSD-3 from `myokit/models`
  (`ohara-cipa-v1-2017.mmt`) — that is the vendored artifact.
- **D4 — physiology tables are Willmann-allometric, not OSPSuite/PK-Sim.**
  The PK-Sim/OSP population database was evaluated and is not used here:
  the OSP runtime is .NET and Linux/Windows-only, with no macOS wheel.  The
  shipped physiology (`drugos.pk.physiology.build_human`) implements the same
  published allometric power laws that PK-Sim's population module embeds —
  Willmann et al. (2007) — with `case_willmann_allometric_physiology` pinning
  reference self-consistency against the Ye (2016) tables, per-organ
  exponents, sex factors, BSA scaling and impairment modifiers.
- **D5 — ADMET-AI is the production-validated ADME/T frontend.** BSD-3 code +
  Zenodo weights, used widely in pharma AI pipelines; it *is* the primary
  "production-validated, downloadable, open-source" concentration-stage model
  (row 1).  Its hERG/TdP heads feed the Stage-2 off-target prior; the
  structure-derived lane for novel chemistry is the ChEMBL kNN resolver of
  `drugos.target.dti` (row 2b, D10), calibrated leave-one-molecule-out.
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
  the Breitwieser DILI-QSP repository is CC-BY-NC-ND, so neither is ported;
  the openly licensed (CC BY 4.0) bile-acid PBK of de Bruijn & Rietjens
  reproduces the *validated* human behavior — clinical cholestasis
  incidence of ~18 marketed drugs — via competitive BSEP inhibition.  Its
  companion GitHub repository (Veronique-de-Bruijn/PBK-model-cholestasis) is
  **CC-BY-NC-ND**: we do *not* vendor or re-derive its R code; the equations
  ported into `organ.liver` are the paper's open-access ones, and only the
  published numeric BSEP IC50 anchors are vendored
  (`data/models/bsep_shh_ic50_reference.json`, CC BY 4.0 attribution in the
  manifest).  R-7 pins the cholestatic-vs-benign *ranking* — the model's key
  falsifiable output — plus Ki=IC50/2 and the organ→submodel wiring.  The
  mitochondrial ETC / redox / hepatocyte-death axes flanking the cholestasis
  anchor are the shipped liver lane (`organ.liver.mitochondrial_block`,
  `organ.liver.redox_state`), pinned by `case_mito_redox_dili`.
- **D8 — CKD-EPI 2021 race-free as the kidney GFR baseline.** The prior
  age/sex default (125 mL/min) is replaced, whenever a measured serum
  creatinine is present on the profile, by the CKD-EPI 2021 creatinine
  equation (Levey et al., NEJM 2021) — the accepted clinical standard, openly
  published, BSA-scaled per subject.  The depth-3 kidney lane is the nephron
  tubular-transport model of `organ.nephron` (`nephron_handling`): proximal
  reabsorption and transporter-mediated secretion resolved from biliary/renal
  kinetics, pinned by `case_nephron_tubular_transport`.  GFR baseline and AKI
  grading do not silently fall back (R-6 asserts the benchmark invariance when
  Scr is absent).
- **D9 — Stage-1 extended clearance/absorption realism, off by default.**
  The validated linear-PBPK baseline is the shipped default; the realism
  extensions in `pk/pbpk_build.py` — active tubular secretion
  (`cl_sec_l_h`), saturable Michaelis–Menten hepatic metabolism
  (`hepatic_vmax_mg_h`/`hepatic_km_mg_l`), per-CYP MM/Hill kinetics
  (`cyp_terms`, abundance-scaled Vmax from the `CYP_ABUNDANCE_PMOL_MG`
  table via `cyp_vmax_mg_h`), biliary drug excretion with
  enterohepatic recirculation (`cl_bil_l_h`/`bile_emptying_1h`), and first-
  pass gut-wall extraction (`gut_extraction_eg`) — default to 0/off so every
  benchmark, L3 PK case and mass-balance case reproduces the linear baseline
  exactly.  They are threaded through `RunSpec`/`spec_from_admet` and pinned
  by the dedicated L2 cases in G4 (`case_clearance_mechanisms`: analytic
  single-pool secretion urine fraction, `F = Fa·(1−Eh)·(1−Eg)`, low-dose
  `Vmax/Km` slope, mass-conservative EHC; `case_cyp_kinetics`: abundance-
  scaled Vmax from the physiology table, isoform additivity to the linear
  twin, doubling-content slope scaling, Hill sigmoid, parameter rejection).

- **D10 — Off-target resolution is two shipped lanes: corpus-calibrated
  specialist heads plus a structure-based ChEMBL kNN resolver.**  The
  bare binary hERG sieve is replaced by a single continuous, strictly
  monotone probability→KD curve (`drugos.target.resolver.kd_from_score`) that
  is *never more potent than* the panel class prior: `P=1` reproduces the
  2 nM hERG prior exactly, `P=0` the 1 mM weak floor, and the classifier
  threshold `P=0.5` lands at ~1.4 µM — the corpus-typical weak potency, within
  10x of the hERG Central median on the conservative side.  The head lanes
  re-score hERG and the Veith CYP2D6/3A4/2C9 inhibition heads through this
  curve.  `drugos.target.dti` adds the second lane: a fingerprint kNN over the
  vendored 16-target snapshot (`predict_kd_nm`, `top_tanimoto`,
  `MIN_NEIGHBOR_TANIMOTO`), which re-binds a mapped panel site to a
  structure-derived KD *only* when chemotype support is present.  Every panel
  site without structure support stays on its disclosed class prior and is
  named in `no_public_data_sites`.  P-gp remains a *substrate* classifier
  (host to inhibition potency); transporters, mitochondrial/endocrine sites
  and unsupported primary-target sites keep class priors.  R-8
  (`case_herg_calibration`, doc/08) pins the curve invariants, the corpus
  anchor window and the head's own concordance with the corpus inhibition
  labels; `case_dti_resolver_calibration` pins the kNN resolver's real
  leave-one-molecule-out accuracy on the snapshot.

- **D11 — Native target-mediated drug disposition as a Stage-1
  mass-balance coupling, off by default.**  The sequential pipeline only
  *approximated* TMDD via the opt-in `feedback_loop` driver (monolithic
  feedback without Stage-2 → Stage-1 flux coupling).  Now
  `PBPKModel.target_binding` (with `mw_g_per_mol`) couples one binding site's
  turnover ODEs directly into a tissue's mass balance, using the exact
  occupancy equations of `drugos.target.occupancy`: the net bound flux
  `kon·D·R − koff·DR` (nmol/h, the volume cancels through the tissue volume)
  leaves the tissue pool, and `kint·DR` internalization drains into a cleared
  sink that is part of `state_total_mass` (receptor is protein and excluded),
  so administered drug mass is conserved-tallied even with a sink.  This is
  the resource the pipeline had asked for on high-affinity, high-abundance
  targets.  Pinned in G4 by `case_tmdd_drug_disposition`: reversible mass
  closure, irreversible sink, dose-disproportional retention (the TMDD
  hallmark), quasi-steady KD recovery `DR/R = D/Kd`, AUC reduction vs the
  free-dispersion twin, off-by-default state count, and degenerate-site
  rejection.  Off by default so every benchmark/L3 case reproduces the linear
  baseline exactly.

- **D12 — Finite-dose multi-layer skin permeation as the transdermal
  absorption path, off by default.**  The transdermal design goal asked for a
  real skin-permeation membrane instead of the generic transdermal depot.
  `AbsorptionParams.skin_layers` (`SkinLayers`) adds four compartments
  (vehicle surface -> stratum corneum -> viable epidermis -> dermis) coupled
  by reversible diffusion-limited links `J_k = D_k*A/L_k*(C_up - C_down/K_k)`
  with interface partition ratios and first-order dermal capillary removal
  into venous blood (weighted by `depot_bioavailability`, complement tallied
  in the skin unabsorbed sink).  The serial 3-resistance membrane reproduces
  the Fick steady state `J_ss = P_eff*A*C_surface` and, at zero flux, the
  interface partitions; all five layer states stay inside `state_total_mass`
  so transdermal mass closes against the dose.  Pinned in G4 by
  `case_transdermal_multi_layer`: mass closure, composite-permeability flux,
  partition recovery, stratum-corneum barrier and diffusivity responsiveness,
  off-by-default state count, degenerate-membrane rejection.  Off by default
  so the validated depot semantics remain the shipped baseline.

- **D13 — Sympathetic suppression as a saturable Emax branch on HR and SV,
  off by default.**  The doc/05 4.3 design note observed that the ERK readout
  cannot encode sympathetic *suppression* (it is a stimulatory downstream-of-β-agonist
  proxy), so a beta-like bradycardia/negative-inotropy drug was not
  representable.  `CardiacParams.sympathetic_tone` (default 1.0) now encodes
  the branch: `sympathetic_tone_from_emax()` maps a free heart exposure to
  residual tone `IC50/(IC50 + C)`, wired from the PBPK heart free compartment
  when `RunSpec.beta_block_ic50_nm` is set, and `simulate_hemodynamics()`
  multiplies it into both heart rate and stroke volume (CO ~ tone²).  Systemic
  resistance stays pinned to the *intact-tone* reference output, so the 
  suppression reads out as hypotension on the two-compartment Windkessel:
  `pa − pv = tone²·(MAP_ref − CVP)` at steady state, feeding the loop
  `co_fraction` and organ perfusion scaling.  Pinned in G4 by
  `case_cardiac_sympathetic_suppression`: IC50 exposure quarters CO exactly,
  zero-exposure/off-by-default baseline, saturating exposure falls to the CVP
  floor, monotone exposure-response, and degenerate-input rejection.  ERK tone
  scales `inotropy/chronotropy` and the sympathetic tone is a *separate*
  multiplicative axis, so the validated stimulatory coupling is unchanged.

- **D14 — Immune-mediated DILI as a saturable adaptive-response hazard,
  off by default.**  Doc/05 4.2 item 5 asked for immune-mediated
  hepatocyte killing.  ``LiverParams.immune_ic50_nm`` now
  supplies a Hill-sigmoid hapten/danger hazard (``immune_hazard()``) from
  free hepatic exposure; the hazard drives a recruitment/decay adaptive
  immune-response ODE (``dI = k_recruit·hazard·(1−I) − k_decay·I``, second
  liver state) whose output loads the fourth ``combined_stress()`` axis via
  ``immune_weight`` (default 0).  The steady-state ``I_ss =
  k_recruit·hazard/(k_recruit + k_decay)`` is physically intuitive (fast
  recruit, slow decay → persistent signal) and pinned analytically by the
  L2 case.  Off by default (``immune_ic50_nm=None``, ``immune_weight=0``)
  so every validated cholestatic/mitochondrial/redox composition is
  reproduced exactly.  Wired from ``RunSpec.dili_immune_ic50_nm`` /
  ``dili_immune_weight`` with a sensible default ``immune_weight=0.5``
  when the IC50 anchor is provided; ``ExposureProfile.dili_immune_ic50_nm``
  surfaces the anchor in the contract.

- **D15 — ACAT-lite multi-segment SI dissolution, off by default.**
  Doc/05 §1.6 noted "no multi-segment small-intestinal dissolution
  resolution."  ``AbsorptionParams.si_segments`` now splits the single SI
  compartment into N sequential equal-volume sub-compartments when set
  (``si_segments`` defaults to ``None``, single lump).  Each sub-compartment
  has its own dissolution cap (``solubility_mg_ml * gi_volume_ml / N`` mg)
  and first-order ``k_si_absorption``; per-segment transit is scaled to
  ``k_si_transit * N`` so total SI transit time (``1/k_si_transit``) is
  preserved.  Bile secretion enters segment 0 (proximal SI).  Off by default
  so the validated single-compartment ACAT-lite is reproduced exactly.
  Wired from ``AbsorptionParams.si_segments`` (no ``RunSpec`` passthrough;
  direct model-construction parameter).

## 5. Gate-keeping notes

- Every data/ file above is sha256-pinned (`data/manifest.json`); adding any
  new file without a manifest entry fails the checksum gate in
  `scripts/release.py`.
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
  baseline), **R-7** (de Bruijn & Rietjens bile-acid PBK cholestasis ranking),
  **D9-L2** (`case_clearance_mechanisms` — extended Stage-1 clearance/
  absorption mechanics), **R-8** (`case_herg_calibration` — calibrated hERG
  head vs corpus, doc/12 D10), **D-L2** (`case_cyp_kinetics` — per-CYP
  MM/Hill clearance), **D-L2** (`case_cheng_prusoff_conversion` — IC50→Ki),
  **TMDD-L2** (`case_tmdd_drug_disposition` — native target-mediated drug
  disposition mass-balance coupling), **SKIN-L2** (`case_transdermal_multi_layer`
  — finite-dose multi-layer skin-permeation membrane), **SYMP-L2**
  (`case_cardiac_sympathetic_suppression` — saturable Emax beta-like
  sympathetic-suppression branch on HR/SV with CO ~ tone², hypotension on the
  fixed-resistance Windkessel, monotone exposure-response, saturating CVP
  floor).  **DILI-L2** (`case_dili_immune_activation` — adaptive immune-response
  DILI QST: saturable hapten hazard driving a recruitment/decay ODE,
  steady-state I_ss analytic, immune_weight inert at 0, monotone
  dose-response, degenerate-input rejection).  **ACAT-L2**
  (`case_acat_multisegment_si` — ACAT-lite multi-segment SI dissolution/
  absorption: N sequential equal-volume SI sub-compartments with per-segment
  dissolution caps, bile enters segment 0, per-segment transit scaled to
  ``k_si_transit * N`` (total transit time preserved); model-layer
  default-off (auto-engaged in full fidelity),
  mass conservation, solubility cap, segment-count sensitivity, degenerate
  rejection).  **DTI-L2** (`case_dti_resolver_calibration` —
  leave-one-molecule-out calibration of the ChEMBL kNN resolver on the
  vendored snapshot: per-target positive correlation, honest scatter, band-
  clamped non-extrapolation; plus the gated re-bind / CNS-anchor wiring claims
  on haloperidol / dofetilide).  **PHYS-L2**
  (`case_willmann_allometric_physiology` — 70 kg reference self-consistency
  vs the Ye 2016 tables, Willmann 2007 allometric power laws, sex factors,
  Mosteller BSA, cardiac-output scaling, impairment modifiers).
  **NEPH-L2** (`case_nephron_tubular_transport` — tubular reabsorption /
  secretion handling on the nephron, feeder of the engaged renal lane).
  **ORD-L19** (`case_ord_multi_ionic_qt` — ORd multi-ionic QT block from the
  CiPA-v1 2017 retune).  **MITO-L2** (`case_mito_redox_dili` — mitochondrial
  ETC / redox / hepatocyte-death axes flanking the cholestasis anchor).
  **PK-GMFE-L25** (`case_pk_gmfe_endpoints` — endpoint-band GMFE/APE over
  published PK parameters: 5 compounds, 13 pairs, pooled GMFE 1.296, APE
  58.0 %, 12/13 within 2×).  All 40 analytic cases (with the 5 Tier-1
  benchmarks, 45 total) run in the G4 gate.
- **E2E** (`case_full_chain_admet_to_report` — end-to-end full-chain
  verification: structure → ADMET-AI → `spec_from_admet` → `run_pipeline` →
  `to_contract` → `render_all`/`write_report`, two unrelated molecules;
  asserts contract determinism, the eight-section contract shape, Cmax > 0,
  `Fa ∈ (0, 1]`, the per-molecule trust record (anchors, G5 policy, class
  priors, mechanism terms engaged in full fidelity), and that all report
  artifacts are written.  Real ADMET-AI runtime: missing runtime = hard FAIL,
  never a silent fallback.  Fast (fake-predictor) twin runs in the G3 unit
  suite as `tests/test_full_chain.py`).

## 6. Shipped-Item Register

Every shipped baseline default and pipeline anchor, with its current state and
the evidence that certified it.  Status legend for this register:

- **DONE (DEFAULT-ENGAGED)** — implemented, validated by an L2 case, and
  engaged in every `fidelity="full"` run (the mandated default) with a
  disclosed auto-anchor in the trust record's `estimates`; the linear lane is
  reachable only by an explicit `fidelity="baseline"` opt-out, and a full run
  that cannot engage a term raises `RealismError` instead of degrading
  silently (G5).
- **WIRED (live)** — active in the default path or in the validation gate and
  disclosed in every contract's `trust` record.
- **WIRED** — vendored/installed, actively computed in the validation suite
  and gated (sha256 + case).
- **BY DESIGN** — a deliberate shipped-baseline choice, enforced by its case.

> **Full-fidelity rollup note.** The decision records (D9–D15) below describe
> each feature at the point it shipped as a linear-lane opt-in.  Since the
> full-fidelity rollup (`RunSpec.fidelity`, default `"full"`), the ledger rows
> L1–L11 are **DEFAULT-ENGAGED** with disclosed auto-anchors; the validated
> linear lane survives as the explicit `fidelity="baseline"` opt-out, and a
> full run that cannot engage a term raises `RealismError` (G5, no silent
> degradation).  Case titles/metrics that still say "off by default" refer to
> the *model-layer* default (term absent in the `PBPKModel`), which the
> pipeline auto-engagement passes over in full runs.

| # | Baseline / shipped item | Documented | Status | Admission evidence |
|---|---|---|---|---|
| L1 | Active tubular secretion (`cl_sec_l_h`) | doc/05 §1.4, doc/12 D9 | **DONE (DEFAULT-ENGAGED)** — `pk/pbpk_build.py`; `case_clearance_mechanisms` (analytic single-pool secretion urine fraction); full runs auto-anchor `cl_sec = CL_renal·(0.5 + 0.5·Pgp)` | extended by the nephron tubular-transport lane (L17) |
| L2 | Saturable MM hepatic clearance (`hepatic_vmax_mg_h`/`hepatic_km_mg_l`) | doc/05 §1.3, doc/12 D9 | **DONE (DEFAULT-ENGAGED)** — low-dose `Vmax/Km` slope pinned by `case_clearance_mechanisms`; full runs auto-anchor `Vmax = CL_h·Km` (Km = 1 mg/L prior) | — |
| L3 | Per-CYP abundance-scaled MM/Hill kinetics (`cyp_terms`) | doc/05 §1.3 + Barter table note, doc/12 D9 | **DONE (DEFAULT-ENGAGED)** — `case_cyp_kinetics` (abundance-scaled Vmax, isoform additivity, doubling slope, Hill, parameter rejection); the MM auto-anchor yields to explicit per-CYP terms | — |
| L4 | First-pass gut-wall extraction (`gut_extraction_eg`) | doc/05 §1.6, doc/12 D9 | **DONE (DEFAULT-ENGAGED)** — `F = Fa·(1−Eh)·(1−Eg)` pinned; full runs auto-anchor `Eg = 0.10 + 0.30·P(CYP3A4 inhibition)` | — |
| L5 | Biliary excretion + enterohepatic recirculation (`cl_bil_l_h`) | doc/05 §1.3, doc/12 D9 | **DONE (DEFAULT-ENGAGED)** — mass-conservative EHC pinned; full runs auto-anchor `cl_bil = CL_h·fraction(logP)` | — |
| L6 | Native TMDD mass-balance coupling (`target_binding`) | doc/05 §2.4, doc/12 D11 | **DONE (DEFAULT-ENGAGED)** — `case_tmdd_drug_disposition` (mass closure, irreversible sink, dose-disproportional retention, KD recovery, AUC reduction, degenerate rejection); full runs auto-engage the primary-affinity site (heart for hERG-name sites, else liver) | — |
| L7 | Finite-dose multi-layer skin permeation (`skin_layers`) | doc/05 §1.4, doc/12 D12 | **DONE (DEFAULT-ENGAGED)** — `case_transdermal_multi_layer` (Fick steady flux, partitions, barrier/diffusivity, mass closure); full runs auto-engage 4-layer membrane on the transdermal route | — |
| L8 | Sympathetic-suppression cardiac branch (`beta_block_ic50_nm`) | doc/05 4.3, doc/12 D13 | **DONE (DEFAULT-ENGAGED)** — `case_cardiac_sympathetic_suppression` (CO ~ tone², hypotension, monotone, CVP floor); full runs auto-anchor a disclosed null-effect IC50 for non-blockers | — |
| L9 | Immune-mediated DILI axis (`dili_immune_ic50_nm`/`dili_immune_weight`) | doc/05 4.2 (former design item), doc/12 D14 | **DONE (DEFAULT-ENGAGED)** — `case_dili_immune_activation` (hapten hazard ODE, I_ss analytic, monotone, degenerate rejection); full runs auto-anchor `0.30·DILI IC50` (or disclosed null) | — |
| L10 | ACAT multi-segment SI dissolution (`si_segments`) | doc/05 §1.6, doc/12 D15 | **DONE (DEFAULT-ENGAGED)** — `case_acat_multisegment_si` (per-segment caps, bile at segment 0, transit scaling, mass conservation); full runs auto-engage 3 segments | — |
| L11 | Organ-feedback outer loop (clearance coupling) | doc/08 risk 5, `RunSpec.feedback_loop` | **DONE (DEFAULT-ENGAGED)** — coupled-clearance driver; full runs set one coupled re-run (pathway stays single-pass; the loop uses a fast liver/kidney/cardiac pass) | — |
| L12 | CNS partitioning & grading anchoring | doc/07, R-4 | **WIRED (live)** — `case_admet_bbb_cns` pins the BBB_Martins kpu decision; the CNS *endpoint* line is structure-anchored only where chemotype support exists (`drugos.target.dti.cns_ic50_nm`), else the disclosed 0.20 class prior, masked from grading; `cns_grading_anchored` disclosed per run | haloperidol honored sub-µM (DRD2 top-match 0.67); desert queries stay at the prior |
| L13 | Non-hERG off-target resolution (panel) | doc/10, doc/12 D10 | **WIRED (live)** — head lanes re-score hERG + CYP2D6/3A4/2C9 inhibitors through the corpus-calibrated P→KD curve; the ChEMBL kNN resolver (row 2b) re-binds mapped sites with chemotype support; without support the site keeps its class prior, named in `no_public_data_sites` | LOO calibration (`case_dti_resolver_calibration`): per-target positive correlation, tight scatter, band-clamped |
| L14 | Willmann allometric physiology | doc/10, doc/12 row 1c / D4 | **WIRED** — `drugos.pk.physiology.build_human`/`HumanPhysiology`; `case_willmann_allometric_physiology` pins 70 kg reference self-consistency (Ye 2016), per-organ Willmann (2007) exponents, sex factors, BSA scaling, CO scaling, impairment modifiers | PHYS-L2 (gate) |
| L15 | ChEMBL kNN general DTI resolver | doc/12 row 2b, D10 | **WIRED** — `predict_kd_nm`/`top_tanimoto` with the `MIN_NEIGHBOR_TANIMOTO` evidence gate; `resolve_offtarget_panel` re-binds supported sites, `no_public_data_sites` names the rest; never invents potency beyond the recorded band | leave-one-molecule-out calibration aboard the vendored snapshot (DTI-L2, gate) |
| L16 | Mitochondrial/redox DILI axes | doc/12 row 4c / D7 | **WIRED** — `organ/liver.py` mitochondrial ETC (adaptive mitogenesis), redox/GSH, hepatocyte-death ODEs flanking the cholestasis anchor; `case_mito_redox_dili` pins ETC block ↔ redox depletion ↔ ALT-graded death | MITO-L2 (gate) |
| L17 | Nephron tubular transport | doc/05 §1.5, doc/12 row 4b | **WIRED** — `organ.nephron.nephron_handling`: proximal reabsorption / transporter-mediated secretion kinetics, feeder of the engaged renal-elimination lane; `case_nephron_tubular_transport`; GFR baseline anchored via CKD-EPI 2021 (R-6) | NEPH-L2 (gate) |
| L19 | ORd main-path multi-ionic QT block (QTw) | doc/12 D2, doc/11 | **WIRED** — CiPA-v1 2017 retune vendored (`data/models/ohara-cipa-v1-2017.mmt`) as `organ.cardiac_ap.cipa_apd90`; the ORd lane covers ICaL/INaL/IKs/IK1 block; `case_ord_multi_ionic_qt` spans the R-3 axis and the CiPA rank corollary | ORD-L19 (gate) |
| L25 | Published-PK endpoint GMFE benchmark | doc/12 §5 PK-GMFE-L25; `data/benchmarks/published_pk.json` | **WIRED (live)** — `case_pk_gmfe_endpoints` certifies endpoint geometry over published PK parameter bands (CL, t½, Fa, fe, Vss): 5 compounds, 13 compound–metric pairs across the L3 Tier-1 set; pooled GMFE 1.296, APE 58.0 %, 12/13 within 2× | band check per pair; single disclosed outlier acetaminophen t½ (2.5×) carried in the case verdict |
| L26 | GFR-only renal elimination | doc/05 §1.5, DISCLAIMER §2 | **BY DESIGN (baseline lane)** — `fidelity="baseline"` keeps linear-renal; full fidelity auto-engages tubular secretion (`cl_sec = CL_renal·(0.5 + 0.5·Pgp)`, disclosed) | the nephron lane (L17) adds depth to the engaged secretion term |
| L27 | Predictive-regime reliability + measured/empirical overrides | DISCLAIMER §2, doc/07 | **DONE (DEFAULT)** — `drugos/reliability.py`: six regimes (weakest axis dominates) disclosed as `trust.reliability` (regime/label/reliability/basis/band_cv/disclaimer) on every run; regime CV feeds the D21 ensemble (`EnsembleConfig.cv`); `--measured`/`Measurements` override true parameters into the same full-fidelity pipeline (full PK → `measured_in_range_on_label`; measured zeros stay honest, engaged-and-disclosed terms); `--empirical`/`EmpiricalObservations` report observed-vs-predicted fold-error / within-2x rows in `trust.empirical_agreement` (disagreement = expected state, not a bug); `case_predictive_regime` (L1) | measured `fa`/`logP`/organ-flux scalars widen the mid-evidence surface; full-PK measured + on-label requires a scaffold anchor on dose/route |

## 7. Full-Chain Trust & End-to-End Verification Mechanism

Goal (doc/05 design intent, made executable): *a run on an arbitrary drug is
simulated end-to-end — structure → ADME/T → spec → PBPK → binding → pathway →
organ → clinical — and every admitted source is auditable in the report.*  A
full-chain run must not silently borrow a class prior, a default, or a
head-prediction without saying so.

### 7.1 Per-run trust record (`fidelity_provenance`)

`drugos.pipeline.fidelity_provenance(spec, exposure)` is a pure function whose
return value is emitted as the **`trust`** section of the JSON contract
(`RunResult.to_contract()`) and rendered into the report.  Fields:

| Field | Meaning |
|---|---|
| `fidelity` | `"full"` (default: every realism term engaged with disclosed auto-anchors) or `"baseline"` (explicit disclosed opt-out that keeps the validated linear lane) |
| `policy` | G5 no-silent-fallback statement: every admitted source is named in this record |
| `anchors_wired` | Production anchors actually bound into this run: CKD-EPI 2021 GFR (R-6) and the GCDCA bile-acid cholestasis PBK (R-7) are always wired; the ADMET-AI BBB_Martins brain-partition head (R-4) is disclosed whenever an ADMET prediction is present; the corpus-calibrated hERG P→KD sieve (R-8) whenever the hERG head is present |
| `cns_grading_anchored` | Whether the CNS *endpoint* line is anchored by a chemotype-supported `cns_ic50_nm` (else 0.20 class prior, masked from grading) |
| `mechanism_terms_engaged` | Which realism terms are active for this run (tubular secretion, biliary/EHC, gut-wall extraction, saturable MM hepatic, per-CYP kinetics, TMDD, skin layers, sympathetic branch, immune DILI, feedback loop, ACAT SI, SC/IM depot — see §6 L1–L11); non-empty in every `fidelity="full"` run |
| `mechanism_terms_degraded` | Which terms could **not** be engaged and why (route-inapplicable or missing-route terms in a baseline run). Empty in a full run — a full run whose terms would be missing instead raises `RealismError` (no silent degradation, G5) |
| `estimates` | Per-term auto-anchor basis lines (sorted, term → basis) for the values `_engage_full_fidelity` synthesized (MM Vmax/Km, tubular secretion, biliary fraction, gut-wall Eg, ACAT segments, feedback loop, SC/IM ka, skin layers, immune-DILI IC50, sympathetic null IC50, TMDD site) |
| `no_public_data_sites` | Target sites resolved by class priors because the chemotype-desert / corpus-coverage gap gives them no structure support, sorted |
| `admet_ml_estimates` | Which ADMET-AI heads contributed (all non-`raw` populated fields), or `None` if the run is ADMET-free |
| `reliability` | Predictive-regime disclosure (DISCLAIMER §2, §7.3): `regime`, `label`, `reliability`, `basis`, `band_cv` (recommended parameter-ensemble CV), `disclaimer` |
| `empirical_agreement` | Present when observed clinical data are supplied (`EmpiricalObservations`): per-endpoint `fold_error` / `within_2x` rows plus a written policy that disagreement with empirical data is the expected state of the mechanistic model, not a bug |

### 7.2 End-to-end gate

- **G4 case `case_full_chain_admet_to_report`** (L1) drives the real chain on
  two unrelated molecules (acetaminophen, caffeine): `parse_structure` →
  `predict_admet` (real ADMET-AI; missing runtime = hard FAIL) →
  `spec_from_admet` → `run_pipeline` → `to_contract` → `render_all` +
  `write_report`.  Metrics: contract determinism on a re-run, eight-section
  contract shape, `Cmax > 0`, `Fa ∈ (0,1]`, anchors disclosed, G5 policy,
  class priors disclosed, mechanism terms engaged in full fidelity
  (`fidelity == "full"`, `mechanism_terms_degraded == []`), report artifacts
  written.
- **G3 twin `tests/test_full_chain.py`** runs the same chain against a
  deterministic fake predictor so every fidelity branch (default panel,
  full-realism terms, mixed-confidence priors, ADMET-free run, BBB-anchored vs
  masked CNS) is branch-covered inside the unit suite.
- The trust record is branch-covered in both lanes; any anchor admitted to
  `anchors_wired` arrives with its own L2 equivalence case (G4 factory rule).

### 7.3 Predictive-regime reliability & true-parameter overrides (DISCLAIMER §2)

The clause — *"predictions for novel molecules, extrapolated doses, or
off-label routes are the least reliable; disagreement with empirical data is
the expected state of a mechanistic model, not a bug"* — is made executable by
`drugos/reliability.py`.  Every run is classified on three evidence axes and
the **weakest axis dominates** (a run never reports more confidence than its
least supported axis):

| Axis | Evidence |
|---|---|
| Chemistry | likely. Measured full PK (fup + hepatic + renal clearance) = strongest; some measured PK/absorption = middle; a name-anchored benchmark scaffold solidifies only if corroborated by measurements; otherwise machine-predicted chemistry = weakest |
| Route | On-label = the scaffold's validated route; off-label = any other route |
| Dose | Inside the scaffold's validated reference window `[0.5×, 2×]` vs outside it |

Six regimes, weakest → strongest, each with a recommended parameter-ensemble
CV (`band_cv`) that widens when evidence weakens: `novel_molecule` (0.60),
the mid-evidence chemistry regime (0.45), `validated_offlabel_route` (0.50),
`validated_extrapolated_dose` (0.35), `validated_in_range_on_label` (0.20),
`measured_in_range_on_label` (0.15).  The regime CV flows into the D21
parameter-ensemble uncertainty stage when `EnsembleConfig.cv` is left at
its `None` default, keeping disclosed reliability and the reported
uncertainty band coupled.

**Measured true parameters.** `--measured {…}` (or `/api/run` `measured`
payload, or the `Measurements` dataclass) overrides any subset of `fup`,
`cl_hep_l_h`, `cl_renal_l_h`, `cl_sec_l_h`, `cl_bil_l_h`, `fa`,
`hepatic_vmax_mg_h`, `hepatic_km_mg_l`, `qt_ic50_nm`, `dili_ic50_nm`,
`dili_immune_ic50_nm`, `cns_ic50_nm`, `beta_block_ic50_nm`.  Overrides are
applied at the top of `run_pipeline` (before the full-fidelity gate) and the
run uses them directly — never a synthesized stand-in — with each override
credited in `estimates` as a "measured true-parameter override".  A measured
zero (e.g. `cl_sec_l_h = 0`, no active tubular secretion) is an honest
determination: the term stays engaged and is disclosed "(measured)" rather
than skipped.  Full measured PK upgrades the run to the most reliable
`measured_in_range_on_label` regime; the DISCLAIMER target is then literally
achieved: given all true parameters of a drug and a human body, the model
predicts the real effects at full fidelity.

**Empirical agreement.** `--empirical {…}` (or `/api/run` `empirical`
payload, or `EmpiricalObservations`) supplies observed `plasma_cmax_mg_l`,
`auc_last_mg_h_l`, `peak_delta_qtc_ms`, `peak_alt_uln`.  The run reports each
one as an observed-vs-predicted `fold_error` + `within_2x` row in
`trust.empirical_agreement`, alongside the written policy that disagreement
with empirical data is the expected state of the mechanistic model — it is
reported, never synthesized away.

Validation: `case_predictive_regime` (L1) certifies disclosure determinism,
the band-CV ordering, the measured-override upgrade being honoured by the run,
the empirical-agreement rows, and the regime CV reaching the uncertainty
stage.  Clinical accuracy itself remains the job of the L3 Tier-1 benchmark
cases.

## 8. References (additions beyond doc/02)

- O'Hara T, Virág L, Varró A, Rudy Y. *Simulation of the undiseased human
  cardiac ventricular action potential...* PLoS Comput Biol 2011;7(5):e1002061.
- Dutta S, Chang KC, et al. (Li Z). *Optimization of an in-silico cardiac cell
  model for proarrhythmia risk assessment.* Front Physiol 2017;8:616.
- Swanson K, Walther P, et al. *ADMET-AI: a machine learning ADMET platform.*
  Nat Mach Intell 2024.
- Du F, et al. *hERG Central.* J Chem Inf Model 2022 (Harvard Dataverse
  doi:10.7910/DVN/7BVDG8).
- Wilhelms M, et al. / CiPA in-silico guidance for multi-channel ORd use.
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
- Willmann S, et al. *Physiological model for simulating gastrointestinal
  transit and drug absorption in rats and humans.* J Pharmacokinet
  Pharmacodyn 2003;30:109-24 (allometric power laws embedded by the
  PK-Sim population module — row 1c, D4).