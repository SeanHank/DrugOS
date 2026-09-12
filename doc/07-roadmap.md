# System Realization & Shipped Subsystems

The five-stage pipeline — `drug ──▶ in-vivo concentration ──▶ target binding
──▶ signaling pathway ──▶ organ function ──▶ clinical phenotype` (doc/05) —
is shipped end to end. Every stage is wired to a production-validated,
open-source anchor and certified by the validation suite (doc/12). The G4
gate runs the 40-case analytic suite plus the five-benchmark GMFE corpus
(45 registered steps), all green, with the generated report in
`validation/report.md`. Each shipped item below is disclosed in the per-run
trust record (`doc/12 §6`/`§7`).

## Stage 1 — In-vivo Concentration & ADME/T

- **ADMET-AI frontend (WIRED).** Stage-1 ADME/T primitives
  (`drugos.pk.admet.ADMETPredictor`) run the Zenodo-archived ADMET-AI
  platform (Swanson et al. 2024, *Nat Mach Intell*; BSD-3 code, CC-BY-4.0
  weights) for the physchem and ADME/T priors.
- **R literature PK cross-check (WIRED + REQUIRED, R-1).** Every pipeline run
  executes the literature PK estimator `rbridge/literature_pk.R` on R ≥ 4.5
  (rpy2 ≥ 3.6) and gates the CL/AUC/t½/Vss agreement ≤ 2 %
  (`case_r_bridge`).
- **Willmann allometric physiology (L14).** `drugos.pk.physiology.build_human`
  resolves organ volumes, flows and tissue composition from the Willmann et
  al. (2007) population allometric power laws against the Ye (2016) reference
  tables; `case_willmann_allometric_physiology` pins reference
  self-consistency, per-organ exponents, sex factors, Mosteller BSA,
  cardiac-output scaling and impairment modifiers.

## Stage 2 — Target Binding

- **Measured hERG + hERG Central corpus (R-2).**
  `data/benchmarks/herg_measured_nm.json` (ChEMBL target CHEMBL240) and the
  306,893-row `data/corpora/herg_central.tsv.gz` (Du et al. 2022, Dataverse
  CC0) drive the corpus calibration and the per-compound measured hERG
  override (`case_corpus_calibration`).
- **Corpus-calibrated hERG sieve (R-8).**
  `drugos.target.resolver.kd_from_score` maps the classifier score onto a
  strictly monotone probability→KD curve that is never more potent than the
  2 nM hERG class prior; `case_herg_calibration` pins the curve invariants,
  the anchor window and the head's concordance with the corpus inhibition
  labels.
- **ChEMBL kNN general DTI resolver (L15).** `drugos.target.dti.predict_kd_nm`
  queries the vendored 16-target ChEMBL snapshot
  (`data/chembl/offtarget_snapshot.json`, 3035 rows) and re-binds a panel
  site to a structure-derived KD only when best-pair Tanimoto ≥
  `MIN_NEIGHBOR_TANIMOTO`; sites without chemotype support keep disclosed
  class priors named in `no_public_data_sites`.
  `case_dti_resolver_calibration` certifies the leave-one-molecule-out
  accuracy and the band-clamped non-extrapolation.

## Stage 3 — Signaling Pathway

- **Huang & Ferrell 1996 SBML MAPK cascade (R-5).** The vendored
  `data/models/huang1996-mapk-cascade.xml` (BioModels BIOMD0000000009, CC0)
  is parsed by python-libsbml through `drugos.pathway.sbml_pathway` and is
  the pipeline-default signal lane with doubly-phosphorylated ERK (PP_K)
  readout; occupancy suppresses the readout monotonically.
  `case_sbml_mapk_validation` pins parse, steady-state and monotone
  inhibition. Non-MAPK scaffolds (e.g. the Rohwer 2000 PTS SBML) certify the
  parser equivalence gate (`case_sbml_scaffold_equivalence`, L18).

## Stage 4 — Organ Function

- **Cardiac repolarisation / QT.** The O'Hara-Rudy 2011 human ventricular AP
  model (`data/models/ohara-2011.mmt`, BSD-3 Myokit) runs through
  `drugos.organ.cardiac_ap` (R-3): baseline APD90 266 ms, dofetilide at its
  measured IC50 ΔAPD90 ≈ 115 ms, monotone block ramp, warfarin control 0.
  The CiPA-v1 2017 retune (`data/models/ohara-cipa-v1-2017.mmt`) is the
  multi-ionic ORd lane (ICaL/INaL/IKs/IK1) pinned by
  `case_ord_multi_ionic_qt`.
- **Sympathetic-suppression branch.** A saturable Emax beta-adrenergic block
  on heart rate and stroke volume (CO ~ tone²), with hypotension read out on
  the fixed-resistance Windkessel; `case_cardiac_sympathetic_suppression`
  pins the quartering, the zero-exposure baseline and the CVP floor.
- **Kidney GFR / nephron.** The CKD-EPI 2021 race-free creatinine equation
  sets baseline GFR whenever serum creatinine is carried
  (`drugos.organ.kidney.ckdepi_2021_egfr`, Mosteller BSA; R-6); the nephron
  tubular-transport lane (`drugos.organ.nephron.nephron_handling`) resolves
  proximal reabsorption and transporter-mediated secretion
  (`case_nephron_tubular_transport`).
- **Liver cholestasis + mitochondrial/redox DILI.** The GCDCA bile-acid PBK
  (`drugos.organ.liver.simulate_gcdca_pbk`, BSEP Ki = IC50/2; R-7) anchors
  the cholestasis axis from the vendored BSEP IC50 anchors
  (`data/models/bsep_shh_ic50_reference.json`). The mitochondrial ETC block,
  redox/GSH and hepatocyte-death ODEs flank the anchor
  (`case_mito_redox_dili`), and the adaptive immune-response hazard loads the
  fourth `combined_stress` axis (`LiverParams.immune_ic50_nm`,
  `case_dili_immune_activation`).

## Stage 5 — Clinical Phenotype

- **ICH E14 / CTCAE ladders + Redfern TdP bands.** `drugos/clinical/` grades
  QTc/TdP, DILI, AKI and CNS against the published phenotype ground truth,
  with the ORd anchor supplying the independent ionic confirmation
  (`case_clinical_grading`, `case_risk_ordering`).

## Full-Fidelity Engagement & Reliability

- **Full-fidelity defaults.** Every realism term (tubular secretion, biliary
  EHC, gut-wall extraction, saturable MM hepatic, per-CYP kinetics, TMDD
  coupling, multi-layer skin, sympathetic branch, immune DILI, feedback loop,
  ACAT multi-segment SI, SC/IM depot) is DEFAULT-ENGAGED in `fidelity="full"`
  runs with a disclosed auto-anchor in `trust.estimates`; a full run that
  cannot engage a term raises `RealismError` rather than degrading silently
  (G5). The validated linear lane remains reachable only through the explicit
  `fidelity="baseline"` opt-out (doc/12 §6 L1–L11).
- **Robustness stages.** `drugos/robustness/` runs the log-normal
  multiplicative parameter ensemble with 90% bands (`uncertainty.py`),
  virtual-patient anthropometric sampling with grade incidence
  (`population.py`), and OAT + Saltelli Sobol driver attribution
  (`sensitivity.py`).
- **Predictive-regime reliability (L27).** `drugos/reliability.py` classifies
  every run into six regimes with the weakest evidence axis dominating,
  disclosed as `trust.reliability` (regime/label/reliability/basis/band_cv/
  disclaimer); the regime CV feeds the parameter-ensemble uncertainty stage.
  Measured true-parameter overrides (`--measured`/`Measurements`) and
  observed-data comparison (`--empirical`/`EmpiricalObservations`) enter the
  same full-fidelity pipeline, with `trust.empirical_agreement` reporting
  observed-vs-predicted fold-error and within-2x rows;
  `case_predictive_regime` certifies the disclosure determinism, band-CV
  ordering, override handling and CV routing.

## Web Service & CNS Panel

- **Flask playground (WIRED).** `src/drugos/web/app.py` (+
  `src/drugos/web/templates/index.html`, deep-purple dark theme) exposes
  `/api/run` with `with_uncertainty`, `with_population` and `with_sensitivity`
  toggles, a CNS brain-exposure plot, and robustness panels (90% bands,
  cohort incidence, local/Sobol driver tables).
- **CNS organ panel.** `drugos.organ.cns.py` computes passive brain
  free-exposure; the CNS risk line is drawn only when a chemotype-supported
  `cns_ic50_nm` is present, else the disclosed 0.20 class prior (masked from
  grading); `case_admet_bbb_cns` pins the BBB_Martins kpu decision.

## Full-Chain Verification & Data Registry

- **End-to-end full-chain verification.** `case_full_chain_admet_to_report`
  (G4, real ADMET-AI, two molecules) plus `tests/test_full_chain.py` (G3,
  deterministic fake predictor) run the chain parse → ADME/T → spec →
  pipeline → contract → report with contract determinism and report-artifact
  checks (doc/12 §7.2).
- **Per-run trust record.** `drugos.pipeline.fidelity_provenance` emits the
  fidelity lane, the G5 no-silent-fallback policy, wired production anchors,
  engaged/degraded mechanism terms, auto-anchor estimates, class-prior sites
  and ADMET-AI head contributions into every contract's `trust` section and
  the report (doc/12 §7).
- **Data registry.** Every vendored data/model file is sha256-pinned in
  `data/manifest.json`; regeneration and fetch scripts live under
  `scripts/data/` (doc/11).