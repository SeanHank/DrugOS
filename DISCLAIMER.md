# Disclaimer

**DrugOS is a research-grade mechanistic modeling platform. It is not a
medical device, not a clinical-decision-support product, and not validated for
regulatory or diagnostic use.**

By using DrugOS — whether as source code, an installed package, a CLI tool, a
web application, or any report it produces — you acknowledge and agree to the
following:

## 1. No medical, clinical, or regulatory product

- DrugOS output (risk verdicts, grades, trajectories, uncertainty bands,
  virtual-cohort incidence, or any derived quantity) is a **simulation of model
  assumptions**, generated from a baseline parameterization. It is not a
  reproduced clinical trial, not an individual patient assessment, and not a
  substitute for physician judgment.
- Nothing in DrugOS is FDA-, EMA-, or otherwise regulatory-cleared. It must not
  be used to support dosing decisions, patient management, label claims,
  submissions, or any activity where a wrong prediction could harm a person or
  cost money on the strength of the prediction.

## 2. Model composition and fidelity

- Where production-validated, open-source models exist, DrugOS integrates them
  as first-class components rather than re-deriving them: the
  **O'Hara–Rudy 2011** human ventricular action-potential model (BSD-3, vendored
  and sha256-pinned), the **Huang & Ferrell 1996 ultrasensitive MAPK cascade**
  (CC0 SBML, BioModels BIOMD0000000009, vendored and sha256-pinned as the
  Stage-3 pathway default), the **CKD-EPI 2021 race-free creatinine equation**
  (open-published kidney GFR baseline, R-6), the **de Bruijn & Rietjens 2024
  GCDCA bile-acid PBK** (open-published cholestasis anchor — equations ported
  from the CC BY 4.0 paper; the authors' CC-BY-NC-ND repository code is not
  copied, R-7), **ADMET-AI** structural ADME/T predictions (BSD-3
  code, CC-BY weights), **measured hERG dose-response corpora** (ChEMBL / hERG
  Central, CC0/CC), and a **required R implementation** of the literature PK
  estimators (Wagner 1976, Gibaldi & Perrier 1982, Rowland & Tozer 2010; see
  `doc/12-production-models.md`). These assets keep their upstream licenses and
  are pinned by checksum in `data/manifest.json`.
- Every run emits a machine-readable **trust record** (`fidelity_provenance`,
  `doc/12 §7`) naming the run's `fidelity` lane (the default `full` engages
  every realism term with a disclosed auto-anchor, `baseline` is the explicit
  disclosed opt-out), the terms that could not engage
  (`mechanism_terms_degraded`, empty in full — full runs that cannot engage a
  term raise `RealismError`), the per-term auto-anchor basis lines
  (`estimates`), exactly which production anchors are wired into that run
  (CKD-EPI GFR R-6, de Bruijn & Rietjens bile-acid cholestasis PBK R-7,
  ADMET-AI BBB_Martins R-4, corpus-calibrated hERG P→KD sieve R-8), which
  ADMET-AI heads contributed, which target sites still fall back to class
  priors, and that no term outside the record contributes to any readout
  (P5–P9 are planned release tracks in `doc/07`; none is a runtime entry). A
  run never reports more confidence than its disclosed composition supports;
  the G5 no-silent-fallback rule is part of that contract.
- An end-to-end validation case drives the real full chain per run —
  structure → ADMET-AI → `spec_from_admet` → pipeline → JSON contract →
  rendered report — on multiple unrelated molecules, and fails the gate if
  any contract section, disclosure or report artifact is missing or
  non-deterministic (`doc/12 §7.2`, `case_full_chain_admet_to_report`).
- The pipeline runs **full fidelity by default** (`fidelity="full"`): every
  realism extension is engaged and auto-anchored with disclosed estimates
  in the trust record — active tubular secretion, saturable Michaelis–
  Menten hepatic metabolism, per-CYP abundance-scaled kinetics, first-pass
  gut-wall extraction, biliary excretion with enterohepatic recirculation,
  native TMDD target binding, finite-dose multi-layer skin permeation,
  sympathetic-suppression cardiac branch, immune-mediated DILI axis, ACAT
  multi-segment SI dissolution, and the organ-feedback loop — each validated
  by its own G4 case.  The validated linear lane is available via an explicit
  `fidelity="baseline"` opt-out; a full-fidelity run that cannot engage a
  term raises an explicit `RealismError` rather than degrading silently (G5).
  Off-by-default means every shipped default run reproduces the benchmarked
  baseline exactly (`doc/12 §6` L1–L11).
- The shipped lane realizes the validated baselines directly: class-median
  priors for non-hERG panel sites, GFR-based renal elimination as the renal
  default (the deeper nephron model is a planned release in `doc/07` P7),
  lumped-compartments ACAT-family intestinal absorption as the absorption
  default, and the in-house Mito/redox/cell-death liver axes that flank the
  validated cholestasis anchor — each implemented, engaged and covered by its
  own G4 case, and none a retained stand-in. The G6 gate
  (`scripts/release.py gates` / `marker-audit`) fails any marker in the set at
  any count in any scanned code or doc file, keeping no allowance file, no
  whitelist and no ledger of permitted items. The audit is read-only over the
  citation/license surface — the attributed sources named throughout these
  documents and in `doc/*` are never added to or rewritten by it.
- Kidney, liver and CNS subsystems are the in-house QST implementations today,
  with the validated anchors wired into their organ readouts; CNS
  brain-partition verdicts are decided by the ADMET-AI BBB_Martins head
  whenever a prediction is present (R-4), while CNS *endpoint* grading is
  anchored only by an explicit `cns_ic50_nm` (otherwise the 0.20 class prior,
  disclosed). The license-restricted external engines (DILIsym, CMR nephron,
  Physiome/Reactome scaffolds, PK-Sim physiology) are never linked into the
  runtime; their workflows are addressed by the in-house implementations
  documented in `doc/12`, and each planned release-track engine in `doc/07`
  P5–P9 is gated by an L2 equivalence case before admission (G4).
- Real drugs interact with biology in ways any finite model — including the
  integrated production-validated components — does not represent. Integrating
  a trusted model does not make the integrated prediction trustworthy.
- Values reproduced from `admet_ai`, R, or the vendored corpora are returned
  as-is; their accuracy is governed by the upstream artifacts and licenses,
  not by this project.
- Predictions for novel molecules, extrapolated doses, or off-label routes are
  the least reliable. Disagreement with empirical data is the expected state of
  a mechanistic model, not a bug.
- Every run is classified into a **predictive regime** (weakest evidence axis
  dominates) and the trust record carries that disclosure
  (`trust.reliability`): `novel_molecule` < `partial_evidence` <
  `validated_offlabel_route` / `validated_extrapolated_dose` <
  `validated_in_range_on_label` < `measured_in_range_on_label`. The regime
  sets a recommended parameter-ensemble CV, so a weakly-evidenced run reports
  a correspondingly wider uncertainty band, never a tighter one (`doc/12
  §7.3`, `reliability.py`).
- Measured **true parameters** can be supplied as overrides (`--measured`
  JSON, `Measurements`): a measured parameter set is used by the exact same
  full-fidelity pipeline — never a synthesized stand-in — and a full measured
  PK set upgrades the run to the most reliable measured regime. A measured
  zero (e.g. no active tubular secretion) is an honest determination, engaged
  and disclosed as measured rather than skipped (G5).
- Observed clinical data can additionally be supplied (`--empirical`,
  `EmpiricalObservations`); disagreement with those observations is reported
  per-endpoint as fold-error / within-2x rows in the trust record
  (`trust.empirical_agreement`) with a written policy that disagreement is the
  expected state of the mechanistic model — it is never synthesized away.
  This makes the clause above executable rather than aspirational: the
  goal of the platform is, given all true parameters of a drug and a human
  body, to predict the real effects of that drug, with reliability and
  disagreement always disclosed.

## 3. No warranty

- DrugOS is provided **AS IS**, **WITHOUT WARRANTY OF ANY KIND**, express or
  implied, including but not limited to the warranties of merchantability,
  fitness for a particular purpose, and non-infringement.

## 4. Limitation of liability

- In no event shall the project authors, contributors, or copyright holders be
  liable for any claim, damages, or other liability arising from, out of, or in
  connection with the use of DrugOS or any output it produces — whether in an
  action of contract, tort, or otherwise.

## 5. License

- DrugOS is licensed under the **GNU Affero General Public License v3.0
  (AGPL-3.0)**. Nothing in this disclaimer grants you additional rights beyond
  that license; everything in this disclaimer survives the license terms.

**If your use case requires validated, regulator-approved prediction, do not
use this software.**