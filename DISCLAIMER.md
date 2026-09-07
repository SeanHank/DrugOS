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
- The pipeline still deliberately simplifies biology where no downloadable
  production model is wired yet: lumped hepatic clearance, ACAT-lite absorption,
  class-typical priors for non-hERG targets, GFR-only renal elimination, and the
  liver Mito/redox/cell-death axes that flank the validated cholestasis anchor
  (see `doc/05`, `doc/08`, `doc/12`).
- Kidney, liver and CNS subsystems are the in-house QST implementations today,
  with the validated anchors (CKD-EPI GFR baseline, bile-acid cholestasis PBK)
  wired into their organ readouts; CNS brain-partition verdicts are decided by
  the ADMET-AI BBB_Martins head when a prediction is present (R-4). The named
  production models catalogued in `doc/12` (DILIsym, CMR nephron,
  Physiome/Reactome scaffolds) are the remaining P5–P9 integration targets.
- Real drugs interact with biology in ways any finite model — including the
  integrated production-validated components — does not represent. Integrating
  a trusted model does not make the integrated prediction trustworthy.
- Values reproduced from `admet_ai`, R, or the vendored corpora are returned
  as-is; their accuracy is governed by the upstream artifacts and licenses,
  not by this project.
- Predictions for novel molecules, extrapolated doses, or off-label routes are
  the least reliable. Disagreement with empirical data is the expected state of
  a mechanistic model, not a bug.

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