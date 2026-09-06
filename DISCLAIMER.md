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

## 2. Limited scope and fidelity

- The baseline (2026.9.0) deliberately simplifies biology: lumped hepatic
  clearance, ACAT-lite absorption, a single MAPK-cascade pathway kernel, a
  static off-target panel with class-typical priors, and GFR-only renal
  elimination (see `doc/05`, `doc/08`). Real drugs interact with biology in
  ways this baseline does not represent.
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