# Literature Review

This document reviews the peer-reviewed and community literature that underpins each stage of the DrugOS pipeline. It is organized along the translational chain `drug -> concentration -> target -> pathway -> organ -> phenotype`.

---

## 1. Stage 1: Drug -> In-Vivo Concentration (Pharmacokinetics)

### 1.1 Whole-Body Physiologically Based Pharmacokinetics (PBPK)

**Core reference concept.** Whole-body PBPK models explicitly represent organs/tissues that matter for absorption, distribution, metabolism and excretion (ADME), and they are parameterized from anatomical and physiological information (blood flow rates, volumes of vascular/interstitial/cellular spaces) plus substance-specific physico-chemical and in-vitro properties. This is the methodological foundation for predicting tissue exposure without any prior clinical data (Nestorov 2007; Jones et al. 2009; Rowland, Peck & Tucker 2011).

**Reference software.** The Open Systems Pharmacology (OSP) Suite provides PK-Sim (whole-body PBPK for small molecules and proteins) and MoBi (expert-level model building and QSP). The suite ships:
- A physiological database for humans and common lab animals (mouse, rat, minipig, dog, monkey, beagle, rabbit).
- Generic passive processes (blood-flow distribution) and active processes (enzyme-mediated metabolism) that are automatically handled.
- The OSP PBPK Model Library: validated, evaluation-report-backed models for ~40+ compounds (antibiotics, antivirals, sedatives, cardiovascular agents, monoclonal antibodies, etc.).
- An R interface (`ospsuite` R package) and population simulation (1000+ virtual subjects with inter-individual variability).

**Population modeling.** Willmann et al. (2007) describe a physiology-based whole-body population model generating individual virtual subjects from age/sex/gender/body-composition covariates, used to propagate inter-individual variability into PK predictions. This is the natural mechanism for the "full human parameters" input of DrugOS.

### 1.2 Structure-to-ADME Machine Learning

**ADMET-AI** (Swanson et al. 2024) is a graph neural network (Chemprop-RDKit) trained on 41 ADMET datasets from the Therapeutics Data Commons (TDC). It holds the highest average rank on the TDC ADMET Leaderboard and predicts, from a SMILES string alone:
- Solubility, LogP, logD, plasma-protein binding, half-life, clearance
- BBB penetration, oral bioavailability
- hERG blockade, AMES mutagenicity, hepatotoxicity, DILI, clinical toxicity probability

It reports each prediction as a percentile relative to a DrugBank reference set, enabling contextualized risk interpretation. It is open source and runs batch predictions (1M molecules in ~3.1 h). ADMET-AI provides the "bottom-up" parameterization path when no measured ADME data exists.

**Review basis.** Obrezanova et al. (2022) demonstrated prediction of in-vivo PK parameters and time-exposure curves directly from chemical structure, confirming that structure -> PK-parameter inference is tractable and already used in discovery.

### 1.3 High-Throughput / Mechanistic Absorption Modeling

Simulations Plus ADMET Predictor + GastroPlus use the ACAT (Advanced Compartmental Absorption and Transit) model to estimate fraction absorbed, volume of distribution and systemic clearance from structure, then feed high-throughput PBPK simulation. The ACAT concept (Yu et al. 1996; Agoram et al. 2001) forms the recommended absorption sub-model for the oral route: gastrointestinal tract divided into physiologically meaningful compartments with transit, dissolution and permeability kinetics.

### 1.4 Python-Based PK/PD Tooling

- **SciPy-native ODE stack (selected over PySB/pysb-pkpd)**: PySB provides rule-based construction of PK/PD and QSP/QST ODE models, but its `sympy<1.12` dependency pin conflicts with the modern torch/admet-ai line (`sympy>=1.13.3`). Per the resolution recorded in `06-technology-stack.md`, the baseline implements all ODE systems (compartmental PK, occupancy, pathway) directly with `scipy.integrate.solve_ivp` (LSODA/Radau); any number of compartments is supported.
- **OpenPKPD**: open-source population PK/PD toolkit with NONMEM-style control streams, SAEM/NUTS estimation, VPC, SBML import, PBPK, TMDD and DDE support — a candidate for the downstream pharmacometric/estimation layer.

---

## 2. Stage 2: Concentration -> Target Binding

### 2.1 Target Identification and Pharmacology Data

- **DrugBank**: authoritative target records (UniProt IDs, gene symbols), mechanism of action, binding constants where available; also serves as the reference set for ADMET-AI percentiles.
- **ChEMBL**: large-scale bioactivity data (IC50/EC50/Ki/Kd where available).
- **PDB / AlphaFold**: experimental and predicted structures for the primary and off-target proteins.

### 2.2 AI Binding-Affinity Prediction

Recent reviews (Zhou et al. 2025; Hsu 2025; comprehensive DTI reviews covering 2016-2025: 180 methods analyzed) identify three families:

1. **Structure-based ML** — use docking poses / complex structures as input; deep learning (3D CNN/point clouds/graphs) then scores affinity. Notable for generalizability: CORDIAL (interaction-only deep learning, PNAS 2025) focuses on physico-chemically intuitive interface properties to maintain predictive performance on protein families unseen during training; IPBind (geometric deep learning using interatomic potential between bound/unbound states, arXiv 2025) improves generalization and gives atom-level insight.
2. **Sequence/descriptor-based ML** — e.g., AttentionDTA (attention-based CNN on protein sequence + SMILES), MINDG (integrated learning), usable even without a 3D protein structure.
3. **Physics-based baselines** — molecular docking scoring functions and MD/FEP for high-accuracy, small-batch refinement.

**Recommendation.** For DrugOS, sequence/descriptor-based models (or DrugBank *in-vivo* targets when the drug is approved) initialize the target set; structure-based ML affinity estimates then parameterize binding kinetics. Where literature kinetics (kon, koff) exist, they override ML estimates.

### 2.3 Binding Kinetics and Occupancy

- **Target occupancy (TO) models** explicitly track drug-target complex formation with kon (second-order) / koff (first-order) and target turnover rate ρ (Daryaee & Tonge 2019 review; classical PK/PD literature). This is the preferred Stage-2 mechanism because it captures non-equilibrium behavior of the human body.
- **Target-Mediated Drug Disposition (TMDD)** and its quasi-steady-state (QSS) approximation (Mager & Jusko 2001; Gibiansky et al.) account for saturable loss of drug through high-affinity binding — essential when target abundance is comparable to drug amount.

---

## 3. Stage 3: Target Binding -> Signaling Pathway (QSP)

### 3.1 Quantitative Systems Pharmacology

QSP is defined (Sorger et al. 2011; Vicini & van der Graaf 2013) as the integration of quantitative approaches from systems biology, engineering and pharmaceutical sciences to link drug-interaction kinetics to cellular/organismal response in the context of physiology. QSP models are multiscale: they connect protein- or drug-interaction kinetics to cellular response, then to clinical biomarkers.

**Representative mechanistic example:** a first-of-its-kind MET-pathway QSP model (Front. Pharmacol. 2025, combining mechanistic QSP) built ~130 molecular species and ~69 ODEs covering receptor-ligand binding, receptor phosphorylation, downstream AKT/MAPK signaling, and regulatory crosstalk among MET/EGFR/ALK/ROS1; it was calibrated on cell-line, Animal, and clinical data for 16 drugs. This is exactly the "signaling network" grain size DrugOS needs, and demonstrates the standard formalism: **mass-action kinetics + Hill-type reaction laws**.

### 3.2 Pathway Databases as Model Sources

- **KEGG, Reactome, PANTHER, WikiPathways** provide hand-curated pathway networks (nodes, edges, gene/protein annotations). The QSP review (PMC4917453) treats them as the scaffold for building quantitative signaling models: upstream/downstream interaction identification from network topology, then conversion to ODE systems.

### 3.3 Enzyme/Kinase and Reaction Kinetics

Kinase enzyme kinetics (MM), cooperative binding (Hill), transcription factor activation (Hill-type), and protein turnover (synthesis/degradation) are the standard building blocks. Receptor-occupancy-to-effect transduction follows stimulus-response cascades (Emax with EC50 < Kd under amplification; Hill coefficient emergent from cascade depth).

---

## 4. Stage 4: Signaling Pathway -> Organ Function (QST and Organ Models)

### 4.1 Quantitative Systems Toxicology (QST)

QST extends QSP to safety endpoints. DILIsym — the most widely used QST platform — mechanistically represents drug-induced liver injury with sub-models for **(a) bile-acid transporter inhibition, (b) mitochondrial dysfunction/ETC inhibition, (c) oxidative stress, (d) immune-mediated pathways**, plus adaptive mechanisms (mitochondrial biogenesis). It takes in-vitro assay parameterization and PBPK-derived liver exposure as inputs and simulates virtual populations for liver-safety endpoints (e.g., ALT > 3x ULN, bilirubin > 2x ULN, Hy's Law cases). It has supported regulatory submissions: e.g., differentiating hepatotoxic vs. non-hepatotoxic CGRP antagonists (predicted before clinical confirmation) and guiding fezolinetant Phase-3 dose selection (DILIsym v7A, CPT: Pharmacometrics & Systems Pharmacology 2026).

**On-target vs off-target toxicity** (Beattie & Sher 2025 review): on-target toxicity = exaggerated pharmacology of the intended target (modeled as QSP with safety biomarkers out-of-range); off-target toxicity = adverse events via additional targets (modeled as universally-applicable mechanism panels, often organ-specific).

### 4.2 Organ-Level Physiological Models (Physiome Project)

The Physiome Project / Virtual Physiological Human provides validated organ and organ-system models (cardiovascular circulation with ventricular interaction and valve dynamics; neural nephron system for kidney; cardiac electrophysiology). These provide the physiological substrate to translate pathway perturbations into organ-output changes (cardiac output, blood pressure, glomerular filtration, ECG/QT).

### 4.3 Adverse Outcome Pathway (AOP) Framework

AOPs formally encode Molecular Initiating Event (MIE) -> Key Events (KEs) -> Adverse Outcome. They are the canonical organizational structure to bridge Stage 2/3 molecular events with Stage 5 toxicological endpoints — and are explicitly used to structure QST models (Beattie & Sher 2025).

---

## 5. Stage 5: Organ Function -> Clinical Phenotype

### 5.1 Translational PK/PD to Biomarkers and Endpoints

- The translational QSP literature (e.g., Translational QSP review, AAPS J 2019) shows QSP models quantitatively link drug-induced signaling to clinical biomarkers and proposed target-disease endpoints.
- Emax / sigmoid Emax theory (Hill equation; Goutelle et al. 2008) provides the canonical dose-concentration-effect mapping; it is the minimum viable mechanism for any biomarker or physiological indicator with saturable response.
- **Exposure-based DILI risk:** a retrospective analysis of 241 drugs (Wang et al. 2023 / similar integrative mechanistic studies) showed that the ratio of predicted Cmax to in-vitro toxicity IC50 correlates with clinical DILI risk with ROC AUC up to 0.91-0.96 prospectively — providing an explicitly validated, implementable toxicity-score formula.

### 5.2 Physiologically Plausible Biomarker Translation

For each organ model, output physiological indicators are mapped to clinical-grade measurements (ALT/AST/bilirubin for liver; QTc/HR/BP for cardiovascular; GFR/serum creatinine for kidney), with reference ranges and severity tiers (CTCAE-like grading) for phenotype output.

---

## 6. Validation Paradigm Presented in Literature

The literature supports the following validation ladder that DrugOS will follow:
1. **Benchmark compounds** with extensive PK/PD data (from OSP PBPK Model Library + literature profiles).
2. **Prospective prediction tests**: e.g., predicting DILI or pharmacodynamic readouts before comparing to observed clinical data (as in DILIsym publications and CORDIAL's prospective screening design).
3. **Quantitative goodness-of-fit**: predicted-vs-observed profiles, geometric-mean fold error in population simulations, AUC under ROC for toxicity classifiers.

---

## 7. Key Literature References (Selected)

1. Nestorov, I. Whole-body PBPK models. Expert Opin. Drug Metab. Toxicol. 2007.
2. Rowland, M.; Peck, C.; Tucker, G. PBPK in drug development and regulatory science. Annu. Rev. Pharmacol. Toxicol. 2011.
3. Willmann, S. et al. Physiology-based whole-body population model. J. Pharmacokinet. Pharmacodyn. 2007.
4. Swanson, K. et al. ADMET-AI: a machine learning ADMET platform. Bioinformatics 2024. (doi:10.1093/bioinformatics/btae416)
5. Obrezanova, O. et al. Prediction of in-vivo PK parameters from chemical structure. Mol. Pharmaceutics 2022.
6. PK-Sim / Open Systems Pharmacology documentation and OSP PBPK Model Library (open-systems-pharmacology.org).
7. Daryaee, F.; Tonge, P. J. PK/PD models that incorporate drug-target binding kinetics. Curr. Opin. Chem. Biol. 2019.
8. Mager, D. E.; Jusko, W. J. Target-mediated drug disposition. J. Pharmacokinet. Pharmacodyn. 2001.
9. Zhou, H. et al. Recent advances in ML predictions of protein-ligand binding affinities. Curr. Opin. Struct. Biol. 2026.
10. IPBind: geometric deep learning binding affinity prediction, arXiv:2504.16261 (2025).
11. CORDIAL: interaction-only deep learning affinity ranking, PNAS 2025 (10.1073/pnas.2508998122).
12. Zia, A. et al. Drug-target interaction/affinity prediction: deep learning models and advances review. Comput. Biol. Med. 2025.
13. MET-pathway multiscale QSP model. Front. Pharmacol. 2025 (10.3389/fphar.2025.1685468).
14. Sorger, P. K. et al. Quantitative and Systems Pharmacology in the Post-genomic Era (NIH white paper) 2011.
15. Physiome Project / Virtual Physiological Human (physiomeproject.org; Interface Focus 2011).
16. DILIsym: QST impacting drug development (Watkins / DILI-sim Initiative; CPT:PSP 2019; Toxicol. Sci. 2020).
17. DILIsym fezolinetant Phase-3 dose-selection study. CPT 2026 (10.1002/cpt.70194).
18. Beattie, K. A.; Sher, A. Application of mechanistic mathematical modeling to toxicology: QST. Handb. Exp. Pharmacol. 2025.
19. Goutelle, S. et al. The Hill equation: a review of its capabilities in pharmacological modelling. Fundam. Clin. Pharmacol. 2008.
20. Mager/Bloomingdale et al. Quantitative systems toxicology. Curr. Opin. Toxicol. 2017.