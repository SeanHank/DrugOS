# Credits

DrugOS is an independent, from-scratch implementation of the design in
`doc/*.md`. The platform and its surrounding research build on the public work,
datasets, and tools of others; this page records the ones the project
explicitly draws from and stands on.

## Core engineering dependencies

- **RDKit** (RDKit Open-Source Cheminformatics Software) — molecular
  descriptors, pKa / LogP parameterization, stereochemistry.
- **SciPy / NumPy** — sparse ODE solvers (`solve_ivp`, LSODA), Sobol/Latin
  hypercube sampling, statistical utilities.
- **Pydantic** — typed configuration and input validation.
- **Flask** — the web playground.
- **Myokit** (BSD-3-Clause, Mirams et al.) — cardiac cell-modelling framework
  that compiles and runs the vendored human ventricular AP models; C backend
  uses **SUNDIALS** (BSD-3).
- **R language runtime** (GPL, via `rpy2`) — required to execute the literature
  PK estimator (`src/drugos/rbridge/literature_pk.R`).
- **Matplotlib / Plotly / Jinja2** — visualization and report rendering for the
  web and document outputs.
- The project is **Python 3.11+** and enforces `mypy --strict` and 100 %
  branch coverage on every commit.

## Data and knowledge resources referenced by the design

- **DrugBank** — target sets for approved/reference drugs (principal
  off-targets: hERG, CYP enzymes, BSEP, transporters, nuclear receptors).
- **ChEMBL / BindingDB** — measured binding kinetics and affinities for
  target-parameterization priors.
- **KEGG / Reactome / PANTHER** — pathway-graph assembly for the QSP layer
  (targeted for post-baseline releases; the shipped pipeline default is the
  Huang & Ferrell 1996 ultrasensitive MAPK cascade in SBML, vendored from
  BioModels BIOMD0000000009 (CC0) and validated as R-5).
- Tissue expression / quantitative proteomics tables — target abundance
  priors per tissue.
- **ADMET-AI** (Swanson et al., structural ADMET prediction) — machine-learning
  ADMET/physchem estimates feeding the PBPK parameterization (BSD-3 code,
  Zenodo weights).
- **hERG Central** (Du et al., *J Chem Inf Model* 2022) — dose-response patch-
  clamp corpus (Harvard Dataverse doi:10.7910/DVN/7BVDG8) calibration anchor.
- **ChEMBL** (EMBL-EBI) — measured hERG IC50 rows for the benchmark targets.

## Modeling prior art the methodology builds on

- **PBPK / ACAT-style absorption** — Rodgers & Rowland and Poulin–Theil
  tissue-partition methods; physiologically based scaling of clearance and
  distribution.
- **Target occupancy kinetics** — classic target-mediated drug disposition and
  the Daryaee & Tonge occupancy/turnover formulations.
- **QSP pathway models** — the MET-QSP lineage (Erk/PI3K module conventions)
  for reaction-geometry and conservation-based steady-state initialization.
- **DILIsym** — bile-acid transport inhibition and mitochondrial dysfunction
  mechanisms for the liver (DILI) organ sub-model; DILIsym itself is
  proprietary, so the liver cholestasis axis is anchored on the openly
  published de Bruijn & Rietjens bile-acid PBK (below).
- **de Bruijn & Rietjens 2024** (*Arch. Toxicol.* 98:3077-3095,
  doi:10.1007/s00204-024-03775-6, **CC BY 4.0**) — the GCDCA intrahepatic
  bile-acid PBK validated against the clinical cholestasis incidence of ~18
  marketed drugs; its equations are implemented in `drugos.organ.liver` as
  the cholestasis axis anchor (R-7).  The authors' companion GitHub
  repository (Veronique-de-Bruijn/PBK-model-cholestasis) is **CC-BY-NC-ND**
  and its R code is intentionally *not* vendored or ported; only the
  open-access paper's equations and the published numeric BSEP SHH IC50
  anchors (vendored in `data/models/bsep_shh_ic50_reference.json`) are used.
  Underlying SHH efflux-inhibition measurements: Morgan et al. 2010;
  Marchant et al. 2019.
- **CKD-EPI 2021 race-free** (Inker, Levey, Eneanya et al., *N Engl J Med*
  2021;385:1737-49) — the published standard creatinine-based eGFR equation,
  implemented in `drugos.organ.kidney` and applied as the kidney GFR baseline
  from a measured serum creatinine (R-6); BSA scaling via Mosteller 1987.
- **Organ-systems / quantitative systems toxicology** — the lumped-circulation
  cardiac-electrical QST and nephron/GFR kidney lineage used for the
  cardiovascular and renal panels.
- **O'Hara, Virág, Varró & Rudy 2011** (*PLoS Comput Biol* e1002061) — the
  undiseased human ventricular action-potential model, vendored in Myokit
  `.mmt` encoding as the cardiac validation anchor (R-3); the **ORd–CiPA-v1
  (2017)** retune (Dutta, Chang et al., *Front Physiol*) is vendored for the
  multi-ion-channel upgrade path. Myokit-encoded files originate from
  `github.com/myokit/models` (BSD-3) and reproduce the published equations;
  per the project's terms, the original journal papers should be cited, not
  the repository.
- **Open Systems Pharmacology (OSP)** — an upstream cross-validation
  evaluation target for model-library compounds.

## The authors' note

This project is written, tested and maintained by **Sean Hank**. The baseline
numbering and quality-gate discipline ("four gates, 100 % coverage, no blanket
suppressions, validation-before-feature") follow the repository's own
governance in `doc/09-quality-gate.md`.

If a resource you contributed is represented inaccurately here, please open an
issue — corrections are welcome and credited.