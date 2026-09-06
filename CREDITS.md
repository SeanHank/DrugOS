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
  (targeted for post-baseline releases; the shipped baseline compiles a MAPK
  cascade kernel).
- Tissue expression / quantitative proteomics tables — target abundance
  priors per tissue.
- **ADMET-AI** (Swanson et al., structural ADMET prediction) — machine-learning
  ADMET/physchem estimates feeding the PBPK parameterization.

## Modeling prior art the methodology builds on

- **PBPK / ACAT-style absorption** — Rodgers & Rowland and Poulin–Theil
  tissue-partition methods; physiologically based scaling of clearance and
  distribution.
- **Target occupancy kinetics** — classic target-mediated drug disposition and
  the Daryaee & Tonge occupancy/turnover formulations.
- **QSP pathway models** — the MET-QSP lineage (Erk/PI3K module conventions)
  for reaction-geometry and conservation-based steady-state initialization.
- **DILIsym** — bile-acid transport inhibition and mitochondrial dysfunction
  mechanisms for the liver (DILI) organ sub-model.
- **Organ-systems / quantitative systems toxicology** — the lumped-circulation
  cardiac-electrical QST and nephron/GFR kidney lineage used for the
  cardiovascular and renal panels.
- **Open Systems Pharmacology (OSP)** — an upstream cross-validation
  evaluation target for model-library compounds.

## The authors' note

This project is written, tested and maintained by **Sean Hank**. The baseline
numbering and quality-gate discipline ("four gates, 100 % coverage, no blanket
suppressions, validation-before-feature") follow the repository's own
governance in `doc/09-quality-gate.md`.

If a resource you contributed is represented inaccurately here, please open an
issue — corrections are welcome and credited.