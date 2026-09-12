<div align="center">

# DrugOS

**Multiscale, mechanism-based modeling of drug response in the human body.**

_All models are wrong, but some are useful. — George E. P. Box_

From a SMILES string to a graded, evidence-attributed toxicity verdict — with
the whole biology on the way rendered as equations you can read.

</div>

<p align="center">
  <a href="https://www.python.org/downloads/"><img alt="Python 3.11+" src="https://img.shields.io/badge/Python-3.11%2B-4b7bec?logo=python&logoColor=fff&labelColor=232946"></a>
  <a href="#validation--quality-gates"><img alt="Validation 38/38" src="https://img.shields.io/badge/validation-38%2F38%20green-2acc74?style=flat"></a>
  <a href="#validation--quality-gates"><img alt="Coverage 100%" src="https://img.shields.io/badge/Coverage-100%25-2acc74?style=flat"></a>
  <a href="https://github.com/"><img alt="AGPL-3.0" src="https://img.shields.io/badge/License-AGPL--3.0-5865f2?style=flat"></a>
  <a href="doc/09-quality-gate.md"><img alt="Lint" src="https://img.shields.io/badge/Lint-ruff%20+%20mypy%20--strict-9855e2?style=flat"></a>
</p>

---

## What is DrugOS?

DrugOS is a **multiscale human drug-response simulator**. Give it a molecule, a
dose, a route and a human profile and it simulates the compound end-to-end —
whole-body exposure down to organ-level clinical biomarkers — returning a
transparent, mechanism-attributed risk verdict.

It is a **release-grade research platform**: every layer ships with documented
equations (`doc/05`), a callable CLI and a web playground, and an enforceable
validation suite whose **cases must all stay green before anything is released**
(`doc/09`).

```
 SMILES + dose + route + human profile
   │   RDKit descriptors + ADMET parameterization
   ▼
 Exposure   whole-body PBPK → tissue concentration (83 compartments, lumped circulation)
   │
   ▼
 Target     receptor binding kinetics → occupancy & time-at-target
   │
   ▼
Pathway    signal transduction (QSP ODEs — Huang/Levchenko SBML MAPK cascade)
    │
    ▼
  Organ      liver (DILI, bile-acid cholestasis PBK) · cardiac (QTc/TdP) · kidney (AKI/GFR, CKD-EPI) · CNS brain exposure
   │
   ▼
 Clinical   graded biomarkers + composite toxicity verdict (DILI / QT / AKI / CNS)
   │
   ▼
 Decision   uncertainty bands · virtual cohort · sensitivity · prospective limits (D21–D24)
```

## Features

| Layer | What it simulates | Mechanism |
|---|---|---|
| **Exposure** | C(t) in 83 tissue/plasma compartments; permeability/Fa-gated and logS solubility-limited oral absorption; full-fidelity (default-on, auto-anchored) tubular secretion, Michaelis–Menten hepatic clearance, biliary excretion with enterohepatic recirculation and first-pass gut-wall extraction | PBPK + ADMET + RDKit molecular descriptors |
| **Target** | receptor binding & occupancy kinetics | affinity-driven occupancy model |
| **Pathway** | signal transduction after exposure | QSP ODEs (Huang/Levchenko SBML MAPK cascade, vendored) |
| **Organ — liver** | DILI trajectory, ALT/AST/bilirubin, Hy's Law | QST dose–response + mechanistic grade; cholestasis axis anchored to the de Bruijn & Rietjens 2024 bile-acid PBK (R-7) |
| **Organ — cardiac** | QTc/TdP risk, ΔQTc, MAP | ion-channel QST + Fridericia correction; ORd-2011 ionic cross-check (R-3) |
| **Organ — kidney** | AKI grade, serum creatinine, GFR | GFR/creatinine QST; GFR baseline = CKD-EPI 2021 race-free from serum creatinine when present (R-6) |
| **Organ — CNS** | brain free exposure trajectory | passive blood–brain barrier model |
| **Clinical** | graded biomarkers + composite verdict | three-line fusion (mechanistic / exposure / structural) |
| **Decision** | 90 % uncertainty bands · cohort incidence · Sobol/OAT sensitivity · prospective anchoring | D21–D24 ensembles over the full pipeline |
| **Reliability** | predictive-regime disclosure on every run (`trust.reliability`: novel < partial < off-label/extrapolated dose < validated in-range < measured), regime-driven ensemble breadth, measured true-parameter overrides, observed-vs-predicted agreement rows (`trust.empirical_agreement`) | `reliability.py` (DISCLAIMER §2 made executable — doc/12 §7.3) |

## Try it in under a minute

```bash
# development install
python -m pip install -r requirements.txt -r requirements-dev.txt
python -m pip install -e ".[dev]"

# full quality gate: lint + types + 100 % branch coverage + validation suite
python scripts/release.py gates

# or just fire off a run right now
python -m drugos run --benchmark dofetilide --dose 0.5 --route oral
```

Or install a published artifact (built automatically by CI — see
[Release Engineering](#release-engineering)):

```bash
python -m pip install dist/drugos-2026.9.1-py3-none-any.whl
```

## CLI

```bash
drugos --version                    # drugos 2026.9.1
drugos benchmarks                   # acetaminophen warfarin midazolam ciprofloxacin dofetilide

# full pipeline report (markdown on stdout; json/html via --out DIR)
drugos run --smiles "CC(=O)Nc1ccc(O)cc1" --dose 1000 --route oral
drugos run --benchmark dofetilide --dose 0.5 --format json --out reports/

# decision study: uncertainty ensemble + virtual cohort + global sensitivity
drugos study --benchmark dofetilide --n-unc 17 --n-pop 16 --n-sobol 8 --out studies/

# reliability: measured true-parameter overrides + observed-data comparison
drugos run --benchmark acetaminophen --measured '{"fup": 0.75, "cl_hep_l_h": 22.0, "cl_renal_l_h": 2.0}' --empirical '{"plasma_cmax_mg_l": 12.0}' --format json --out reports/

# interactive web playground
drugos serve --port 8080            # open http://127.0.0.1:8080
```

Example output (`drugos run --benchmark dofetilide --dose 0.5 --route oral`):

```markdown
# DrugOS pipeline report — dofetilide

- Dose 0.5 mg (oral) · Cmax 0.002052 mg/L · AUC0-t 0.02953 mg·h/L · tmax 3.06 h

## Target engagement
| Target | Peak occupancy | Time-at-target (h) |
|---|---|---|
| hERG (Kv11.1) | 0.006375 | 0.06864 |
| Estrogen receptor | 0.001381 | 0.007194 |
| … 15 more off-target rows from the safety panel |

## Pathway signaling
- **erk_active** peak fold-change 37.29 vs drug-free baseline.

## Organ trajectories
- **Liver (DILI):** grade 0, ALT 1 xULN, bilirubin 1 xULN, Hy's Law not met
- **Cardiac:** ΔQTc 20.86 ms, peak QTc 435.9 ms (none band), MAP 93 mmHg
- **Kidney:** AKI grade 0, peak Scr ratio 1, min GFR 125 mL/min

## Composite toxicity
| Endpoint          | Risk   | Driver       | Reason                       |
|-------------------|--------|--------------|------------------------------|
| Drug-induced liver injury | 0.031 | mechanistic | organ QST grade 0 (normal) |
| QT prolongation / TdP | 0.813 | mechanistic | organ QST grade 1 (mild) |
| Acute kidney injury  | 0.220 | mechanistic | organ QST grade 0 (normal) |
| CNS liability        | 0.200 | mechanistic | class prior only (unanchored) |

**Verdict:** High composite risk (81%, driver qt)
```

## Web Playground

`drugos serve` launches a Flask single-page app — deep-purple dark theme —
that runs the live pipeline in your browser:

- molecule input via SMILES **or** one of the five benchmarks, dose, route, profile;
- **Target engagement** (occupancy table) and **Pathway signaling** panels;
- **CNS · brain free exposure** trajectory in the organ views;
- decision toggles on any run:
  - **Uncertainty (D21)** — 90 % credible bands around every readout;
  - **Population (D22)** — a virtual cohort with incidence of grade ≥ 1 outcomes;
  - **Sensitivity (D23)** — local drivers + first/total Sobol indices;
- everything backed by the same `POST /api/run` contract the CLI uses.

## Validation & Quality Gates

**Validation tier** — the pipeline's stages are anchored to
**production-validated, open-source models** (`doc/12-production-models.md`):
wheel-strength ADMET/T priors from **ADMET-AI**, measured hERG + hERG Central
corpora for target binding (R-2), the **O'Hara–Rudy 2011** human ventricular
AP model (BSD-3 Myokit encoding) as the cardiac cross-check (R-3), the
**CKD-EPI 2021 race-free** creatinine equation as the kidney GFR baseline
(R-6), the **de Bruijn & Rietjens 2024 GCDCA bile-acid PBK** (CC BY 4.0) as
the liver cholestasis anchor (R-7), and a **required-R** literature-PK
estimator that re-derives every clearance/AUC result in the R runtime (R-1,
agreement gated ≤ 2 %).

`validation/` holds the tier ladder; `python validation/run_validation.py`
regenerates `validation/report.md` (currently **38/38 cases green**):

| Tier | Case | Checks |
|---|---|---|
| L1 analytic | robustness sanity (D21–D23 behavior) | 90 %-band sanity, non-negativity, Sobol bounds, OAT sign |
| L1 analytic | prospective fidelity (D24) | deterministic re-runs, held-out high-risk profile retention |
| L2 limit | clinical grading | exposure line vs prior-only fusion, mechanistic grade ladder |
| L2 limit | risk ordering | dofetilide QT ≫ warfarin QT, dose–DILI separability |
| L3 empirical | benchmark suite | mass balance, dose proportionality, occupancy, pathway, organ PK |
| L3 empirical | required-R (R-1) | R-literature CL/AUC/t½ agreement ≤ 2 % on all benchmarks, `r:agree` |
| L3 empirical | corpus calibration (R-2) | measured hERG IC50 geomean ≈ 26 nM; corpus %-inhibition long tail → per-compound hERG override |
| L3 empirical | ORd cardiac (R-3) | ORd 2011 baseline APD90 266 ms ∈ [200,350]; IKr block prolongs ΔAPD90 ≥ 30 ms (measured ≈ 115 ms); warfarin control = 0 |
| L2 limit | CKD-EPI 2021 (R-6) | male 60 y Scr 1.0 → 86.2 mL/min/1.73 m²; female lower at same Scr; BSA-scaled absolute GFR; Scr-profile drives `gfr_ml_min`; no-Scr keeps default |
| L2 limit | bile-acid cholestasis PBK (R-7) | ritonavir-class IC50 0.2 µM → ~10.3× intrahepatic GCDCA (cholestatic); itraconazole-class 10 mM → 1.0× (benign); rank order in Ki; Ki=IC50/2 |
| L2 limit | calibrated hERG head (R-8) | continuous monotone P→KD (floor 1 mM / ceil 2 nM / threshold ~1.4 µM within 10x of corpus median); ADMET hERG head ranks corpus `hERG_inhib` actives above inactives |
| L2 limit | clearance realism (secretion / gut-wall / MM / EHC) | single-pool secretion urine fraction analytic; F = Fa·(1−Eh)·(1−Eg); MM low-dose Vmax/Km slope + saturation; mass-conservative EHC |
| L2 limit | Cheng–Prusoff IC50→Ki | Ki = IC50/(1+[S]/Km); default [S]/Km=1 reproduces Ki=IC50/2; BSEP/CYP3A4/hERG anchors; invalid input rejected |
| L2 limit | per-CYP hepatic kinetics | abundance-scaled Vmax from the `CYP_ABUNDANCE_PMOL_MG` table; isoform additivity = lumped twin; doubling content → double slope; Hill sigmoid; parameter rejection |
| L2 limit | native TMDD drug disposition | reversible mass closure; internalized sink clears drug; dose-disproportional retention; quasi-steady DR/R = D/Kd; AUC falls vs free twin; default-on (auto-engaged at the primary-affinity panel site in full runs) |
| L2 limit | multi-layer transdermal permeation | Fick steady flux = composite P_eff·A·C_surf (±2%); interface partitions recovered; SC barrier + diffusivity responsiveness; mass closure; default-on (auto-engaged for the transdermal route in full runs) |
| L2 limit | sympathetic-suppression branch | IC50 free exposure halves HR and SV → CO = Q/4 (pa−pv = tone²·(MAP−CVP)); zero-exposure baseline preserved; saturating exposure → CVP floor; monotone exposure-response; degenerate inputs rejected; default-on (null-effect anchor auto-assigned in full runs) |
| L2 limit | immune-mediated DILI QST | steady-state I_ss = k_recruit·hazard/(k_recruit+k_decay); immune_weight=0 inert in the baseline lane; immune_weight>0 raises dead_frac; monotone exposure-response; degenerate inputs rejected; default-on (auto-anchored in full runs) |
| L2 limit | ACAT-lite multi-segment SI | mass conservation (zero clearance: state_total = dose); per-segment solubility cap produces dissolution-limited feces; more segments change dissolution dynamics; si_segments=0 rejected; baseline lane identical to single-SI; default-on (si_segments=3 in full runs) |
| L1 self-consistency | **E2E full-chain + trust** | real chain — structure → ADMET-AI → spec → pipeline → contract → rendered report — on two unrelated molecules: determinism, 8-section contract, Cmax > 0, Fa ∈ (0,1], per-molecule `trust` record (fidelity lane / anchors / G5 policy / class priors / engaged terms / estimates), artifacts written |

**Quality gate** — `python scripts/release.py gates` is the single release gate:
`scripts/release.py` (the merged successor of the former `quality_gate.sh`)
runs G1–G4 and a **no-silent-fallback audit**:

- **G1** ruff (E/F/W/I/UP/B, line length 100) — clean
- **G2** mypy `--strict` across `src/drugos` — no errors
- **G3** pytest with **100 % branch coverage** of `src/drugos`, no
  exclusions
- **G4** validation suite — 38/38 green, report regenerated
- **G5** fallback audit — every `except` handler in the package must surface an
  explicit error; silent swallowing is a hard failure (inventory pinned in
  `scripts/fallback_allowlist.json`)

The ORd lane needs `myokit>=1.39` (pip) plus system **SUNDIALS** headers for
myokit's C codegen (`conda install -c conda-forge sundials`) — see
`doc/06 §1`.

See [`doc/09-quality-gate.md`](doc/09-quality-gate.md) for the full contract.

## Release Engineering

CI (`.github/workflows/ci.yml`) automates releases on top of the quality gate:

1. **quality** — runs G1–G5 (via `python scripts/release.py gates`) on every
   push, PR and tag; uploads the gate log and `validation/report.md` as
   artifacts;
2. **build** (needs `quality`) — reads the release version from
   `pyproject.toml`, `python -m build` produces the **wheel**
   and **sdist**, validates them with `twine check`,
   smoke-tests the wheel in a clean virtualenv, and opens a **GitHub Release** with
   the artifacts attached.

Locally, `scripts/release.py` does the same bookkeeping in one shot: runs the
six gates (G1-G6, including the fallback audit and the marker audit), syncs
the project-wide version (`YYYY.M.V`) from an
explicit `--version` (never auto-bumped), rewrites the version / status numbers
across `README.md`, `doc/*.md` and the web banner, and writes
`build/release_status.json`.

Release → push to `main` or
`git tag 2026.9.1 && git push --tags`. Version scheme `YYYY.M.V`.

## Package layout

```
src/drugos/
  inputs/       # structure / route / dose / human-profile parsing
  pk/           # ADMET parameterization + whole-body PBPK
  target/       # target identification + binding occupancy
  pathway/      # signal-transduction QSP (vendored Huang/Levchenko SBML MAPK default)
  organ/        # liver · cardiac (QTc) · kidney · CNS QST
  clinical/     # biomarkers + composite toxicity fusion
  robustness/   # D21 uncertainty · D22 population · D23 sensitivity
  pipeline.py   # RunSpec → RunResult orchestration contract
  report/       # JSON / markdown / HTML report rendering
  web/          # Flask playground + /api/run (shipped in the wheel)
  cli.py        # run · benchmarks · study · serve · version subcommands
  version.py    # YYYY.M.V
```

## Design documents

Start with `doc/01-project-overview.md`, then `doc/05-methodology-pipeline.md`:

`doc/01` overview · `doc/02` data & datasets · `doc/03` system architecture ·
`doc/04` module design · `doc/05` methodology & pipeline · `doc/06` tech stack ·
`doc/07` roadmap · `doc/08` scope & limits · `doc/09` quality gate ·
`doc/10` dataset comparability · `doc/11` download status ·
`doc/12` production-validated model matrix

## Community

- [CONTRIBUTING.md](CONTRIBUTING.md) — how to add a validation case, open an
  issue, or propose a gate change.
- [CREDITS.md](CREDITS.md) — tooling, datasets and prior art the platform builds on.

## License

This project is licensed under the **GNU Affero General Public License v3.0** (AGPLv3).

See [DISCLAIMER.md](DISCLAIMER.md) for important legal notices and limitations.

Copyright © 2026 Sean Hank.