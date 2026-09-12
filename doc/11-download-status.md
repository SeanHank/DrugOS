# Vendored Data & Model Registry

Every shipped data artifact is vendored under `data/` and pinned by sha256 in
`data/manifest.json`; the checksum gate in `scripts/release.py` fails the
release on any mismatch or shadow edit. Dependencies distributed as packages
(ADMET-AI weights, RDKit, myokit, python-libsbml) are installed pin-exact via
pip/conda and are never vendored. Published equations that anchor a stage are
implemented directly in `src/drugos` from their open-access sources (listed
below). Companion docs: `doc/10` (dataset comparability), `doc/12` (model
integration matrix), `data/README.md` (vendoring policy).

## Vendored artifacts (sha256-pinned)

| # | Data / model | Source | License | Pinned artifact | Wired use / regeneration |
|---|---|---|---|---|---|
| 1 | Published PK reference bands (CL, t½, Fa, fe, Vss) for the five benchmarks | Goodman & Gilman 13th ed.; host USPIs; Bergan 1986; Greenblatt 1984; Smith 1992; Rodgers & Rowland 2006 / Ye et al. 2016 partitions | Published numeric ranges, cited academic use | `data/benchmarks/published_pk.json` (`11855a66…2115`) | GMFE endpoint bands (L25, `case_pk_gmfe_endpoints`) |
| 2 | Measured hERG IC50 rows (five benchmarks) | ChEMBL REST API, target CHEMBL240 | CC0 / public-domain summary | `data/benchmarks/herg_measured_nm.json` (`11ade841…4726`) | measured-hERG override (R-2); rebuilt by `scripts/data/fetch_chembl_herg.py` |
| 3 | Off-target bioactivity snapshot — 16 targets, 3035 rows | ChEMBL REST API (measured `=` IC50/Ki rows) | CC0 / public-domain summary | `data/chembl/offtarget_snapshot.json` (`b577438d…f1ac`) | kNN DTI resolver (L15, `drugos.target.dti`); rebuilt by `scripts/data/fetch_chembl_offtarget.py` |
| 4 | hERG / TdP assay corpus — 306,893 rows | hERG Central (Du et al. 2022), Dataverse doi:10.7910/DVN/7BVDG8, datafile 5724875; TDC `Herg` provenance | CC0 | `data/corpora/herg_central.tsv.gz` (`d44a4e14…86e`) | corpus calibration (R-2/R-8) |
| 5 | ORd 2011 human ventricular AP model | O'Hara-Rudy 2011 (PLoS Comput Biol e1002061); Myokit encoding from `github.com/myokit/models` | BSD-3 (Myokit); CC-BY publication | `data/models/ohara-2011.mmt` (`0ba58e14…c27`) | R-3 cardiac cross-check via `drugos.organ.cardiac_ap` (needs `myokit` + SUNDIALS headers — doc/12 §5) |
| 6 | ORd-CiPA-v1 2017 retune | Dutta et al. 2017 (Front Physiol 8:616); Myokit encoding from `github.com/myokit/models` | BSD-3 (Myokit); CC-BY publication | `data/models/ohara-cipa-v1-2017.mmt` (`bac881a5…410`) | multi-ionic ORd lane (L19, `case_ord_multi_ionic_qt`) |
| 7 | Huang & Ferrell 1996 MAPK cascade SBML | BioModels BIOMD0000000009 | CC0 | `data/models/huang1996-mapk-cascade.xml` (`1f95d793…ab98`) | Stage-3 pathway lane (R-5), readout PP_K, via `drugos.pathway.sbml_pathway` (python-libsbml, pinned in `pyproject.toml`) |
| 8 | Rohwer 2000 E. coli PTS SBML | BioModels BIOMD0000000038 | CC0 | `data/models/rohwer2000-pts.xml` (`a7ef9080…701a`) | non-MAPK scaffold for the SBML equivalence gate (L18, `case_sbml_scaffold_equivalence`) |
| 9 | BSEP (ABCB11) SHH efflux-inhibition IC50 anchors | de Bruijn & Rietjens 2024 (*Arch. Toxicol.* 98:3077, CC BY 4.0); primary BSEP measurements (Morgan 2010; Marchant 2019) | CC BY 4.0 (paper); transcribed numeric table | `data/models/bsep_shh_ic50_reference.json` (`26639e08…65b`) | R-7 cholestasis lane (`drugos.organ.liver.simulate_gcdca_pbk`, Ki = IC50/2). Companion model repo is CC-BY-NC-ND; its code is not ported — only the paper's open-access equations and these numeric anchors |

## Code-implemented published anchors (no file)

- **Willmann et al. (2007) allometric physiology** — the population power
  laws implemented in `drugos.pk.physiology.build_human` against the Ye
  (2016) reference tables; `case_willmann_allometric_physiology` (L14).
- **CKD-EPI 2021 race-free GFR equation** — Levey et al., *N Engl J Med*
  2021;385:1737, implemented as `drugos.organ.kidney.ckdepi_2021_egfr`
  (Mosteller BSA scaling; R-6), applied whenever serum creatinine is on the
  profile.
- **Safety-panel class priors, organ IC50 priors, CTCAE ladders, toxicity
  fusion priors** — literature-derived class medians in
  `drugos/target/targets.py`, `drugos/organ/*` and `drugos/clinical/*`, with
  URL-annotated references in the code.

## Package-distributed dependencies (installed, never vendored)

- **ADMET-AI** GNN weights + DrugBank reference percentiles — `pip install
  admet-ai` (pinned in `pyproject.toml`); Zenodo-archived model, BSD-3 code,
  CC-BY-4.0 weights.
- **RDKit** (BSD-3) — physchem, SMILES parsing.
- **myokit** (BSD-3) — ORd model execution; requires system SUNDIALS headers.
- **python-libsbml** (BSD) — SBML parsing, pinned in `pyproject.toml`.
- **R ≥ 4.5 + rpy2 ≥ 3.6** — executes the literature PK estimator
  `rbridge/literature_pk.R` on every run (doc/12 row 1b).

## Reproducibility

- Every vendored file is recorded in `data/manifest.json` with `path`,
  `sha256`, `source`, `license`, `status` and `retrieved_utc`; the `gates`
  checksum step recomputes each sha256 and fails on any mismatch.
- Regeneration is scripted: `scripts/data/fetch_chembl_herg.py` rebuilds the
  measured hERG rows and `scripts/data/fetch_chembl_offtarget.py` rebuilds
  the off-target snapshot from the ChEMBL REST API; each script is named in
  the manifest `source` note for that artifact.

## Summary

- 9 vendored sha256-pinned files across `benchmarks/`, `corpora/`, `chembl/`
  and `models/`.
- 3 code-implemented published anchors (Willmann physiology, CKD-EPI 2021,
  priors/ladders) and 5 package-distributed runtime dependencies.