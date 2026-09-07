# Dataset Download-Status Checklist

Status of every dataset the pipeline needs (`doc/10` details the science;
`doc/04` the catalog; `data/README.md` the vendoring policy).

**Legend**
- `VENDORED` — on disk under `data/`, pinned by sha256 in `data/manifest.json`,
  verified by `scripts/release.py gates`.
- `VIA-PKG` — installed as a dependency (pip/conda); never vendored.
- `NOT-DOWNLOADED` — not on disk; required before the corresponding realism
  upgrade is enabled.
- `LICENSE-GATED` — availability requires an academic-license application or
  rate-limited API fetch; not downloaded.
- `N/A` — no file; quantity currently a hardcoded prior.

| # | Dataset / quantity | Source | Status | Where / how to obtain |
|---|---|---|---|---|
| 1 | Published PK bands (CL, t½, Fa, fe, Vss) for the 5-benchmark corpus | Goodman & Gilman 13th ed.; USPIs; Bergan 1986; Greenblatt 1984; Smith 1992 | **VENDORED** | `data/benchmarks/published_pk.json` (`sha256 11855a66…2115`) |
| 2 | ADMET-AI GNN weights + DrugBank reference percentiles | Zenodo-archived model, distributed via PyPI | **VIA-PKG** | `pip install admet-ai` (pinned in `pyproject.toml`) |
| 3 | RDKit cheminformatics (physchem, SMILES parsing) | RDKit (BSD-3) | **VIA-PKG** | `pip install rdkit` |
| 4 | OSP / PK-Sim physiology database (volumes, flows, tissue composition) | Open Systems Pharmacology (Apache-2.0) | **NOT-DOWNLOADED** | `github.com/Open-Systems-Pharmacology`; import via `ospsuite` |
| 5 | Willmann 2007 population equations | Literature (published equations) | N/A (implemented as power-laws in `physiology.py`) | — |
| 6 | ChEMBL bioactivity (IC50/EC50/Ki/Kd per target) | EMBL-EBI (CC BY-SA 4.0) | **PARTIAL — measured subset VENDORED** | bulk `chembl_webresource_client` for the 5-benchmark hERG IC50; `scripts/data/fetch_chembl_herg.py` → `data/benchmarks/herg_measured_nm.json` (`sha256 11ade841…4726`) |
| 7 | DrugBank full distribution (approved targets, pharmacology) | DrugBank (academic licence) | **NOT-DOWNLOADED / LICENSE-GATED** | drugbank.com academic application |
| 8 | hERG / TdP assay corpus (patch-clamp IC50) | hERG Central (Du et al. 2022), Dataverse doi:10.7910/DVN/7BVDG8; TDC `Herg` provenance | **VENDORED** | `data/corpora/herg_central.tsv.gz` (`sha256 d44a4e14…86e`); original .tab id 5724875 |
| 9 | DILI labels (LiverTox; TDC DILI) | NIH/NCATS; TDC (open) | **NOT-DOWNLOADED / LICENSE-GATED (LiverTox bulk)** | LiverTox query/NCATS terms; TDC `tdcommons` loader blocked in this env (py3.12 `pkg_resources` removal); fetch raw dataverse/CSV directly instead |
| 10 | In-vitro toxicity assay parameters (BSEP inh., ETC, oxidative stress) | DILIsym literature consortia | **PARTIAL — BSEP SHH IC50 anchors VENDORED** (row 19b); ETC/oxidative-stress tables still NOT-DOWNLOADED | cited papers (numerical tables) |
| 11 | KEGG pathway graphs (hsa:) | KEGG (academic licence) | **NOT-DOWNLOADED / LICENSE-GATED** | KEGG FTP mirror |
| 12 | Reactome pathway graphs | Reactome (CC BY 4.0) | **NOT-DOWNLOADED** | Reactome graph DB / SBML download |
| 13 | PANTHER pathway terms | PANTHER (free) | **NOT-DOWNLOADED** | pantherdb.org |
| 14 | Physiome organ models (heart/nephron/liver SBML) | Physiome Model Repository (CC BY-SA) | **NOT-DOWNLOADED** | physiomeproject.org |
| 14b | Human ventricular action-potential model (ORd 2011) + CiPA-v1 2017 retune | O'Hara-Rudy 2011 (PLoS CB); ORd-CiPA-v1 (Dutta et al. 2017) via Myokit `.mmt` | **VENDORED** | `data/models/ohara-2011.mmt` (`sha256 0ba58e14…c27`) + `data/models/ohara-cipa-v1-2017.mmt` (`sha256 bac881a5…410`); from `github.com/myokit/models` (BSD-3); drives the R-3 cardiac cross-check lane via `drugos.organ.cardiac_ap` (needs `myokit` pip + SUNDIALS headers, doc/12/06) |
| 14c | CKD-EPI 2021 race-free GFR equation (kidney baseline, R-6) | Levey/Inker et al., *N Engl J Med* 2021;385:1737 (open publication) | N/A (published equation, implemented in `organ.kidney`; no file) | — |
| 19b | BSEP (ABCB11) SHH efflux-inhibition IC50 anchors (cholestasis lane, R-7) | de Bruijn & Rietjens 2024 paper, *Arch. Toxicol.* 98:3077 (CC BY 4.0); underlying SHH IC50 in primary BSEP literature (Morgan 2010; Marchant 2019) | **VENDORED** | `data/models/bsep_shh_ic50_reference.json` (`sha256 26639e08…65b`).  Note: the author's model repository is CC-BY-NC-ND — its R code is *not* ported; only the paper's open-access equations (in `organ.liver`) and these published numeric anchors are used |
| 15 | Population anthropometry + covariate correlations | NHANES (US CDC, public) | **NOT-DOWNLOADED** | CDC NHANES bulk XPT/SAS |
| 16 | QTc / ΔQTc clinical reference distribution | Published cardio-safety aggregates | **NOT-DOWNLOADED** | ICH E14 / literature summary tables (no raw XML redistribution) |
| 17 | Endpoint co-occurrence (SIDER / Offsides / FAERS) | SIDER (CC BY-NC-SA); Offsides; FAERS (public) | **NOT-DOWNLOADED** | SIDER/Offsides dumps; FAERS quarterly files (needs dedup script) |
| 18 | Pathway kinetic rate constants (Kd/Vmax/kcat) | BRENDA, SABIO-RK | **NOT-DOWNLOADED** | via curated extraction, rate-limited API |
| 18b | Ultrasensitive MAPK cascade SBML (pathway Stage-3 production anchor, R-5) | Huang & Ferrell 1996 via BioModels BIOMD0000000009 (CC0) | **VENDORED** | `data/models/huang1996-mapk-cascade.xml` (`sha256 1f95d793…ab98`); parsed by `drugos.pathway.sbml_pathway` via `python-libsbml` (pip, pinned `pyproject.toml`); drives the pipeline-default signal lane (readout PP_K). Stub `stubs/libsbml/` + doc/06 for G2 |
| 19 | Safety-panel class priors (16 sites), organ IC50 priors, CTCAE ladders, toxicity fusion priors | Literature-derived class medians | N/A (hardcoded in `targets.py`/`organ/*`/`clinical/*`; URL-annotated references in code) | — |

## Summary

- **VENDORED (checksum-pinned): 7** files (`published_pk.json`, measured hERG
  submation `herg_measured_nm.json`, `corpora/herg_central.tsv.gz`, ORd 2011
  + ORd-CiPA-v1 2017 ventricular models, the Huang/Levchenko MAPK SBML
  `models/huang1996-mapk-cascade.xml`, and the BSEP IC50 reference anchors
  `models/bsep_shh_ic50_reference.json`).
- **VIA-PKG (always available with the install): 4** (`admet-ai`, `rdkit`,
  `myokit`, `python-libsbml`; myokit additionally needs system SUNDIALS
  headers for its C codegen).
- **On-disk effective datasets/models: 8 units** serving 5 benchmarks + custom
  SMILES (published PK bands, per-compound measured hERG, 306 k-row hERG
  corpus, fully mechanistically runnable human ventricular AP models, the
  MAPK cascade SBML for the pathway lane, and the cholestasis IC50 anchors).
- **NOT-DOWNLOADED: 12** bulk corpora (rows 4, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18).
- **LICENSE-GATED subset: 3** (DrugBank full distribution, LiverTox bulk, KEGG).
- **N/A (implemented priors / published equations, code-is-the-truth): 5** (rows 5, 14c, 18 fronts, 19).

Any file added under `data/` MUST be recorded in `data/manifest.json`
(`path`, `sha256`, `source`, `license`, `status`, `retrieved_utc`) or a
`gates` run fails on the checksum mismatch — see `data/README.md`.

The acquisition targets above are the ones that unlock the P1–P9 realism
increments ranked in `doc/10 §4`.