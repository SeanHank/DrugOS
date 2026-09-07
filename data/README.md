# DrugOS `data/` — vendored datasets

Everything under `data/` is pinned by sha256 in `data/manifest.json`, and
`python scripts/release.py gates` verifies every checksum before a release.
**Do not edit or add files here without updating the manifest on purpose.**

## Vendored now

| Table | File | Content | Consumers | Checksum source |
|---|---|---|---|---|
| `pk` | `benchmarks/published_pk.json` | Original experimental/PK pass-bands (CL, t½, Fa, fe, Vss) for the 5-benchmark compound corpus (midazolam, acetaminophen, warfarin, ciprofloxacin, dofetilide) | `validation/benchmarks/base.py` (`load_published`/`load_cites`/`load_about`) → every benchmark compound module | `data/manifest.json` |
| organ | `models/bsep_shh_ic50_reference.json` | Representative BSEP (ABCB11) SHH efflux-inhibition IC50 anchors (µM) transcribed from the open-access de Bruijn & Rietjens 2024 paper (CC BY 4.0), used to pin the liver cholestasis anchor (R-7, `simulate_gcdca_pbk`) | `validation/cases/case_liver_cholestasis_pbk.py` (R-7) and the organ cholestasis lane | `data/manifest.json` |

The bands are the published measured ranges used as validation pass criteria;
they are loaded from this file (never inlined), so the reference and the code
cannot drift apart. Full citation strings are in
`validation/benchmarks/citations.py` (Goodman & Gilman 13th ed.; host USPIs;
Bergan 1986; Greenblatt 1984; Smith 1992; partition model Rodgers & Rowland
2006, Ye et al. 2016).

## Why the rest is *not* vendored

doc/04 catalogs the full research corpus. The bulk sources are either
distributed as installable Python packages (ADMET-AI weights via pip,
RDKit cheminformatics, the physiology tables already compiled into
`src/drugos/pk/physiology.py` / `partitions.py`) or are license-restricted
bulk downloads that must be obtained under academic license with
rate-limit-aware fetch (full DrugBank distribution, LiverTox). Vendoring them
would either duplicate the package supply chain (checksum drift risk) or ship
content we are not licensed to redistribute. Each such source is documented in
doc/04 with its intended acquisition route.

## Adding a dataset

1. Drop the file under `data/<category>/`.
2. Record it in `data/manifest.json`: `path`, `sha256` (compute with
   `shasum -a 256`), `source`, `license`, `status`, `retrieved_utc`.
3. Point the consuming code at the file and fail loudly (never a clever
   default) when it is missing or malformed.
4. The gate will verify the checksum; a release with a tampered/unpinned file
   is blocked.