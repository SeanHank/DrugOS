# DrugOS — Design Documents

Multiscale, mechanism-based modeling of drug response in the human body, following the chain:

```
drug -> in-vivo concentration -> target binding -> signaling pathway
     -> organ function -> clinical phenotype
```

| | |
|---|---|
| Author | Sean Hank |
| License | AGPLv3 |
| Version | 2026.9.1 (scheme: `YYYY.M.V` — year, month, intra-month revision from 0) |
| Python env | `/opt/anaconda3/envs/drug_os/bin/python` |

All inputs (chemical structure, route of administration, dose, complete human parameters) are **arbitrarily configurable** — see `01-project-overview.md` section 3.1.

## Document Index

| Doc | Content |
|---|---|
| [01-project-overview.md](01-project-overview.md) | Vision, causal chain, inputs/outputs, scope, principles |
| [02-literature-review.md](02-literature-review.md) | Literature basis per pipeline stage (PBPK, DTI-ML, QSP, QST, organ models, AOP) |
| [03-system-architecture.md](03-system-architecture.md) | Six-stage architecture, component breakdown, data contract, execution modes |
| [04-data-sources.md](04-data-sources.md) | External data register, ingestion/versioning plan, quality gates |
| [05-methodology-pipeline.md](05-methodology-pipeline.md) | Mathematical formulation and implementation plan per stage |
| [06-technology-stack.md](06-technology-stack.md) | Libraries, environment, file layout, reproducibility |
| doc/07 | Realization retrospective and milestones (M1-M6) |
| [08-validation-and-risk.md](08-validation-and-risk.md) | Three-tier validation ladder, benchmark set, risk register |

## Recommended Reading Order

1. `01-project-overview.md` — what and why
2. `05-methodology-pipeline.md` — how (the core)
3. `03-system-architecture.md` — how the modules fit together
4. `02-literature-review.md` — why (evidence base)
5. `doc/07` — realization retrospective (phases and milestones)
6. `04-data-sources.md`, `06-technology-stack.md`, `08-validation-and-risk.md` — supporting detail

## Quick Start for Contributors

- `doc/07`: start with the foundations phase, then the Stage-1 PK phase.
- The data contract (doc 03, SS3) is the integration API between stage modules.
- Stage 1 must be validated before Stage 2 coupling (see doc 07 M1 gate).