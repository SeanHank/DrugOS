# Methodology: End-to-End Prediction Pipeline

This is the core design document. It specifies the mathematical formulation, sub-model choices, and implementation plan for each pipeline stage. Every stage is grounded in the literature reviewed in `02-literature-review.md`.

---

## Stage 1 — Drug -> In-Vivo Concentration (PK)

### 1.1 Goal

Given a structure, route, dose, and human profile, produce `plasma` and per-organ concentration-time profiles.

### 1.2 Step 1A: Structure -> ADME/PhysChem Parameterization

Inputs: canonical SMILES. Outputs: the parameter set required by PBPK.

| Parameter | Derivation |
|---|---|
| Molecular weight, LogP, LogD(pH), pKa, TPSA, HBD/HBA | RDKit descriptors |
| Aqueous solubility (logS) | ADMET-AI / RDKit-ESOL |
| Fraction unbound in plasma (fu) | ADMET-AI (plasma protein binding) |
| Intrinsic clearance (Clint) | ADMET-AI metabolic stability; refined via CYP phenotyping (CYP3A4/2D6/2C9/2C19/1A2) from ML or literature |
| Permeability class / fraction absorbed (Fa) | ADMET-AI bioavailability + PAMPA/Caco-2 models; ACAT absorption for oral |
| Efflux/uptake transporter flags (P-gp, OATP1B1/1B3, BCRP, BSEP) | ChEMBL/ADMET-AI transporter inhibition + literature |
| T1/2, Vss priors | ADMET-AI half-life / volume; overridden by PBPK-computed values |

**Accuracy note.** ADMET-AI ranks first on the TDC leaderboard and is suitable for first-pass parameterization; when measured in-vitro ADME data or validated PBPK models exist (e.g., OSP Model Library), they take precedence.

**Measured true-parameter override (DISCLAIMER §2, doc/12 §7.3).** Any of the
ADME/potency scalars can be supplied as *measured* true parameters
(`--measured` JSON, `/api/run` `measured` payload, or the `Measurements`
dataclass: `fup`, `cl_hep_l_h`, `cl_renal_l_h`, `cl_sec_l_h`, `cl_bil_l_h`,
`fa`, `hepatic_vmax_mg_h`, `hepatic_km_mg_l`, `qt_ic50_nm`, `dili_ic50_nm`,
`dili_immune_ic50_nm`, `cns_ic50_nm`, `beta_block_ic50_nm`). Overrides are
applied at the very top of `run_pipeline`, credited in the trust record's
`estimates` as measured, and used directly by the same full-fidelity chain —
never routed through a synthesized stand-in.  A measured zero is an honest
determination (e.g. `cl_sec_l_h = 0` means "no active tubular secretion"):
the term stays engaged and is disclosed "(measured)".  Full measured PK
upgrades the run to the most reliable `measured_in_range_on_label` regime.
Observed clinical data can additionally be supplied (`--empirical` /
`EmpiricalObservations`: `plasma_cmax_mg_l`, `auc_last_mg_h_l`,
`peak_delta_qtc_ms`, `peak_alt_uln`); the run reports each as an
observed-vs-predicted `fold_error` / `within_2x` row in
`trust.empirical_agreement` — disagreement with empirical data is the expected
state of the mechanistic model, reported not suppressed.

### 1.3 Step 1B: Human Parameter Resolution

The user provides a sparse profile: age, sex, height, weight, (optional) organ status, disease state. The resolver computes:

- Organ/tissue volumes (vascular, interstitial), blood-flow fractions, cardiac output -> allometry/covariate equations (Willmann et al. 2007; ICBP reference tables).
- Plasma protein concentrations (albumin, AGP) — age/sex adjusted.
- Enzyme abundances (CYP450 per gram liver) — population references with covariance.
- Target expression levels per tissue (from proteomics/expression databases where available).

Output: a full PBPK parameter vector for this individual (used both for point simulation and as the center of a virtual-population distribution for uncertainty quantification).

### 1.4 Step 1C: Whole-Body PBPK Model

**Structure.** A system of compartments including at minimum: lung, heart, liver (with gut), kidney, brain, muscle, adipose, skin, bone, remainder (spleen, pancreas, fat-free/lean tissue lumping). Each compartment consists of sub-spaces where appropriate: vascular, interstitial, and cellular (intracellular) for small molecules, with passive permeability/partition exchange plus active processes.

**Governing equation (per compartment, unbound mass balance).**

```
dA_comp/dt = Q_comp * (C_in - C_out) + net passive flux + active flux - metabolism - excretion
```

where `C_in` is the inflow (arterial) concentration, `C_out = C_free_tissue / Kpu` (partition-corrected), `Q_comp` is organ blood flow obtained from the resolved physiology. Tissue partition coefficients use the Rodgers & Rowland (or Poulin-Theil) methods parameterized by LogP/pKa/LogD, tissue composition, and fu.

**Absorption (oral).** ACAT-style first-order transit: stomach -> small-intestine segments -> colon, with pH-dependent solubility/dissolution, permeability-driven absorption into the portal vein, first-pass hepatic extraction.

**Absorption (SC/IM/transdermal).** Subcutaneous, intramuscular and
transdermal doses enter a first-order `depot` compartment that feeds venous
blood with rate constant `k_depot_absorption` (default ~0.15/h) and routed
bioavailability `F`; the unabsorbed fraction (`1 - F`) is routed to feces, so
mass is conserved. Users may override the absorption rate per run
(e.g. `--sc-im-ka`), giving these routes a blunted, delayed Cmax relative to
the equivalent IV bolus (see §1.6 for the transdermal baseline note).

**Clearance.** Hepatic metabolism via MM/Hill kinetics on unbound liver concentration for each CYP (with abundance-scaled Vmax and literature/ML Km); renal excretion via glomerular filtration of unbound drug (fu * GFR * fraction unbound-driven), plus tubular secretion terms if transporter data exist.

**Tubular transport (renal).** Beyond glomerular filtration of the free
fraction (`fup*GFR`), an active tubular-secretion term (`cl_sec`, OAT/OCT
transporters) adds a first-order flux on the unbound kidney drug; it is off
by default (pure filtration, `RunSpec.cl_sec_l_h = 0`).

**Biliary excretion & enterohepatic recirculation.** A bile pool
(`PBPKModel` ``bile`` compartment) receives free parent from the liver at a
biliary-secretion rate `cl_bil` and empties by `k_bile_emptying` into the
small-intestine lumen; the emptied mass is then subject to the ordinary SI
reabsorption/transit kinetics, so a fraction recirculates via the portal vein
and the rest reaches the feces sink. Off by default (`cl_bil = 0`).

**First-pass gut-wall extraction.** Enterocyte (intestinal) first-pass
metabolism removes a fraction `gut_extraction_eg` of the SI-absorbed flux
before the portal blood, so oral `F = Fa·(1−Eh)·(1−Eg)` (§1.6). Off by
default (`Eg = 0`).

**Saturable hepatic clearance.** Hepatic metabolism can be switched from the
linear `cl_hep` term to a Michaelis-Menten `Vmax/Km` form on the unbound
liver concentration (`RunSpec.hepatic_vmax_mg_h`/`hepatic_km_mg_l`, set
together); below Km the slope equals `Vmax/Km`, so the linear default is the
low-dose limit. Off by default (linear clearance preserved).

**Per-CYP hepatic kinetics.** The lumped hepatic term can instead be replaced
by a *sum of per-isoform* Michaelis-Menten/Hill fluxes on the unbound liver
concentration (`RunSpec.cyp_terms`, a tuple of `CypTerm`). Each term carries
an isoform, `Km` (mg/L unbound liver), abundance-scaled `Vmax` (mg/h) and an
optional Hill coefficient `n` (default 1 = Michaelis-Menten). The abundance
scaling is the resource this stage had asked for: `cyp_vmax_mg_h` turns a
per-isoform liver content from the physiology CYP-abundance table
(`HumanPhysiology.hepatic_cyp_content_nmol`, Barter et al. 2013) into
`Vmax = kcat · content · MW`, so a run supplies kcat/Km and the table supplies
the capacity. Hill `n > 1` gives a sigmoidal (cooperative) flux — sharper
saturation below Km and faster approach to Vmax above it. Off by default,
and mutually exclusive with the lumped `hepatic_vmax_mg_h`/`hepatic_km_mg_l`
term.

**TMDD coupling (implemented).** If target engagement consumes significant drug (high-affinity, high-abundance target), Stage 2 occupancy fluxes (drug + receptor -> complex -> internalization/clearance) feed back into the tissue mass balances. A native coupling (`PBPKModel.target_binding`, DEFAULT-ENGAGED in full runs — auto-bound to the primary-affinity panel site, heart for hERG-name sites else liver) adds the binding site's turnover ODEs directly to the tissue ODEs, so internalization clears drug within the mass balance (§2.4); the sequential approximation via `feedback_loop` remains available for unbound-exposure-only runs.

**PK metrics.** Cmax, tmax, AUCinf, t1/2, Vss, Cl are computed from the simulated profiles and reported. **Systemic bioavailability `F` is also reported per route** (doc/12 row PK-f): IV = 1; depot routes (SC/IM/transdermal) = the depot bioavailability; oral = `Fa·(1 − Eh)·(1 − Eg)` where `Fa = 1 − feces_fraction` (the colon-transit feces sink), `Eh = CL_h/(Q_h + CL_h)` is the well-stirred first-pass hepatic extraction and `Eg` is the first-pass intestinal extraction (auto-anchored from the CYP3A4-inhibition probability in full runs, 0 in the baseline lane) — the model already routes absorbed oral drug into the liver compartment, so this is the first-pass-corrected value, `PkMetrics.f_abs`. `Fa` itself follows the permeability-gated absorption constant tuned to `RunSpec.fa` (§1.6), so reported `F` reflects both the drug's absorbability and its extraction.

### 1.5 Implementation

- Formulate as sparse coupled-ODE system; solve with `scipy.integrate.solve_ivp` (LSODA) with tight tolerances for stiffness.
- Optional integration: run the equivalent model in the OSP suite for model-library compounds as a cross-validation; keep equations transparent in Python.
- Population mode: sample virtual individuals (n>=100) from physiological covariance, producing 5-95th percentile bands.

### 1.6 Baseline Implementation Scope (2026.9.1)

The shipped engine implements a faithful *baseline* of the mechanics above.
The baseline's narrower span, itemized in doc/01 §5 (out of scope for the
2026.9.1 baseline; each reduction is a documented default, never a hidden
stand-in), is:

- **Clearance** defaults to a lumped single hepatic-metabolic `cl_hep` (from
  the resolved physiology and structure-derived partition) plus renal GFR of
  unbound drug.  Extra realism terms are implemented and engaged in every
  `fidelity="full"` run (the default) with disclosed auto-anchors:
  active tubular secretion (`RunSpec.cl_sec_l_h`,
  auto-anchored `cl_sec = CL_renal·(0.5 + 0.5·Pgp)` when absent),
  a saturable Michaelis-Menten hepatic term
  (`hepatic_vmax_mg_h`/`hepatic_km_mg_l`, auto-anchored `Vmax = CL_h·Km`,
  Km = 1 mg/L prior when absent),
  biliary drug excretion with enterohepatic recirculation
  (`cl_bil_l_h`/`bile_emptying_1h`, auto-anchored from logP) and
  first-pass gut-wall extraction (`gut_extraction_eg`,
  auto-anchored from CYP3A4 inhibition probability, §1.4).
  The validated linear lane remains reachable by an explicit
  `fidelity="baseline"` opt-out.  Per-CYP MM/Hill richness is implemented
  and engaged in full runs (`PBPKModel.cyp_terms` /
  `RunSpec.cyp_terms`); explicit per-CYP terms yield to the lumped MM
  auto-anchor when both are absent.  Abundance-scaled Vmax is built from the
  `CYP_ABUNDANCE_PMOL_MG` table (`cyp_vmax_mg_h`).
- **Oral absorption** is ACAT-lite: first-order stomach -> small-intestine ->
  colon transit with pH-dependent dissolution.  **ACAT-lite multi-segment SI (DEFAULT-ENGAGED in full runs, doc/05 §1.6,
  2026.9.1):** ``AbsorptionParams.si_segments``
  splits the small intestine into N sequential equal-volume sub-compartments,
  each with its own dissolution cap (``solubility_mg_ml * gi_volume_ml / N``)
  and first-order absorption/transit; per-segment transit is scaled to
  ``k_si_transit * N`` so total SI transit time is preserved, and bile
  secretion enters segment 0 (proximal SI).  When ``si_segments`` is None
  (baseline lane), the model is the single-compartment ACAT-lite baseline;
  full runs auto-engage ``si_segments = 3``.  **The
  small-intestine absorption rate is permeability/Fa-gated**: `_absorption_params` tunes
  `k_si_absorption` to the compound's fraction-absorbed target
  (`absorption_rate_from_fa`, `RunSpec.fa`) — from the benchmark published
  `f_abs` band (or `fa_override`), from the predicted HIA head for ADMET-AI
  SMILES runs (`_admet_fa`), and from the default endogenous rate otherwise —
  so a low-permeability molecule peaks lower/later and loses more to the
  colon-feces sink while a fast one absorbs in the proximal gut.  **Solubility-
  limited dissolution** is gated on the ADMET-AI `logS` for novel molecules:
  the dissolved small-intestine pool is capped at `S * V_gi`
  (`AbsorptionParams.solubility_mg_ml`), so an excess low-solubility dose
  spills undissolved mass forward to the colon/feces sink.  Benchmark runs
  (no ADMET call) are unaffected.
- **SC/IM absorption** is a first-order depot (mass-conserving, user-tunable
  `k_depot_absorption` / `F`); no microcirculation/lymphatic resolution in the
  subcutaneous tissue.
- **Transdermal absorption** defaults to the same first-order depot, but a
  finite-dose **multi-layer skin-permeation membrane** is auto-engaged in
  full runs on the transdermal route via
  `AbsorptionParams.skin_layers` / `RunSpec.skin_layers` (`SkinLayers`, §1.4),
  explicit values overriding it when set: vehicle surface -> stratum corneum -> viable
  epidermis -> dermis, each interface a reversible diffusion-limited link
  `J_k = D_k*A/L_k*(C_up - C_down/K_k)`, with first-order dermal capillary
  removal into venous blood weighted by `depot_bioavailability`.  The serial
  3-resistance membrane reaches the Fick steady state
  `J_ss = P_eff * A * C_surface` with
  `P_eff = 1/(L_sc/D_sc + L_ve/(D_ve*K_sc) + L_der/(D_der*K_ve*K_sc))`, and a
  no-flux limit returns the interface partitions; transdermal mass closes
  against the dose across surface + membranes + systemic + unabsorbed
  (`case_transdermal_multi_layer`).
- **ADMET -> PK auto-wiring** is first-class: `pipeline.spec_from_admet` builds
  a full `RunSpec` for a novel molecule from `fup_plasma`,
  `cl_int_hep_ml_min_kg` (scaled by body weight), GFR-derived renal clearance
  and the HIA absorption gate — shared verbatim by the CLI
  (`cli._smiles_spec`) and direct `run_pipeline` consumers, so the chain
  "new chemistry -> in-vivo concentration" has exactly one implementation.
- **TMDD** has a native mass-balance coupling (`PBPKModel.target_binding`,
  §1.4/§2.4): binding-site turnover ODEs live inside the tissue ODEs and
  internalization clears drug directly, with the bound/internalized masses
  part of the drug tally.  The DEFAULT-ENGAGED `feedback_loop` driver (§4.5)
  stays enabled for full runs (one coupled re-run through a fast
  liver/kidney/cardiac pass) and is disabled only in the baseline lane.

---

## Stage 2 — Concentration -> Target Binding

### 2.1 Goal

Given free drug concentrations at target sites, predict fractional target occupancy and bound-complex kinetics over time, plus off-target engagement.

### 2.2 Step 2A: Target Identification

- **Approved/reference drugs:** DrugBank targets for the given drug (if present in DB) — primary targets + known off-targets (hERG, CYP enzymes, BSEP, transporters, nuclear receptors).
- **Novel molecules:** DTI/DTA machine learning (sequence-based AttentionDTA/MINDG-style or structure-based when a structure is available) ranks candidate targets from the safety-critical panel and mechanism hypotheses.
- **Recommended default off-target safety panel** (baseline 2026.9.1): hERG (Kv11.1), CYP3A4/2D6/2C9 (inhibition), BSEP, MRP3/MRP4, OATP1B1, P-glycoprotein, mitochondrial complex I/II/III/IV/MCT, and the glucocorticoid/sex-hormone receptors (endocrine effects).
- **Baseline scope (2026.9.1).** The shipped engine resolves targets for built-in
  benchmark compounds only and defaults to the off-target safety panel with
  class-typical IC50 priors; the DrugBank resolution kernel and the DTI/DTA-ML
  target ranker for novel molecules are planned for later releases (doc/01 §5).
- **Provenance disclosure (2026.9.1).** Every contract carries a machine-readable
  `trust` record (`fidelity_provenance`, doc/12 §7): the wired production
  anchors (R-4/R-6/R-7/R-8 as applicable), the engaged vs off-but-available
  realism terms, the class-prior target sites, the ADMET-AI heads that actually
  contributed, and the set of planned release tracks (P5–P9) — so the resolution
  path used for *this* molecule (measured / re-scored head / class prior) is
  auditable in the report. The end-to-end case
  `case_full_chain_admet_to_report` verifies the full chain per run (doc/08,
  doc/12 §7.2).

### 2.3 Step 2B: Binding Kinetics Parameterization

Priority order for Kd / kon / koff of drug-target pairs:

1. **Measured kinetics/affinity** from ChEMBL, BindingDB, or literature (kon, koff, Kd) — highest priority.
2. **Structure-based AI affinity** (CORDIAL/IPBind-class) mapped from predicted binding energy to Kd when a 3D complex is obtainable (docking or AlphaFold target + RDKit pose via docking).
3. **Sequence/descriptor-based DTA** (AttentionDTA-class) otherwise.
4. **Fallback default** for safety targets: use class-typical IC50/EC50 medians for pharmacology-informed priors and flag low confidence.

> **Off-target resolution — two paths (2026.9.1, doc/12 D10).** For novel
> molecules scored by ADMET-AI without a measured potency, the inhibitor-class
> heads that name a panel site re-score that site through a single
> corpus-calibrated, monotone probability→KD curve
> (`drugos.target.resolver.kd_from_score`):
>
> - *Path B (hERG):* the ADMET-AI hERG-head probability maps onto a continuous
>   KD — never more potent than the 2 nM panel prior (the flagship conservative
>   anchor, reproduced exactly at `P=1`), degrading exponentially toward the
>   1 mM weak floor as confidence falls (`P=0`), with the classifier threshold
>   (`P=0.5`) landing at the corpus-typical weak potency (~1.4 µM, within 10x of
>   the hERG Central median IC50 on the conservative side; doc/08 R-8).
> - *Path A (resolver seam):* the same curve re-binds the Veith CYP2D6/3A4/2C9
>   inhibition heads onto their panel priors; P-gp, the transporters,
>   mitochondrial and endocrine sites — and any primary-target site — stay on
>   class priors until the sequence-DTI resolver (doc/12 §1 row 2b, P5) ships a
>   corpus-calibration + equivalence case under G4/G5.
>
> The effective KD (resolved KD, or a per-compound `qt_ic50_nm` override)
> replaces the hERG site in the safety panel once, so Stage-2 occupancy, the
> QT drive and the exposure anchor all stay mutually consistent
> (`_effective_herg_kd_nm` / `_herg_sieved_panel`). Benchmarks (no ADMET-AI
> call) fall back to the panel prior unchanged.

Target abundance (Rtot per tissue) comes from tissue expression/quantitative proteomics tables; store per tissue.

### 2.4 Step 2C: Occupancy Model

Reversible binding with target turnover (Daryaee & Tonge 2019; classic TO model):

```
dR/dt   = ksyn - ρ * R - kon * D_free * R + koff * DR
dD_free/dt = (delivery from tissue PK) - (clearance) - kon*D_free*R + koff*DR
dDR/dt  = kon * D_free * R - koff * DR - kint * DR
```

- `ρ` = target turnover rate; `ksyn = ρ * R0` sets baseline; `kint` = complex internalization rate.
- Fractional occupancy = DR / Rtot.
- **TMDD (implemented)**: when [drug]/[Rtot] is not large, apply the full TMDD equations or the quasi-steady-state approximation. The native mass-balance coupling (`PBPKModel.target_binding`, doc/05 §1.4) implements exactly the `dD_free` leg above — the net bound flux leaves the tissue pool and `kint·DR` drains into a cleared sink, so cumulative internalized drug mass is conserved-tallied with the administered dose (`case_tmdd_drug_disposition`).

### 2.5 Outputs

- Per-target, per-tissue occupancy-time profiles.
- Time-at-target (occupancy-integral) used by Stage 3 as the "exposure signal".

---

## Stage 3 — Target Binding -> Signaling Pathway (QSP)

### 3.1 Goal

Translate occupancy into downstream pathway activity and relative-to-baseline physiological signals.

### 3.2 Step 3A: Pathway Graph Assembly

For each engaged target, assemble the relevant sub-network from KEGG / Reactome / PANTHER (e.g., RTK -> RAS -> RAF -> MEK -> ERK; PI3K -> AKT -> mTOR; apoptosis/caspase; bile-acid/enterohepatic signaling; calcium/ion-channel signaling; mitochondrial respiration coupling). Prune to a **tractable graph** (20-200 reactions) that covers:
- The intended therapeutic route (efficacy pathway)
- The principal toxicity-relevant routes (see Stage 4 organ panels)

> **Baseline scope (2026.9.1).** The shipped kernel is the Huang/Levchenko
> MAPK cascade (R-5): 20 irreversible mass-action reactions over 22 species
> (RTK -> RAS -> RAF -> MEK -> ERK, dual-phosphorylation readout), compiled
> from the in-repo SBML and executed by ``simulate_sbml_pathway``.  A
> ~12-reaction / ~10-species in-repo DSL compiler is retained as the fast
> path, but the pipeline default is the SBML cascade (the 20-200-reaction
> target).  KEGG/Reactome/PANTHER ingestion and larger toxicity-route graphs
> are planned for later releases (doc/01 §5); the compiler contract and ODE export support
> them unchanged.

### 3.3 Step 3B: ODE Generation (direct SciPy implementation)

Encode reactions with:
- **Mass-action kinetics** for receptor-ligand binding, phosphorylation, PP2A/p65 dephosphorylation, degradation, complex formation.
- **Hill-type laws** for transcription-factor activation, gene-regulatory and feed-back terms, and for enzyme-catalyzed steps (Michaelis-Menten).
- **Conservation/steady-state** initialization from physiological baselines (basal receptor density, kinase levels, substrate pools).

```
dx_i/dt = Σ reactions ± (Hill activation/inhibition) - degradation
```

Model geometry mirrors the MET-QSP precedent (~100 species / ~70 ODEs) and is exported to the shared ODE contract.

### 3.4 Step 3C: Perturbation Simulation

- Input: occupancy-time profile from Stage 2 (drug as perturbation term in receptor mass balance).
- Outputs: activity-time traces of key nodes (e.g., pERK fraction, pAKT, p53, RIPK1/RIPK3, mitochondrial membrane potential proxy, bile-acid pool interference index).
- Baseline normalization: report fold-change relative to drug-free steady state.

### 3.5 Validation Points

- Dose-response shapes (Hill/Emax) must emerge naturally: verify EC50 and Hill slope against literature where available (amplification yields EC50 < Kd).
- Known pathway readouts (e.g., phospho-protein biomarker data) compared if published.

---

## Stage 4 — Signaling Pathway -> Organ Function (QST / Physiology)

### 4.1 Goal

Map pathway/toxicity signals to organ-level functional endpoints. The baseline (2026.9.1) covers liver, cardiovascular, and kidney.

> **CYP450 abundance table (Stage 4 upstream data, implemented):** `HumanPhysiology` exposes per-isoform hepatic microsomal abundances (`CYP_ABUNDANCE_PMOL_MG`, Barter et al. 2013 central immunoquantified values for CYP1A2/2A6/2B6/2C8/2C9/2C19/2D6/2E1/3A4) and the derived per-isoform liver content (`hepatic_cyp_content_nmol`). The per-CYP enzyme-kinetics stage now consumes this table: `cyp_vmax_mg_h` turns a per-isoform liver content into an abundance-scaled capacity (kcat·content·MW), and `PBPKModel.cyp_terms` sums Michaelis-Menten/Hill fluxes per isoform on the unbound liver concentration. Like every realism term it is **DEFAULT-ENGAGED in full runs** (explicit per-CYP terms yield to the lumped MM auto-anchor when absent) and reachable in the `fidelity="baseline"` lane by setting `cyp_terms` explicitly.

### 4.2 Liver (DILI) — QST Sub-Model (DILIsym-style)

Mechanistic sub-models:

1. **Bile-acid transport inhibition** — BSEP/NTCP/MRP3/MRP4 inhibition from Stage 2 safety-target data (IC50); enterocyte/hepatocyte bile-acid pools; basolateral and canalicular efflux; predicts cholestasis and biliary biomarkers. **The cholestasis axis is anchored to the production-validated GCDCA bile-acid PBK of de Bruijn & Rietjens 2024** (paper CC BY 4.0; validated vs clinical cholestasis incidence of ~18 marketed drugs, R-7, doc/12 row 4c): free-hepatic drug competitively inhibits BSEP efflux (`Km_app = Km_BSEP·(1 + C_free/Ki)`, `Ki = IC50/2`) and the intrahepatic pool fold-ratio above the 1.5× risk threshold drives cholestasis.
2. **Mitochondrial dysfunction** — ETC complex inhibition (in-vitro IC50), ATP production shortfall, decrease in mitochondrial membrane potential; with adaptive mitogenesis term.
3. **Oxidative stress** — reactive-oxygen-species generation vs. glutathione buffer; GSH depletion; protein/lipid damage.
4. **Hepatocyte death** — apoptosis (caspase-driven via TNF/ligand pathway) + necrosis (ATP/oxidative-threshold-driven); includes regeneration dynamics.
5. **Immune-mediated component** (DEFAULT-ENGAGED in full runs):
   ``LiverParams.immune_ic50_nm``
   supplies a saturable hapten/danger hazard (``immune_hazard()``, Hill
   sigmoid) driven by free hepatic exposure; the hazard loads an adaptive
   immune-response ODE ``dI/dt = k_recruit·hazard·(1−I) − k_decay·I`` whose
   output level fills the fourth ``combined_stress()`` axis via
   ``immune_weight`` (0.5 in full runs) — auto-anchored
   ``IC50_immune = 0.30·DILI IC50`` (or a disclosed null-effect anchor for
   non-immune-triggering molecules) when absent.  The adaptive-recruitment
   kinetics produce a delayed, persistent immune signal (physiologically
   distinct from the direct cholestatic/mitochondrial/redox axes) whose
   steady-state level ``I_ss = k_recruit·hazard/(k_recruit + k_decay)`` is
   pinned by the analytic-L2 case.

**Inputs.** PBPK liver tissue exposure (`C_liver(t)`, free), plus in-vitro IC50/assay parameters (BSEP inhibition, ETC inhibition, oxidative stress) from literature/ChEMBL/in-vitro consortia data.
**Outputs.** ALT/AST/total-bilirubin serum trajectories; dead-cell fraction; DILI-grade classification (ALT > 3x ULN; Hy's Law criteria). Follows the fezolinetant DILIsym precedent where PBPK exposure + in-vitro toxicity parameters predict liver signal.

> **Implemented (2026.9.1):** `src/drugos/organ/liver.py` — cholestasis (de Bruijn & Rietjens bile-acid PBK, R-7), mitochondrial (ETC, adaptive mitogenesis), redox/GSH, hepatocyte death (sigmoid kill + regeneration) ODEs; `LiverParams` carries compound-specific IC50s, `liver_params_from_panel()` wires the Stage-2 safety panel. **Pathway→organ closure:** the Stage-3 ERK/MAPK readout fold-change is the hepatocyte `proliferation_signal` that scales `regeneration_1h` via `regeneration_scale()` (damped `regen_pathway_gain=0.5`, clamped to `[regen_pathway_min=0.5, regen_pathway_max=1.5]`, so an off-target pathway suppression can attenuate but never ablate regeneration); the reported total-bilirubin rise is capped at `bile_rise_max_fold`×ULN while the un-clamped excess stays in the mechanical `stress`. **Immune-mediated DILI axis (DEFAULT-ENGAGED in full runs):** `immune_hazard()` (Hill sigmoid) converts free hepatic exposure to a hapten/danger hazard driving a recruitment/decay adaptive immune-response ODE (2nd liver state, `immune_recruit_1h`/`immune_decay_1h`) whose output loads the fourth `combined_stress()` axis through `immune_weight` (0.5 in full runs; the auto-anchor `IC50_immune = 0.30·DILI IC50`, or a disclosed null-effect anchor, is set whenever the term is absent).  The adaptive response produces a delayed, persistent signal whose steady state `I_ss = k_recruit·hazard/(k_recruit + k_decay)` is pinned by the analytic-L2 case, and the 2-state ODE feeds the validated hepatocyte-death dynamics with the immune contribution active only when `immune_ic50_nm` and `immune_weight > 0` are both set. Validated: acetaminophen overdose reproduces ALT > 3x ULN with Hy's Law while therapeutic dosing stays grade 0/1, and the regeneration-coupling + bilirubin ceiling are pinned by L2 cases.

### 4.3 Cardiovascular — Lumped Circulation + Electrical Axis

- **Lumped circulation model** (Physiome-style minimal hemodynamic system with ventricular interaction and valve dynamics): cardiac output, stroke volume, mean arterial pressure, central venous pressure as state variables.
- **hERG/QT axis**: from hERG blockade fraction (Stage 2 off-target), estimate IKr reduction and a QTc-prolongation model (multiplicative / Emax on hERG channel current; literature QTc-hERG relationships). Maps to torsades-de-pointes risk band.
- **Inotropy/chronotropy modulators**: beta/catecholamine pathway effects feed into contractility and rate.

> **Implemented (2026.9.1):** `src/drugos/organ/cardiac.py` — `predict_qtc()` Emax hERG->IKr->QTc axis with TdP banding (450/480/500 ms) and a two-compartment Windkessel (`simulate_hemodynamics()`) for MAP/CVP/CO/SV. The hERG blockade fraction is driven by the **PBPK cardiac (heart) free exposure** (`unbound_tissues["heart"]`), not the hepatic one. **Inotropy/chronotropy coupling (pathway→organ):** the Stage-3 ERK/MAPK amplification ratio gates `CardiacParams.inotropy = chronotropy` via `_cardiac_tone_scale()` (damped `gain=0.2`, capped at 1.3×); the ERK proxy encodes only the *stimulatory* branch (ERK1/2 is downstream of beta-adrenergic E-C coupling), so a baseline or suppressed readout keeps tone exactly 1.0 and the loop `co_fraction` stays neutral for benign drugs, while a genuinely amplified readout raises CO/HR. **Sympathetic suppression branch (DEFAULT-ENGAGED in full runs):** `CardiacParams.sympathetic_tone` (1.0 in the baseline lane) realizes the roadmap suppression axis as an explicit, saturable Emax (beta-adrenergic site blockade) via `sympathetic_tone_from_emax()` — `tone = IC50/(IC50 + C_free_heart)`, wired from the PBPK heart free exposure when `RunSpec.beta_block_ic50_nm` is set (a disclosed null-effect anchor is auto-assigned in full runs when absent). The tone multiplies both heart rate and stroke volume, so cardiac output scales with tone²; systemic resistance is pinned to the intact-tone reference output, so a suppression that lowers CO reads out as hypotension (arterial-venous pressure drop falls with tone², `pa - pv = tone²·(MAP_ref − CVP)` in the Windkessel steady state), which feeds the loop `co_fraction` and perfusion scaling. Validated: 0.5 mg dofetilide peak Delta-QTc 20 ms lies in the published prolongation band vs a negligible-hERG control (L3 case); the sympathetic axis is pinned by the IC50-exposure quarters-CO analytic L2 case.

### 4.4 Kidney — Nephron + GFR Model

- **Nephron-level model** (Physiome neural-nephron lineage, planned release P7): glomerular filtration, tubular reabsorption/secretion; clearance coupling with Stage 1 renal elimination.
- **Nephrotoxicity endpoints**: acute kidney injury proxies (GFR decline, tubular injury biomarker (KIM-1 heteromer in practice)) — the baseline keeps serum creatinine + GFR from renal function sub-model.

> **Implemented (2026.9.1):** `src/drugos/organ/kidney.py` — nephron injury sigmoid drives a floored GFR; serum creatinine from the closed-form balance Scr = P/GFR; KDIGO AKI stage from Scr ratio/GFR drop. **The baseline GFR is anchored to the CKD-EPI 2021 race-free creatinine equation** (R-6, doc/12 row 4b) whenever a measured serum creatinine is carried on the profile (`ckdepi_2021_egfr`, BSA-scaled via Mosteller). **KIM-1 tubular biomarker:** the proximal-tubule injury signal is translated to a graded urinary KIM-1 xUNL row (`summarize_clinical_kidney`: `1 + 8·injury`) with CTCAE-style thresholds (1.5/3/5/10 xUNL), so nephrotoxicity is also surfaced by a tubule-specific marker, not only by GFR/creatinine. Validated: exact Scr=P/GFR at zero exposure, monotonic KDIGO escalation, and CKD-EPI reference points / pipeline wiring (L2 cases).

### 4.5 Feedback to PK

Organ dysfunction feeds back to PK: reduced hepatic clearance (liver), reduced GFR (kidney), reduced cardiac output (perfusion-limited distribution). Implement in monolithic mode directly; in sequential mode via a slow-timescale outer loop (recompute PBPK with updated organ parameters once or twice).

> **Implemented (2026.9.1):** `src/drugos/organ/feedback.py` — `feedback_from_results()` maps dead fraction / GFR/CO factors to `OrganFeedback`; `apply_pk_scaling()` returns a PBPK model copy with scaled hepatic/renal clearances and perfusion (hepatic/renal/CO), leaving the original untouched. The cardiac-output factor is **derived from the simulated hemodynamics**: `co_fraction = CO_simulated / CO_reference` (from `CardiacResult.co_l_min`), so a reduced-cardiac-output physiology scales perfusion-limited distribution on the re-run (neutral 1.0 with intact cardiac panels).

---

## Stage 5 — Organ Function -> Clinical Phenotype + Toxicity

### 5.1 Biomarker Translation & Grading

- Map organ outputs to clinical-grade biomarkers with units and reference ranges (ALT, AST, bilirubin, GFR, creatinine, **urinary KIM-1**, HR, BP, QTc).
- Grade each (CTCAE-style Grade 0-4) with severity thresholds; produce time-to-onset and duration.

> **Kidney row set.** `summarize_clinical_kidney` grades GFR (mL/min), serum
> creatinine (xULN) and urinary KIM-1 (xUNL, derived from the nephron injury
> fraction) so AKI is surfaced both by the functional markers and by a
> proximal-tubule injury biomarker (`BIOMARKERS["KIM_1"]`).

> **Spec-pinned grading.** Biomarker grading threads the *spec-driven*
> constants (e.g. `cns_ic50_nm`, the effective hERG KD) into the organ
> summary instead of re-instantiating model defaults (`_grade_organ` carries
> the `CnsParams` used by the simulation), so the reported margins always
> match the simulated exposure-ratio anchors.

### 5.2 Composite Toxicity Scoring

Three independent evidence lines combined with uncertainty:

1. **Mechanistic QST outputs** (Stage 4) — mechanistic DILI/QT/AKI signals, driven by PBPK exposure + in-vitro parameters.
2. **Exposure-ratio score** — validated formula from literature: `Cmax_unbound / in-vitro toxicity_IC50` (e.g., lowest lethal/toxicity concentration across assays); flag if ratio exceeds established thresholds (ROC AUC ~0.91-0.96 in 241-drug DILI analyses).
3. **Structural ADMET flags** — ADMET-AI predictions (hERG blockade probability, AMES, hepatotoxicity probability, DILI probability) as independent structural-molecule evidence.

**Fusion.** Shipped implementation (`drugos/clinical/toxicity.py`): each endpoint derives a probability from its evidence lines —
`P_line = sigmoid(slope * (log10(Cmax_free/IC50) - ratio_center))` for the exposure line, `RULES[].grade_probs[grade]` for the mechanistic line — and the lines are combined in log-odds against a per-endpoint class prior,
`logit(P) = logit(prior) + Σ w_i (logit(P_i) - logit(base_i))`, with a Beta CI shrunk by the sum of evidence weights. Per-endpoint priors DILI 0.25 / QT 0.20 / AKI 0.22 / CNS 0.20; weights MECH 1.2 / EXP 1.0 / STRUCT 0.7.

**Four endpoints (Stage 4 + CNS).**

| Endpoint | Evidence enabled by | Note |
|---|---|---|
| DILI | liver organ + `dili_ic50_nm` (in-vitro override or class prior via BSEP) | acetaminophen class prior gives the DILI anchor |
| QT / TdP | cardiac organ + `qt_ic50_nm` (hERG; dofetilide override vs negligible-hERG defaults) | warfarin/hERG-default compounds stay flat |
| AKI | kidney organ + pharmacophore flag | KDIGO-style GFR/creatinine escalation |
| CNS | `brain_free_exposure` from `organ/cns.py` + `cns_ic50_nm` | **only fused when the molecule is CNS-anchored**; otherwise the 0.20 class prior holds (default 100 µM would wrongly flag APAP) |

Grading: CTCAE-like Grade 0-4 ladders per biomarker (`grade_value`), time-to-onset + duration from threshold crossing (`crossing_interval`). The overall composite risk is resolved to a **low / elevated / High verdict** with the dominant driver named via `overall_driver()`.

### 5.3 Clinical Phenotype Report

Structured output:
- Predicted concentration-time and PK metrics
- Target occupancy & pathway activity profiles
- Predicted physiological indicator trajectories (each graded + CNS brain-free exposure curve)
- Adverse-event risk table per endpoint (DILI, QT prolongation/TdP, nephrotoxicity, CNS effects) with mechanism attribution and confidence
- Sensitivity/driver analysis (which parameters drive the risk margins)

---

## Phase 6 — Uncertainty & Decision Layers (D21-D24)

Shipped implementations in `drugos/robustness/` (all deterministic via fixed seeds):

1. **D21 Parameter ensemble** (`uncertainty.py`) — multiplicative log-normal perturbation (CV: regime-driven by default, see doc/12 §7.3 — `EnsembleConfig.cv = None` resolves to the run's predictive-regime `band_cv`, so a weaker-evidence run automatically sweeps wider; an explicit `cv` overrides; the fall-back scalar default is 0.30) of the PK/potency scalars (`cl_hep`, `cl_renal`, `fup`, `qt/dili/cns_ic50`); each member is a full pipeline run; 90% percentile bands (`band()`, q5/q50/q95) for plasma/organ curves, endpoint-risk quantiles, and per-member verdict counts.
2. **D22 Virtual cohort** (`population.py`) — anthropometric sampling (sex/age/height via BMI, deterministic seed); one full pipeline run per individual at the same dose; `risk_quantiles()`, grade ≥ 1 / ≥ 2 `incidence()`, verdict counts.
3. **D23 Sensitivity & drivers** (`sensitivity.py`) — `local_sensitivity` one-at-a-time ±10% normalized log-sensitivities d ln y / d ln p with a ranked `drivers` list; `run_sobol_sensitivity` Saltelli first/total indices (scrambled Sobol design, independent evaluator injectable for CI speed).
4. **D24 Prospective rerun** — reproducibility runbook (identical fixed-seed runs bit-equal) + held-out profile/dose stability; enforced by registration as validation cases (doc/08; see the D21-D24 case in `validation/cases/`).

The Playground exposes D21-D23 as UI toggles (`with_uncertainty` / `with_population`
/ `with_sensitivity`); the CNS brain-free-exposure curve is always plotted in the
organ-panel view and drives the CNS composite risk line only when the molecule is
CNS-anchored (un-gated structural bands are excluded, §5.2).

---

## Monolithic vs Sequential: Recommended Default

- **Default execution = sequential** for the baseline's rigor: run Stage 1 fully, ingest profiles into Stage 2, then Stage 3 and 4, then Stage 5. Simpler to debug and validate per stage.
- **Monolithic execution** optional for small pathway models (single ODE solve); enables full feedback. Enable when the pathway graph <= ~30 nodes.

This staged structure also cleanly maps to test/validation plans in `08-validation-and-risk.md`.