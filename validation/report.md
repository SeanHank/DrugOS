# DrugOS Validation Report

- Status: **PASS** (25/25 cases passed)
- Generated: 2026-09-07 05:02 UTC
- DrugOS version: 2026.9.0
- Interpreter: /opt/anaconda3/envs/drug_os/bin/python
- Fold-error allowance: within 2x of the published band centre (doc/08 Tier 2, GMFE <= 2); Fa bands additionally clamp to 1.

## Results

| Case | Level | Metric | Predicted [allowed band] | Criterion |
|---|---|---|---|---|
| midazolam | L3 | cl_plasma_l_h | 26.86 [12.73, 50.91] L/h | pass |
| acetaminophen | L3 | cl_plasma_l_h | 24.35 [10.06, 40.25] L/h | pass |
|  |  | t_half_h | 2.462 [1.146, 4.583] h | pass |
|  |  | f_abs | 0.867 [0.4427, 1] fraction | pass |
| warfarin | L3 | cl_plasma_l_h | 0.232 [0.0866, 0.3464] L/h | pass |
|  |  | vss_l | 7.035 [4, 16] L | pass |
|  |  | t_half_h | 23.77 [18.37, 73.48] h | pass |
| ciprofloxacin | L3 | cl_plasma_l_h | 23.49 [10.06, 40.25] L/h | pass |
|  |  | t_half_h | 2.562 [2.121, 8.485] h | pass |
|  |  | f_abs | 0.8611 [0.3775, 1] fraction | pass |
|  |  | urine_fraction | 0.4042 [0.2646, 1] fraction | pass |
| dofetilide | L3 | cl_plasma_l_h | 16.92 [7.746, 30.98] L/h | pass |
|  |  | t_half_h | 5.508 [3.969, 15.87] h | pass |
|  |  | f_abs | 0.8611 [0.4472, 1] fraction | pass |
|  |  | urine_fraction | 0.5666 [0.3571, 1] fraction | pass |
| mass-balance (midazolam IV, no elimination) | L1 | mass_budget_rel_error | 1.066e-15 [0, 1e-06] fraction | pass |
| dose-proportionality (midazolam IV) | L1 | auc_40mg_over_10mg | 4 [3.6, 4.4] ratio | pass |
| one-compartment analytic limit | L2 | auc_inf_mgh_l | 10.02 [9.5, 10.5] mg.h/L | pass |
|  |  | t_half_h | 46.85 [41.44, 50.65] h | pass |
| occupancy equilibrium (drug=Kd) | L2 | occupancy_at_free_equals_Kd | 0.5 [0.48, 0.52] fraction | pass |
| pathway amplification (MAPK EC50<Kd equivalent) | L2 | emax_over_baseline_fold | 714.1 [2, 1e+04] ratio | pass |
|  |  | ec50_signal | 0.0008154 [0, 0.5] fraction | pass |
|  |  | drug_ec50_nm_at_kd_1 | 0.0008161 [0, 0.5] nM | pass |
| liver DILI dose-response (acetaminophen) | L2 | therapeutic_peak_alt_uln | 1.064 [0, 2] xULN | pass |
|  |  | overdose_peak_alt_uln | 3.313 [3, 20] xULN | pass |
|  |  | therapeutic_to_overdose_escalation | 3.113 [2, 1000] fold | pass |
|  |  | overdose_hy_law | 1 [1, 1] flag | pass |
| cardiac QTc prolongation (dofetilide hERG) | L3 | dofetilide_delta_qtc_ms | 20.23 [5, 55] ms | pass |
|  |  | warfarin_delta_qtc_ms_control | 0.01702 [0, 1] ms | pass |
| kidney GFR/AKI escalation (KDIGO) | L2 | scr_at_zero_fn | 1 [1, 1] ratio | pass |
|  |  | low_exposure_scr_ratio | 1.001 [1, 1.05] ratio | pass |
|  |  | high_exposure_scr_ratio | 4.542 [2, 8] ratio | pass |
|  |  | aki_grade_escalation | 3 [1, 3] grade | pass |
| Stage-5 clinical grading: analytic point-matches | L2 | exposure_line_probability | 0.8264 [0.8264, 0.8264] P(risk) | pass |
|  |  | prior_only_posterior | 0.25 [0.25, 0.25] P(risk) | pass |
|  |  | mechanistic_grade3_probability | 0.45 [0.45, 0.45] P(risk) | pass |
|  |  | grade_ladder_higher | 2 [2, 2] grade | pass |
|  |  | grade_ladder_lower | 2 [2, 2] grade | pass |
|  |  | crossing_onset_duration | 1 [1, 1] h | pass |
| Corpus calibration cross-check (R-2) | L3 | dofetilide_herg_measured_geomean_nm | 26.4 [0, 100] nM | pass |
|  |  | model_prior_to_measured_ratio | 0.07575 [0.05, 20] fold | pass |
|  |  | herg_corpus_median_pct_inh_at_1uM | 8.254 [0, 30] % inhibition | pass |
|  |  | herg_corpus_p99_9_pct_inh_at_1uM | 31.64 [25, 100] % inhibition | pass |
| Stage-5 composite risk ordering vs clinical anchors | L3 | dofetilide_qt_risk | 0.5259 [0.45, 0.95] P(risk) | pass |
|  |  | warfarin_qt_risk | 0.0053 [0, 0.35] P(risk) | pass |
|  |  | acetaminophen_20g_dili_risk | 0.97 [0.6, 1] P(risk) | pass |
|  |  | acetaminophen_1g_dili_risk | 0.1177 [0, 0.4] P(risk) | pass |
|  |  | unanchored_cns_prior | 0.2 [0.2, 0.2] P(risk) | pass |
|  |  | held_out_dose_profile_stable | 0 [0, 0] disagreements | pass |
| Phase-6 robustness engines: D21-D24 self-consistency | L1 | ensemble_reproducible_max_band_diff | 0 [0, 0] au | pass |
|  |  | band_monotonic_violations | 0 [0, 0] points | pass |
|  |  | population_min_risk_nonneg | 0.1132 [0, 0] P(risk) | pass |
|  |  | sobol_first_total_in_range | 0 [0, 0] flags | pass |
|  |  | dili_ic50_sensitivity_sign | -0.4624 [-1, 0] dlnR/dlnIC50 | pass |
| SC/IM depot analytic (Bateman single pool) | L2 | cmax_mg_l | 0.123 [0.1171, 0.1294] mg/L | pass |
|  |  | tmax_h | 12 [11.34, 13.86] h | pass |
|  |  | auc_inf_mgh_l | 18.1 [17.1, 18.9] mg.h/L | pass |
|  |  | unabsorbed_feces_mg | 1 [0.9, 1.1] mg | pass |
| D24 prospective rerun fidelity (dofetilide QTc) | L1 | test_retest_max_risk_diff | 0 [0, 0] P(risk) | pass |
|  |  | held_out_subject_qt_risk | 0.8345 [0.7, 0.99] P(risk) | pass |
|  |  | held_out_verdict_disagreements | 0 [0, 0] count | pass |
| R literature-PK cross-check (R-1) | L3 | r_literature_cl_agreement_max | 1.931e-15 [0, 0.02] fraction | pass |
|  |  | r_verdicts_not_agree | 0 [0, 0] count | pass |
|  |  | r_two_comp_fits | 0 [0, 5] count | pass |
| cardiac AP cross-check (ORd/IKr) | L3 | ord_apd90_base_ms | 266.3 [200, 350] ms | pass |
|  |  | dofetilide_at_ic50_delta_apd90_ms | 114.8 [30, 1000] ms | pass |
|  |  | apd90_monotone_25_to_50_ms | 67.7 [0, 1000] ms | pass |
|  |  | warfarin_control_delta_apd90_ms | 0 [-1e-06, 1e-06] ms | pass |
| ADMET-AI BBB_Martins -> CNS partition (R-4) | L2 | cns_kpu_brain_penetrant | 1 [1, 1] ratio | pass |
|  |  | cns_kpu_brain_non_penetrant | 0.2 [0.2, 0.2] ratio | pass |
|  |  | benchmark_no_admet_kpu_brain | 1 [1, 1] ratio | pass |
|  |  | non_penetrant_peak_brain_free_ratio_of_penetrant | 0.2 [0.2, 0.2] ratio | pass |
|  |  | non_penetrant_cns_grade_lte_penetrant | 0 [-1e+09, 0] grade | pass |
| Huang/Levchenko SBML MAPK cascade integration (R-5) | L2 | reactions_count | 20 [20, 20] count | pass |
|  |  | species_contain_ppk | 1 [1, 1] flag | pass |
|  |  | drug_free_ppk_steady_state | 0.9818 [0.5, 2] conc | pass |
|  |  | full_occupancy_ppk_suppression | 0 [0, 0.05] conc | pass |
|  |  | monotone_inhibition_in_occupancy | 1 [1, 1] flag | pass |
|  |  | inhibited_fold_change_lte_1 | 0 [-1, 1] fold | pass |
|  |  | readout_is_dual_phospho_erk | 1 [1, 1] flag | pass |
| CKD-EPI 2021 race-free GFR baseline (R-6) | L2 | egfr_male_60_scr1_pct2_5 | 86.16 [85, 92] mL/min/1.73 m2 | pass |
|  |  | egfr_female_lower_same_scr | 64.5 [0, 86.16] mL/min/1.73 m2 | pass |
|  |  | egfr_bsa_scaled_absolute | 98.22 [86.16, 112] mL/min | pass |
|  |  | profile_scr_gfr_uses_ckdepi | 98.22 [98.12, 98.31] mL/min | pass |
|  |  | profile_without_scr_keeps_default | 125 [112.5, 125] mL/min | pass |
| Bile-acid cholestasis PBK (R-7) | L2 | baseline_no_drug_fold | 1 [0.995, 1.005] fold | pass |
|  |  | cholestatic_ritonavir_stress | 0.9975 [0.5, 1] 0..1 (1.5x threshold = 0.5) | pass |
|  |  | benign_itraconazole_fold | 1 [1, 1.15] fold | pass |
|  |  | ranking_ki_monotone | 9.321 [9.321, 9.321] fold delta | pass |
|  |  | ki_conversion_ic50_over_2 | 0.1 [0.1, 0.1] umol/L | pass |
|  |  | organ_wiring_chol_consistency | 0.9975 [0.9475, 1.048] 0..1 | pass |

## Evidence levels

Results are graded by how much epistemic weight they carry (doc/08 §1.1-1.4 tier ladder, strongest first):

| Level | Meaning | Basis | What it certifies |
|---|---|---|---|
| **L3** | Empirically anchored (Tier 1) | Predicted vs published human clinical ranges (USPI / literature) under the 2x GMFE allowance. | The pipeline reproduces clinically observed human PK within the fold allowance — the strongest evidence in this suite. (doc/08 §1.1) |
| **L2** | Analytic / mechanistic limit (Tier 2) | Closed-form or mechanistically forced point-match derived from the stage's own equations. | The stage ODEs solve their intended dynamics correctly; it does not by itself certify human predictivity (needs L3). (doc/08 §1.2) |
| **L1** | Internal consistency / CI (Tier 3) | ODE / compiler self-consistency: mass conservation and linearity axioms, no external data. | Numerical correctness of the solver and mass bookkeeping; weakest in epistemic weight — 'just computes right'. (doc/08 §1.3-1.4) |

Per-level status:

- **L3** (Empirically anchored (Tier 1)): 10/10 cases green.
- **L2** (Analytic / mechanistic limit (Tier 2)): 11/11 cases green.
- **L1** (Internal consistency / CI (Tier 3)): 4/4 cases green.

## Notes & limitations

- **midazolam**: reported: CL=26.9 L/h, Vss(MRT)=192 L, t1/2=5.11 h Vss/t1/2 reported but not asserted: lumped R&R partition overpredicts the apparent Vss of low-fup lipophilic bases (doc/08 risk #4); clearance axis is the validated output.
- **acetaminophen**: reported: CL=24.3 L/h, Vss(MRT)=118 L, t1/2=2.46 h
- **warfarin**: reported: CL=0.232 L/h, Vss(MRT)=7.04 L, t1/2=23.8 h
- **ciprofloxacin**: reported: CL=23.5 L/h, Vss(MRT)=119 L, t1/2=2.56 h
- **dofetilide**: reported: CL=16.9 L/h, Vss(MRT)=152 L, t1/2=5.51 h
- **mass-balance (midazolam IV, no elimination)**: no-elimination IV bolus; end-state body mass 5 mg vs dose 5 mg (rel. err 1.07e-15)
- **dose-proportionality (midazolam IV)**: AUC10=0.372, AUC40=1.489 mg.h/L (linearity ~4.00)
- **one-compartment analytic limit**: V=66.4 L; analytic AUC=10.00 mg.h/L, t1/2=46.0 h
- **occupancy equilibrium (drug=Kd)**: steady-state occupancy 0.500; analytic D/(D+Kd) = 0.5
- **pathway amplification (MAPK EC50<Kd equivalent)**: r0=0.01322, emax=9.437, ec50=0.0008154, hill=2.15; EC50 signal 1e-3 << half-saturation 0.5
- **liver DILI dose-response (acetaminophen)**: therapeutic ALT 1.06xULN (safe), 20 g overdose ALT 3.31xULN with total bilirubin 2.62xULN -> Hy's Law met; DILI grade 0->3. Acetaminophen overdose injury is characterized by centrilobular necrosis, transaminase >3x ULN and mixed cholestasis (DILI severity scales, e.g. Maria & Victorino / Hy's Law criteria).
- **cardiac QTc prolongation (dofetilide hERG)**: 0.5 mg dofetilide peak Delta-QTc 20.2 ms vs published ~20-60 ms/QTc-prolonging clinical band (low-nM hERG block); warfarin control 0.02 ms.  Published: dofetilide (Tikosyn) USPI lists QT/QTc prolongation; peak Delta-QTc in the 0.5 mg single-dose range is around 10-35 ms and TdP aggregates in QTc > 500 ms.
- **kidney GFR/AKI escalation (KDIGO)**: closed-form Scr=P/GFR exact at zero exposure; 1.0 mg/L free kidney exposure -> Scr ratio 4.54 (KDIGO stage 3) vs 1.001 (stage 0); GFR floor 26 mL/min.  KDIGO criteria: Scr x2 -> stage 2, x3 -> stage 3 (or GFR drop).
- **Stage-5 clinical grading: analytic point-matches**: Exposure ROC line reproduces the closed form sigmoid(1.2*(0.0 - (-1.3))) = 0.8264; empty-evidence fusion returns the DILI prior 0.25 exactly; grade ladder and crossing windows match the CTCAE conventions of doc/05 5.1-5.2.
- **Corpus calibration cross-check (R-2)**: dofetilide ChEMBL hERG IC50 geomean 26.4 nM (core rows, outlier >=10 uM excluded; regen by scripts/data/fetch_chembl_herg.py); model class prior 2.0 nM = 0.1x of measured (conservative direction, within the 20x envelope). hERG Central corpus: 306893 PMID-anchored rows; median %inhibition at 1 uM = 8.3, P99.9 = 31.6 (blockade is the exception, so a per-compound hERG override is the honest modelling choice).
- **Stage-5 composite risk ordering vs clinical anchors**: dofetilide QT 0.526 (qt-driven) > warfarin QT 0.005; APAP 20 g DILI 0.970 > 1 g DILI 0.118 (dili-driven); unanchored CNS sits on the 0.20 class prior.  Published anchors: dofetilide (Tikosyn) is a QT-prolonging hERG blocker and is contraindicated with renal/QT risk; massive acetaminophen overdose causes centrilobular hepatic necrosis (DILI), while warfarin is not a QT liability.
- **Phase-6 robustness engines: D21-D24 self-consistency**: Fixed-seed D21 ensemble reproduces itself exactly (max median-band diff 0); 90% band monotone with 0 violations; D22 cohort incidence non-negative; D23 first/total indices inside [-1,1]/[0,1]; DILI risk strictly decreases with a rising IC50 (-0.4624 per +10% IC50).
- **SC/IM depot analytic (Bateman single pool)**: V=66.4 L, ka=0.3/h, F=0.9; analytic Cmax=0.123 mg/L, Tmax=12.6 h, AUC=18.0 mg.h/L
- **D24 prospective rerun fidelity (dofetilide QTc)**: Repeated identical runs agree to 0.0e+00 in risk and keep verdict 'High composite risk (53%, driver qt)'; an independent female-70 profile also sustains the high-QT regime (dofetilide QT 0.835, driver qt).  Basis: reproducibility is the precondition of the runbook; the QTc band itself is anchored by the L3 dofetilide Tier-1 case (see case_cardiac_qtc).
- **R literature-PK cross-check (R-1)**: worst |CL_r - CL_py|/CL_py over midazolam / acetaminophen / warfarin / ciprofloxacin / dofetilide: 1.93e-15 (acetaminophen); midazolam=r:agree / acetaminophen=r:agree / warfarin=r:agree / ciprofloxacin=r:agree / dofetilide=r:agree; method-of-residuals two-comp fits: 0/5
- **cardiac AP cross-check (ORd/IKr)**: ORd 2011 (myokit, endo, 50 pre-paces @1 Hz): baseline APD90=266.3 ms; delta-APD90 @25% block=47.1 ms, @50% block (measured-IC50 concentration)=114.8 ms, control (0% block)=0.0000 ms; monotone +67.7 ms between block levels; warfarin control confirmed zero prolongation (matches encoder ordering warfarin 0.017 ms << dofetilide 20.2 ms)
- **ADMET-AI BBB_Martins -> CNS partition (R-4)**: BBB_Martins P=0.9 -> kpu_brain 1.00, P=0.1 -> 0.20 (source: ADMET-AI BBB_Martins head); restricted brain peak = 20% of penetrant; CNS grades 0 <= 0; baseline (no ADMET-AI) kpu=1.00 untouched — benchmark anchors unchanged
- **Huang/Levchenko SBML MAPK cascade integration (R-5)**: parsed 22 species / 20 reactions from BIOMD0000000009 (volume 4.0e-12 L); drug-free PP_K steady state 0.982; occupancy monotone 0.6->0.917, 0.9->0.005, 1.0->0.000; full-signal fold-change 0.000 (inhibition).
- **CKD-EPI 2021 race-free GFR baseline (R-6)**: CKD-EPI 2021 race-free: male 60 y SCR 1.0 -> 86.2 mL/min/1.73 m2 (CKD-2 band); female same Scr 64.5; BSA-scaled 98.2 mL/min; Scr-carrying profile gfr=98.2 mL/min; no-Scr default untouched (125.0 mL/min).
- **Bile-acid cholestasis PBK (R-7)**: de Bruijn & Rietjens (2024) bile-acid PBK reproduced at 1 uM free-hepatic exposure: ritonavir-class (IC50 0.2 uM) fold 10.32x / stress 1.00 (cholestatic), itraconazole-class (IC50 10 mM) fold 1.00x (benign); Ki=IC50/2 pinned; organ cholestasis 1.00 matches the submodel.

## Tier-1 geometric-mean fold error (L3 asserted metrics)

| Benchmark | GMFE |
|---|---|
| midazolam | 1.06x |
| acetaminophen | 1.19x |
| warfarin | 1.33x |
| ciprofloxacin | 1.36x |
| dofetilide | 1.21x |

## Tier coverage (doc/08)

- Stage 1 (PK): benchmark compounds + analytic limit + mass budget + dose-proportionality — **green**.
- Stage 2 (occupancy): target-turnover equilibrium ODE vs analytic D/(D+Kd) point-wise match — **green**.
- Stage 3 (pathway): 3-tier MAPK amplifier — steady-state EC50 below the receptor-Kd-equivalent signal (Emax/Hill fit, EC50<0.5) and >2x baseline amplification — **green**.
- Stage 4 (organ): liver DILI dose-response (ALT/bilirubin/Hy's Law at overdose), cardiac QTc prolongation vs the published dofetilide Delta-QTc band, and kidney GFR/AKI KDIGO escalation — **green**.
- Stage 5 clinical / report / CLI tiers: scheduled with their stage modules (see roadmap); each new model feature must add a validation case before merge (G4 rule).

## Citations

- **midazolam**: IV PK summary: Greenblatt D.J. et al. (Clin Pharmacokinet 1984); Goodman & Gilman's Pharmacological Basis of Therapeutics, 13th ed.; Versed (midazolam) Injection USPI PK table. Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph.
- **acetaminophen**: Goodman & Gilman 13th ed.; Tylenol (acetaminophen) USPI: CL ~0.24 L/h/kg, Vss ~0.7-0.9 L/kg, t1/2 ~1.5-3 h, Fa ~0.85-0.95. Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph.
- **warfarin**: Goodman & Gilman 13th ed.; Coumadin (warfarin) USPI: high plasma protein binding (fup ~1%), small Vd (8-14 L), long t1/2. Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph.
- **ciprofloxacin**: Cipro (ciprofloxacin) USPI; Bergan T. (J Antimicrob Chemother 1986): CL_oral ~21.6 L/h, ~60-70% renal, Fa ~0.7, t1/2 ~3-6 h. Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph.
- **dofetilide**: Tikosyn (dofetilide) USPI: fup 0.36, CL ~16.3 L/h (~80% renal), Vss 3.3-4.7 L/kg, t1/2 ~7-9 h, Fa ~0.9; Smith D.A. et al. (Br J Clin Pharmacol 1992). Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph.
- Partition & whole-body graph: Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph..
- Occupancy model: Daryaee F., Tonge P.J. (Annu Rev Pharmacol Toxicol 2019;59:507-529) pharmacometric target-turnover occupancy model; see doc/05 2.4..
- Pathway amplification: Huang C.-Y., Ferrell J.E. Jr. (Proc Natl Acad Sci USA 1996;93:10078-83) ultrasensitivity and signal amplification in the MAPK cascade; DrugOS 3-tier MAPK with basal+signal arms, doc/05 3.2-3.5..

This report is machine-generated by `validation/run_validation.py`; hand edits are overwritten.