# PKPD-02: numerical verification for a one- versus two-compartment amikacin forecasting study

**Status: research in progress. Solver and event-handling verification has been done on synthetic data. No population model has been fitted, no patient-held-out forecast has been made, and nothing here is a clinical or dosing result.**

**Verification is not validation.** What is tested here is whether the code solves the model equations and handles dosing events correctly, judged against an independent calculation. Whether either model describes amikacin pharmacokinetics or forecasts concentrations well has not been investigated.

PKPD-02 is a pre-specified, exploratory plan to ask whether the structure of a population pharmacokinetic (PopPK) model, one compartment or two, changes how reliable the **90 % posterior predictive interval** is for a *future* amikacin concentration in therapeutic drug monitoring (TDM). Before any real forecasting is allowed, the numerical machinery (ODE solving, infusion and covariate handling, patient-level event construction) has to be shown correct. This repository holds that verification work, the protocol it implements, and the evidence behind each result, including a gate that **failed** and was fixed in the open.

Independent research project by Don Aborah ([@pkaysantana](https://github.com/pkaysantana)), developed with substantial assistance from AI coding assistants (Claude Code and Devin); see [Authorship, AI assistance and provenance](#authorship-ai-assistance-and-provenance).

| | |
|---|---|
| **Question** | Does one- vs two-compartment structure change the accuracy, calibration and sharpness of future-concentration predictive intervals? |
| **Done** | Gate S-1 (**failed**), remediation S-1r (**conditional pass**), S-2 preparation (event representation, synthetic batteries, structure-only audit, static fit-interface inspection). |
| **Not done** | Population estimation, posterior forecasting, scoring, leakage audit on real folds, any use of measured concentrations. |
| **Stack** | R only (R 4.5.2; `rxode2` / `nlmixr2`, pinned in `renv.lock`). |
| **Data** | Public Strathclyde amikacin dataset, CC BY 4.0. **Not included here** (see [Provenance](#provenance-references-and-responsible-use)). |

**A reviewer's first five minutes:** read §3 below, then `reports/S1_verification_report.md` (the failure and its diagnosis), `reports/S1r_verification_report.md` (the remediation), and `R/pkpd02_reference_matexp.R` (the independent reference solver). `Rscript tools/verify_public.R` reproduces the four synthetic gates without data (§4).

## 1. Research question and motivation

Bayesian PopPK models turn one or two TDM samples into an individual prediction. A point forecast can be accurate while its interval is too narrow (overconfident) or uselessly wide, so the interval is what matters for decisions. The original Strathclyde analysis kept a one-compartment structure for amikacin; a later exploratory reuse of the same data preferred two compartments on objective function value, for a different purpose. That disagreement motivates a controlled comparison; it does not show which structure forecasts better.

The proposed comparison (`protocol/PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md`) is **paired and patient-held-out**: each model is fitted in five patient-level folds that exclude the held-out patient; one earlier sample updates the held-out patient's Bayesian estimate; the next later sample is the withheld target; the primary endpoint is the paired difference in 90 % interval score, with coverage and width reported separately. It is an *internal, retrospective, dosing-conditional* evaluation, not external validation and not a claim that either structure is clinically superior. A null result is a valid outcome.

**What this repository shows and does not show.** It shows verification of the software and numerics the comparison depends on. It does not contain the comparison.

## 2. Models and endpoint

Both structures use intravenous zero-order infusion input $R(t)$ (mg/h), clearance driven by renal function (`CRCA`, mL/min) and volume by adjusted body weight (`AJBW`, kg), with covariates piecewise-constant in time (addendum v1.1, §2).

**Model A, one compartment**

$$
\mathrm{CL}_i(t)=\theta_{CL}\frac{\mathrm{CRCA}_i(t)}{100}e^{\eta_{CL,i}},\qquad V_i(t)=\theta_V\frac{\mathrm{AJBW}_i(t)}{70}e^{\eta_{V,i}}
$$

$$
\frac{dA_c}{dt}=R(t)-\frac{\mathrm{CL}}{V}A_c,\qquad f=\frac{A_c}{V}
$$

**Model B, two compartments**

$$
\frac{dA_c}{dt}=R(t)-\frac{\mathrm{CL}}{V_1}A_c-Q\Big(\frac{A_c}{V_1}-\frac{A_p}{V_2}\Big),\qquad f=\frac{A_c}{V_1}
$$

$$
\frac{dA_p}{dt}=Q\Big(\frac{A_c}{V_1}-\frac{A_p}{V_2}\Big)
$$

with $V_1,V_2\propto \mathrm{AJBW}/70$ and between-patient variability on CL and the central volume only. Amounts are zero at each episode reset. The residual model is the variance-additive `combined2`, $v(f)=\sigma_{add}^2+\sigma_{prop}^2 f^2$, identical for both arms.

**Endpoint.** For observed target $y$ and central 90 % limits $L,U$,

$$
\mathrm{IS}_{0.10}(L,U;y)=(U-L)+20(L-y)\mathbb{1}(y \lt L)+20(y-U)\mathbb{1}(y \gt U),
$$

a proper scoring rule (Gneiting & Raftery, 2007). The primary contrast is $`\Delta_i=\mathrm{IS}_{B,i}-\mathrm{IS}_{A,i}`$ per patient.

**What the verification covers.** An amount-based ODE (explicit states `Ac`, `Ap`) is compared with an independent base-R solver, a per-interval exact solution by Padé scaling-and-squaring matrix exponential (`R/pkpd02_reference_matexp.R`), itself checked against closed-form 1C/2C solutions and `expm::expm`. The rule, fixed before the runs, is $`|C_{test}-C_{ref}|\le\max(10^{-5}\,\mathrm{mg/L},\,10^{-6}|C_{ref}|)`$ at every output row, with amounts checked the same way. Amounts, not concentrations, are carried across covariate changes. Two input patterns the unguarded ODE mishandles (a covariate change carried only by an observation row; a reset during or exactly at the end of an infusion) are canonicalised or rejected by a guard. **Identifiability and fit behaviour are not investigated**: no model has been fitted.

## 3. Verification status

Numbers are from the cited report and evidence file. "Synthetic" = generated by seeded code in this repository. Key: **FAIL** = gate failed and is recorded as failed; **Conditional pass** = passes only under stated conditions; **Pass for the classes tested** = pass within the stated synthetic scope only; **Not performed** = no result exists.

| Stage | What was tested | Data | Outcome | Source |
|---|---|---|---|---|
| **S-1** | `rxode2` `linCmt()` vs an independent amount-based ODE, 1C and 2C, time-varying covariates; 102 scenarios | Synthetic | **FAIL (fail-closed).** 5,063 of 6,354 points breach the gate; max abs error 4.0e-2 / 4.6e-2 mg/L against a 1e-5 allowance. The ODE agreed with the independent reference (0 failures, ≤ 1.4e-9 mg/L); the discrepancy agrees with `linCmt()` applying the new covariates over the interval *ending* at a record (a next-observation-carried-backward reading) rather than the intended last-observation-carried-forward, for the installed `rxode2` 5.1.8 (non-gating diagnostic in the report). Nothing was relaxed to obtain agreement. | `reports/S1_verification_report.md`, `evidence/s1_results.json`; preserved in `archive/s1_failed_2026-10-30/` |
| **S-1r** | Explicit amount-based ODE through three solve paths; guarded wrapper; 102 archived scenarios + 21 boundary cases | Synthetic | **Conditional pass** (`PASS_CONDITIONAL_ON_INPUT_CONTRACT`). 0 failing of 19,998 rows; max abs error ≤ 1.24e-9 mg/L. Conditional on an input contract (R-COV, R-RESET) still to be ratified; the unguarded path fails both hazards (errors up to 9.5e-2 and 14 mg/L). Library-default tolerances fail (3,021 / 3,039 points); a candidate `atol = rtol = 1e-10` passes. | `reports/S1r_verification_report.md`, `evidence/s1r_results.json` |
| **S-2 preparation** | Canonical event representation and episode timelines vs the reference: batteries A (150 patients), B (250), C (200, episode-relative clock, gap-invariance 0.001–100 h) | Synthetic | **Pass for the classes tested**: 0 failures; max abs error 4.2e-10, 1.6e-9, 4.9e-10. | `reports/S2prep_report.md`, `evidence/s2prep_canonical_equivalence.json` |
| S-2 prep, real-structure check | Recorded event structure of all 124 patients / 148 episodes (event and covariate metadata only) through the pipeline with synthetic parameters | Real metadata, no concentrations | 2,310 rows, 0 failures (max 3.9e-10 / 2.2e-8 / 7.6e-10). **Not reproducible from this repository**: needs the source workbook and a protocol preflight file not distributed here. | `evidence/s2prep_real_structure_domain_check.json` |
| S-2 prep, metadata audit | Counts only: 2,308 rows, 124 patients, 148 episodes, 24 resets; eligibility ceilings 120 episodes / 101 patients and 93 / 79; 2 observation-only covariate changes; 0 reset-tied samples | Real metadata, no concentrations | Descriptive upper bounds, not the analysed sample. **Not reproducible from this repository** (needs the source workbook). Patient-level tables, and two records that named a patient, episode and time, are deliberately not published; the JSON is a marked, redacted copy (`PUBLIC_RELEASE_NOTES.md` §3). | `evidence/s2prep_metadata_audit.json` |
| S-2 prep, static inspection | Compilation of explicit-ODE models; fitting-control interface; which acceptance criteria are checkable without a fit | None | Several items marked **unverifiable without a fit** (covariance coordinates, success-code field, bound tolerance, fallback behaviour). | `evidence/s2_static_schema_inspection.json` |
| S-2 synthetic fit tests; S-2a (≥ 1,000 random histories); sensitivity-state tolerances and run time | | | **Not performed.** Test design written, unapproved (`reports/S2_synthetic_test_design.md`). | |
| Estimation, posterior forecasting, scoring, leakage audit, coverage / interval-score results | | | **Not performed. No empirical result exists.** | |

## 4. Run it

Verified on one machine only: R 4.5.2, Ubuntu under WSL2. **Not verified:** a from-scratch `renv::restore()`, macOS, native Windows, other R or `rxode2` versions.

```sh
git clone https://github.com/pkaysantana/pkpd02-amikacin-solver-verification.git
cd pkpd02-amikacin-solver-verification
Rscript -e 'renv::restore()'          # needs internet and a C/C++ toolchain (rxode2)
Rscript tools/verify_public.R         # all four synthetic gates, about 2.5 minutes, no data needed (exit 0 = all behaved as recorded)
Rscript tools/verify_public.R s1r     # or a single gate: s1 | s1r | s2prep_canonical | s2_static
```

`tools/verify_public.R` runs each gate in its **own temporary copy** of the repository, so committed evidence and the preserved failure record are never modified, and gates cannot affect each other. It prints, per gate, whether the exit status and verdict text match the recorded outcome, which evidence files differ from the committed copies, and whether the repository was left unchanged. **`s1` is expected to exit 1**: its recorded result is FAIL.

Why not run the scripts directly? They write evidence in place. Running `tests/s1_synthetic_solver_equivalence.R` rewrites `evidence/s1_results.json`, after which `tests/s1r_run_all.R` correctly reports that the in-place S-1 files no longer match their archive and fails. The S-1 script is part of the hash-manifested failure record, so it was left untouched.

**Last run (3 Nov 2026, on this release's code, tests and evidence; the runner reads no documentation):** all four gates behaved as recorded (`s1` exit 1 with `S-1 VERDICT: FAIL`; `s1r`, `s2prep_canonical`, `s2_static` exit 0 with their recorded verdict text) and the repository was left unchanged. The run used a pre-restored `renv` library with the same `renv.lock` (SHA-256 `c2fc9a55…`), not a fresh `renv::restore()`. Every deterministic evidence file the gates rewrite matched the committed copy byte-for-byte, except fields recording a timestamp, a solve time or an environment listing (`s1_results.json`, `s1_sessionInfo.txt`, `s1r_results.json`, `s1r_equivalence_results.json`, `s2_pilot_sensitivity_accuracy_nonfit.csv`, `s2_static_schema_inspection.json`). Without a restored library the first command bootstraps only the `renv` package itself and needs internet access.

**File integrity.** `sha256sum -c RELEASE_MANIFEST.sha256` (Linux; `shasum -a 256 -c` on macOS) checks every released file against the hashes recorded at release. `.gitattributes` stops Git converting line endings, because `tests/s1r_run_all.R` checks archived files by hash and a CRLF conversion would break that check for a reason unrelated to the science.

**Cannot be run from a clone:** `tests/s2prep_metadata_audit.R`, `tests/s2prep_real_structure_domain_check.R` and `tests/s2prep_run_all.R` read the source workbook's event and covariate columns (the concentration column is skipped at read time) and a protocol preflight file that is not distributed. Their paths are taken from `PKPD02_SOURCE_XLSX` and `PKPD02_PREFLIGHT_JSON`. They are included for inspection; their recorded outputs are described in `PUBLIC_RELEASE_NOTES.md`.

## 5. Repository structure

| Path | Contents |
|---|---|
| `R/` | ODE models, guarded solver wrapper, independent matrix-exponential reference, event canonicalisation, metadata-only source reader |
| `tests/` | Gate scripts: `s1_*` (original), `s1r_*` (remediation), `s2prep_*` and `s2_static_schema_inspection.R` (S-2 preparation) |
| `tools/` | `verify_public.R`, the non-destructive runner |
| `protocol/` | Protocol v0.9, addendum v1.1, **draft, unfrozen** amendment v1.2 |
| `reports/` | One report per stage, with decisions still open |
| `evidence/` | Machine-readable results and point-level comparison tables written by the gates |
| `archive/` | Hash-manifested snapshots of the failed S-1 and the S-1r state |
| `PUBLIC_RELEASE_NOTES.md` | How this public copy differs from the original working repository, with before/after hashes; what was withheld or redacted; data-rights status |
| `RELEASE_MANIFEST.sha256`, `.gitattributes` | SHA-256 of every released file; byte-exact checkout settings |
| `LICENSE`, `THIRD_PARTY_NOTICES.md` | MIT Licence for the original work; notices for third-party files and dependencies |

## 6. Known limitations

* **No empirical result.** Nothing here says which structure forecasts better or whether either is adequate.
* **Verification is not validation.** Agreement with an independent numerical solution shows the code implements the stated equations and event rules. It says nothing about whether the model is a good description of amikacin pharmacokinetics.
* **Some documents the protocol cites are not distributed**: earlier internal notes (a dataset feasibility audit, a critical appraisal, a novelty check, a technical review) and two preflight JSON files. They were never part of this repository and cannot be inspected from it.
* **Real-structure results are historical records.** The three that used the source dataset's event metadata cannot be re-run from this repository (workbook not distributed, per-patient detail withheld). Everything else is reproducible without data.
* **S-1r is conditional.** The pass depends on the input contract and on the unguarded hazards staying behind the wrapper, for `rxode2` 5.1.8 as installed. No cross-version comparison was made.
* **Same-author reference.** The independent solver was written with the same event conventions as the code it checks, so agreement validates those conventions, not their clinical correctness. Not covered: `addl`, lag times, steady state, durations over 24 h, extreme parameter excursions during estimation.
* **Representation assumptions awaiting ratification.** `TIME` is episode-relative in the source; a single-ID patient table needs an explicit episode offset (prediction-invariant on synthetic tests, still an assumption). Whether any reset interrupts an infusion cannot be determined from recorded times; one patient is flagged as possible.
* **Fit and uncertainty untested.** Covariance extraction, the two-compartment acceptance criteria, sensitivity-state tolerances and ODE run time.
* **Study design limits.** One institution, retrospective, conditional on realised doses, sample times standing in for result-availability times, and the same dataset motivated the choice of structures.
* **Not clinical software.** No clinical utility, dosing safety or fitness for patient care is shown.

## 7. Roadmap (planned, not completed)

**Decisions needed from the author:** ratify the input contract, the episode-relative clock and the reset-tie rule; decide handling of an infusion interrupted by a reset; approve synthetic-data fitting (currently blocked); specify DV-validity and LLOQ handling; freeze the modelling addendum, version pins and a pre-fit manifest.

**Verification gates still to pass:** S-2 synthetic fit tests and S-2a on at least 1,000 random histories; production tolerances for sensitivity states and ODE run time; an exact-sensitivity reference; leakage audit (protocol Gate D) and uncertainty-engine validation (Gate C) before any score is computed.

**Implementation:** consolidate duplicated event preprocessing; add a guarded `foceiControl` constructor and fit-object extraction that fails a stage when a field is missing.

**Empirical stages (not started):** patient-held-out fitting, posterior forecasting, interval scoring and protocol reporting, only after the gates above.

## Provenance, references and responsible use

**Data source and attribution (CC BY 4.0).** Thomson A, Robb F, Siebinga H (2021). *Amikacin data for population pharmacokinetic modelling.* University of Strathclyde. https://doi.org/10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851. The data file is licensed **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/) and described by the publisher as anonymised patient data extracted from routine clinical files. **Changes:** no copy of the data is distributed. This repository contains no rows, concentration values or per-patient values from the dataset; it contains software, synthetic data, and aggregate counts and structure statistics computed from the event and covariate columns. The dataset's creators and the University of Strathclyde do not endorse this work.

**Additional rights wording, and what is unresolved.** The publisher's page footer reads (checked 3 Nov 2026): "Copyright © 2026 University of Strathclyde, its licensors, and contributors. All rights are reserved, including those for text and data mining, AI training, and similar technologies. For all open access content, the relevant licensing terms apply." This release treats that as a site-wide notice that defers to the CC BY 4.0 licence on the open-access data file. The publisher has not been asked to confirm this, so the point is unresolved; `PUBLIC_RELEASE_NOTES.md` §5 states exactly which files would be removed if the notice were read more strictly. If you reuse this repository, cite the dataset and follow CC BY 4.0 for any data you obtain from it.

**Literature** (bibliographic details checked against Crossref)

1. Siebinga H, Robb F, Thomson AH. Population pharmacokinetic evaluation and optimization of amikacin dosage regimens for the management of mycobacterial infections. *J Antimicrob Chemother* 2020;75(10):2933-2940. https://doi.org/10.1093/jac/dkaa277
2. Ideno Y, Kasai H, Tanigawara Y. SCOUT: An exploratory approach to scouting dose-relevant covariates. *CPT Pharmacometrics Syst Pharmacol* 2026;15(4):e70235. https://doi.org/10.1002/psp4.70235
3. Uster DW, Stocker SL, Carland JE, et al. A model averaging/selection approach improves the predictive performance of model-informed precision dosing: vancomycin as a case study. *Clin Pharmacol Ther* 109(1):175-183 (published online 2020). https://doi.org/10.1002/cpt.2065. Erratum: *Clin Pharmacol Ther* 110(5):1401 (2021). https://doi.org/10.1002/cpt.2353
4. Witta S, Wicha SG. Evaluation of a model averaging algorithm for model-informed precision dosing in the context of parameter misspecifications. *J Pharmacokinet Pharmacodyn* 2026;53(6):54. https://doi.org/10.1007/s10928-026-10064-5
5. Gneiting T, Raftery AE. Strictly proper scoring rules, prediction, and estimation. *J Am Stat Assoc* 2007;102(477):359-378. https://doi.org/10.1198/016214506000001437
6. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ* 2024;385:e078378. https://doi.org/10.1136/bmj-2023-078378

### Authorship, AI assistance and provenance

This is an independently directed research project by Don Aborah, developed with substantial assistance from AI coding assistants: Anthropic's Claude Code, and Devin (an AI software-engineering agent, used for audit and reconciliation work). Three things are kept distinct:

| | What this release says |
|---|---|
| **Scientific direction** | The author directed the project. This release does not itemise which specific designs, specifications, decisions or reviews were the author's own, and none should be inferred. |
| **Implementation** | Code, tests, reports and protocol drafts were produced with substantial AI assistance. The reports are written in an implementer's first person because they were produced that way; read them as assistant-drafted documents the author directed, not as hand-written ones. This release does not claim that the author manually wrote, verified or independently derived every implementation. |
| **Reproducible evidence** | The numerical results come from seeded, executable tests that write their outputs to `evidence/`; anyone can re-run them (§4). The independent reference solver (`R/pkpd02_reference_matexp.R`, 74 lines) can be read and checked without trusting either the author or the assistants. It shares event conventions with the code it checks (§6). |

**Provenance.** This repository is a history-free public copy of the author's working repository, whose history includes commits by the author and by an AI agent; none of that history, nor any branch or pull request, is carried over. `PUBLIC_RELEASE_NOTES.md` records every difference from the working repository.

### Responsible use, affiliation and licence

* **Non-clinical use.** A research and software-verification project. It must not be used to make or inform patient-care or dosing decisions, is not a validated medical device, and claims no clinical performance.
* **Affiliation.** Independent student project; not affiliated with, sponsored by or endorsed by any company or institution named in the references.
* **Licence.** The original code, tests, reports and protocol written for this repository are released under the MIT Licence (`LICENSE`). Third-party material keeps its own licence (`THIRD_PARTY_NOTICES.md`): the generated `renv/` bootstrap files are renv's (MIT, Posit Software, PBC), and R packages such as `rxode2` and `nlmixr2` (GPL-licensed) are not redistributed here. The source data are governed by their own licence (CC BY 4.0), above; no data are included.
