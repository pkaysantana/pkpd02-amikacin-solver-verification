# PKPD-02 — Corrected Model-Specification and Leakage-Control Addendum (v1.1)

**Prepared:** 30 October 2026  
**Parent:** `PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md`  
**Supersedes:** proposed v1.0 addendum, as supplied in `Pasted text(20261030-133529).txt`  
**Status:** **Technical specification ready for implementation verification; NOT yet frozen, publicly preregistered, or authorised for PK model fitting.**

## Executive decision

The scientific comparison remains **one-compartment versus two-compartment population amikacin PK** for **retrospective, patient-held-out, future-assay posterior predictive intervals**, *conditional on the dosing subsequently recorded*. The principal paired endpoint remains the **90% equal-tailed interval score** (lower is better), with coverage and width reported separately. The inference is exploratory; the dataset does not establish a new dosing algorithm, prospective benefit, or external transportability.

Version 1.1 closes the **written-specification defects** identified in the independent v1.0 review, but a written decision is not the same as a passed software or leakage test. **Gate A metadata preflight and an independent synthetic mathematical cross-check have been performed; the `nlmixr2`/`rxode2` version-specific model tests and all patient-held-out fits remain unperformed.** No repository edits, Git actions, population-model fits or patient-concentration forecasts have been undertaken.

### Decision log: what changed from v1.0

| ID | v1.0 issue | v1.1 disposition | Consequence |
|---|---|---|---|
| V11-01 | `linCmt()` claimed to guarantee LOCF semantics | Require pinned-version event-table **equivalence test** against state-continuous ODE; read-back of setting alone is inadequate | Solver test **blocking before any real fit** |
| V11-02 | Covariance labelled exactly `r,s`; Omega assumed always absent | Accept only verified *sandwich R,S method* from actual fit metadata; inspect named covariance, scale and Omega handling | Covariance audit **blocking for U2 and estimation precision checks** |
| V11-03 | Two-compartment OFV constrained as if strictly nested | Reclassify OFV comparison as diagnostic **C3-D**, never a pass/fail gate | Eight active criteria: C1, C2, C4–C9 |
| V11-04 | `combined2` called an exact replication of source | **Pre-specified variance-additive modelling choice**, not demonstrated reproduction of source NONMEM `$ERROR` | Keep both arms on `combined2`; source-equivalence claim removed |
| V11-05 | “Exact posterior” and “unmask” terminology | Numerical posterior integration; **engineered information isolation**, not prospective blinded data | Claim and audit semantics corrected |
| V11-06 | One- versus two-prior forecasts potentially have different targets | One-prior uses **most recent** of two preceding samples; two-prior uses both; **same target** | Clean secondary paired estimand |
| V11-07 | Negative interval limits clipped without distributional clarification | **Raw equal-tailed 90% interval is primary**; clipping at zero is sensitivity only, with model content disclosed | Primary score corresponds to underlying likelihood |
| V11-08 | All dose records implicitly assumed 30-minute infusions | Preserve positive `AMT/RATE` values as recorded; separately flag off-30-min entries; six RATE=0 cases use predeclared 0.5 h approximation and bolus sensitivity | Dose-history provenance limitation, not silently fixed |
| V11-09 | Uniform serum-creatinine floor inferred | Use recorded CRCA as supplied; acknowledge historically **inconsistent** creatinine flooring | No invented uniform adjustment |
| V11-10 | R numerical implementation already considered “locked” | Freeze *only after* version pinning, synthetic solver and fit-object extraction tests; pre-fit manifest then Git tag if explicitly authorised | Draft can be evaluated before fitting |

## 1. Data, unit of analysis, and independent information

**Original workbook:** `Amikacin_Data_File.xlsx`  
**SHA-256:** `521daef9e7d753419006259ed15e3ac4f26a1345c290c79c1cacf42e11c767cb`  
**Fields:** `#ID, TIME, TAD, EVID, AMT, RATE, DV, SEX, AGE, WT, IBW, AJBW, HT, SCREAT, CRCL, CRCA`.

Source: Thomson, Robb and Siebinga, Strathclyde data deposit, DOI https://doi.org/10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851; Siebinga et al., *J Antimicrob Chemother* 2020, DOI https://doi.org/10.1093/jac/dkaa277.

There are **2,308 event rows**, **124 patients**, **684 observation records**, **1,600 ordinary dose records** and **24 reset-and-dose records**. The preliminary episode parser yields **148 treatment episodes**. A strict historical observation-pair rule yields **120 eligible episodes in 101 patients**, and a two-prior rule yields **93 eligible episodes in 79 patients**. Only the first eligible episode per patient contributes to a given primary analysis; repeated episodes from that patient stay in the same population-model training or evaluation fold. These are **preflight counts**, not fitted-model results, and may change only through documented Gate A corrections.

**Primary forecast:** at the first valid prior sample with at least one strictly earlier administered dose, predict the first later measured assay with at least one dose whose **start time is strictly between** origin and target. Forecasts condition on all **realised** intervening doses, not necessarily doses known at the origin. The observed target value is excluded from the individual's update and the population-training fold. No post-origin clinical covariate measurement is available to the evaluation patient's predictor.

**Secondary two-prior forecast:** select the first target for which two preceding usable observations exist and a dose starts strictly between the more recent prior observation and target. For the same patient and **identical target and dose path**, compare (i) updating with the **most recent** prior concentration only and (ii) updating with both prior concentrations. Changing this target rule after observing prediction performance is prohibited.

**Patient-level fold partition:** five folds stratified using *timing-derived eligibility class only* (0/1/2), with the saved random-number-generator specification in v1.0 and seed 20261030. All patient episodes belong to exactly one fold. No source-published population parameter fitted on all 124 patients may be used as an evaluation-fold prior or as a starting value.

## 2. PK model equations and observation model

The initial amount in every state is zero at each episode reset. Both structures receive the same IV input record, elimination from the central compartment, prescribed patient covariates, two Gaussian between-patient random effects, and identical residual likelihood.

For patient *i*, at time *t*, let `R_i(t)` denote the sum of infusion input rates from active recorded doses in mg/h. Recorded `CRCA` is in mL/min; `AJBW` is in kg. Within a training episode the covariates are piecewise constant from the record timestamp forward; in the evaluation episode they are fixed to their first-record baseline values for the full update and forecast.

**Model A — one compartment**

\[
\mathrm{CL}_i(t)=\theta_{\mathrm{CL}}\,\frac{\mathrm{CRCA}_i(t)}{100}e^{\eta_{\mathrm{CL},i}},\quad
V_i(t)=\theta_V\,\frac{\mathrm{AJBW}_i(t)}{70}e^{\eta_{V,i}}.
\]
\[
\frac{d A_{c,i}}{dt}=R_i(t)-\frac{\mathrm{CL}_i(t)}{V_i(t)}A_{c,i}(t),\quad
f_i(t)=\frac{A_{c,i}(t)}{V_i(t)}.
\]

**Model B — two compartments**

\[
\mathrm{CL}_i(t)=\theta_{\mathrm{CL}}\frac{\mathrm{CRCA}_i(t)}{100}e^{\eta_{\mathrm{CL},i}},\quad
V_{1,i}(t)=\theta_{V1}\frac{\mathrm{AJBW}_i(t)}{70}e^{\eta_{V1,i}},
\]
\[
Q=\theta_Q,\quad V_{2,i}(t)=\theta_{V2}\frac{\mathrm{AJBW}_i(t)}{70}.
\]
\[
\dot A_{c,i}=R_i(t)-\frac{\mathrm{CL}_i(t)}{V_{1,i}(t)}A_{c,i}
-Q\left(\frac{A_{c,i}}{V_{1,i}(t)}-\frac{A_{p,i}}{V_{2,i}(t)}\right),
\]
\[
\dot A_{p,i}=Q\left(\frac{A_{c,i}}{V_{1,i}(t)}-\frac{A_{p,i}}{V_{2,i}(t)}\right),\quad
f_i(t)=\frac{A_{c,i}(t)}{V_{1,i}(t)}.
\]

**IIV:** A: `(eta_CL,eta_V)` and B: `(eta_CL,eta_V1)` are independently normal with diagonal variance matrix `diag(omega_CL², omega_centralV²)`; *no* IIV on Q or V2, no Omega cross-correlation, and no interoccasion variability. These are **design constraints**, not assertions that excluded variability does not exist.

**Primary residual likelihood (`combined2`):**

\[
y_{ij}\mid \eta_i,\psi \sim N\left(f_{ij},\;v(f_{ij})\right),\quad
v(f)=\sigma_{\mathrm{add}}^2+\sigma_{\mathrm{prop}}^2 f^2.
\]

This is the variance-additive interpretation documented by `rxode2`, **not a verified reproduction of the original paper's undocumented NONMEM error expression**. The same formulation applies to A and B. `combined1` and proportional-only errors may be evaluated as *separately labelled pre-specified sensitivities* only and must not replace the primary likelihood based on scores.

**Parameter bounds retained from v1.0:** CL 0.5–50 L/h; central V 2–200 L; Q 0.01–100 L/h; peripheral V 0.5–500 L; additive SD 0.001–20 mg/L; proportional SD 0.001–2; positive IIV. They are estimation safeguards rather than established physiological truncation limits. Bound transforms and their exact scales must be inspected in the pinned package before fit; a run hitting a bound is not accepted.

**Covariate provenance:** use the supplied `CRCA` and `AJBW`; do not recompute or apply a new universal creatinine floor. The source records include creatinine measurements both set to 60 μmol/L and below 60 μmol/L. Where height is −99, retain the source's corresponding adjusted/actual-weight covariates, with a consistency flag. The provenance and actual historical measurement availability remain limitations.

## 3. Solver versions and pre-fit solver-equivalence test (blocking)

**Do not assert version-agnostic equivalence.** The official `rxode2` changelog documents a correction to the time-varying-covariate semantics of `linCmt()`; earlier releases could effectively apply the end-of-interval value despite `covsInterpolation='locf'`.

Before fitting, record **exact version, package source and resolved commit/CRAN build** for `R`, `nlmixr2`, `nlmixr2est`, `rxode2`, `lotri`, `dparser` and other relevant dependencies in `renv.lock`. Record the requested controls and their values from actual fit objects. Do not specify a package version that has not been installed and tested. 

**Required executable solver acceptance test S-1:** construct one synthetic patient event table with (a) two infusions, (b) a covariate step at a time not coinciding with either infusion, (c) another step during an active infusion, (d) an observation just before and just after each change, (e) a simultaneous sample/start-dose tie, (f) an episode reset, and (g) nonzero prior central/peripheral amounts when a step occurs. Solve **the identical event table** with `linCmt()` and an independent state-continuous ODE with explicit LOCF (`covsInterpolation='locf'`). Compare both 1C and 2C.

**Numerical threshold:** for every compared value, require `|C_linCmt - C_ODE| <= max(1e-5 mg/L, 1e-6*|C_ODE|)`, and check the states immediately before/after step changes. Report solver tolerance, maximum absolute error, model version and exact events. If it fails, stop before fitting; revise the engine or package selection by a dated amendment. Merely verifying the string `locf` is insufficient.

**Separate S-2:** confirm the forecast closed-form infusion-superposition implementation (with fixed evaluation covariates) against an independent ODE for at least 1,000 randomized histories including overlapping or interrupted infusions, resets and observation/dose ties. Use the same mixed absolute/relative threshold. The preflight Python mathematical check in §9 does **not** substitute for this production-path test.

## 4. Model fitting and two-compartment identifiability

**Estimator:** five patient-held-out training rotations, `nlmixr2` FOCE with interaction (`est='focei'`), same residual likelihood and covariates for A/B; requested `covMethod='r,s'`. Retain the previously specified R0 (9 starts), R1 (27 starts) and R2 (27 starts with `bobyqa`) numerical ladder. Report the actual optimizer success code and boundary warnings. A multi-start run must not be selected by a lower OFV if it failed convergence.

**Active acceptance criteria** — applied to training data **before scoring any held-out target**:

| Rule | Acceptance requirement |
|---|---|
| C1 | Optimiser success code from pinned package, no boundary warning; record exact code and extraction path. |
| C2 | Among *C1-converged runs only*, at least 3/9 in R0 or 5/27 in R1/R2 have OFV ≤ best converged OFV + 0.5 and agree on specified log parameters (≤0.10 for CL, central V, additive/proportional SD; ≤0.20 for Q, V2, and IIV variances). Preserve original v1.0 thresholds. |
| C3-D | **Diagnostic only:** record OFV(B)–OFV(A) and optimizer traces; **no automatic rejection**. With Q strictly bounded above zero, B is not guaranteed to achieve the A OFV by exact nesting. |
| C4 | Successfully computed, positive-definite **named** covariance with verified sandwich R/S method and **no undocumented fallback**. Accept labels such as `r,s (full)` only when method and coordinates are verified—not by string equality. |
| C5 | Condition number of the *structural log-parameter correlation block* <1,000. Confirm coordinate conversion first. |
| C6 | All absolute pairwise correlations of structural log-parameters ≤0.95. |
| C7 | Log-scale SE ≤0.30 for CL and central V; ≤0.50 for Q and V2; RSE ≤50% for both residual SDs. Confirm which parameters appear in the returned covariance before calculation. |
| C8 | Structural parameters farther than 0.05 on **estimation log scale** from hard bounds; additive SD ≥0.05 mg/L and proportional SD ≥0.02. |
| C9 | Both IIV variances ≥0.001. |

Run stages stop at the first satisfying all **eight active** conditions (C1,C2,C4–C9). The accepted fit is the best **converged** OFV candidate passing every active rule in the stage, and all attempts/OFVs/warnings are logged. If the named covariance or numerical bound data cannot be extracted reliably from the pinned fit object, that stage **fails**, rather than silently omitting a check.

**Fold consequences unchanged:** if either model fails in ≥2 of 5 folds, do not perform the full pre-specified paired analysis. With ≤1 failure in each arm, restrict to folds where both pass and report effective n; if fewer than 80% of Gate A eligible patients remain, downgrade to descriptive pilot. No new parameter fixing, covariate change, or outcome-informed bounds are permitted.

## 5. Posterior predictive uncertainty: U1 and U2

Let `psi_hat` denote a population-model fit estimated exclusively on other patient folds. Let `y0` be the permitted earlier measured concentrations, and `D` the allowed recorded dose path before the target. The conditional individual posterior is:

\[
\pi(\eta\mid y_0,D,\hat\psi) \propto
\phi_2(\eta;0,\hat\Omega)
\prod_{j\in\mathrm{allowed\ priors}}
\phi\{y_{0j};f(t_{0j},D;\eta),v(f(t_{0j},D;\eta))\}.
\]

The **primary U1 predictive distribution** is:

\[
F_{1}(y)=\int \Phi\!\left(\frac{y-f(t_1,D;\eta)}{\sqrt{v(f(t_1,D;\eta))}}\right)
\pi(\eta\mid y_0,D,\hat\psi)\,d\eta.
\]

This is **numerically integrated**, not exactly computed. Include both parameter uncertainty about a patient's random effects and the model's residual assay variability. Retain the proposed importance-sampling scheme (bivariate t with 4 df, MAP-centered Hessian proposal, N=100,000 doubling to 400,000, ESS≥20,000 and PSIS k<0.7) and an independently validated deterministic quadrature fallback; log proposal, ESS, k, tail stability, used N and failures. Cross-check raw 5%/95% quantiles under independent seeds to within 2% of interval width. No source-target DV may enter either distribution before scoring.

**Primary reported interval:** `[q0.05(F1), q0.95(F1)]`, **unclipped** and on the same continuous concentration scale as the specified `combined2` likelihood. Its primary 90% interval score is

\[
IS_{0.1}(L,U;y)=(U-L)+20(L-y){\mathbf{1}}(y<L)+20(y-U){\mathbf{1}}(y>U).
\]

Empirical 90% coverage, average raw width, raw median prediction error and weighted interval score (if justified independently) are descriptive secondary endpoints. A *nonnegative clipped* interval `[max(0,L),U]` is a separate sensitivity; report `F1(0)`, altered interval width, raw-versus-clipped scores and the actual implied distributional content. Do not refer to the clipped interval as exactly 90% model probability when it is not.

**U2 sensitivity:** examine approximate sampling uncertainty of **estimated fixed typical effects and residual SDs**, with `Omega` held at its accepted fold estimate *by design*. The returned covariance may contain a full block including Omega; **do not** assume it is absent.

1. Require a positive-definite, **named** covariance matrix, with documented estimated-parameter coordinates and estimator scale. Compare parameter names against the model's `theta_CL`, central-volume term, plus `theta_Q` and `theta_V2` for B, `sigma_add`, `sigma_prop`, and any other explicitly included fixed effects. Select a named block `Sigma_z` rather than relying on raw matrix position.
2. Define transformation from fitted coordinates `z` to the positive physical parameter **log vector** `xi=log(physical theta and sigma)`; derive the Jacobian `J_{ab}=d xi_a/d z_b` from the verified fit-object parameterisation. Then `Sigma_xi = J Sigma_z J^T`. Do **not** assume `J=diag(1,1/sigma,...)` unless that matches the inspected stored coordinates.
3. Sample *M=200* asymptotic normal draws on `xi` (as previously proposed), map back to physical parameters, and verify positivity, valid covariate scaling, no hard-bound violation and consistency of the reconstructed model. **Do not silently discard invalid draws**: if any invalid draw or covariance mapping occurs, U2 is marked **not computable** for that fold; log rate and reason. A prior validated constrained reparameterisation would require a pre-fit amendment.
4. Numerically compute conditional individual posteriors and average their predictive CDFs for the accepted draws; preserve the 20,000-proposal U2 target and the ESS/PSIS thresholds described in v1.0, separately labelled. Omega uncertainty is **excluded intentionally**. This is an **approximate hybrid sensitivity**, not a full Bayesian joint posterior over all population parameters.

If U2 cannot pass its parameter-coordinate or engine tests, **U1 remains primary**, and U2 is reported unavailable; it must not be repaired after inspection of held-out target scores.

## 6. Event parsing and dosing-time conventions (Gate A)

The workbook remains unchanged. A canonical metadata-only event table carries original file row, patient ID, episode ID, original `TIME`, episode-relative time, `EVID`, dose amount/rate, `TAD` and a *valid-observation flag* (numeric concentration >0), but the pre-specified forecast index has **no measured concentration values**. All perturbations of outcome values should preserve the validity flag.

**T1** `EVID=4` starts a new episode and also delivers its dose if `AMT>0`. Initial states zero, with patient episodes never mixed.  
**T2** chronological order is by (patient, episode, numeric `TIME`, type-dependent tie rank, original file row); a time decrease within source episode is a predeclared exclusion requiring logging.  
**T3** for simultaneous sample and dose, `TAD>0` implies sample predose; `TAD=0` implies sample after dose-start for ordering only. Dose input is zero-order, so an instantaneous dose start does not jump concentration. Log these events.  
**T4** dose at origin time, ranked after origin, may enter the **conditional realised dose path** but does **not** satisfy the strictly intervening-dose requirement.  
**T5** a dose co-timed with the target, ranked after sample, is not included in that prediction and does not qualify as intervening.  
**T6** keep coincident dose rows as recorded for primary analysis. Flag cases with duplicate-like, split-like or unclassified amounts. The predeclared sensitivity removes the second coincident dose when that patient's forecast uses the event; it is not assumed to be a true duplicate.  
**T7** for `RATE=0` with positive `AMT`, assign **0.5-h infusion as a transparent working approximation**, consistent with the published clinical description; explicitly reanalyse with bolus interpretation if such a patient supplies a primary forecast.  
**T8** for `RATE>0`, retain the *recorded* `AMT/RATE` duration, including the substantial number not equal to 0.5 h. Flag `|AMT/RATE−0.5 h|>0.05 h`; do not silently force 0.5 h. For invalid negative rates or duration >24 h, use 0.5 h only under a logged predeclared fallback. A **whole-cohort 0.5-h duration sensitivity** may be fitted and scored separately if explicitly configured before pre-fit freeze; it never replaces the as-coded primary results.  
**T9** if `TAD` disagrees by >0.25 h with the previous dose-start time difference, flag but retain `TIME`; shifting time by `TAD` is a forecast-specific sensitivity only.  
**T10** keep any sample during infusion and flag; verify against the source's claimed exclusion.  
**T11** a valid primary origin must strictly follow an episode dose start; do not count pretreatment baseline observations as informative Bayesian origin concentrations.  
**T12** if `HT=-99`, log whether supplied `AJBW` and `CRCA` approximately equal actual-weight fields; retain supplied covariates.  
**T13** retain positive observations below the unknown assay LLOQ as coded; describe this measurement limitation, but do not inspect their numerical values to select forecasters.  
**T14** invalid or nonnumeric observation concentrations do not enter training or prediction updates; validity flag is computed in a one-way preprocessing step.  
**T15** keep fully identical file rows for the primary analysis, log duplicate counts, and apply T6 sensitivity only for coincident dose rows in a used forecast.  
**T16** accept `TAD>24 h`; not automatically erroneous when dosing is intermittent/three times weekly.

**Preflight observations verified on the uploaded file:** 148 episodes; 0 decreasing-`TIME` episodes; 6 positive-dose `RATE=0` rows (patient 46); two exact sample/dose co-timed timestamps (patients 50 and 73); 34 observations with `TAD>24`; **579** positive-rate doses with duration outside 0.45–0.55 h; **4 duplicate additional event rows** across repeated full-record entries. No pharmacokinetic fitting or concentration-conditioned selection was used in these checks. Counts and cases are preserved in `PKPD02_v1_1_Metadata_Only_Preflight.json`.

**Critical provenance note:** the published source reports **30-minute administration**, yet the actual `AMT/RATE` coding frequently differs. Both facts must be stated; the cause is unknown. This contradiction is a limit on the meaning of the as-coded primary analysis, not a licence to rewrite historical rates without corroboration.

## 7. Leakage-control and release gates

The original published source and workbook have already been examined; the relevant goal is **computational isolation**, not genuinely unseen prospective outcome blinding. Preserve L1–L13 from v1.0, with the following enforceable interpretations:

- **L1–L2:** five patient folds, no patient/episode overlaps in held-out training.  
- **L3:** replacing every held-out patient's concentration with other **positive** synthetic values must leave the fit for that fold unchanged; re-fit is a later implementation gate.  
- **L4:** forecast input contains **only** permitted prior observed values, baseline covariates, doses before target and training-fold parameter hash, and no future target.  
- **L5:** replacing target and post-origin measured concentrations with positive canary values must leave the **timing-derived index and forecasts** identical (provided validity flags unchanged).  
- **L6:** perturbing post-baseline covariates cannot change held-out forecasts; this tests the baseline-frozen evaluation design.  
- **L7:** evaluation forecast reads no other treatment episode from the same patient.  
- **L8:** no published all-data Siebinga or SCOUT parameter estimate appears in starting values.  
- **L9:** requested `covsInterpolation` metadata **and** S-1 synthetic numerical behaviour both pass.  
- **L10:** synthetic and empirical input paths are separated.  
- **L11–L12:** scoring is the **only** join of forecasts to held-out targets, after manifest gates.  
- **L13:** each forecast's population-parameter hash corresponds to the rotation excluding that patient.

**Order of operations (no premature fitting):**

1. Finalise v1.1 and an analysis config; use Gate A timing-validity-only outputs and initial fold allocation.  
2. Install/record and pin R, nlmixr2, nlmixr2est, rxode2 and dependencies; run the **synthetic solver-equivalence test S-1** and **fit-object-schema inspection** using fully synthetic data. Confirm that C1–C9 and U2 can be evaluated from that version. The fitting of a *synthetic* model for compatibility is not a Strathclyde study fit, but is still a future implementation activity requiring its own approval.  
3. Construct **pre-fit manifest** with SHA-256 for original workbook, parent protocol, v1.1, exact analysis config, canonical event table, audit log, fold map, `renv.lock`, code commit and test log. Use a local immutable snapshot initially; publish a Git tag/remote link **only with explicit permission**.  
4. Only then permit five-rotation training fits, fit acceptance, posterior-engine synthetic validation and target-free forecast files.  
5. Construct a **final pre-scoring manifest** containing accepted-fit hashes, S-2/Gate C test reports, leakage tests, and `forecasts.csv`. Do not join targets before this.  
6. Join scored targets once, evaluate the locked 90% interval score and patient-level BCa bootstrap (10,000 resamples, seed 20261031), report coverage and sharpness, and retain failures/deviations.

**Gate thresholds retained:** model failure in ≥2 folds means no full paired head-to-head; losing >20% of eligible patient forecasts downgrades to descriptive pilot. All decisions and amendments have dated reasons, and explicitly state whether any held-out outcome evaluation preceded the change.

## 8. Scientific interpretation limits

- Population-model **form and covariates were historically informed by the full source**, so fold-wise parameter refitting does not erase the legacy advantage of one compartment. Label this **within-source exploratory cross-validation**, not independent external validation.
- Structural fit, numerical stability, future-assay error, interval coverage and interval width are **different outcomes**. A lower training OFV cannot decide a winner in the primary held-out interval score.
- Predictive calibration does not establish correct antimicrobial target attainment, nephrotoxicity protection, bacterial killing, AMR reduction, or clinical deployment.
- With approximately 101 independent patients, the design is primarily sensitive to *large* departures from nominal coverage. Report uncertainty; do not infer equivalence from nonsignificance.
- If the two-compartment fit fails the predeclared identifiability rules, a reportable negative result is preferable to post-hoc parameter fixing.

## 9. Completed independent prechecks and what they do **not** prove

| Check | Result | What it proves | What it does not prove |
|---|---|---|---|
| Original file checksum | `521daef9…767cb`, matches quoted hash | File identity | Accuracy of every historical clinical field |
| Metadata-only event audit | 2,308 rows, 148 episodes, 101/79 eligible patients | Reproducibility of event-level preliminary screening | Correctness of every dosing interpretation |
| Synthetic analytic vs independent numerical ODE, 1C | 50 synthetic histories, max absolute concentration discrepancy **5.64×10⁻¹¹ mg/L** | Mathematical 1C infusion formula cross-check | `nlmixr2` or production R solver behaviour |
| Synthetic analytic vs independent numerical ODE, 2C | Same 50 synthetic histories, max discrepancy **5.22×10⁻¹¹ mg/L** | Mathematical 2C infusion formula cross-check | Posterior quantile computation or real patient forecasts |
| R runtime/version-specific `linCmt()` test | **Not run** (no R runtime available in this environment) | None | Cannot close S-1 or numerical-fitting readiness |
| Named fit covariance/parameter-transform schema | **Not run** | None | Cannot lock C4–C7 or U2 execution compatibility |
| U1/U2 posterior and leakage tests | **Not run** | None | Cannot claim forecasting engine is validated |

Precheck output files: `PKPD02_v1_1_Metadata_Only_Preflight.json`, `PKPD02_v1_1_Synthetic_Math_Precheck.json`.

## 10. Gate decision — 30 October 2026

**Scientific addendum: conditionally specified.** v1.1 repairs the identified written mathematical and methodological defects. **Pre-fit freeze: NOT APPROVED YET.** Exact R environment and implementation-specific numerical tests (S-1, fit-object covariance extraction, S-2, Gate C, L1–L13) remain outstanding. Crucially, no patient-fold population parameter fit should be started before S-1 and the fit-schema test pass and the pre-fit manifest exists. Revisions made to repair implementation incompatibilities before any held-out scoring must be dated and transparent.

**Sources and relevant documentation**

- Siebinga H, Robb F, Thomson AH (2020). *Population pharmacokinetic evaluation and optimization of amikacin dosage regimens for the management of mycobacterial infections.* DOI https://doi.org/10.1093/jac/dkaa277.
- Thomson A, Robb F, Siebinga H (2021). *Amikacin data for population pharmacokinetic modelling.* https://doi.org/10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851.
- `rxode2` release notes, especially time-varying `linCmt()` correction: https://nlmixr2.github.io/rxode2/news/index.html.
- `rxode2` variance-additive residual documentation: https://nlmixr2.github.io/rxode2/reference/rxode2.html.
- `nlmixr2est` covariance controls and full named covariance: https://cran.r-universe.dev/nlmixr2est/doc/manual.html.
- Independent source audit: `PKPD02_Addendum_v1_0_Technical_Review.md`; parent protocol `PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md`.
