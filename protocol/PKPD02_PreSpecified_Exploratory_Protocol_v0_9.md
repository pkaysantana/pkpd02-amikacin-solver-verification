# PKPD-02 | Pre-specified exploratory evaluation protocol

**Working title:** Structural-model sensitivity of Bayesian predictive-interval reliability for longitudinal amikacin therapeutic drug monitoring  
**Version:** 0.9 (protocol-design stage, 30 October 2026)  
**Status:** Draft **analysis plan**, not formally preregistered; not an ethics approval, clinical trial protocol or completed validation  
**Data:** Strathclyde public amikacin workbook, `Amikacin_Data_File.xlsx`; DOI https://doi.org/10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851  
**Provenance:** Siebinga, Robb & Thomson (2020), with Ideno et al. (SCOUT, 2026) as an alternative structural-analysis precedent  
**Research phase:** specification only; **no population models fitted or GitHub changes made**  

> **Scientific positioning.** This is an exploratory *internal* evaluation of an existing population and data source. It is not a claim that Bayesian forecasting, interval coverage, structural misspecification or model averaging are new. We will not claim a new drug-dosing algorithm or a clinical safety benefit. The main research value, if any, must come from a rigorously characterised sensitivity or failure mechanism under independent patient-level assessment.

## 1. Research problem and rationale

Bayesian population pharmacokinetic (PopPK) models are used to update individual drug exposure predictions after therapeutic drug monitoring (TDM). A point forecast can be reasonably accurate even when the associated predictive interval is too narrow, or unhelpfully wide. For a future observed concentration, a valid predictive distribution must incorporate the uncertainty remaining about the patient's PK parameters and the assay/observation process; for some models, interoccasion effects and other uncertainties also matter.

The original Strathclyde study retained a **one-compartment** amikacin structure after considering numerical identifiability and diagnostics, whereas the subsequent SCOUT exploratory reuse selected a **two-compartment** structure using objective function value, with a different analysis objective and without covariates in its initial amikacin fit. This disagreement motivates a controlled comparison, but does not show which structure forecasts future assays better. Prior literature already evaluates future-concentration forecasting and 90% interval coverage (e.g., Jung et al., 2026, tacrolimus); Uster et al. (2021), Hughes and Keizer (2021) and Witta and Wicha (2026) constrain methodological novelty.

**Target population:** patients represented in the Strathclyde amikacin monitoring records (treatment for mycobacterial infection). **Generalisability target:** similar patients and care pathways, not all amikacin indications, intensive care, paediatrics or unseen hospitals.

## 2. Study aims and hypotheses

**Primary objective.** Quantify the *paired* difference in accuracy–sharpness performance of **central 90% posterior predictive intervals for a later measured amikacin concentration** from a pre-specified one-compartment model and a structurally extended two-compartment model, using population parameters fitted exclusively on other patients.

**Primary estimand:** population mean of the patient-level difference in **90% interval score**, `IS_2comp - IS_1comp`, evaluated at the **first eligible future target in the first eligible episode of each held-out patient**; negative means the two-compartment model yields a better score. The comparison is about *forecast reliability*, not relative fit to the same concentrations.

**Secondary objectives:** describe empirical 90% coverage (and its uncertainty), median/mean interval width, point-error metrics, error direction, sensitivity to number of prior concentrations, forecast horizon and time since dose. Assess whether discrepancies indicate model-structure sensitivity or are primarily driven by uncertainty about parameters, renal function, sampling phase or a few influential observations.

**Hypotheses (falsifiable):**

- **H1 (structural sensitivity):** interval scores and/or empirical interval coverage materially differ between the candidate structures on withheld future assays. Null: no practically consequential difference, with uncertainty explicitly reported. This is a question, not an expected positive result.
- **H2 (failure mechanism, exploratory):** any degradation is associated with regions plausibly informative for distribution kinetics (especially early post-infusion samples) or with long between-measurement horizons. These are descriptive and underpowered unless the relevant subsets are large.
- **H3 (sampling history, exploratory):** using two prior TDM concentrations instead of one may alter sharpness/calibration, but this improvement is not assumed and the samples are not randomised.

No hypothesis claims novelty of predictive intervals or proposes a new method.

## 3. Source data and present feasibility audit

Public workbook: `AMIKACIN  DATA` worksheet, 16 columns; 2,308 event rows after two header rows; 124 patient identifiers; 684 concentration (`EVID=0`) records; 1,600 conventional dose (`EVID=1`) records and 24 reset-and-dose (`EVID=4`) records, i.e. 1,624 dosing events in total. Splitting a person's records at each reset yields **148 treatment episodes**. There are 120 episodes from **101 distinct patients** with an observed concentration followed by an administered dose strictly between that origin and a later measured concentration; for a two-prior concentration scenario, 93 episodes across **79 patients** pass a related rule. These are *mechanical eligibility ceilings*, not necessarily the final analysed sample.

Variables: `#ID`, `TIME`, `TAD`, `EVID`, `AMT`, `RATE`, `DV`, `SEX`, `AGE`, `WT`, `IBW`, `AJBW`, `HT`, `SCREAT`, `CRCL`, `CRCA`. The literal height sentinel `-99` denotes missing height; do not treat it as a physical measurement. The workbook supplies no MICs, treatment outcomes, toxicity events, actual laboratory result-availability timestamps or independent hospital cohort.

**Audit before estimation:** verify each episode's timestamp order; resolve simultaneous dose and sampling events using `TAD` and documented pre-dose semantics; verify `RATE` and `AMT` describe physically possible infusion durations; flag nonpositive or nonnumeric measured concentrations, patient covariate anomalies, potentially repeated rows and measurement/dosing ambiguities. Maintain an exclusion log at the *event, episode and person* levels. Do not exclude a patient based on whether their held-out concentration is hard to predict.

## 4. Design, information flow and estimand population

**Design:** retrospective observational, out-of-patient evaluation with internally cross-fitted population models. The primary analysis will be **one eligible forecast per patient**. The dataset was collected and modelled before this protocol: this is a pre-specified *reanalysis*, not a prospective study.

### 4.1 Episode and forecast construction — locked rules

1. Sort by `#ID` and the original event sequence, preserving the order of tied timestamps; define each `EVID=4` as a new episode boundary and simultaneously an administered dose if its amount is valid. **Never** carry a compartment state across episode resets.
2. Within each episode identify the earliest `EVID=0` measurement as the candidate origin, at timestamp `t0`. The concentration at `t0` is the one individual TDM sample available for the update. Its recorded time is treated as the idealised result-availability time because actual assay turnaround is absent.
3. Select the **earliest** subsequent `EVID=0` target at `t1` for which at least one administered dose occurs **strictly** within `t0 < t_dose < t1`. The rule deliberately does not count a dose timestamped at the target sampling time as preceding the measured target (e.g. the prior identified patient-73 same-time ambiguity).
4. Pick the **chronologically first eligible episode per patient** for primary scoring. Patient identifiers are unit-of-analysis and unit-of-fold assignment; no patient contributes twice to the primary endpoint.
5. All pre-origin dose events in the same episode are permitted as known treatment history. For prediction from `t0` to `t1`, include the *eventually realised* doses strictly after origin and before target (and any justified target-time dose ordering) in both model calculations. This is therefore a **conditional retrospective future-assay forecast given realised dosing**, **not** a dose-plan forecast a clinician necessarily could have issued at `t0`.
6. Do not use the target concentration, other later concentrations, future creatinine readings or later model-updating information in the posterior calculation for that target. At-origin covariates must be derived from the first recorded episode information or another measurement documented as no later than origin. If creatinine timestamps cannot be independently established, use a **baseline-only covariate snapshot** as the primary choice and explicitly label its availability as assumed, not proven.
7. Primary eligibility has **no post-hoc concentration cutoffs**, no arbitrary upper forecast-horizon exclusion and no requirement for a certain pharmacokinetic sampling phase. Record all exclusions; report descriptive horizons to reveal heterogeneity.

**Important temporal interpretation:** `t0` and `t1` are sample collection times, *not guaranteed decision or result-availability times*. The study's claimed validity is limited to retrospective, dosing-conditional forecasts.

### 4.2 Two-prior-concentration analysis — secondary

A separate, explicitly selected cohort may use the earliest **two** within-episode measured concentrations preceding a later target, requiring at least one dosing event strictly between the **second** origin and target. The target must not be used to update individual parameters. Report the changed eligibility denominator (at most 79 unique patients by the initial audit); **do not compare its scores naively with the entire one-prior cohort**, because selection and target times differ. For a fair paired within-cohort contrast, use the *same* patients and future targets, updating once with the first previous assay and once with both previous assays; document dose history in each case. This secondary evaluation remains observational.

## 5. Candidate models and fitting firewall

**Model A — structural reference:** intravenous one-compartment elimination/distribution model, with prespecified population covariates motivated by the original Siebinga paper: renal function for clearance and body size for volume. Use only variables actually present and available in the training records; no outcome-driven screening of covariates. Individual random effects and the residual observation model must be explicit.

**Model B — structural alternative:** intravenous two-compartment disposition (`CL`, `V1`, `Q`, `V2`) with the **same pre-specified clearance and weight covariate hypotheses** as Model A where biologically and statistically meaningful. This is a *controlled structural comparison*, **not** a reproduction of SCOUT's covariate-free exploratory fit. State the deviation from SCOUT explicitly. The same primary concentration error family is required for both models to avoid attributing a changed error model to the number of compartments.

**Fixed prior structural assumptions before target inspection:** intravenous dosing from `AMT`/`RATE`; positive typical parameters; positive between-patient variability components where identifiable; a proportional residual error family for the primary comparison, with a combined additive/proportional alternative reserved as a sensitivity analysis only if identifiable. Nonnegative, plausible covariate scaling must be specified mathematically before fitting. The published original equations may serve as a *specification reference*, but **published parameters trained on the entire 124-patient source must not become the evaluation-fold population prior**.

**Structural comparability guardrail:** both models get the same patient-level training folds, predictor definitions and residual error convention; estimation methods and convergence criteria are documented. Do not select between structures on evaluation outcomes. If a two-compartment fit is numerically unstable or its parameter covariance cannot support valid posterior uncertainty under predeclared conditions, label that fold/model **non-estimable**, rather than silently fixing parameters or reporting narrower intervals. A comprehensive fallback/abort rule must be frozen before numerical fitting (see section 11).

### 5.1 Population estimation and evaluation firewall

- Assign all **124 patients** to **five patient-level folds**, with a documented fixed seed, prior to any fitting. Every episode for a patient remains in their fold. In each rotation fit population parameters using the other four folds, retaining all qualifying training events. The evaluation fold contributes *only* its allowed pre-origin samples for its individual's Bayesian update and its withheld target for final scoring.
- Fit Models A and B **separately** in each training set. Fix their population estimates during each held-out person's posterior update; no use of the evaluation patient to estimate population fixed effects, random-effect variances, model-selection thresholds or residual error.
- Record all folds, optimiser versions, starting-value rules, parameter bounds, estimation failures, covariance failures and changes in effective evaluable n **before** examining predictions by arm. Training parameters need not recover published full-source values exactly. No tuning on held-out target results.
- Because the source literature originally used this same dataset to motivate our model comparison, cross-fitting removes **direct population-parameter overlap**, but does not make the entire research question historically independent. Call the resulting assessment **internal cross-fitted evaluation**, never external validation.

**Uncertainty distribution:** sample/predict future *assay observations* using the training-fold population uncertainty assumptions, the held-out patient's posterior parameter uncertainty given the allowed prior concentration(s), and residual assay error. Report explicitly whether population-estimation parameter uncertainty and time-varying random effects were integrated or treated as fixed; if excluded, the resulting distribution is **conditional on fitted population parameters** and may understate total uncertainty. A plug-in MAP interval obtained by varying residual noise alone does **not** constitute a full Bayesian posterior predictive distribution. Numerical sampling error and convergence are independently checked on training/synthetic cases before outcomes are unlocked.

## 6. Endpoints and statistical analysis plan

### 6.1 Primary endpoint — proper interval score

For observed target concentration `y` and pre-specified central 90% predictive limits `L` and `U` (the 0.05 and 0.95 predictive quantiles), define

`IS_0.10(L,U;y) = (U-L) + 20*(L-y)*I(y<L) + 20*(y-U)*I(y>U)`.

Smaller is better: this **proper interval scoring rule** penalises wide intervals and observations outside them, avoiding the trivial solution of predicting infinitely wide intervals. The patient-specific primary contrast is `Δ_i = IS_B,i - IS_A,i` in mg/L. Report the paired mean `mean(Δ_i)` with a **95% patient-level resampling interval** and the mean score for each arm. Also present the median and distribution to show influential observations. For clinical comparability, present the relative mean-score difference as a secondary contextual metric, not a new primary outcome.

**Inference:** two-sided exploratory testing of no mean difference only if both models are estimable under the locked plan. Use patient-resampled bootstrap (predeclare ≥5,000 resamples and a fixed seed at implementation) for uncertainty, with the caveat that fixed cross-fitted population training overlaps mean such intervals are **conditional/approximate** and do not fully propagate the uncertainty from population re-estimation. If feasible, nested refitting bootstrap is a sensitivity check, not a prerequisite for the descriptive primary result. Do not claim superiority if the interval includes zero, and do not claim equivalence merely because it does.

**Materiality lens, not a formal clinical threshold:** a ≥10% relative reduction in mean interval score is designated as a *working practical-interest flag*, not a clinically validated minimal important difference and not a substitute for its uncertainty interval. Do not choose a more favourable threshold after results.

### 6.2 Secondary reporting — calibration, sharpness and point accuracy

- **Empirical central 90% coverage** separately for A/B: `C = n^{-1}Σ I(L_i≤y_i≤U_i)`; report counts, percent and **95% Wilson binomial interval** for one-forecast-per-patient analysis. Also report paired model coverage discordance and uncertainty for the **paired difference** in coverage; do not infer paired-model power from the single-model binomial calculation.
- **Sharpness:** mean and median `U-L` in mg/L, alongside interval score; broad intervals must not count as success simply because they cover observations.
- **Point accuracy:** median prediction error (signed), median absolute error, mean absolute error and root mean squared error for posterior predictive median; display observed-vs-predicted and residual-by-horizon plots.
- **Calibration across nominal levels:** 50%, 80%, 90% and 95% predictive coverage curves *descriptively*; small n means these plots must not be overinterpreted.
- **Sensitivity analyses:** target within versus beyond 168 hours from origin; early post-dose samples `<2 h` versus later; severe outliers/influential observations; baseline versus at-origin last-known renal covariates **only if provenance can be established**; one versus two prior concentrations in a paired 79-patient-maximum subset; descriptive impact of absolute concentration level and renal function. Do not fit numerous subgroup-specific models.
- **Missing forecasts:** publish by-model/by-fold failures and patient characteristics. Primary complete-pair inference requires predictions from **both** arms; report exclusions and ascertain whether missingness differs between arms. No convenient post-hoc replacement of failed intervals.

### 6.3 Precision and interpretation

At the data-audit ceiling, **101 independently represented patients** yield approximate 95% Wilson bounds of **82.7–94.5%** when 91/101 future concentrations fall inside a nominal 90% interval. Prior exact binomial calculations suggested ~82% power to detect 80% rather than 90% coverage but ~34% for 85% rather than 90%, *assuming independent patient forecasts and no further exclusions*. Those are design calculations **not observed coverage**. The true paired model comparison's precision depends on paired differences and discordant coverage indicators, which are unknown before forecasts. This is **exploratory**, not powered equivalence/non-inferiority or small-effect superiority.

## 7. Leakage prevention and audit trail

**Hard prohibitions:** no evaluation patient's later target assay in fold-level parameter estimation, posterior updating, covariate derivation or method tuning; no mixing of episodes across resets; no reusing future creatinine as a known-at-origin covariate; no model ranking by held-out observed concentrations before protocols are frozen; no adjustment of forecast horizon or target based on goodness of fit; no interpretation of retrospectively realised dosing as prospectively known dosing. No mixing synthetic observations with the real clinical outcome denominator.

**Reproducible outputs at implementation stage:** immutable raw data checksum; canonical event-cleaning log; patient/episode forecast index; frozen training-fold assignment; held-out target mask; model specifications; full numerical diagnostics; per-person predictions from each arm; pre-specified scoring script; exclusion table; versioned results and limitations. The coding and repository actions remain **outside this protocol-authoring stage**.

## 8. Bias, ethics and reporting limitations

This is a secondary analysis of a public anonymised research dataset. Check attribution and licence conditions and any institution-specific governance requirements before dissemination. No clinical decisions will be generated or implemented.

Major foreseeable risks: non-random monitoring/follow-up; predominantly early-phase target concentrations; uncertainty about assay timing and LLOQ; time-varying kidney function with incomplete timestamp provenance; structural nonidentifiability, especially for the two-compartment model under sparse observations; differences in illness and care pathways from other amikacin settings; lack of an independent centre; training-fold correlation; legacy source-structure selection from this same dataset; use of actual future doses in the conditional forecast; no clinical endpoint or MIC information. Follow the transparency principles of TRIPOD+AI (2024) where applicable to prediction-model evaluation; it is **reporting guidance**, not a proof of methodological validity.

## 9. Additional controlled simulation — separate optional arm

If the empirical results suggest a specific failure mechanism, a **separately declared simulation extension** can generate data under a known one- or two-compartment truth and alter one factor at a time (structural form, renal-function misspecification, residual error magnitude, sampling schedule, prior dispersion). Evaluate the same interval score and coverage, then assess whether the real-data pattern is mechanistically consistent. Simulations must be designed without calibrating scenarios to make Model B win and should not be used to manufacture evidence of real-world clinical benefit. Because this extension is not currently fully specified, it is **not part of the locked primary study** and must receive its own versioned supplement before implementation.

## 10. Pre-specified result interpretation

| Scenario | Allowed conclusion | Prohibited conclusion |
|---|---|---|
| Similar interval scores and coverage; CIs broad | No persuasive difference detected; study imprecise | Structural equivalence proved |
| B has a lower interval score with interpretable uncertainty | B offers better scores *in this internal, conditional evaluation* | Two compartments clinically superior everywhere |
| B lower score but substantial undercoverage | Accuracy/sharpness trade-off; possible overconfidence | B unequivocally trustworthy |
| A superior despite B fitting better | Within-source structural fit does not guarantee future-assay performance | Two-compartment PK is false generally |
| 2-comp fails identifiability in multiple folds | Model complexity is not robustly estimable in the planned design | Post hoc parameter fixes justify a clean head-to-head win |
| Neither arm calibrates reliably | Reliability problem exposed; method improvement remains future research | Novel calibration algorithm demonstrated |
| Simulation reproduces a failure mechanism | Mechanistic plausibility in controlled conditions | Prospective clinical proof |

## 11. Protocol readiness gates and amendments

**Gate A — data audit:** finalise event-time ordering, dose input semantics, covariate provenance and episode-reset parsing. If >20% of the nominal 101 patients are excluded for ambiguity, downgrade to a **descriptive pilot** and transparently reassess precision *without rewriting the original primary design to chase significance*.

**Gate B — independent model fit:** define exact parameterisation, random effects, residual error and convergence tolerances in a frozen modelling addendum, **before target concentrations are unmasked for scoring**. A 2-compartment model requiring extensive post-hoc constraints in more than one fold will be designated **non-estimable for the pre-specified comparison**. The one-compartment descriptive baseline may still be reported. Thresholds and allowable remediation must be fully specified in that addendum, not improvised after scores are seen.

**Gate C — uncertainty engine:** posterior predictive draws must include observation error; validate against analytically checkable simple cases and ensure stable estimated quantiles. An engine delivering only MAP point estimates or residual-only intervals cannot answer the primary question.

**Gate D — leakage audit:** independently verify patient-fold exclusions and no target or later information in updates. If audit fails, invalidate/re-run the affected analysis under unchanged rules and document the correction.

**Gate E — final interpretation:** no claims of external transportability, improved AMR outcomes, dose safety or a new algorithm; no statistically definitive small-effect claim with the available sample. Report all deviations, regardless of outcome.

**Amendment discipline:** every change after this draft obtains a dated version, specific reason and explicit statement whether it was adopted **before or after** outcome-based model evaluation. The user has not formally registered the plan; use “pre-specified analysis protocol” rather than asserting public preregistration.

## 12. Decisions locked now vs still needed

| Status | Item | Resolution |
|---|---|---|
| **Specified** | Source population and file | Public Strathclyde mycobacterial amikacin records |
| **Specified** | Primary unit of analysis | One first-eligible forecast per patient |
| **Specified** | Target selection | First later assay after at least one strictly intervening dose |
| **Specified** | Forecast interpretation | Retrospective, **conditional on realised intervening doses** |
| **Specified** | Patient independence | 5-fold patient-held-out population estimation |
| **Specified** | Structural comparison | 1-compartment vs 2-compartment, with matched primary covariate/error assumptions |
| **Specified** | Primary endpoint | Paired mean difference in 90% interval score |
| **Specified** | Secondary endpoints | Empirical coverage/Wilson CI; width, errors, stratified descriptions |
| **Open — required before fitting** | Exact covariate equations and base weights | Model specification addendum referencing original parameters |
| **Open — required before fitting** | Two-compartment stability/fallback | Freeze numeric identification thresholds and valid fixed effects |
| **Open — required before fitting** | Bayesian uncertainty computation | Full posterior approximation/sampling and simulation diagnostics |
| **Open — required before fitting** | Timing conventions and ambiguous events | Resolve with data dictionary/source record; document exclusion rules |
| **Open — recommended** | Formal public registration | Only after all analytical choices are frozen and before outcome scoring |

## 13. References (anchor sources)

1. Siebinga H, Robb F, Thomson AH. *Population pharmacokinetic evaluation and optimization of amikacin dosage regimens for the management of mycobacterial infections.* J Antimicrob Chemother. 2020. https://doi.org/10.1093/jac/dkaa277. Dataset: https://doi.org/10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851.
2. Ideno S, Kasai H, Tanigawara Y. *SCOUT: An Exploratory Approach to Scouting Dose-Relevant Covariates.* CPT Pharmacometrics Syst Pharmacol. 2026. https://doi.org/10.1002/psp4.70235.
3. Uster DW et al. *A Model Averaging/Selection Approach Improves Predictive Performance of Model-Informed Precision Dosing.* Clin Pharmacol Ther. https://doi.org/10.1002/cpt.2065. Published mathematical correction: https://doi.org/10.1002/cpt.2353.
4. Hughes JH, Keizer RJ. *A hybrid machine-learning pharmacokinetic approach* and related prior-adaptation work, 2021. See `PKPD02_Core_Paper_Critical_Appraisal_v1.md` for exact provenance and limitations.
5. Witta M, Wicha SG. Study of model averaging under misspecified population PK parameters, 2026. https://doi.org/10.1007/s10928-026-10064-5.
6. Jung et al. Paediatric tacrolimus future-occasion forecasting and empirical 90% interval coverage, 2026. https://doi.org/10.3389/fphar.2026.1896778.
7. Gneiting T, Raftery AE. *Strictly proper scoring rules, prediction, and estimation.* J Am Stat Assoc. 2007;102:359–378. Interval score implementation reference: https://search.r-project.org/CRAN/refmans/scoringutils/html/interval_score.html.
8. Collins GS et al. *TRIPOD+AI statement: updated guidance for reporting clinical prediction models.* BMJ. 2024;385:e078378. https://doi.org/10.1136/bmj-2023-078378.
9. Supporting internal documents: `PKPD02_Strathclyde_Dataset_Feasibility_Audit_v1.md`, `PKPD02_Core_Paper_Critical_Appraisal_v1.md`, `PKPD02_Targeted_Novelty_Kill_Test_and_Precision_v1.md`.

---

**Protocol design decision, 30 October 2026:** Proceed to a *model-specification and leakage-check addendum*, not directly to implementation. Primary test is exploratory and conditional; detailed numerical choices must be locked before any future concentration is used to assess candidate forecasts. A negative or inconclusive result remains a valid outcome.
