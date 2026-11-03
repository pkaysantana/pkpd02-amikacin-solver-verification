# PKPD-02 — Input-contract adjudication and S-2 preparation: report

**Date:** 30 October 2026 · **Scope:** no fitting (synthetic or real), no concentration (DV) read, no forecast targets, no pre-fit manifest, no freeze, no GitHub push. Re-run everything with `Rscript tests/s2prep_run_all.R`; machine-readable results in `evidence/s2prep_results.json`.

## 1. Summary

| Decision | Outcome |
|---|---|
| 1 R-COV fail-closed + canonical representation | Done. Canonicalisation re-carries an observation-only change on an `EVID=2` row at the same time with the same values; conflicts rejected; invariant to same-time row order; equivalent to the independent reference. |
| 2 R-RESET fail-closed + episode timelines | Done. Episode-specific timelines (no residual stop events) match the reference for resets inside, at the end of, and after infusions. |
| 3 Patient identity | Done. `PATIENT_ID`/`EPISODE` kept; η per patient; per-episode η refused without an explicit decision. |
| 4 Reset-tied sample quarantine | Done by file order + `TAD` only. **0 such samples exist in the data.** |
| 5 Metadata-only audit | Done; DV never loaded. Counts in §3. |
| 6 Amendment | Revision 2 written, still an unfrozen draft; revision 1 and all earlier artefacts preserved. |
| 7 Static S-2 inspection | Done; several coordinate and fallback questions cannot be settled without a fit (§5). |
| 8 S-2 test design | Written, unapproved, not run (`reports/S2_synthetic_test_design.md`). |
| 9 Domain verification | **Not "entire"**; see §6. |

All runner components pass; the S-1r regression passes with the same results as before.

## 2. A finding that changed the work: `TIME` is episode-relative

All 24 resets and every patient's first row have `TIME = 0`. Consequences, in order of importance:
- My **first-pass audit was wrong** (it grouped rows across episodes: spurious tied-covariate, coincident-event and "no active infusion at reset" counts). It was detected because 24 time decreases matched the 24 resets, and was **withdrawn and redone within episodes**. None of the first-pass numbers appears in any deliverable.
- "An infusion active at the reset" cannot be assessed from recorded times (unknown wall-clock gap). I report what is determinable instead (§3).
- A single-ID patient table (needed for shared η) requires an **explicit episode offset**; I tested that predictions are invariant to it, but it is a representation assumption for you to ratify.

## 3. Metadata audit results (`evidence/s2prep_metadata_audit.json`)

Reproduces the v1.1 preflight (2,308 rows, 124 patients, 148 episodes; eligibility 120/101 and 93/79; both sample/dose ties; all five example rows after the Excel-row offset).

| Check | Count |
|---|---|
| Observation-only covariate changes (R-COV) | **2** (2 patients, 2 training-relevant episodes); 22 other changes already on restart rows; 0 tied rows with conflicting covariates |
| Resets with all prior infusions demonstrably finished by the previous episode's last record | 23 of 24 |
| Resets where interruption is possible but not evidenced (previous episode ends with a dose record; patient 36) | **1** |
| Reset-tied samples (post-reset / indeterminate / quarantined) | 0 / 0 / **0** |
| Coincident events within an episode | 4 groups (2 sample+dose, 2 dose+dose; 2 exact-duplicate groups) |
| Real structure through the pipeline | 148 timelines and the 2,310-row patient-level table pass the validator with 0 messages |

**Exclusions and eligibility effect.** Samples excluded by the tie rule: **0**. Episodes excluded for R-RESET: **0**. Episodes excluded for R-COV: **0** with canonicalisation (2 if canonicalisation were refused: eligibility would fall from 120/101 to 118/100 one-prior and 93/79 to 91/78 two-prior). Primary forecasts touching a quarantined sample or a flagged reset: **0**. Gate A's 20% ambiguity threshold is not approached. *Limitation:* validity of `DV` (numeric, positive) was not assessed, so counts are upper bounds on that dimension.

## 4. Verification evidence (all against the independent reference; rule `max(1e-5 mg/L, 1e-6·|C|)`)

| Test | Rows per model | Failures | Max abs error |
|---|---|---|---|
| A: 150 patients, 390 observation-only covariate changes, short intervals | 2,539 | 0 | 4.2e-10 |
| B: 250 patients, 154 multi-episode, resets inside/at end of infusions (absolute clock) | 7,451 | 0 | 1.6e-9 |
| C: 200 patients, episode-relative clock; timelines + patient-level table at gaps 0.001/1/100 h | 5,748 / 5,927 ×3 | 0 | 4.9e-10; gap-invariant to 5.6e-10 |
| Real event structure, synthetic parameters: timelines (1e-12), production prediction model (1e-10), patient-level table | 2,310 / 2,310 / 2,334 | 0 | 3.9e-10 / 2.2e-8 / 7.6e-10 |
| S-1r regression (102 scenarios + boundary cases) | unchanged | 0 | unchanged |

Hazard power: unguarded single-table solves of battery B reach 244 mg/L error and −4,134 mg amounts; the replayed R-COV hazard reaches 9.5e-2 mg/L (1C) unguarded and 0 after canonicalisation. **Caveat:** the generic random observation-only battery did *not* trigger the unguarded R-COV failure by itself, because it is configuration-dependent; power for R-COV rests on the replayed hazard tables.

Defects I found in my own work and fixed (all before the final passing runs): a wrong tie-placement shift for pre-reset samples; a test-generator labelling error for pure-reset episodes; a fail-open gap (`ui$foceiModel` is NULL unless `nlmixr2est` is attached, and `rxSolve(NULL, …)` silently returned a control object), now an error in the guard; and a relative-error summary that blew up at near-zero concentrations (now reported only where |C_ref| ≥ 1e-3 mg/L). `R/pkpd02_solve.R` therefore differs from the S-1r-time version (hash recorded in the S-1r results; the archive holds the later file).

## 5. Static S-2 inspection (`evidence/s2_static_schema_inspection.json`)

Verified: explicit ODE compiles with amount and sensitivity states and no `linCmt()`; tolerances only via `rxControl=`; nlmixr2 defaults (1e-6/1e-3) give 3e-3 sensitivity errors versus ~2e-7 at ≥ 1e-8; `covFallback = list()` is accepted; bounds encode on the log scale. **Unverifiable without a fit (not approved):** covariance coordinates (Omega optimised on a `sqrt` coordinate), the C1 success-code field, `boundTol` for C8, whether fallback is really disabled. Fit-object accessors exist for `atol`, `rtol`, `covsInterpolation`, `control`, `seed`, `ui`; `fit$cov` and `fit$objf` are not accessor methods and need inspection of a real object.

## 6. Is the entire admissible input domain verified?

**No, and I do not claim it.** Verified: (i) the classes of event table exercised by S-1r and the three S-2-prep batteries; (ii) the **actual event structure of all 124 patients/148 episodes** with synthetic parameters, through both representations and the production prediction model at the proposed tolerance. Not verified: structures outside those classes (e.g. `addl`, lags, steady state, bolus handling beyond the T7 preprocessing, durations > 24 h, EVID=3 as it never occurs in the data); parameter excursions during estimation (bounds, extreme η); the FOCEI inner sensitivity model beyond a 32-patient pilot; any estimation behaviour; and the assumption that episode offsets and unknown wall-clock gaps are immaterial. The guard makes unverified structures fail rather than run, which is the intended safety property.

## 7. Remaining blockers

1. **Ratify** F-1/F-2 (episode-relative clock; offset representation) and the reset-tie rule.
2. **Approval to fit synthetic data** to read covariance coordinates, C1 field names, fallback behaviour: blocks C4–C7 and U2.
3. Production **sensitivity tolerances** must be set and tested (defaults fail).
4. **Exact-sensitivity reference** (augmented matrix exponential) for gradient tests is not yet written.
5. **Runtime** of ODE-based FOCEI is unknown (budget proposed in the design).
6. **DV validity/LLOQ** assessment needs an explicit permission and a one-way flag design.
7. Patient 36's possible interruption, dose+dose ties (T6) and duplicate rows (T15): decisions to keep as documented.
8. S-2a (v1.1 forecast-engine test, ≥ 1,000 histories) not yet run on the explicit-ODE engine.
9. rxode2 5.1.8 hazards persist; any package change requires re-running these tests and updating `renv.lock`.

## 8. Proposed code changes (none beyond the tested libraries are implemented)

Move the T7/T8 preprocessing now duplicated in two scripts into one library function; add a guarded `foceiControl` constructor (rxControl only, `covsInterpolation="locf"`, `covFallback=list()`); encode v1.1 bounds in `ini()`; implement exact-sensitivity reference and fit-object pathway assertions; add extraction functions for C1/C4–C9/U2 that fail the stage when a field is missing; keep the guard in front of every solve.

## 9. Files

`R/pkpd02_events.R`, `R/pkpd02_source_reader.R`, `R/pkpd02_solve.R` (guard + NULL check) · `tests/s2prep_canonical_equivalence.R`, `s2prep_metadata_audit.R`, `s2prep_real_structure_domain_check.R`, `s2_static_schema_inspection.R`, `s2prep_run_all.R` · `evidence/s2prep_*`, `s2_static_schema_inspection.json`, `s2_pilot_sensitivity_accuracy_nonfit.csv`, `s2_extraction_map.csv` · `protocol/PKPD02_v12_amendment_DRAFT.md` (rev 2) · `reports/S2_synthetic_test_design.md` · archives: `archive/s1_failed_2026-10-30/`, `archive/s1r_pass_conditional_2026-10-30/` (includes amendment rev 1).
