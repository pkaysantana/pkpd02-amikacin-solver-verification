# PKPD-02 — Technical Protocol Amendment v1.2 (DRAFT, revision 2): explicit amount-based ODE engine and canonical event representation

**Date:** 30 October 2026
**Parents:** `PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md` (SHA-256 `41153a5007547ae00938229d1b3acd9c140b6d82973ba6b62b74b41a68536920`), `PKPD02_v11_addendum.md` (SHA-256 `10afc03f2321990cd24541d792517ad83072d9b476bec76caba57c707fc720e8`)
**Status:** **DRAFT PROPOSAL. NOT frozen, NOT preregistered. No pre-fit manifest exists. No PK model fitting (synthetic or real), forecasting, forecast-target construction or concentration (DV) inspection has been performed.** Changes are adopted, if at all, before any outcome-based model evaluation.

## Revision history

| Rev | Content | Preserved at |
|---|---|---|
| 1 | Failure of `linCmt()` S-1; proposal C-1…C-5 (explicit ODE, pathway guard, input contract R-COV/R-RESET, tolerances, L9 wording) | `archive/s1r_pass_conditional_2026-10-30/protocol/PKPD02_v12_amendment_DRAFT.md` (manifest SHA-256 `a32a3d1d354f643874231f792f42a473807f6e7c8b61837c30c04387d6a7ae52`) |
| 2 | Reviewer adjudication of R-COV and R-RESET, canonical event representation, patient-identity rule, reset-tie quarantine, metadata-only audit results, S-2 preparation | this file |

The original failed S-1 artefacts remain in `archive/s1_failed_2026-10-30/` (manifest SHA-256 `67f940ad0af5595fc29f284cbcf7464253c648427005a459bf5110de54926cf8`).

## 1. Standing proposal (from revision 1, unchanged)

C-1: explicit amount-based ODEs replace `linCmt()` in population fitting and forecasting (amounts continuous across covariate changes). C-2: all solves pass through a guard that forces `useLinCmt=FALSE`, LOCF and explicit tolerances and asserts the executed pathway. C-4: production tolerances must be set explicitly (proposal: `atol = rtol = 1e-10`; nlmixr2 defaults 1e-6/1e-3 fail the S-1 rule). C-5: L9 is read as S-1r plus the C-2 assertions. All v1.1 model equations, covariates, bounds, `combined2` likelihood, C1–C9, folds, seeds and leakage gates are unchanged.

## 2. Reviewer decisions incorporated (30 October 2026)

1. **R-COV is a fail-closed guard.** Observation-only covariate changes are never silently accepted by the solver. A sanctioned *canonicalisation* re-carries such a change on an explicit `EVID=2` row at the **same time with the same values**; no value or effective time is invented; tied rows with conflicting covariates are rejected, not resolved.
2. **R-RESET is a fail-closed guard.** A genuine reset zeroes all amounts and cancels earlier infusions; this is implemented by **independently constructed episode-specific timelines**, which contain no residual infusion-stop events. It is not assumed that any reset interrupts a clinically ongoing infusion.
3. **Patient identity is preserved.** Episodes carry `PATIENT_ID` and `EPISODE`; `SOLVE_ID` is only a solver key. Random effects are **per patient** (v1.1 §2: no IOV). Independent per-episode η requires an explicit model decision and is refused otherwise.
4. **Reset-tied samples are quarantined when ambiguous**, by source event metadata only (file order and `TAD`), never by a measured value.

## 3. Findings that change the representation

**F-1: `TIME` is episode-relative in the source.** All 24 reset rows and the first row of all 124 patients have `TIME = 0`; there are 0 decreases within an episode and 24 across a patient's rows. A patient-level table with recorded `TIME` is therefore non-monotone and cannot be given to rxode2. Consequences:
- the wall-clock gap between episodes is not recorded, so "an infusion active at the reset" **cannot be assessed from recorded times**; an earlier first-pass audit that reported zero such cases was an artefact of comparing different clocks and is withdrawn;
- an earlier infusion can never overlap the next episode's clock in the timeline representation, so R-RESET cannot arise there; it can arise only in a patient-level table built with an absolute clock.

**F-2: single-ID patient tables need an explicit episode offset.** nlmixr2 attaches η to `ID`, so one patient-level ID is required to share η across episodes. Episode *k>1* is shifted by `offset_k = max(last record, last nominal infusion end of episode k−1) + gap` and its first row restored to `EVID=4/3`. The offset is a **representation device with no clinical meaning**: because states are zeroed at the reset, predictions must be invariant to `gap`. This invariance is tested (gaps 0.001, 1, 100 h; max difference ≤ 5.6e-10 mg/L), but it is a representation assumption for the reviewer to ratify.

## 4. Specification of the canonical event representation (draft)

1. Episode = rows between resets (`EVID` 3 or 4) in **file order**; `EVID=4` starts an episode and also delivers its dose.
2. Effective infusion duration per v1.1 T7/T8 (RATE=0 with AMT>0 → 0.5 h; negative rate or > 24 h → 0.5 h fallback; recorded `AMT/RATE` otherwise); every dose then has `RATE>0`.
3. Within an episode, same-time ordering: covariate (`EVID=2`) row, dose rows, then samples. A zero-order input produces no state jump at a dose start, so the T3 sample/dose order does not change predictions; covariate carriage at the same time is governed by the conflict rule.
4. **Reset-tied sample rule (episode-relative clock).** An observation in the same episode with `TIME` equal to the reset row's `TIME` follows the reset row in file order; it is `post_reset` iff `TAD == 0`, otherwise (including `TAD` missing) `indeterminate`. **Indeterminate observations are excluded from primary fitting and forecasting and reported.** (Absolute-clock variant, used only in synthetic stress tests: `pre_reset` iff file-order before the reset row and `TAD>0`; `post_reset` iff after and `TAD==0`; else indeterminate.)
5. Two routes: **episode timelines** (one solver timeline per episode; η looked up by `PATIENT_ID`) and the **patient-level offset table** (for FOCEI). Both pass the validator and both are checked against an independent reference.
6. Solver guard unchanged: R-COV and R-RESET inputs are rejected with an error; canonicalisation/episode splitting is the only sanctioned route.

## 5. Metadata-only audit of the Strathclyde event records (`evidence/s2prep_metadata_audit.json`)

Target-free: only `#ID, TIME, TAD, EVID, AMT, RATE, AJBW, CRCA` were loaded; `DV` and all other columns were skipped at read time (audit trail in the JSON; workbook SHA-256 verified `521daef9…767cb`). Reproduces the v1.1 preflight: 2,308 rows, 124 patients, 148 episodes, 684/1,600/24 events (EVID 0/1/4); eligibility 120 episodes / 101 patients (one-prior) and 93 / 79 (two-prior), all five preflight example rows (as Excel row numbers).

| Item | Count |
|---|---|
| Time decreases inside an episode | 0 (24 across rows of a patient = the resets) |
| Dose rows with RATE=0 (T7) | 6 (patient 46); 579 recorded durations outside 0.45–0.55 h |
| **R-COV:** covariate-change groups (within episode) | 24; 22 carried by a restart row; **2 observation-only** (2 patients, 2 episodes) |
| Tied rows with conflicting covariates | 0 |
| Coincident events (within episode) | 4 time-groups, 4 patients: 2 sample+dose (patients 50, 73; `TAD` 24, sample first), 2 dose+dose; 2 exact-duplicate-row groups |
| **R-RESET:** reset rows | 24, all EVID=4, none EVID=3 |
| — previous episode demonstrably finished all infusions by its last record | 23 |
| — indeterminate (previous episode's last record is a dose whose infusion extends beyond it; patient 36) | 1 (interruption possible but **not evidenced**) |
| **Reset-tied samples** | 0 (0 post-reset, 0 indeterminate, 0 quarantined) |
| Canonicalisation inserted rows (real structure) | 2 (the two obs-only changes) |
| Real structure through the sanctioned pipeline | 148 episode timelines and the patient-level offset table (2,310 rows) both pass the validator with 0 messages; the naive recorded-TIME patient table gives 42 messages (19 non-monotone time; the 6 R-COV and 17 tied-covariate messages are largely cross-episode `TIME=0` collisions, not real events) |

**Effect on eligibility (counts of episodes/patients):** as recorded 120/101 (one-prior), 93/79 (two-prior). After quarantine: unchanged (0 rows). If resets with an indeterminate interruption were excluded: unchanged (the one episode is not a primary forecast). If observation-only covariate changes were instead excluded (no canonicalisation): 118/100 and 91/78 (2 episodes). Primary forecasts touching a quarantined sample: 0; touching a flagged reset: 0; episodes with an obs-only covariate change: 2 (evaluation covariates are baseline-fixed, so this affects training episodes only). Gate A's 20% ambiguity rule is not approached. *Limitation:* observation validity (numeric positive `DV`) was deliberately not assessed, so these counts are upper bounds on that dimension.

## 6. Verification evidence (synthetic parameters; no concentrations)

- `tests/s2prep_canonical_equivalence.R` (seed 20261103): batteries against the independent reference, all rows within `max(1e-5 mg/L, 1e-6·|C|)`: A) 150 patients with 390 observation-only covariate changes (short intervals) → 2,539 rows per model, 0 failures, max error 4.2e-10; B) 250 patients, 154 multi-episode, 69 patients with a reset inside and 64 with a reset exactly at the end of an infusion (absolute clock) → 7,451 rows per model, 0 failures, max 1.6e-9; C) episode-relative clock (as in the data), 200 patients → episode timelines and patient-level offset table (three gaps), 0 failures, max 4.9e-10, invariant to the gap.
- Inputs the solver guard rejects: unguarded single-table solves of battery B give errors up to 244 mg/L and amounts down to −4,134 mg; the S-1r R-COV hazard replay shows unguarded errors up to 9.5e-2 mg/L (1C) and 1.1e-3 (2C), and 0 after canonicalisation. The *generic* random obs-only battery did **not** trigger the unguarded failure by itself (it is configuration-dependent); detection power for R-COV rests on the replayed hazard tables.
- Hand-built reset-tie cases (7 absolute-clock, 3 relative-clock) classify, place and value correctly; per-episode η are refused without an explicit decision.
- `tests/s2prep_real_structure_domain_check.R` (seed 20261104): the real event structure (124 patients, 148 episodes) with synthetic parameters/etas passes on episode timelines (ODE 1e-12), the production prediction model at 1e-10 (max error 2.2e-8 mg/L), and the patient-level offset table, for 1C and 2C. Aggregate output only; no concentrations.

## 7. S-2 preparation (`evidence/s2_static_schema_inspection.json`, `reports/S2_synthetic_test_design.md`)

Static inspection only. Verified: the compiled FOCEI models for the explicit ODE specification contain explicit amount states plus sensitivity states and no `linCmt()`; tolerances reach the solver only via `foceiControl(rxControl=…)`; `covFallback = list()` is accepted; v1.1 bounds encode as log-scale `ini()` bounds; nlmixr2 defaults (atol 1e-6, rtol 1e-3, sensitivities likewise) give 3e-3 sensitivity errors in a pilot solve versus ~2e-7 at ≥ 1e-8. **Unverified** (needs a synthetic fit, not yet approved): the coordinates of the returned covariance (Omega `sqrt` coordinate, bounded and log parameters), the C1 success-code field, `boundTol` mapping for C8, fallback behaviour. The test design (T-P, T-G, T-R, T-C, T-T, T-V, T-RT, T-H) is proposed with pre-declared thresholds, unapproved and not run.

## 8. Unresolved ambiguities and limitations

1. Wall-clock gaps between episodes are unknown; an interrupted infusion at a reset is possible for one episode (patient 36) but not evidenced. By construction it cannot affect any observation (states are zero after the reset), but the offset representation assumes it is immaterial.
2. Episode offsets (F-2) are a representation assumption requiring ratification.
3. Observation validity (DV numeric/positive) was not assessed; LLOQ handling (T13) is untouched.
4. Dose+dose ties (2 groups) and exact-duplicate rows (2 groups) remain governed by T6/T15 sensitivities; they were not altered.
5. Installed rxode2 5.1.8 retains the R-COV and R-RESET hazards; safety depends on the guard.
6. Per-episode η (IOV) is excluded by v1.1 §2; any change is an explicit model decision.
7. Whether the whole admissible input domain is verified: see `reports/S2prep_report.md` §6 (verified for the tested classes and the real event structure, not for unseen structures or for estimation behaviour).

## 9. Status statement

Unfrozen draft. Nothing here freezes the study, creates a pre-fit manifest, authorises fitting, reveals forecast targets, or has been pushed to any repository.
