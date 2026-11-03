# PKPD-02 Gate S-1r — Explicit-ODE Pathway Verification (remediation of failed S-1)

**Date:** 30 October 2026 · **Data:** synthetic only · **Fitting / forecasting / Strathclyde data / targets:** none · **Pushes, freeze, pre-fit manifest:** none

## Verdict

**Gate S-1 passes for the corrected explicit amount-based ODE pathway, conditional on an input contract that you must ratify.** The unguarded rxode2 ODE path has two further failure modes (§4); the guarded wrapper rejects those inputs rather than solving them. This is not a clean unconditional pass, and I have not treated it as one. Machine-readable verdict: `evidence/s1r_results.json` → `overall_verdict: PASS_CONDITIONAL_ON_INPUT_CONTRACT`. Exit status of `Rscript tests/s1r_run_all.R` is 0 only when every component passes.

## 1. What was preserved (item 1)

The failed S-1 script, helper, report, comparison points, package versions, session info and `renv.lock` are copied to `archive/s1_failed_2026-10-30/` (20 files, `MANIFEST.sha256`, SHA-256 `67f940ad…26cf8`, read-only). The originals remain in place, read-only, and `s1r_run_all.R` re-verifies on every run that both the archive and the in-place originals still match the manifest (**verified: TRUE / TRUE**). The S-1r scripts read the archived event tables (SHA-verified), so the 102 scenarios are exactly those of the failed run.

## 2. Pathway inspection (item 5; `evidence/s1r_pathway_inspection.json`)

- `useLinCmt=FALSE` is supported: `rxSolve.rxUi` has a `useLinCmt` argument defaulting to `getOption("rxode2.useLinCmt", FALSE)`; an explicit `FALSE` overrides a global `TRUE`.
- **Silent conversion is real:** with `useLinCmt=TRUE` (or the global option set) the ODE model is converted (model flag `linCmt` = 1, state `central` instead of `Ac`) and reproduces the S-1 error (8.28e-3 mg/L, canonical case). Default and explicit `FALSE` stay explicit ODE (flag −100, state `Ac`, error 2e-10).
- `nlmixr2est` references `useLinCmt` nowhere. The FOCEI inner, pred-only and pred-no-lhs models compiled from the ODE specification (compiled only, **not fitted**) contain the explicit `Ac` (and `Ap`) states and no `linCmt`.
- The wrapper `pkpd02_solve_ode()` asserts after every solve: no `linCmt` flag, expected states, recorded `useLinCmt=FALSE`, LOCF, tolerances and method. It is immune to the global option, and rejects a `linCmt()` model, `nocb`, RATE=0 doses and malformed event tables (all tested, `evidence/s1r_guard_checks.csv`).
- `combined2` is unchanged: in the compiled production model `rx_r_ = add.sd² + prop.sd²·f²` to 4e-15.
- Not verified: an actual `nlmixr2` fit call (no fitting permitted). The later fit-schema stage must re-assert the pathway from a fit object.

## 3. Equivalence results (items 3, 4, 6, 7, 9)

Independent reference: `R/pkpd02_reference_matexp.R`, per-interval exact solution via Padé scaling-and-squaring (a different algorithm from the archived Taylor reference), handling infusion starts/stops, covariate changes (LOCF), EVID 3/4 resets, ties. Self-test (`evidence/s1r_reference_selftest.csv`, 6/6 pass): vs closed-form 1C (3e-14) and 2C (4e-14), vs `expm::expm` (7.6e-13 relative), mass balance with CL=0 (8.6e-12 mg), no-op row invariance, reset semantics.

ODE: three paths, all with `atol = rtol = 1e-12`, liblsoda, LOCF: **P1** plain rxode2 ODE; **P2** nlmixr2 UI model; **P3** the FOCEI `predOnly` model that nlmixr2 would use. Acceptance rule unchanged, applied against the reference: `|C_ODE − C_ref| ≤ max(1e-5 mg/L, 1e-6·|C_ref|)`; amounts (`Ac`, `Ap`) are compared with the same rule scaled by V₁.

| Suite (per model, per path) | Points | Failing | Max abs error (mg/L) | Max rel error | Max amount error `Ac` / `Ap` (mg) |
|---|---|---|---|---|---|
| 102 archived scenarios | 3,177 | **0** | ≤ 1.24e-9 | ≤ 1.3e-11 | ≤ 2.1e-8 / ≤ 4.5e-9 |
| 21 boundary cases (valid inputs) | 156 | **0** | ≤ 2.8e-10 | ≤ 6.9e-11 | ≤ 3.4e-9 / ≤ 6.2e-10 |

19,998 rows compared in total (3,333 × 2 models × 3 paths); **no failing point**, no NA. The three ODE paths agree with each other to ≤ 5.3e-10 mg/L. Covariate steps: 449 per model, 227 during an active infusion; the ODE central and peripheral amounts at each change equal the reference left-limit amounts (max 2.0e-8 / 3.6e-9 mg; 2,694 checks, 0 failures). Concentration continuity would have been wrong at those steps (median 27% difference), confirming amounts, not concentrations, are carried. Scenario coverage items (a)–(g) are asserted programmatically for both canonical patients, both models, including nonzero peripheral amounts at every step.

Boundary cases (valid): step at infusion start/end, sample at infusion start/end and at time 0 (both row orders), two infusions with same start/end, back-to-back infusions, 1e-6 h and 1e-9 h step spacing, covariate on dose rows, 48 h and 0.001 h infusions, sample/step/infusion-end coincident, EVID=3/4 resets after infusions end, short (0.1 h, 0.001 h) covariate blips carried on EVID=2 rows.

## 4. Unresolved solver hazards in the unguarded ODE path (`evidence/s1r_hazard_cases.csv`)

These were found by the boundary tests. I did not alter the reference or the rule; the wrapper now **fails closed** on both.

| Hazard | Unguarded ODE vs reference | Evidence |
|---|---|---|
| **R-COV:** covariate change carried only by an observation row (EVID=0) and reverting at the next row | Disagrees: up to 9.6e-4 mg/L (0.001 h) and 9.5e-2 mg/L (0.1 h, 1C). Not reproduced for 2C at 0.1 h; absent at ≥ 0.25 h and with an EVID=2 carrying row | Disappears with `hmax` below the interval or with `dop853`; consistent with the integrator stepping past the change, but I did not prove that mechanism |
| **R-RESET:** EVID 3/4 while an infusion is active, or exactly at its end | Phantom negative infusion: concentrations off by up to 14 mg/L, amounts down to −226 mg (1C), −188 mg (2C) | State traces show the cancelled infusion's stop event still applied after the reset; same with `dop853` |

**Non-gating tie-order observation** (`s1r_probe_sample_vs_reset_tie_nongating.csv`): a sample at the same TIME as an EVID=4 reset is evaluated on the post-reset state (0) whichever row comes first, whereas the pre-reset concentration was 21.5 mg/L. v1.1 T3 orders sample vs dose ties but does not cover resets, so this needs a stated convention before real data are used; it was not tested for the gate.

Because EVID=2 carriage of covariate changes was exact in every test, R-COV is cheap to satisfy. **R-RESET is a protocol question, not a coding one:** what should happen to an infusion interrupted by a reset (truncate, split, exclude)? I made no choice and inspected no patient data.

## 5. Numerical tolerances for later production fitting (item 8; `evidence/s1r_tolerance_sweep.csv`)

Swept on the production `predOnly` model over all 3,333 points per model; the candidate `atol = rtol = 1e-10`, liblsoda, was declared before the run, with a requirement of ≤ 10% of the allowance.

| Setting | Failing points (1C / 2C) | Max error as multiple of allowance |
|---|---|---|
| nlmixr2 `foceiControl()` defaults (1e-6 / 1e-3) | 3,021 / 3,039 | ≈ 6,000 / 11,000 |
| rxode2 defaults (1e-8 / 1e-6) | 1,567 / 837 | 11.2 / 6.4 |
| 1e-8 / 1e-8 | 0 / 0 | 0.135 / 0.086 (1C misses the 10% margin) |
| **1e-10 / 1e-10 (candidate)** | **0 / 0** | **1.1e-3 / 8.8e-4: adequate** |
| 1e-12 / 1e-12 | 0 / 0 | 1.2e-5 / 1.3e-5 |
| lsoda and dop853 at 1e-10 and 1e-12 | 0 | ≤ 1.1e-3 |

Library defaults would silently violate the acceptance rule, so production fits must set tolerances explicitly. **Not tested:** the FOCEI inner model's sensitivity states (`atolSens`/`rtolSens`) and ODE run time for the multi-start ladder.

## 6. Environment and reproducibility

`evidence/s1r_sessionInfo.txt`, `s1r_package_versions.csv`, `s1r_results.json` (versions, `renv.lock` entries, SHA-256 of every script, protocol, amendment and lockfile). R 4.5.2, Ubuntu 26.04.1 WSL2; rxode2 5.1.8; nlmixr2est 7.2.1; nlmixr2 7.0.1; rxode2lincmt 0.1.0; lotri 1.0.5; dparser 1.3.1.14. Seeds: failed-S-1 scenario generator 20261101 (archived), reference self-test 20261102; S-1r scripts generate no random numbers. Version-specific behaviour observed: `linCmt()` NOCB-like timing and ignoring `covsInterpolation` (S-1); automatic ODE→`linCmt()` conversion on `useLinCmt=TRUE`; R-COV and R-RESET above; ODE tolerance defaults inadequate.

## 7. Limitations

- The reference is independent of rxode2 and cross-checked against closed forms and `expm`, but it was written by the same implementer with the same event conventions; the 20,000-row agreement checks those conventions, not their clinical correctness.
- A passing S-1 does not validate fitting, covariance extraction, U1/U2, S-2 or leakage gates; none were run.
- The verdict is conditional on the input contract (R-COV, R-RESET) and on the unguarded hazards persisting in rxode2 5.1.8 only being reachable via the wrapper.
- Parameter values and etas are arbitrary synthetic values.

## 8. Decisions for you (nothing is frozen)

See `protocol/PKPD02_v12_amendment_DRAFT.md` §5: ratify or amend the input contract; decide reset-during-infusion semantics; confirm the ODE engine for S-2/forecasting; decide whether sensitivity-state tolerances and run time need a test before fit-schema inspection. I have not proceeded to FOCEI fitting.

## 9. Files

`R/pkpd02_ode_models.R`, `R/pkpd02_solve.R`, `R/pkpd02_reference_matexp.R` · `tests/s1r_reference_selftest.R`, `s1r_pathway_inspection.R`, `s1r_ode_equivalence.R`, `s1r_run_all.R` · `evidence/s1r_results.json`, `s1r_equivalence_results.json`, `s1r_pathway_inspection.json`, `s1r_comparison_points.csv` (every row, all paths), `s1r_failing_points.csv` (empty), `s1r_covariate_step_checks.csv`, `s1r_tolerance_sweep.csv`, `s1r_hazard_cases.csv`, `s1r_guard_checks.csv`, `s1r_ode_path_concordance.csv`, `s1r_summary_by_model_path_suite.csv`, `s1r_scenario_coverage.csv`, `s1r_boundary_events.csv`, `s1r_hazard_events.csv`, `s1r_probe_sample_vs_reset_tie_nongating.csv`, `s1r_reference_selftest.csv`, `s1r_sessionInfo.txt`, `s1r_package_versions.csv` · `protocol/PKPD02_v12_amendment_DRAFT.md` · `archive/s1_failed_2026-10-30/`.
