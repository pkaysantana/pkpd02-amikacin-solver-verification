# PKPD-02 Gate S-1 — Synthetic Solver-Equivalence Verification Report

**Verdict: FAIL (fail closed).** `linCmt()` and the independent amount-based ODE disagree beyond the pre-specified threshold for both the one- and two-compartment models. No scientific assumption, tolerance or model was altered to obtain agreement. Gate S-1 is **not passed**; per addendum v1.1 §3 work stops before any fitting pending reviewer decision and a dated amendment.

- Run: see `evidence/s1_results.json` (`timestamp_utc`). Re-run: `Rscript tests/s1_synthetic_solver_equivalence.R` from the project root (exit status 1 = FAIL, 0 = PASS). Runtime ≈ 8 s.
- Specification: protocol v0.9 and addendum v1.1 §3 (S-1), copied read-only into `protocol/` with SHA-256 recorded in the results JSON. Note: the documents were not in the project folder; they were taken from `<local PKPD folder>/`.
- Data: **synthetic only.** No Strathclyde data, no model fitting, no forecasting, no manifest, no freeze claim, no Git actions.

## 1. What was tested

| Component | Implementation |
|---|---|
| Analytic engine | `rxode2` `linCmt()`, CL and V(s) written from time-varying `CRCA`/`AJBW` per addendum §2 (1C: CL, V; 2C: CL, V1, Q, V2) |
| ODE reference | Separate rxode2 ODE with explicit amount states `Ac` (and `Ap`), `d/dt(Ac) = R − CL/V·Ac [− Q(Ac/V1 − Ap/V2)]`, `cp = Ac/V`; built from the equations, **not** from `linCmt()`. `covsInterpolation="locf"`, liblsoda, atol = rtol = 1e-12 |
| Adjudicating reference | Base-R piecewise-exact matrix-exponential solver (`tests/s1_helpers.R`), no rxode2 code; self-checked against the closed-form 1C solution. Used only to attribute a disagreement; also gating (ODE must agree with it) |
| Event table | One identical table per case is passed to both rxode2 models (`evidence/s1_events.csv`) |

**Scenarios (102 cases × 2 models; 3,177 comparison points per model).**
Two hand-built canonical patients (identical except for row order of the sample/dose tie) and 100 seeded random patients. Canonical: infusion 1 (0–2 h) → covariate step at 3.25 h outside any infusion → infusion 2 (7–11 h) with a coincident sample at t = 7 → covariate step at 9 h *during* infusion 2 → episode reset + dose (EVID=4) at 14 h → step at 16 h with residual amounts. Samples at t−0.001 h, t, t+0.001 h around every step. Random patients: 2–4 overlapping/irregular infusions, 2–7 steps (≥1 during an infusion), optional reset, one sample/dose tie, and η ≈ N(0, 0.3).
Items (a)–(g) of S-1 are asserted programmatically (`evidence/s1_scenario_coverage.csv`): all met for the canonical cases, including nonzero central amounts (min 317.7 mg) and nonzero peripheral amounts (2C) at every step. Overall, 422 steps occur during an active infusion.

**Tolerance:** `|C_linCmt − C_ODE| ≤ max(1e-5 mg/L, 1e-6·|C_ODE|)` at every output row (observations, dose rows, EVID=2 covariate rows).
**Continuity:** amounts, not concentrations, are carried across covariate changes (no concentration-continuity assumption is made anywhere). Step checks compare `C·V_new` at the change time with the pre-change central amount from the independent reference.

## 2. Results

| | 1C | 2C |
|---|---|---|
| Points compared | 3,177 | 3,177 |
| **Points failing gate (linCmt vs ODE)** | **2,530** | **2,533** |
| Max absolute error (mg/L) | 3.998e-2 | 4.568e-2 |
| Max relative error (vs ODE) | 6.37e-4 | 7.00e-4 |
| Points failing, ODE vs independent reference | 0 | 0 |
| Max \|ODE − reference\| (mg/L) | 1.4e-9 | 4.7e-10 |
| Canonical cases only: failing / total points | 42 / 64 | 42 / 64 |
| Canonical max abs error (mg/L) | 8.28e-3 | 1.30e-2 |
| Covariate-step amount checks failing | 845 of 852 (both models pooled) | |
| Max central-amount error at a step: linCmt / ODE (mg) | 0.529 / 2.2e-8 | |

**Failing cases: all 102 of 102 in both models** (`evidence/s1_case_summary.csv`; point-level list in `evidence/s1_failing_points.csv`; no NA comparisons). The 7 step checks that "pass" (all random cases, central-amount differences 1.5e-5–3.5e-4 mg, i.e. within the amount tolerance) still show non-zero deviation; their volume changes are not small (|ln V_new/V_old| 0.04–0.58), so they do not indicate agreement in general. Rows with `time` before a case's first covariate change: **0 failures out of 398 per model**; rows at/after it: 2,530 and 2,533 failures out of 2,779.

## 3. Diagnosis (non-gating; evidence in `s1_interpolation_semantics_diagnostic.csv`)

- The ODE with `locf` agrees with the independent LOCF reference (≤ 1.4e-9 mg/L). The ODE with `nocb` agrees with an independent NOCB reference (≤ 5.7e-10). The ODE is therefore a trustworthy reference.
- `linCmt()` agrees with the independent **NOCB-in-interval** reference (3.0e-11 for 1C, 1.5e-11 for 2C) and not with LOCF (up to 4.0e-2 / 4.6e-2). That is, over the interval ending at a record carrying new covariates, `linCmt()` applied the new CL/V, even though the previous record's covariates should apply under LOCF. This is the end-of-interval behaviour that addendum V11-01 anticipated.
- `linCmt()` output is **bit-identical** under `covsInterpolation="locf"` and `"nocb"` (max difference 0), and identical with `nonmem=TRUE`. So in this install it does not respond to the setting; reading back the string `locf` would have been misleading, as the addendum cautioned.
- Amount continuity is respected by the ODE (step-amount error ≤ 2.2e-8 mg) but not by `linCmt()` (up to 0.53 mg), consistent with the wrong-interval coefficients.
- The addendum states that the official rxode2 changelog documents a correction of this behaviour. I did **not** locate such an entry in the installed `NEWS.md` and did not verify why the installed 5.1.8 build still shows it; this is an observation for the installed stack, not a statement about other versions.

## 4. Version-specific solver behaviour observed

1. `linCmt()` ignores `covsInterpolation` for time-varying covariates in rxode2 5.1.8 / rxode2lincmt 0.1.0 (CRAN builds per `renv.lock`), as above. `covsInterpolation` code mapping: locf = 1, nocb = 2; default = 1.
2. ODE tolerance matters: at rxode2 default tolerances (atol 1e-8, rtol 1e-6, liblsoda) the ODE misses the gate against the exact reference (max 4.5e-4 mg/L 1C, 3.1e-4 2C; 1,520 and 816 points). At atol = rtol = 1e-12, liblsoda, lsoda and dop853 all agree with the exact reference within gate (≤ 1.4e-9). Any ODE-based production path would need tolerances verified against this.
3. In rxode2 5.1.8 an EVID=2 row is emitted in the output with `addDosing=TRUE`; a reset-dose (EVID=4) row appears as an extra EVID=3 row (cp = 0) and is included in the comparison.

## 5. Recorded environment and reproducibility

Details: `evidence/s1_sessionInfo.txt`, `evidence/s1_package_versions.csv`, `evidence/s1_results.json`, `renv.lock` (SHA-256 in JSON).
R 4.5.2 (x86_64, Ubuntu 26.04.1 under WSL2, reference BLAS/LAPACK 3.12); rxode2 5.1.8; nlmixr2 7.0.1; nlmixr2est 7.2.1; lotri 1.0.5; dparser 1.3.1.14; rxode2lincmt 0.1.0; symengine 0.2.14; renv 1.3.1. Single-threaded solves (`cores = 1`).
**Seed:** 20261101 (Mersenne-Twister / Inversion / Rejection). This was chosen by the implementer for S-1 scenario generation; the protocol assigns no S-1 seed, and it is distinct from the fold seed 20261030 and bootstrap seed 20261031.
**Tolerances:** 1e-5 mg/L absolute, 1e-6 relative to |C_ODE| (gate); ODE atol = rtol = 1e-12; ±0.001 h around each change.
**Synthetic parameters (arbitrary, not fitted, not derived from patient data):** 1C CL 4.0 L/h at CRCA 100, V 18 L at 70 kg; 2C CL 4.0, V1 14, Q 3.0, V2 9 (L/h, L, L/h, L). Individual η (CL, V/V1) in `evidence/s1_individual_etas.csv`; none on Q/V2.

## 6. Limitations and what this does not establish

- Peripheral amounts are not exposed by `linCmt()`, so peripheral-state continuity is verified only indirectly via subsequent central concentrations (and directly between the ODE and the independent reference).
- Not covered: an episode reset during an active infusion (the random generator places resets after infusions end; semantics would need a protocol decision); bolus (RATE=0) doses; steady-state dosing; `addl`; lag times.
- Parameters, covariate ranges and η variability are synthetic and arbitrary; a single installed stack was tested; no cross-version comparison was made.
- The fixed-baseline evaluation covariates (addendum §2) would avoid covariate steps in evaluation-patient forecasts, whereas training episodes use time-varying covariates, so the mismatch would bear on training-fold fits. I did not test this, and that reading of scope is mine, not a protocol statement.
- S-1 passing was the requirement; S-2, fit-object schema inspection, leakage tests and all fits remain unperformed and unreviewed.

## 7. Decisions needed from the reviewer (nothing has been changed)

Per §3, options for a dated amendment *before any fitting*; I have implemented none of them:
1. Select an engine that honours LOCF for time-varying covariates (e.g. the amount-based ODE with verified tolerances, or per-interval matrix-exponential propagation), and re-run S-1 against it.
2. Keep `linCmt()` only for settings with covariates constant within each episode, and document that restriction.
3. Investigate the installed rxode2/rxode2lincmt behaviour upstream or test another version, with exact versions pinned in `renv.lock`.

I do **not** recommend loosening the 1e-5 / 1e-6 threshold: the observed errors (up to 4.6e-2 mg/L, ≈ 4,600× the absolute tolerance) are systematic rather than numerical noise.

## 8. Files

| Path | Content |
|---|---|
| `tests/s1_synthetic_solver_equivalence.R` | Executable test (generates all evidence) |
| `tests/s1_helpers.R` | Event builder, independent base-R reference solver, tolerance rule |
| `evidence/s1_results.json` | Machine-readable verdict, tolerances, seed, versions, diagnostics |
| `evidence/s1_comparison_points.csv` | Every compared point: C_lin, C_ode, C_ref, absolute/relative error, tolerance, pass |
| `evidence/s1_failing_points.csv`, `s1_case_summary.csv` | Failures and per-case summaries |
| `evidence/s1_covariate_step_checks.csv` | Amount-continuity checks at every step |
| `evidence/s1_scenario_coverage.csv` | Assertions for S-1 items (a)–(g) |
| `evidence/s1_negative_controls_nocb.csv`, `s1_interpolation_semantics_diagnostic.csv`, `s1_ode_solver_sensitivity.csv` | Non-gating diagnostics |
| `evidence/s1_events.csv`, `s1_individual_etas.csv` | Exact synthetic events and η |
| `evidence/s1_sessionInfo.txt`, `s1_package_versions.csv` | Environment |
