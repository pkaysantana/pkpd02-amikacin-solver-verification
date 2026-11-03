# PKPD-02 — Technical Protocol Amendment v1.2 (DRAFT): explicit amount-based ODE engine

**Date:** 30 October 2026
**Parent documents:** `PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md` (SHA-256 `41153a5007547ae00938229d1b3acd9c140b6d82973ba6b62b74b41a68536920`), `PKPD02_v11_addendum.md` (SHA-256 `10afc03f2321990cd24541d792517ad83072d9b476bec76caba57c707fc720e8`)
**Status:** **DRAFT PROPOSAL. NOT frozen, NOT preregistered, NOT approved by a reviewer. No pre-fit manifest exists. No PK model fitting, forecasting, Strathclyde-data analysis or scoring has been performed.** Adopted, if at all, *before* any outcome-based model evaluation.

## 1. Why this amendment is proposed

Addendum v1.1 §3 made the `linCmt()` vs state-continuous ODE test (S-1) **blocking before any fit** and said that, on failure, the engine or package selection must be revised by a dated amendment. S-1 was run on synthetic data and **failed closed** (30 Oct 2026):

| Item | Observed (synthetic data only) |
|---|---|
| Installed stack | R 4.5.2; rxode2 5.1.8; nlmixr2est 7.2.1; nlmixr2 7.0.1; rxode2lincmt 0.1.0 (CRAN builds, `renv.lock`) |
| Result | `linCmt()` vs ODE: 2,530/3,177 (1C) and 2,533/3,177 (2C) points exceed `max(1e-5 mg/L, 1e-6·\|C\|)`; max abs error 0.040 / 0.046 mg/L; all 102 cases fail |
| Attribution | ODE (`locf`) agrees with an independent base-R matrix-exponential reference to ≤ 1.4e-9 mg/L; `linCmt()` agrees with an independent **NOCB-in-interval** reference to ~3e-11 and ignores `covsInterpolation` (bit-identical under `locf`/`nocb`, and under `nonmem=TRUE`) |
| Preserved evidence | `archive/s1_failed_2026-10-30/` (20 files, manifest SHA-256 `67f940ad0af5595fc29f284cbcf7464253c648427005a459bf5110de54926cf8`); originals left in place, read-only, unmodified |

Scope of the failure: time-varying CL/V inside an episode (training data). The evaluation design fixes covariates at baseline (v1.1 §2), so forecasting alone would not exercise this path, but the same engine must serve both, and an engine whose covariate timing cannot be trusted is not acceptable for either.

## 2. Proposed changes

**C-1 (engine).** Replace `linCmt()` by explicit amount-based ODEs in *both* population fitting and forecasting: states are drug **amounts** (1C: `Ac`; 2C: `Ac`, `Ap`); zero-order IV input from the recorded `AMT/RATE`; concentration `f = Ac/V` (central). Amounts stay continuous across covariate changes; concentration continuity is never assumed. Reference implementation: `R/pkpd02_ode_models.R`.

**C-2 (pathway guard).** Every solve is made through `pkpd02_solve_ode()` (`R/pkpd02_solve.R`), which passes `useLinCmt=FALSE` (overriding `options(rxode2.useLinCmt=TRUE)`), forces `covsInterpolation="locf"`, and **asserts after solving** that the executed model has no `linCmt()` (flag −100), has the expected amount states, and ran with the requested controls. Inspection of the installed versions: `useLinCmt` is an `rxSolve.rxUi` argument defaulting to `getOption("rxode2.useLinCmt", FALSE)`; with `useLinCmt=TRUE` (or the option set) an ODE model **is silently converted** to `linCmt()` and reproduces the S-1 error (8.28e-3 mg/L on the canonical case); `nlmixr2est` contains no reference to `useLinCmt`, and the FOCEI inner/pred models it compiles for the ODE specification contain the explicit `Ac` state. The fit-schema stage must re-assert this from an actual fit object (not done here: no fitting).

**C-3 (event-table input contract).** Verification showed two rxode2 ODE-path behaviours that violate the intended semantics and are therefore **rejected as invalid input** rather than tolerated:
- **R-COV:** a covariate change carried only by observation rows (EVID=0) followed by a short interval was applied with the end-of-interval value in the LSODA path (errors up to 8.5e-2 mg/L for a 0.1 h interval; not reproduced in every configuration, e.g. 2C at 0.1 h; absent with `dop853` or with `hmax` below the interval). Every covariate change must be carried by a restart row (EVID 1–4); in all verified cases EVID=2 carriage was exact.
- **R-RESET:** a reset (EVID 3/4) while an infusion is active, or exactly at its end, leaves the cancelled infusion's stop event active (phantom negative input; amounts down to −226 mg). The wrapper rejects it. **The protocol semantics for this case are not decided here** (see §5).

**C-4 (numerical tolerances).** Production ODE tolerances must be set explicitly. Proposed (pre-declared before the S-1r run, not frozen): `liblsoda`, `atol = rtol = 1e-10`. nlmixr2 `foceiControl()` defaults (`atol` 1e-6, `rtol` 1e-3) and rxode2 defaults (`atol` 1e-8, `rtol` 1e-6) **fail** the acceptance rule on the synthetic suite. The tolerance of the FOCEI inner (sensitivity) model was *not* tested.

**C-5 (L9 / S-1 wording).** v1.1 L9 ("requested `covsInterpolation` metadata **and** S-1 numerical behaviour both pass") is to be read as S-1r (explicit ODE vs independent reference) passing, plus the C-2 pathway assertions.

## 3. Explicitly unchanged

All v1.1 model equations, covariate definitions and scaling, IIV structure, `combined2` residual likelihood (verified numerically: `rx_r_ = add.sd² + prop.sd²·f²`, 4e-15), parameter bounds, C1–C9 criteria, fold design and seeds, U1/U2 definitions, the S-1 numerical rule `|ΔC| ≤ max(1e-5 mg/L, 1e-6·|C|)`, and all leakage gates. Typical values are written on the log scale in the nlmixr2 specification (mathematically identical; the estimation-scale and bounds-transform choice remains for the fit-schema stage).

## 4. Evidence for the amendment (S-1r, synthetic only)

See `reports/S1r_verification_report.md` and `evidence/s1r_*`. In brief: explicit ODE vs independent matrix-exponential reference passes at every compared row (19,998 rows over three ODE paths; scenario max abs error ≤ 1.3e-9 mg/L), amounts continuous at all 449 covariate steps per model (227 during active infusions; 2,694 checks over 3 ODE paths × 2 models), boundary regression passes for validated inputs, and the two hazard classes are rejected by the wrapper.

## 5. Open decisions for the reviewer (not made by the implementer)

1. Accept or amend the input contract (R-COV, R-RESET), including whether to preprocess real data to satisfy it.
2. **Semantics of a reset during an active infusion** (truncate, split, or exclude), interacting with T1, T7 and T8 in v1.1 §6. No Strathclyde data were inspected to see whether this occurs.
3. Whether the S-2 forecast engine is also the ODE (this amendment proposes yes) and how S-2 is re-worded.
4. Tolerance for the FOCEI sensitivity states, and ODE run-time acceptability for the 5-fold multi-start ladder.
5. Whether synthetic fit-object schema inspection may now begin (requires its own approval).

## 6. Status statement

This is an **unfrozen** draft. Nothing here freezes the study, creates a pre-fit manifest, authorises fitting, or has been pushed to any repository.
