#!/usr/bin/env Rscript
# PKPD-02 Gate S-1: synthetic solver-equivalence test
#   rxode2 linCmt()  vs  independent amount-based rxode2 ODE (explicit LOCF)
#   (+ third base-R piecewise-exact reference for adjudication)
# Spec: protocol v0.9 + addendum v1.1 s3.  SYNTHETIC DATA ONLY.
# Usage (from project root): Rscript tests/s1_synthetic_solver_equivalence.R
# Exit status 0 = all gating checks pass; 1 = FAIL (fail closed).

args <- commandArgs(FALSE)
f <- sub("^--file=", "", args[grepl("^--file=", args)])
ROOT <- normalizePath(file.path(dirname(f), ".."))
setwd(ROOT)
suppressMessages({ library(rxode2); library(jsonlite) })
source("tests/s1_helpers.R")

S1_SEED   <- 20261101L  # implementer-chosen for S-1 scenario generation (protocol assigns none)
N_RANDOM  <- 100L
ODE_ATOL  <- 1e-12
ODE_RTOL  <- 1e-12
EPS       <- 1e-3       # "just before / just after" a covariate change (h)
set.seed(S1_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

# ---- models --------------------------------------------------------------------
mods <- list(
  lin = list(
    "1c" = rxode2({
      CL <- tCL * CRCA / 100 * exp(eCL)
      V  <- tV * AJBW / 70 * exp(eV)
      cp <- linCmt()
    }),
    "2c" = rxode2({
      CL <- tCL * CRCA / 100 * exp(eCL)
      V1 <- tV1 * AJBW / 70 * exp(eV)
      Q  <- tQ
      V2 <- tV2 * AJBW / 70
      cp <- linCmt()
    })),
  ode = list(
    "1c" = rxode2({
      CL <- tCL * CRCA / 100 * exp(eCL)
      V  <- tV * AJBW / 70 * exp(eV)
      d/dt(Ac) <- -CL / V * Ac
      cp <- Ac / V
    }),
    "2c" = rxode2({
      CL <- tCL * CRCA / 100 * exp(eCL)
      V1 <- tV1 * AJBW / 70 * exp(eV)
      Q  <- tQ
      V2 <- tV2 * AJBW / 70
      d/dt(Ac) <- -CL / V1 * Ac - Q * (Ac / V1 - Ap / V2)
      d/dt(Ap) <-  Q * (Ac / V1 - Ap / V2)
      cp <- Ac / V1
    }))
)

# ---- synthetic scenarios ---------------------------------------------------------
D <- function(...) data.frame(...)
canon_rows <- function(tie_sample_first) {
  tie <- if (tie_sample_first)
    D(time = c(7, 7), evid = c(0L, 1L), amt = c(0, 700), rate = c(0, 175))
  else
    D(time = c(7, 7), evid = c(1L, 0L), amt = c(700, 0), rate = c(175, 0))
  rbind(
    D(time = 0, evid = 1L, amt = 600, rate = 300),                       # infusion 1: 0-2 h
    D(time = c(0.5, 1.5, 2.5, 3.25 - EPS, 3.25, 3.25 + EPS, 4, 6), evid = 0L, amt = 0, rate = 0),
    tie,                                                                  # infusion 2: 7-11 h (+ tie)
    D(time = c(8, 9 - EPS, 9, 9 + EPS, 10.5, 11.5, 13, 14 - EPS), evid = 0L, amt = 0, rate = 0),
    D(time = 14, evid = 4L, amt = 400, rate = 400),                       # episode reset + dose 14-15 h
    D(time = c(14 + EPS, 14.5, 15.5, 16 - EPS, 16, 16 + EPS, 18, 24), evid = 0L, amt = 0, rate = 0)
  )
}
canon_sched <- D(time = c(0, 3.25, 9, 14, 16),
                 CRCA = c(90, 45, 120, 70, 100),
                 AJBW = c(75, 95, 60, 82, 70))

cases <- list()
add_case <- function(label, ev, eCL, eV, canonical = FALSE) {
  cases[[length(cases) + 1L]] <<- list(label = label, ev = ev, eCL = eCL, eV = eV,
                                       canonical = canonical)
}
add_case("canonical_tie_sample_first", s1_events(1L, canon_rows(TRUE),  canon_sched), 0.10, -0.05, TRUE)
add_case("canonical_tie_dose_first",   s1_events(2L, canon_rows(FALSE), canon_sched), 0.10, -0.05, TRUE)

gen_episode <- function(t0, reset) {
  n_inf <- sample(2:4, 1)
  starts <- t0 + c(0, cumsum(runif(n_inf - 1, 2, 10)))
  dur <- runif(n_inf, 0.25, 4); amt <- round(runif(n_inf, 200, 1200))
  ends <- starts + dur
  hor <- max(ends) + runif(1, 2, 8)
  # steps: >=1 during an infusion, >=1 outside any infusion
  k <- sample(n_inf, 1)
  st <- starts[k] + runif(1, 0.1, 0.9) * dur[k]
  repeat {
    cand <- runif(sample(1:3, 1), t0 + 0.5, hor - 0.5)
    inside <- vapply(cand, function(x) any(x > starts & x < ends), NA)
    if (any(!inside)) break
  }
  steps <- sort(c(st, cand))
  if (min(diff(c(t0, steps, starts))) < 0) NULL
  list(starts = starts, ends = ends, amt = amt, dur = dur, hor = hor, steps = steps)
}
for (i in seq_len(N_RANDOM)) {
  reset <- runif(1) < 0.4
  ep1 <- gen_episode(0, FALSE)
  rows <- D(time = ep1$starts, evid = 1L, amt = ep1$amt, rate = ep1$amt / ep1$dur)
  sched <- D(time = c(0, ep1$steps), CRCA = runif(1 + length(ep1$steps), 20, 150),
             AJBW = runif(1 + length(ep1$steps), 40, 120))
  ot <- c(runif(sample(5:8, 1), 0.01, ep1$hor), ep1$steps - EPS, ep1$steps, ep1$steps + EPS)
  tie_t <- ep1$starts[sample(length(ep1$starts), 1)]
  if (tie_t == 0) tie_t <- ep1$starts[length(ep1$starts)]
  tie_first <- runif(1) < 0.5
  ob <- D(time = ot, evid = 0L, amt = 0, rate = 0)
  tie_o <- D(time = tie_t, evid = 0L, amt = 0, rate = 0)
  if (tie_first) rows <- rbind(tie_o, rows) else rows <- rbind(rows, tie_o)
  rows <- rbind(rows, ob)
  tr <- ep1$hor
  if (reset) {
    t_r <- ep1$hor + runif(1, 0.5, 2)
    ep2 <- gen_episode(t_r, TRUE)
    # reset dose is the first dose of episode 2
    rows <- rbind(rows, D(time = ep2$starts[1], evid = 4L, amt = ep2$amt[1], rate = ep2$amt[1] / ep2$dur[1]))
    if (length(ep2$starts) > 1)
      rows <- rbind(rows, D(time = ep2$starts[-1], evid = 1L, amt = ep2$amt[-1], rate = ep2$amt[-1] / ep2$dur[-1]))
    sched <- rbind(sched, D(time = c(t_r, ep2$steps), CRCA = runif(1 + length(ep2$steps), 20, 150),
                            AJBW = runif(1 + length(ep2$steps), 40, 120)))
    rows <- rbind(rows, D(time = c(runif(sample(5:8, 1), t_r + 0.01, ep2$hor), ep2$steps - EPS,
                                   ep2$steps, ep2$steps + EPS), evid = 0L, amt = 0, rate = 0))
  }
  ev <- s1_events(100L + i, rows, sched)
  # unique times except the single intended sample/dose tie
  stopifnot(!any(duplicated(ev$TIME[ev$EVID == 0])))
  eta <- pmax(pmin(rnorm(2, 0, 0.3), 0.6), -0.6)
  add_case(sprintf("random_%03d", i), ev, eta[1], eta[2])
}

events <- do.call(rbind, lapply(cases, `[[`, "ev"))
icov <- data.frame(ID = vapply(cases, function(x) x$ev$ID[1], 1L),
                   eCL = vapply(cases, `[[`, 0, "eCL"), eV = vapply(cases, `[[`, 0, "eV"))
label_of <- setNames(vapply(cases, `[[`, "", "label"), icov$ID)

# ---- solving ------------------------------------------------------------------------
solve_one <- function(mod, model, interp = "locf", atol = ODE_ATOL, rtol = ODE_RTOL,
                      method = "liblsoda") {
  x <- as.data.frame(rxSolve(mods[[mod]][[model]], S1_PARS[[model]], events, iCov = icov,
                             covsInterpolation = interp, atol = atol, rtol = rtol,
                             method = method, addDosing = TRUE, cores = 1))
  names(x) <- sub("^id$", "ID", names(x))
  x
}
maxerr_vs <- function(a, b) max(abs(a - b))

res <- list(); steps_tbl <- list(); neg <- list(); sens <- list(); semantics <- list()
for (model in c("1c", "2c")) {
  lin <- solve_one("lin", model)
  ode <- solve_one("ode", model)
  # row alignment must be exact between solvers
  stopifnot(nrow(lin) == nrow(ode), all(as.character(lin$ID) == as.character(ode$ID)),
            all(lin$time == ode$time), all(lin$evid == ode$evid))
  # independent reference, looked up by (ID, time)
  refs <- lapply(cases, function(cs) {
    r <- ref_solve(cs$ev, model, S1_PARS[[model]], cs$eCL, cs$eV)
    r$ID <- cs$ev$ID[1]; r
  })
  ref <- do.call(rbind, refs)
  key <- function(id, t) paste(id, sprintf("%.9f", t))
  m <- match(key(lin$ID, lin$time), key(ref$ID, ref$time))
  stopifnot(!anyNA(m))
  pts <- data.frame(model = model, case = label_of[as.character(lin$ID)], ID = as.integer(as.character(lin$ID)),
                    time = lin$time, evid = lin$evid,
                    C_lin = lin$cp, C_ode = ode$cp, C_ref = ref$cp[m],
                    CRCA = ode$CRCA, AJBW = ode$AJBW, row.names = NULL)
  pts$abs_err   <- abs(pts$C_lin - pts$C_ode)
  pts$rel_err   <- ifelse(pts$C_ode == 0, NA_real_, pts$abs_err / abs(pts$C_ode))
  pts$tol       <- s1_tol(pts$C_ode)
  pts$pass      <- !is.na(pts$abs_err) & pts$abs_err <= pts$tol
  pts$ode_vs_ref_abs <- abs(pts$C_ode - pts$C_ref)
  pts$ref_pass  <- !is.na(pts$ode_vs_ref_abs) & pts$ode_vs_ref_abs <= s1_tol(pts$C_ref)
  pts$lin_vs_ref_abs <- abs(pts$C_lin - pts$C_ref)
  res[[model]] <- pts

  # ---- step checks: amount continuity across covariate changes ----
  for (cs in cases) {
    ev <- cs$ev; id <- ev$ID[1]
    chg <- which(c(FALSE, diff(ev$CRCA) != 0 | diff(ev$AJBW) != 0))
    for (j in chg) {
      tc <- ev$TIME[j]; if (ev$EVID[j] == 4) next  # reset: states intentionally zeroed
      V_ratio <- NA
      P <- pts[pts$ID == id & abs(pts$time - tc) < 1e-12, ][1, ]
      rr <- ref[ref$ID == id & abs(ref$time - tc) < 1e-12, ]
      Vold <- ref$V1[ref$ID == id & ref$time == max(ref$time[ref$ID == id & ref$time < tc])]
      Pb <- pts[pts$ID == id & abs(pts$time - (tc - EPS)) < 1e-9, ][1, ]
      steps_tbl[[length(steps_tbl) + 1L]] <- data.frame(
        model = model, case = cs$label, ID = id, t_change = tc,
        during_active_infusion = any(ev$EVID %in% c(1, 4) & ev$TIME < tc &
                                     ev$TIME + ev$AMT / pmax(ev$RATE, 1e-300) > tc & ev$AMT > 0),
        V_old = Vold, V_new = rr$V1,
        Ac_left_ref = rr$Ac_pre, Ap_left_ref = rr$Ap_pre,
        Ac_implied_lin_at_t = P$C_lin * rr$V1,
        Ac_implied_ode_at_t = P$C_ode * rr$V1,
        Ac_amount_abs_diff_lin_vs_ref = abs(P$C_lin * rr$V1 - rr$Ac_pre),
        Ac_amount_abs_diff_ode_vs_ref = abs(P$C_ode * rr$V1 - rr$Ac_pre),
        C_lin_before = if (nrow(Pb)) Pb$C_lin else NA_real_, C_lin_at = P$C_lin,
        # what naive concentration-continuity would have predicted at t (non-gating illustration)
        C_if_conc_continuous = if (nrow(Pb)) Pb$C_lin else NA_real_,
        amount_tol_mg = rr$V1 * s1_tol(P$C_ode))
    }
  }

  # ---- non-gating: harness sensitivity (negative controls) + solver sensitivity ----
  nb_l <- solve_one("lin", model, interp = "nocb")
  nb_o <- solve_one("ode", model, interp = "nocb")
  neg[[model]] <- data.frame(model = model,
    max_abs_lin_nocb_vs_locf = maxerr_vs(nb_l$cp, lin$cp),
    max_abs_ode_nocb_vs_locf = maxerr_vs(nb_o$cp, ode$cp),
    n_points_exceeding_tol_lin_nocb = sum(abs(nb_l$cp - ode$cp) > s1_tol(ode$cp)),
    n_points_exceeding_tol_ode_nocb = sum(abs(nb_o$cp - ode$cp) > s1_tol(ode$cp)))
  sv <- list()
  for (cfg in list(list("liblsoda-default-tol", "liblsoda", 1e-8, 1e-6),
                   list("lsoda-1e-12", "lsoda", 1e-12, 1e-12),
                   list("dop853-1e-12", "dop853", 1e-12, 1e-12),
                   list("liblsoda-1e-12", "liblsoda", 1e-12, 1e-12))) {
    o <- tryCatch(solve_one("ode", model, "locf", cfg[[3]], cfg[[4]], cfg[[2]]), error = function(e) NULL)
    sv[[length(sv) + 1L]] <- data.frame(model = model, ode_config = cfg[[1]],
      status = if (is.null(o)) "solver error" else "ok",
      max_abs_ode_vs_ref = if (is.null(o)) NA else max(abs(o$cp - pts$C_ref)),
      max_abs_ode_vs_lin = if (is.null(o)) NA else max(abs(o$cp - lin$cp)),
      n_exceeding_gate_vs_ref = if (is.null(o)) NA else sum(abs(o$cp - pts$C_ref) > s1_tol(pts$C_ref)))
  }
  sens[[model]] <- do.call(rbind, sv)

  # ---- non-gating diagnostic: which interpolation convention does each engine implement? ----
  refs_nb <- do.call(rbind, lapply(cases, function(cs) {
    r <- ref_solve(cs$ev, model, S1_PARS[[model]], cs$eCL, cs$eV, interval_cov = "nocb"); r$ID <- cs$ev$ID[1]; r }))
  C_ref_nb <- refs_nb$cp[match(key(lin$ID, lin$time), key(refs_nb$ID, refs_nb$time))]
  first_chg <- vapply(cases, function(cs) { e <- cs$ev
    i <- which(c(FALSE, diff(e$CRCA) != 0 | diff(e$AJBW) != 0)); if (length(i)) e$TIME[i[1]] else Inf }, 0)
  after <- pts$time >= first_chg[match(pts$ID, icov$ID)]
  lin_nm <- as.data.frame(rxSolve(mods$lin[[model]], S1_PARS[[model]], events, iCov = icov, covsInterpolation = "locf",
                                  method = "liblsoda", addDosing = TRUE, cores = 1, nonmem = TRUE))
  semantics[[model]] <- data.frame(model = model,
    n_points_before_first_cov_change = sum(!after), n_points_at_or_after = sum(after),
    n_gate_fail_before_first_change = sum(!pts$pass[!after]), n_gate_fail_at_or_after = sum(!pts$pass[after]),
    max_abs_lin_vs_ref_LOCF  = max(abs(lin$cp - pts$C_ref)),
    max_abs_lin_vs_ref_NOCB  = max(abs(lin$cp - C_ref_nb)),
    max_abs_ode_locf_vs_ref_LOCF = max(abs(ode$cp - pts$C_ref)),
    max_abs_ode_nocb_vs_ref_NOCB = max(abs(nb_o$cp - C_ref_nb)),
    max_abs_lin_nonmemTRUE_vs_default = max(abs(lin_nm$cp - lin$cp)),
    n_lin_pts_within_gate_of_NOCB_ref = sum(abs(lin$cp - C_ref_nb) <= s1_tol(C_ref_nb)),
    n_points = nrow(pts))
}
pts <- do.call(rbind, res)
steps <- do.call(rbind, steps_tbl)
neg <- do.call(rbind, neg); sens <- do.call(rbind, sens); semantics <- do.call(rbind, semantics)

# ---- scenario-coverage assertions (design integrity of S-1) ------------------------------
cv <- list()
for (model in c("1c", "2c")) for (lab in c("canonical_tie_sample_first", "canonical_tie_dose_first")) {
  cs <- Filter(function(x) x$label == lab, cases)[[1]]; id <- cs$ev$ID[1]
  s <- steps[steps$model == model & steps$case == lab, ]
  ev <- cs$ev
  cv[[length(cv) + 1L]] <- data.frame(model = model, case = lab,
    a_two_infusions = sum(ev$EVID == 1) + sum(ev$EVID == 4) >= 2,
    b_step_outside_infusion = any(!s$during_active_infusion),
    c_step_inside_infusion = any(s$during_active_infusion),
    d_obs_before_after_each_step = all(vapply(s$t_change, function(t)
      any(abs(ev$TIME - (t - EPS)) < 1e-9 & ev$EVID == 0) && any(abs(ev$TIME - (t + EPS)) < 1e-9 & ev$EVID == 0), NA)),
    e_sample_dose_tie = any(duplicated(ev$TIME) & TRUE) &&
      any(vapply(unique(ev$TIME[ev$EVID == 1]), function(t) any(ev$EVID == 0 & ev$TIME == t), NA)),
    f_episode_reset = any(ev$EVID == 4),
    g_nonzero_Ac_at_all_steps = all(s$Ac_left_ref > 0),
    g_nonzero_Ap_at_all_steps = if (model == "2c") all(s$Ap_left_ref > 0) else NA,
    g_min_Ac_left = min(s$Ac_left_ref), g_min_Ap_left = min(s$Ap_left_ref))
}
cov_tbl <- do.call(rbind, cv)

# ---- verdict -------------------------------------------------------------------------------
summ <- function(p) data.frame(
  n_points = nrow(p), n_fail = sum(!p$pass), n_ref_fail = sum(!p$ref_pass),
  max_abs_err = max(p$abs_err, na.rm = TRUE), max_rel_err = max(p$rel_err, na.rm = TRUE),
  max_ode_vs_ref_abs = max(p$ode_vs_ref_abs, na.rm = TRUE),
  max_lin_vs_ref_abs = max(p$lin_vs_ref_abs, na.rm = TRUE), n_na = sum(is.na(p$abs_err)))
by_model <- do.call(rbind, lapply(split(pts, pts$model), function(p) cbind(model = p$model[1], summ(p))))
canon_pts <- pts[grepl("^canonical", pts$case), ]
by_model_canon <- do.call(rbind, lapply(split(canon_pts, canon_pts$model), function(p) cbind(model = p$model[1], summ(p))))
steps$gate_pass <- steps$Ac_amount_abs_diff_lin_vs_ref <= steps$amount_tol_mg &
                   steps$Ac_amount_abs_diff_ode_vs_ref <= steps$amount_tol_mg
steps$amount_continuity_pass <- steps$gate_pass
coverage_ok <- all(unlist(cov_tbl[, c("a_two_infusions","b_step_outside_infusion","c_step_inside_infusion",
  "d_obs_before_after_each_step","e_sample_dose_tie","f_episode_reset","g_nonzero_Ac_at_all_steps")])) &&
  all(cov_tbl$g_nonzero_Ap_at_all_steps[cov_tbl$model == "2c"])
failing <- pts[!pts$pass | !pts$ref_pass, ]
gate <- list(
  linCmt_vs_ODE_all_points = all(pts$pass),
  ODE_vs_independent_reference_all_points = all(pts$ref_pass),
  step_amount_continuity = all(steps$gate_pass),
  scenario_coverage_a_to_g = isTRUE(coverage_ok),
  no_NA_comparisons = !anyNA(pts$abs_err))
VERDICT <- if (all(unlist(gate))) "PASS" else "FAIL"

# ---- evidence ------------------------------------------------------------------------------------
dir.create("evidence", showWarnings = FALSE)
w <- function(x, f) write.csv(x, file.path("evidence", f), row.names = FALSE)
w(events, "s1_events.csv"); w(icov, "s1_individual_etas.csv"); w(pts, "s1_comparison_points.csv")
w(failing, "s1_failing_points.csv"); w(steps, "s1_covariate_step_checks.csv")
w(cov_tbl, "s1_scenario_coverage.csv"); w(neg, "s1_negative_controls_nocb.csv")
w(sens, "s1_ode_solver_sensitivity.csv"); w(semantics, "s1_interpolation_semantics_diagnostic.csv")
case_tbl <- do.call(rbind, lapply(split(pts, list(pts$model, pts$case), drop = TRUE), function(p)
  cbind(model = p$model[1], case = p$case[1], summ(p))))
w(case_tbl, "s1_case_summary.csv")

pk <- c("rxode2","nlmixr2est","nlmixr2","lotri","dparser","rxode2lincmt","symengine","PreciseSums","n1qn1","lbfgsb3c")
pkv <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")
rl <- jsonlite::fromJSON("renv.lock")
renv_src <- lapply(pk, function(p) { e <- rl$Packages[[p]]; if (is.null(e)) NULL else
  list(version = e$Version, source = e$Source, repository = e$Repository,
       hash = e$Hash, remote_sha = e$RemoteSha) }); names(renv_src) <- pk
interp_codes <- c(locf = rxControl(covsInterpolation = "locf")$covsInterpolation,
                  nocb = rxControl(covsInterpolation = "nocb")$covsInterpolation)
meth_codes <- c(liblsoda = rxControl(method = "liblsoda")$method, lsoda = rxControl(method = "lsoda")$method,
                dop853 = rxControl(method = "dop853")$method)
out <- list(
  gate = "S-1", verdict = VERDICT, synthetic_data_only = TRUE,
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec = c("PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md", "PKPD02_v11_addendum.md s3 (S-1)"),
  protocol_sha256 = list(
    protocol = digest::digest("protocol/PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md", "sha256", file = TRUE),
    addendum = digest::digest("protocol/PKPD02_v11_addendum.md", "sha256", file = TRUE)),
  renv_lock_sha256 = digest::digest("renv.lock", "sha256", file = TRUE),
  seed = list(value = S1_SEED, rng = RNGkind(), note = "implementer-chosen S-1 scenario seed; protocol assigns none for S-1"),
  tolerances = list(abs_mg_per_L = S1_TOL_ABS, rel = S1_TOL_REL, rule = "|C_lin - C_ode| <= max(abs, rel*|C_ode|)",
                    ode_atol = ODE_ATOL, ode_rtol = ODE_RTOL, ode_method = "liblsoda",
                    covsInterpolation = "locf", eps_before_after_h = EPS),
  synthetic_parameters = S1_PARS, n_random_patients = N_RANDOM, n_cases = length(cases),
  gating_checks = gate, by_model_all_cases = by_model, by_model_canonical_only = by_model_canon,
  n_failing_points = nrow(failing),
  failing_cases = unique(failing$case),
  step_checks = list(n = nrow(steps), n_fail = sum(!steps$gate_pass),
                     max_amount_abs_diff_lin_vs_ref_mg = max(steps$Ac_amount_abs_diff_lin_vs_ref),
                     max_amount_abs_diff_ode_vs_ref_mg = max(steps$Ac_amount_abs_diff_ode_vs_ref),
                     n_steps_during_active_infusion = sum(steps$during_active_infusion)),
  negative_controls_nocb = neg, ode_solver_sensitivity = sens,
  interpolation_semantics_diagnostic_nongating = semantics,
  package_versions = as.list(pkv), renv_lock_entries = renv_src,
  rxControl_codes = list(covsInterpolation = as.list(interp_codes), method = as.list(meth_codes),
                         default_covsInterpolation = rxControl()$covsInterpolation,
                         default_atol = rxControl()$atol, default_rtol = rxControl()$rtol,
                         default_method = rxControl()$method),
  R_version = R.version.string, platform = R.version$platform
)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "string"), "evidence/s1_results.json")
writeLines(capture.output(print(sessionInfo()), cat("\nBLAS/LAPACK:\n"), print(extSoftVersion()),
                          cat("\nLa_library:", La_library(), "\nLa_version:", La_version(), "\n")),
           "evidence/s1_sessionInfo.txt")
write.csv(data.frame(package = names(pkv), installed_version = pkv), "evidence/s1_package_versions.csv", row.names = FALSE)

cat("\n==== PKPD-02 Gate S-1 ====\n")
print(by_model, row.names = FALSE); cat("\ngating checks:\n"); print(unlist(gate))
cat("\nfailing points:", nrow(failing), "\nS-1 VERDICT:", VERDICT, "\n")
quit(status = if (VERDICT == "PASS") 0L else 1L)
