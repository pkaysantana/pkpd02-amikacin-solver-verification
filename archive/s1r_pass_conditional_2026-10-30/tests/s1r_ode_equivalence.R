#!/usr/bin/env Rscript
# PKPD-02 Gate S-1r (remediation): explicit amount-based ODE vs INDEPENDENT base-R reference.
# Re-runs the 102 archived synthetic scenarios (identical event tables, SHA-verified) plus
# boundary-event regression cases. SYNTHETIC DATA ONLY; no fitting.
# Acceptance rule (unchanged): |C_ODE - C_ref| <= max(1e-5 mg/L, 1e-6*|C_ref|) at every row.
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(normalizePath(file.path(dirname(f), "..")))
options(error = function() { cat("\nCALL STACK (innermost first):\n", paste(rev(vapply(sys.calls(), function(x) substr(paste(deparse(x), collapse = " "), 1, 150), "")), collapse = "\n"), "\n"); quit(status = 2) })
suppressMessages({ library(rxode2); library(nlmixr2est); library(jsonlite) })
source("R/pkpd02_ode_models.R"); source("R/pkpd02_solve.R"); source("R/pkpd02_reference_matexp.R")
options(rxode2.useLinCmt = NULL)

PRIMARY    <- list(atol = 1e-12, rtol = 1e-12, method = "liblsoda")          # high-accuracy reference setting
CANDIDATE  <- list(atol = 1e-10, rtol = 1e-10, method = "liblsoda")          # production candidate, declared BEFORE this run
MARGIN     <- 0.1                                                           # candidate must use <= 10% of the allowance
EPS        <- 1e-3
ETA_B      <- c(0.10, -0.05)                                                # boundary-case etas
PHYS <- list("1c" = c(tCL = 4.0, tV = 18.0), "2c" = c(tCL = 4.0, tV1 = 14.0, tQ = 3.0, tV2 = 9.0))
LOGTH <- list("1c" = c(tcl = log(4), tv = log(18)), "2c" = c(tcl = log(4), tv1 = log(14), tq = log(3), tv2 = log(9)))

# ---- 0. archived inputs, SHA-verified ---------------------------------------------------
ARCH <- "archive/s1_failed_2026-10-30"
man <- read.table(file.path(ARCH, "MANIFEST.sha256"), col.names = c("sha", "path"), stringsAsFactors = FALSE)
verify <- function(rel) { p <- file.path(ARCH, rel); m <- man$sha[man$path == paste0("./", rel)]
  stopifnot(length(m) == 1, digest::digest(p, "sha256", file = TRUE) == m); m }
sha_events <- verify("evidence/s1_events.csv"); sha_etas <- verify("evidence/s1_individual_etas.csv")
scen <- read.csv(file.path(ARCH, "evidence/s1_events.csv")); et <- read.csv(file.path(ARCH, "evidence/s1_individual_etas.csv"))
scen$CMT <- as.integer(scen$CMT); scen$EVID <- as.integer(scen$EVID)
stopifnot(length(unique(scen$ID)) == 102)

# ---- 1. boundary-event regression cases (explicit per-row covariates) --------------------
B <- function(label, rows) { x <- as.data.frame(do.call(rbind, rows)); names(x) <- c("TIME", "EVID", "AMT", "RATE", "CRCA", "AJBW")
  x$EVID <- as.integer(x$EVID); x <- x[order(x$TIME, seq_len(nrow(x))), ]  # stable: tie order as authored
  data.frame(ID = NA_integer_, x[1:4], CMT = 1L, x[5:6], case = label) }
r <- function(t, evid = 0, amt = 0, rate = 0, crca = 90, w = 75) c(t, evid, amt, rate, crca, w)
bnd <- list(
  B("B01_step_at_infusion_start",      list(r(0,1,600,300), r(1.5), r(2.999), r(3,1,700,350,45,95), r(3,0,0,0,45,95), r(3.001,0,0,0,45,95), r(4,0,0,0,45,95), r(8,0,0,0,45,95))),
  B("B02_step_at_infusion_end",        list(r(0,1,600,300), r(1.999), r(2,2,0,0,50,100), r(2,0,0,0,50,100), r(2.001,0,0,0,50,100), r(4,0,0,0,50,100))),
  B("B03_sample_at_infusion_end",      list(r(0,1,600,300), r(1.999), r(2), r(2.001), r(3.5,2,0,0,40,110), r(3.499), r(3.5,0,0,0,40,110), r(5,0,0,0,40,110))),
  B("B04a_sample_before_dose_tie",     list(r(0,1,600,300), r(3), r(3,1,700,175), r(3.001), r(5), r(9))),
  B("B04b_dose_before_sample_tie",     list(r(0,1,600,300), r(3,1,700,175), r(3), r(3.001), r(5), r(9))),
  B("B05_two_infusions_same_end",      list(r(0,1,400,100), r(1,1,300,100), r(1.999), r(2,2,0,0,60,90), r(2.001,0,0,0,60,90), r(4,2,0,0,120,50), r(3.999,0,0,0,60,90), r(4,0,0,0,120,50), r(4.001,0,0,0,120,50), r(7,0,0,0,120,50))),
  B("B06_two_infusions_same_start",    list(r(0,1,300,100), r(0,1,200,50), r(1), r(1.5,2,0,0,70,95), r(1.501,0,0,0,70,95), r(3.5,0,0,0,70,95), r(6,0,0,0,70,95))),
  B("B07_back_to_back_infusions",      list(r(0,1,600,300), r(1.999), r(2,1,600,200,50,100), r(2,0,0,0,50,100), r(2.001,0,0,0,50,100), r(4,0,0,0,50,100), r(5,0,0,0,50,100), r(8,0,0,0,50,100))),
  B("B08_microsecond_step_sequence",   list(r(0,1,600,100), r(2.9), r(2.999999), r(3,2,0,0,60,80), r(3.0000005,0,0,0,60,80), r(3.000001,2,0,0,70,90), r(3.0000015,0,0,0,70,90), r(3.000002,2,0,0,80,100), r(3.0000025,0,0,0,80,100), r(3.1,0,0,0,80,100), r(4,0,0,0,80,100))),
  B("B09_covariate_changes_on_EVID2_plus_sample_rows", list(r(0,1,600,100), r(1), r(2,2,0,0,45,95), r(2,0,0,0,45,95), r(3,0,0,0,45,95), r(4,2,0,0,120,60), r(4,0,0,0,120,60), r(5,0,0,0,120,60), r(8,0,0,0,120,60))),
  B("B10_step_1e-9_after_first_dose",  list(r(0,1,500,250), r(1e-9,2,0,0,45,100), r(1e-6,0,0,0,45,100), r(0.5,0,0,0,45,100), r(3,0,0,0,45,100))),
  B("B11_pure_reset_EVID3_then_dose",  list(r(0,1,600,300), r(3,0,0,0,90,75), r(5,3,0,0,70,82), r(5.001,0,0,0,70,82), r(6,1,400,400,70,82), r(6.5,0,0,0,70,82), r(8,0,0,0,70,82))),
  B("B12_reset_dose_after_infusion_ended", list(r(0,1,600,300), r(1.999), r(3,4,200,200,60,90), r(3.001,0,0,0,60,90), r(4.5,0,0,0,60,90), r(6.5,0,0,0,60,90), r(10,0,0,0,60,90))),
  B("B13_EVID3_after_infusion_ended_then_dose", list(r(0,1,600,300), r(2.999), r(3,3,0,0,60,90), r(3.001,0,0,0,60,90), r(5,1,300,150,60,90), r(8,0,0,0,60,90))),
  B("B14_long_and_tiny_infusions", list(r(0,1,960,20), r(10,2,0,0,50,100), r(19.999,0,0,0,50,100), r(20,1,1,1000,50,100), r(20.0005,0,0,0,50,100), r(20.001,0,0,0,50,100), r(29.999,0,0,0,50,100), r(30,2,0,0,110,65), r(30,0,0,0,110,65), r(30.001,0,0,0,110,65), r(47,0,0,0,110,65), r(60,0,0,0,110,65))),
  B("B19_short_covariate_blip_on_EVID2_rows_0.1h", list(r(0,1,960,20,50,100), r(10,2,0,0,50,100), r(19.9,2,0,0,90,75), r(19.9,0,0,0,90,75), r(20,2,0,0,50,100), r(20,0,0,0,50,100), r(30,0,0,0,50,100))),
  B("B20_short_covariate_blip_on_EVID2_rows_0.001h", list(r(0,1,960,20,50,100), r(10,2,0,0,50,100), r(19.999,2,0,0,90,75), r(19.999,0,0,0,90,75), r(20,2,0,0,50,100), r(20,0,0,0,50,100), r(30,0,0,0,50,100))),
  B("B15_dose_row_carries_step_inside_infusion", list(r(0,1,600,60), r(3.999), r(4,1,180,90,45,100), r(4.001,0,0,0,45,100), r(5,2,0,0,100,70), r(4.999,0,0,0,45,100), r(5.001,0,0,0,100,70), r(7,0,0,0,100,70), r(12,0,0,0,100,70))),
  B("B16_sample_step_and_infusion_end_coincide", list(r(0,1,600,300), r(1.999), r(2,2,0,0,45,95), r(2,0,0,0,45,95), r(2.001,0,0,0,45,95), r(6,0,0,0,45,95))),
  B("B17_sample_at_time0_both_orders", list(r(0), r(0,1,600,300), r(0.5), r(2.5))),
  B("B18_sample_at_step_and_infusion_start", list(r(0,1,600,300), r(3), r(5,0,0,0,50,100), r(5,1,500,250,50,100), r(5.001,0,0,0,50,100), r(8,0,0,0,50,100))))
bnd <- do.call(rbind, Map(function(i, x) { x$ID <- 1000L + i; x }, seq_along(bnd), bnd))
# ---- 1b. HAZARD cases: inputs the guarded wrapper must REJECT (unguarded ODE known to disagree) ----
hz <- list(
  B("H1_covariate_change_on_obs_row_only_0.001h_then_reverts", list(r(0,1,960,20,50,100), r(10,2,0,0,50,100), r(19.999,0,0,0,90,75), r(20,0,0,0,50,100), r(30,0,0,0,50,100))),
  B("H2_covariate_change_on_obs_row_only_0.1h_then_reverts",   list(r(0,1,960,20,50,100), r(10,2,0,0,50,100), r(19.9,0,0,0,90,75), r(20,0,0,0,50,100), r(30,0,0,0,50,100))),
  B("H3_EVID4_reset_dose_inside_active_infusion",              list(r(0,1,600,100), r(2.999), r(3,4,200,200,60,90), r(3.001,0,0,0,60,90), r(4.5,0,0,0,60,90), r(6.5,0,0,0,60,90), r(10,0,0,0,60,90))),
  B("H4_EVID3_reset_inside_active_infusion",                   list(r(0,1,600,100), r(2.999), r(3,3,0,0,60,90), r(3.001,0,0,0,60,90), r(5,0,0,0,60,90), r(8,0,0,0,60,90))),
  B("H5_EVID3_reset_exactly_at_infusion_end",                  list(r(0,1,600,100), r(5.999), r(6,3), r(6.001), r(8))),
  B("H6_EVID4_reset_dose_exactly_at_infusion_end",             list(r(0,1,600,100), r(5.999), r(6,4,200,200), r(6.001), r(8))))
hz <- do.call(rbind, Map(function(i, x) { x$ID <- 2000L + i; x }, seq_along(hz), hz))
hz_cases <- unique(hz[, c("ID", "case")])

case_of <- c(setNames(rep("scenario", 0), character()))
bnd_cases <- unique(bnd[, c("ID", "case")])
events <- rbind(scen, bnd[, names(scen)])
icov_all <- rbind(data.frame(ID = et$ID, eCL = et$eCL, eV = et$eV), data.frame(ID = bnd_cases$ID, eCL = ETA_B[1], eV = ETA_B[2]))
suite_of <- function(id) ifelse(id >= 1000, "boundary", "scenario")
label_of <- c(setNames(paste0("scenario_ID", et$ID), et$ID), setNames(bnd_cases$case, bnd_cases$ID))
label_of[as.character(c(1, 2))] <- c("canonical_tie_sample_first", "canonical_tie_dose_first")
ve <- pkpd02_validate_events(events); if (length(ve)) stop("event validation: ", paste(ve, collapse = " ; "))

# ---- 2. models and solving paths -----------------------------------------------------------
ui  <- list("1c" = suppressMessages(rxode(pkpd02_1c_ode)), "2c" = suppressMessages(rxode(pkpd02_2c_ode)))
plain <- lapply(pkpd02_rx_ode, function(f) f())
po  <- lapply(ui, function(u) u$foceiModel$predOnly)
path_setup <- function(path, s) {
  if (path == "P1_plain_rxode2_ODE") list(model = plain[[s]], params = PHYS[[s]], icov = icov_all)
  else if (path == "P2_nlmixr2_ui_ODE") list(model = ui[[s]], params = LOGTH[[s]],
         icov = setNames(icov_all, c("ID", "eta.cl", if (s == "1c") "eta.v" else "eta.v1")))
  else { n <- length(LOGTH[[s]]); pr <- data.frame(ID = icov_all$ID, check.names = FALSE)
    for (i in seq_len(n)) pr[[sprintf("THETA[%d]", i)]] <- LOGTH[[s]][[i]]
    pr[[sprintf("THETA[%d]", n + 1)]] <- 0.5; pr[[sprintf("THETA[%d]", n + 2)]] <- 0.1
    pr[["ETA[1]"]] <- icov_all$eCL; pr[["ETA[2]"]] <- icov_all$eV
    list(model = po[[s]], params = pr, icov = NULL) } }
PATHS <- c("P1_plain_rxode2_ODE", "P2_nlmixr2_ui_ODE", "P3_nlmixr2_focei_predOnly_ODE")
run_path <- function(path, s, cfg = PRIMARY) { su <- path_setup(path, s)
  o <- pkpd02_solve_ode(su$model, s, su$params, events, su$icov, cfg$atol, cfg$rtol, cfg$method)
  o$C <- if ("cp" %in% names(o)) o$cp else o$rx_pred_; o }
keyf <- function(id, t) paste(id, sprintf("%.12f", t))

# ---- 3. independent reference per case -----------------------------------------------------
refs <- list()
for (s in c("1c", "2c")) {
  rr <- lapply(split(events, events$ID), function(e) { id <- e$ID[1]; ic <- icov_all[icov_all$ID == id, ]
    x <- ref_matexp(e, s, PHYS[[s]], c(ic$eCL, ic$eV)); x$ID <- id; x })
  refs[[s]] <- do.call(rbind, rr) }

# ---- 4. primary comparison: every output row, every path, both models ---------------------------
cmp <- list(); steps <- list(); pw <- list()
for (s in c("1c", "2c")) for (p in PATHS) {
  o <- run_path(p, s); rf <- refs[[s]]; m <- match(keyf(o$ID, o$time), keyf(rf$ID, rf$TIME)); stopifnot(!anyNA(m))
  pw[[paste(s, p)]] <- as.list(attr(o, "pathway_checks"))
  C_ref <- rf$cp[m]; V1 <- rf$V1[m]; at <- V1 * pkpd02_tol(C_ref)
  d <- data.frame(model = s, path = p, suite = suite_of(o$ID), case = unname(label_of[as.character(o$ID)]), ID = o$ID, time = o$time, evid = o$evid,
                  C_ode = o$C, C_ref = C_ref, Ac_ode = o$Ac, Ac_ref = rf$Ac[m], Ap_ode = if (s == "2c") o$Ap else NA_real_, Ap_ref = if (s == "2c") rf$Ap[m] else NA_real_)
  d$abs_err <- abs(d$C_ode - d$C_ref); d$rel_err <- ifelse(d$C_ref == 0, ifelse(d$abs_err == 0, 0, Inf), d$abs_err / abs(d$C_ref))
  d$tol <- pkpd02_tol(d$C_ref); d$pass <- !is.na(d$abs_err) & d$abs_err <= d$tol
  d$amount_tol <- at; d$Ac_abs_err <- abs(d$Ac_ode - d$Ac_ref); d$Ap_abs_err <- abs(d$Ap_ode - d$Ap_ref)
  d$amount_pass <- d$Ac_abs_err <= at & (s == "1c" | d$Ap_abs_err <= at)
  cmp[[paste(s, p)]] <- d
  # step checks: amounts at the covariate-change row vs reference LEFT-LIMIT amounts (continuity)
  for (id in unique(events$ID)) { e <- events[events$ID == id, ]; epi <- cumsum(e$EVID %in% c(3L, 4L))
    chg <- which(c(FALSE, diff(e$CRCA) != 0 | diff(e$AJBW) != 0) & !(e$EVID %in% c(3L, 4L)))
    for (j in chg) { t <- e$TIME[j]; k <- which(o$ID == id & abs(o$time - t) < 1e-12)[1]; rj <- rf[rf$ID == id & rf$TIME == t, ]
      prevV <- rf$V1[rf$ID == id & rf$TIME == max(rf$TIME[rf$ID == id & rf$TIME < t])]
      dose <- e$EVID %in% c(1L, 4L) & e$AMT > 0 & epi == epi[j] & e$TIME < t & e$TIME + e$AMT / pmax(e$RATE, 1e-300) > t
      tolA <- rj$V1 * pkpd02_tol(rj$cp)
      steps[[length(steps) + 1L]] <- data.frame(model = s, path = p, suite = suite_of(id), case = unname(label_of[as.character(id)]), ID = id, t_change = t,
        during_active_infusion = any(dose), V_old = prevV, V_new = rj$V1, abs_ln_V_ratio = abs(log(rj$V1 / prevV)),
        Ac_left_ref = rj$Ac_left, Ap_left_ref = rj$Ap_left, Ac_ode_at_change = o$Ac[k], Ap_ode_at_change = if (s == "2c") o$Ap[k] else NA_real_,
        Ac_abs_diff = abs(o$Ac[k] - rj$Ac_left), Ap_abs_diff = if (s == "2c") abs(o$Ap[k] - rj$Ap_left) else NA_real_, amount_tol = tolA,
        # illustration only: concentration just before the change times V_new would be wrong if concentration were carried
        C_ode_at_change = o$C[k], C_if_concentration_continuous = rj$Ac_left / prevV,
        pass = abs(o$Ac[k] - rj$Ac_left) <= tolA && (s == "1c" || abs(o$Ap[k] - rj$Ap_left) <= tolA)) } }
}
cmp <- do.call(rbind, cmp); steps <- do.call(rbind, steps)

# ---- 5. concordance among the three ODE paths ---------------------------------------------------
conc <- do.call(rbind, lapply(c("1c", "2c"), function(s) { a <- cmp[cmp$model == s & cmp$path == PATHS[1], ]; do.call(rbind, lapply(PATHS[-1], function(p) {
  b <- cmp[cmp$model == s & cmp$path == p, ]; stopifnot(nrow(a) == nrow(b), all(a$ID == b$ID), all(a$time == b$time))
  data.frame(model = s, comparison = paste(PATHS[1], "vs", p), n_points = nrow(a), max_abs_diff = max(abs(a$C_ode - b$C_ode)),
             n_exceeding_gate = sum(abs(a$C_ode - b$C_ode) > pkpd02_tol(a$C_ref))) })) }))

# ---- 6. tolerance sweep (production-tolerance adequacy) on ALL points via the production predOnly model ----
cfgs <- list(list("nlmixr2_foceiControl_default", 1e-6, 1e-3, "liblsoda"), list("rxode2_default", 1e-8, 1e-6, "liblsoda"),
             list("1e-8", 1e-8, 1e-8, "liblsoda"), list("1e-10", 1e-10, 1e-10, "liblsoda"), list("1e-12", 1e-12, 1e-12, "liblsoda"),
             list("lsoda_1e-10", 1e-10, 1e-10, "lsoda"), list("dop853_1e-10", 1e-10, 1e-10, "dop853"), list("lsoda_1e-12", 1e-12, 1e-12, "lsoda"), list("dop853_1e-12", 1e-12, 1e-12, "dop853"))
sweep <- do.call(rbind, lapply(cfgs, function(cf) do.call(rbind, lapply(c("1c", "2c"), function(s) {
  o <- tryCatch(run_path(PATHS[3], s, list(atol = cf[[2]], rtol = cf[[3]], method = cf[[4]])), error = function(e) conditionMessage(e))
  if (is.character(o)) return(data.frame(config = cf[[1]], atol = cf[[2]], rtol = cf[[3]], method = cf[[4]], model = s, status = o, n_points = NA, n_fail = NA, max_abs_err = NA, max_err_over_allowance = NA))
  rf <- refs[[s]]; m <- match(keyf(o$ID, o$time), keyf(rf$ID, rf$TIME)); e <- abs(o$C - rf$cp[m]); al <- pkpd02_tol(rf$cp[m])
  data.frame(config = cf[[1]], atol = cf[[2]], rtol = cf[[3]], method = cf[[4]], model = s, status = "ok", n_points = length(e), n_fail = sum(e > al),
             max_abs_err = max(e), max_err_over_allowance = max(e / al)) }))))

# ---- 7. guard checks (the guards must work AND the harness must be able to detect the old failure) ----
g <- list(); add <- function(n, ok, detail = "") g[[length(g) + 1L]] <<- data.frame(check = n, pass = isTRUE(ok), detail = detail)
ev1 <- events[events$ID == 1, ]; ic1 <- setNames(icov_all[icov_all$ID == 1, ], c("ID", "eta.cl", "eta.v")); rf1 <- refs[["1c"]][refs[["1c"]]$ID == 1, ]
bad <- ev1; bad$CRCA[bad$TIME == 7][1] <- 10
add("wrapper rejects tied rows with different covariates", inherits(try(pkpd02_solve_ode(ui[["1c"]], "1c", LOGTH[["1c"]], bad, ic1), silent = TRUE), "try-error"))
bolus <- ev1; bolus$RATE[1] <- 0
add("wrapper rejects AMT>0 with RATE=0 (T7 is preprocessing)", inherits(try(pkpd02_solve_ode(ui[["1c"]], "1c", LOGTH[["1c"]], bolus, ic1), silent = TRUE), "try-error"))
lin <- rxode2({ CL <- tCL * CRCA / 100 * exp(eCL); V <- tV * AJBW / 70 * exp(eV); cp <- linCmt() })
add("wrapper rejects a linCmt() model", inherits(try(pkpd02_solve_ode(lin, "1c", PHYS[["1c"]], ev1, icov_all[icov_all$ID == 1, ]), silent = TRUE), "try-error"))
add("wrapper refuses covsInterpolation='nocb'", inherits(try(pkpd02_solve_ode(ui[["1c"]], "1c", LOGTH[["1c"]], ev1, ic1, covsInterpolation = "nocb"), silent = TRUE), "try-error"))
options(rxode2.useLinCmt = TRUE); w <- try(pkpd02_solve_ode(ui[["1c"]], "1c", LOGTH[["1c"]], ev1, ic1), silent = TRUE); options(rxode2.useLinCmt = NULL)
add("wrapper unaffected by options(rxode2.useLinCmt=TRUE)", !inherits(w, "try-error") && all(abs(w$cp - rf1$cp[match(keyf(1, w$time), keyf(1, rf1$TIME))]) <= pkpd02_tol(rf1$cp[match(keyf(1, w$time), keyf(1, rf1$TIME))])))
hzlin <- suppressMessages(as.data.frame(rxSolve(ui[["1c"]], LOGTH[["1c"]], ev1, iCov = ic1, covsInterpolation = "locf", atol = 1e-12, rtol = 1e-12, addDosing = TRUE, cores = 1, useLinCmt = TRUE)))
mm <- match(keyf(1, hzlin$time), keyf(1, rf1$TIME)); add("harness detects silent ODE->linCmt conversion (raw useLinCmt=TRUE fails gate)", any(abs(hzlin$cp - rf1$cp[mm]) > pkpd02_tol(rf1$cp[mm])), sprintf("max abs %.3e mg/L", max(abs(hzlin$cp - rf1$cp[mm]))))
nb <- suppressMessages(as.data.frame(rxSolve(plain[["1c"]], PHYS[["1c"]], ev1, iCov = icov_all[icov_all$ID == 1, ], covsInterpolation = "nocb", atol = 1e-12, rtol = 1e-12, addDosing = TRUE, cores = 1)))
mm <- match(keyf(1, nb$time), keyf(1, rf1$TIME)); add("harness detects NOCB (ODE with covsInterpolation='nocb' fails gate)", any(abs(nb$cp - rf1$cp[mm]) > pkpd02_tol(rf1$cp[mm])), sprintf("max abs %.3e mg/L", max(abs(nb$cp - rf1$cp[mm]))))
ll <- suppressMessages(as.data.frame(rxSolve(lin, PHYS[["1c"]], ev1, iCov = icov_all[icov_all$ID == 1, ], covsInterpolation = "locf", addDosing = TRUE, cores = 1)))
mm <- match(keyf(1, ll$time), keyf(1, rf1$TIME)); add("legacy linCmt() S-1 failure reproduces vs independent reference (informational)", any(abs(ll$cp - rf1$cp[mm]) > pkpd02_tol(rf1$cp[mm])), sprintf("max abs %.3e mg/L", max(abs(ll$cp - rf1$cp[mm]))))
guards <- do.call(rbind, g)


# ---- 7b. hazard cases: wrapper must reject; unguarded ODE disagreement is recorded (informational) ----
hzr <- do.call(rbind, lapply(c("1c", "2c"), function(s) do.call(rbind, lapply(hz_cases$ID, function(id) {
  e <- hz[hz$ID == id, names(scen)]; ic <- data.frame(ID = id, eCL = ETA_B[1], eV = ETA_B[2]); ve <- pkpd02_validate_events(e)
  o <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], e, ic, validate = FALSE)       # UNGUARDED, for demonstration only
  rf <- ref_matexp(e, s, PHYS[[s]], c(ic$eCL, ic$eV)); m <- match(keyf(id, o$time), keyf(id, rf$TIME)); er <- abs(o$cp - rf$cp[m])
  data.frame(model = s, ID = id, case = hz_cases$case[hz_cases$ID == id], wrapper_rejects = length(ve) > 0, rejection_reason = paste(ve, collapse = " | "),
             unguarded_n_rows = nrow(o), unguarded_n_exceeding_gate = sum(er > pkpd02_tol(rf$cp[m])), unguarded_max_abs_err_mg_per_L = max(er),
             unguarded_min_Ac_mg = min(o$Ac), unguarded_disagrees_with_reference = any(er > pkpd02_tol(rf$cp[m]))) }))))

# ---- 8. non-gating probe: sample row vs reset row at the same TIME (not specified by protocol) ----
pr <- function(order) { e <- if (order == "sample_row_first") data.frame(TIME = c(0, 3, 3, 3.001), EVID = c(1L, 0L, 4L, 0L), AMT = c(600, 0, 200, 0), RATE = c(300, 0, 200, 0))
                        else data.frame(TIME = c(0, 3, 3, 3.001), EVID = c(1L, 4L, 0L, 0L), AMT = c(600, 200, 0, 0), RATE = c(300, 200, 0, 0))
  e$ID <- 1L; e$CMT <- 1L; e$CRCA <- 90; e$AJBW <- 75
  o <- pkpd02_solve_ode(plain[["1c"]], "1c", PHYS[["1c"]], e, data.frame(ID = 1L, eCL = 0, eV = 0))
  rr <- ref_matexp(e, "1c", PHYS[["1c"]], c(0, 0)); data.frame(order = order, ode_C_at_t3_obs_rows = paste(round(o$cp[o$time == 3 & o$evid == 0], 6), collapse = ","),
    pre_reset_C = round(rr$Ac_left[rr$TIME == 3] / rr$V1[rr$TIME == 3], 6), post_reset_C = round(rr$cp[rr$TIME == 3], 6)) }
probe <- rbind(pr("sample_row_first"), pr("reset_row_first"))

# ---- 9. verdict ---------------------------------------------------------------------------------------
summ <- function(d) data.frame(n_points = nrow(d), n_fail = sum(!d$pass), max_abs_err = max(d$abs_err), max_rel_err = max(d$rel_err[is.finite(d$rel_err)]),
                               n_amount_fail = sum(!d$amount_pass), max_Ac_abs_err_mg = max(d$Ac_abs_err), max_Ap_abs_err_mg = if (all(is.na(d$Ap_abs_err))) NA else max(d$Ap_abs_err, na.rm = TRUE))
tab <- do.call(rbind, lapply(split(cmp, list(cmp$model, cmp$path, cmp$suite), drop = TRUE), function(x) cbind(model = x$model[1], path = x$path[1], suite = x$suite[1], summ(x))))
rownames(tab) <- NULL
scen_pts <- cmp[cmp$suite == "scenario" & cmp$path == PATHS[1] & cmp$model == "1c", ]
n_cases_scen <- length(unique(cmp$ID[cmp$suite == "scenario"])); n_pts_per_model <- vapply(c("1c", "2c"), function(s) sum(cmp$model == s & cmp$path == PATHS[1] & cmp$suite == "scenario"), 1L)
cov_ids <- c(1L, 2L); cov_tbl <- do.call(rbind, lapply(c("1c", "2c"), function(s) do.call(rbind, lapply(cov_ids, function(id) { st <- steps[steps$model == s & steps$path == PATHS[2] & steps$ID == id, ]; e <- events[events$ID == id, ]
  data.frame(model = s, ID = id, a_two_infusions = sum(e$EVID %in% c(1L, 4L)) >= 2, b_step_outside_infusion = any(!st$during_active_infusion), c_step_inside_infusion = any(st$during_active_infusion),
             d_samples_before_after_each_step = all(vapply(st$t_change, function(t) any(abs(e$TIME - (t - EPS)) < 1e-9 & e$EVID == 0) && any(abs(e$TIME - (t + EPS)) < 1e-9 & e$EVID == 0), NA)),
             e_sample_dose_tie = any(vapply(unique(e$TIME[e$EVID == 1L]), function(t) any(e$EVID == 0 & e$TIME == t), NA)), f_episode_reset = any(e$EVID == 4L),
             g_nonzero_Ac_at_steps = all(st$Ac_left_ref > 0), g_nonzero_Ap_at_steps = if (s == "2c") all(st$Ap_left_ref > 0) else NA) }))))
cov_ok <- all(unlist(cov_tbl[, -(1:2)]), na.rm = TRUE)
cand <- sweep[sweep$config == "1e-10" & sweep$method == "liblsoda", ]
gate <- list(
  pathway_inspection_script_passes = file.exists("evidence/s1r_pathway_inspection.json") && isTRUE(fromJSON("evidence/s1r_pathway_inspection.json")$pass),
  scenario_count_102_and_3177_points_per_model = n_cases_scen == 102 && all(n_pts_per_model == 3177),
  all_points_pass_all_paths_primary_tol = all(cmp$pass),
  all_amount_states_pass_Ac_and_Ap = all(cmp$amount_pass),
  covariate_step_amount_continuity = all(steps$pass),
  boundary_regression_all_pass = all(cmp$pass[cmp$suite == "boundary"]),
  scenario_coverage_a_to_g = isTRUE(cov_ok),
  ode_paths_mutually_concordant = all(conc$n_exceeding_gate == 0),
  all_pathway_assertions_true = all(unlist(pw)),
  guard_checks_pass = all(guards$pass[!grepl("informational", guards$check)]),
  hazard_inputs_rejected_by_wrapper = all(hzr$wrapper_rejects))
tolerance_adequacy <- list(candidate = CANDIDATE, declared_before_run = TRUE, margin_requirement = sprintf("max error <= %.0f%% of allowance at every point", 100 * MARGIN),
  candidate_n_fail = sum(cand$n_fail), candidate_max_err_over_allowance = max(cand$max_err_over_allowance),
  ADEQUATE = all(cand$status == "ok") && all(cand$n_fail == 0) && max(cand$max_err_over_allowance) <= MARGIN,
  nlmixr2_default_tolerances_pass_gate = all(sweep$n_fail[sweep$config == "nlmixr2_foceiControl_default"] == 0, na.rm = TRUE))
VERDICT <- if (all(unlist(gate))) "PASS" else "FAIL"
input_contract <- list(
  "R-COV" = "Every covariate change must be carried by a solver-restart row (EVID 1-4); observation-only carriage is rejected.",
  "R-RESET" = "No reset row (EVID 3/4) while an infusion is active or ending at the reset time; protocol semantics for that case are undecided and the wrapper refuses it.",
  "tolerance" = "Production atol/rtol must be set explicitly; nlmixr2 foceiControl defaults (atol 1e-6, rtol 1e-3) fail the gate.")
conditions <- c("Pass is conditional on reviewer ratification of the input contract (R-COV, R-RESET).",
  "FOCEI inner (sensitivity) model accuracy is not tested here (no fitting); predOnly model only.",
  "Unguarded hazards remain in the installed rxode2 5.1.8 ODE pathway; see s1r_hazard_cases.csv.")

# ---- 10. evidence ---------------------------------------------------------------------------------------------
dir.create("evidence", FALSE); w <- function(x, f) write.csv(x, file.path("evidence", f), row.names = FALSE)
w(cmp, "s1r_comparison_points.csv"); w(cmp[!cmp$pass | !cmp$amount_pass, ], "s1r_failing_points.csv"); w(steps, "s1r_covariate_step_checks.csv")
w(tab, "s1r_summary_by_model_path_suite.csv"); w(conc, "s1r_ode_path_concordance.csv"); w(sweep, "s1r_tolerance_sweep.csv"); w(guards, "s1r_guard_checks.csv")
w(cov_tbl, "s1r_scenario_coverage.csv"); w(probe, "s1r_probe_sample_vs_reset_tie_nongating.csv")
w(events[events$ID >= 1000, ], "s1r_boundary_events.csv"); w(hzr, "s1r_hazard_cases.csv"); w(hz, "s1r_hazard_events.csv")
res <- list(gate = "S-1r", verdict = VERDICT, verdict_qualifier = if (VERDICT == "PASS") "PASS_CONDITIONAL_ON_INPUT_CONTRACT" else "FAIL", input_contract = input_contract, conditions = conditions, hazard_cases = hzr, synthetic_data_only = TRUE, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  gating_checks = gate, tolerance_rule = list(abs_mg_per_L = PKPD02_TOL$abs, rel = PKPD02_TOL$rel, reference = "C_ref (independent base-R matrix exponential)"),
  primary_ode_settings = PRIMARY, production_tolerance_assessment = tolerance_adequacy,
  n_scenarios = n_cases_scen, n_boundary_cases = nrow(bnd_cases), points_per_model_scenario = as.list(n_pts_per_model),
  points_per_model_boundary = list("1c" = sum(cmp$model == "1c" & cmp$path == PATHS[1] & cmp$suite == "boundary"), "2c" = sum(cmp$model == "2c" & cmp$path == PATHS[1] & cmp$suite == "boundary")),
  summary = tab, failing_cases = unique(cmp$case[!cmp$pass | !cmp$amount_pass]), step_checks = list(n = nrow(steps), n_fail = sum(!steps$pass),
    n_during_active_infusion = sum(steps$during_active_infusion & steps$path == PATHS[2]) , max_Ac_abs_diff_mg = max(steps$Ac_abs_diff), max_Ap_abs_diff_mg = max(steps$Ap_abs_diff, na.rm = TRUE)),
  path_concordance = conc, guard_checks = guards, archived_inputs_sha256 = list(events = sha_events, individual_etas = sha_etas),
  event_tables_identical_to_failed_S1 = TRUE)
writeLines(toJSON(res, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "string"), "evidence/s1r_equivalence_results.json")
cat("\n==== PKPD-02 Gate S-1r: explicit ODE vs independent reference ====\n"); print(tab, row.names = FALSE)
cat("\ngating checks:\n"); print(unlist(gate)); cat("\ntolerance sweep:\n"); print(sweep[, c("config", "model", "status", "n_fail", "max_abs_err", "max_err_over_allowance")], row.names = FALSE)
cat("\nhazard cases (unguarded ODE vs reference; wrapper rejection):\n"); print(hzr[, c("model", "case", "wrapper_rejects", "unguarded_n_exceeding_gate", "unguarded_max_abs_err_mg_per_L", "unguarded_min_Ac_mg")], row.names = FALSE)
cat("\nproduction candidate ADEQUATE:", tolerance_adequacy$ADEQUATE, "\nS-1r EQUIVALENCE VERDICT:", VERDICT, "\n")
quit(status = if (VERDICT == "PASS") 0L else 1L)
