#!/usr/bin/env Rscript
# Domain check: the REAL Strathclyde event structure (event times, doses/durations per v1.1 T7/T8, CRCA/AJBW timing,
# resets) run through the sanctioned pipeline with SYNTHETIC parameters and etas, compared with the independent
# per-episode reference. Not fitting; no concentration (DV) is read; only aggregate error summaries are written
# (no per-row output), to keep synthetic and empirical paths clearly separated (leakage gate L10).
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)]); setwd(normalizePath(file.path(dirname(f), "..")))
suppressMessages({ library(rxode2); library(nlmixr2est); library(jsonlite) })   # nlmixr2est must be attached for ui$foceiModel
source("R/pkpd02_ode_models.R"); source("R/pkpd02_solve.R"); source("R/pkpd02_reference_matexp.R"); source("R/pkpd02_events.R"); source("R/pkpd02_source_reader.R")
SEED <- 20261104L; set.seed(SEED)
PHYS <- list("1c" = c(tCL = 4.0, tV = 18.0), "2c" = c(tCL = 4.0, tV1 = 14.0, tQ = 3.0, tV2 = 9.0)); LOGTH <- list("1c" = c(tcl = log(4), tv = log(18)), "2c" = c(tcl = log(4), tv1 = log(14), tq = log(3), tv2 = log(9)))
plain <- lapply(pkpd02_rx_ode, function(f) f()); ui <- list("1c" = suppressMessages(rxode(pkpd02_1c_ode)), "2c" = suppressMessages(rxode(pkpd02_2c_ode))); po <- lapply(ui, function(u) u$foceiModel$predOnly)
ev <- pkpd02_read_source_metadata(); trail <- attr(ev, "trail"); ev$EPISODE <- pkpd02_episode_index(ev)
dur <- pkpd02_effective_duration(ev$AMT, ev$RATE, ev$EVID); ev$IS_DOSE <- ev$EVID %in% c(1L, 4L) & !is.na(ev$AMT) & ev$AMT > 0
tab <- ev; tab$DUR_RULE <- dur$rule; tab$AMT[!tab$IS_DOSE] <- 0; tab$RATE <- ifelse(tab$IS_DOSE, ifelse(dur$rule == "as_recorded", ev$RATE, tab$AMT / dur$dur), 0)
ids <- unique(tab$ID); eta <- data.frame(PATIENT_ID = ids, eCL = pmax(pmin(rnorm(length(ids), 0, .3), .6), -.6), eV = pmax(pmin(rnorm(length(ids), 0, .3), .6), -.6))
keyf <- function(p, e, t) paste(p, e, sprintf("%.9f", t))
tl <- pkpd02_episode_timelines(tab, pkpd02_classify_reset_ties(tab, "episode_relative"), "episode_relative"); cz <- pkpd02_canonicalize_events(tl$timelines); stopifnot(nrow(cz$rejected) == 0, length(pkpd02_validate_events(cz$events[, PKPD02_EV_COLS])) == 0)
ic <- pkpd02_icov_by_patient(tl$timelines, eta); pkpd02_assert_shared_eta(ic)
pl <- pkpd02_patient_level_table(tl$timelines, gap = 1); czp <- pkpd02_canonicalize_events(pl); stopifnot(nrow(czp$rejected) == 0, length(pkpd02_validate_events(czp$events[, PKPD02_EV_COLS])) == 0)
offs <- attr(pl, "offsets"); res <- list()
for (s in c("1c", "2c")) {
  rf <- do.call(rbind, lapply(split(tab, paste(tab$ID, tab$EPISODE)), function(e) { et <- eta[eta$PATIENT_ID == e$ID[1], ]; x <- ref_matexp(e[, PKPD02_EV_COLS], s, PHYS[[s]], c(et$eCL, et$eV)); x$PATIENT_ID <- e$ID[1]; x$EPISODE <- e$EPISODE[1]; x }))
  summ <- function(route, cp, ref) { er <- abs(cp - ref); data.frame(model = s, route = route, n_rows = length(er), n_fail = sum(er > pkpd02_tol(ref)), max_abs_err_mg_per_L = max(er), n_ref_below_1e3 = sum(abs(ref) < 1e-3), max_rel_err_where_ref_ge_1e3 = max((er / abs(ref))[abs(ref) >= 1e-3])) }   # relative error is meaningless at near-zero reference concentrations (the acceptance rule uses an absolute floor there)
  runtl <- function(model, params, icov, atol, rtol, col) { o <- pkpd02_solve_ode(model, s, params, cz$events[, PKPD02_EV_COLS], icov, atol, rtol); m <- match(keyf(o$ID %/% 100L, o$ID %% 100L, o$time), keyf(rf$PATIENT_ID, rf$EPISODE, rf$TIME)); stopifnot(!anyNA(m)); list(cp = o[[col]], ref = rf$cp[m]) }
  a <- runtl(plain[[s]], PHYS[[s]], ic[, c("ID", "eCL", "eV")], 1e-12, 1e-12, "cp"); res[[length(res) + 1L]] <- summ("episode_timelines_P1_atol1e-12", a$cp, a$ref)
  # production prediction model (predOnly) at the proposed production tolerance 1e-10
  pr <- data.frame(ID = ic$ID, check.names = FALSE); n <- length(LOGTH[[s]]); for (i in seq_len(n)) pr[[sprintf("THETA[%d]", i)]] <- LOGTH[[s]][[i]]; pr[[sprintf("THETA[%d]", n + 1)]] <- 0.5; pr[[sprintf("THETA[%d]", n + 2)]] <- 0.1; pr[["ETA[1]"]] <- ic$eCL; pr[["ETA[2]"]] <- ic$eV
  b <- runtl(po[[s]], pr, NULL, 1e-10, 1e-10, "rx_pred_"); res[[length(res) + 1L]] <- summ("episode_timelines_P3_predOnly_atol1e-10", b$cp, b$ref)
  # patient-level offset table (single ID per patient; offsets are a representation device)
  icp <- data.frame(ID = eta$PATIENT_ID, eCL = eta$eCL, eV = eta$eV); op <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], czp$events[, PKPD02_EV_COLS], icp, 1e-12, 1e-12)
  epi <- vapply(seq_len(nrow(op)), function(i) { of <- offs[offs$PATIENT_ID == op$ID[i], ]; max(of$EPISODE[of$offset <= op$time[i] + 1e-9]) }, 1L)
  rel <- op$time - vapply(seq_len(nrow(op)), function(i) offs$offset[offs$PATIENT_ID == op$ID[i] & offs$EPISODE == epi[i]], 0)
  m <- match(keyf(op$ID, epi, round(rel, 6)), keyf(rf$PATIENT_ID, rf$EPISODE, round(rf$TIME, 6))); stopifnot(!anyNA(m))
  res[[length(res) + 1L]] <- summ("patient_level_offset_table_P1_atol1e-12", op$cp, rf$cp[m]) }
tab_res <- do.call(rbind, res); rownames(tab_res) <- NULL
VERDICT <- if (all(tab_res$n_fail == 0)) "PASS" else "FAIL"
out <- list(test = "s2prep_real_structure_domain_check", verdict = VERDICT, scope = "real event STRUCTURE with SYNTHETIC parameters/etas; no DV read; no fitting; aggregate output only", seed = SEED, source_audit_trail = trail,
            n_patients = length(ids), n_episodes = length(unique(paste(tab$ID, tab$EPISODE))), n_event_rows = nrow(tab), canonicalisation_inserted_rows = nrow(cz$inserted), results = tab_res,
            pre_validated = list(episode_timelines_pass_validator = TRUE, patient_level_offset_table_pass_validator = TRUE))
dir.create("evidence", FALSE); writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA), "evidence/s2prep_real_structure_domain_check.json")
print(tab_res, row.names = FALSE); cat("DOMAIN CHECK:", VERDICT, "\n"); quit(status = if (VERDICT == "PASS") 0L else 1L)
