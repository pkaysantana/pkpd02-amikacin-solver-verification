#!/usr/bin/env Rscript
# S-2 preparation, decisions 1-4: canonical event representation vs the independent reference.
#  (1) R-COV: canonicalised (explicit EVID=2) tables == reference LOCF on the ORIGINAL table
#  (2) R-RESET: episode-specific timelines == reference reset semantics (zero state, cancelled infusions)
#  (3) patient identity: shared patient-level etas across episodes; per-episode etas refused
#  (4) reset-tied samples: metadata-only classification; indeterminate ones quarantined
# SYNTHETIC DATA ONLY. No DV exists anywhere in this test.
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(normalizePath(file.path(dirname(f), "..")))
suppressMessages({ library(rxode2); library(jsonlite) })
source("R/pkpd02_ode_models.R"); source("R/pkpd02_solve.R"); source("R/pkpd02_reference_matexp.R"); source("R/pkpd02_events.R")
SEED <- 20261103L; set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
PHYS <- list("1c" = c(tCL = 4.0, tV = 18.0), "2c" = c(tCL = 4.0, tV1 = 14.0, tQ = 3.0, tV2 = 9.0))
plain <- lapply(pkpd02_rx_ode, function(f) f())
keyf <- function(id, t) paste(id, sprintf("%.12f", t))
D <- function(...) data.frame(..., stringsAsFactors = FALSE)

# ---- random patient generator (raw, source-like) --------------------------------------------------
gen_patient <- function(id, allow_reset, covmode) {
  n_ep <- if (allow_reset) sample(1:3, 1) else 1L; rows <- list(); sched <- list(); truth <- list(reset_inside = FALSE, reset_at_end = FALSE, obs_only_changes = 0L)
  ts <- 0; pure_prev <- FALSE
  for (e in seq_len(n_ep)) {
    n_inf <- sample(1:3, 1); st <- ts + c(0, cumsum(runif(n_inf - 1, 0.2, 6))); dur <- runif(n_inf, 0.25, 3); rate <- runif(n_inf, 50, 400); amt <- rate * dur
    en <- st + dur; te <- max(en) + runif(1, 1, 6); last <- e == n_ep; tr <- NA_real_; kind <- "after"
    pure_next <- if (!last) runif(1) < 0.4 else FALSE   # whether the NEXT episode starts with a pure EVID=3 reset (its first infusion is then dropped)
    if (!last) { kind <- sample(c("after", "inside", "at_end"), 1, prob = c(.3, .4, .3)); k <- sample(n_inf, 1)
      if (e > 1 && pure_prev && k == 1L && kind != "after") { kind <- "after"; k <- 1L }
      tr <- switch(kind, after = te + runif(1, 0.3, 2), inside = st[k] + runif(1, 0.2, 0.8) * dur[k], at_end = en[k]); te <- tr
      if (kind == "inside") truth$reset_inside <- TRUE; if (kind == "at_end") truth$reset_at_end <- TRUE }
    # episode doses (the first of episode>=2 is carried by the reset row)
    keepd <- st < te | e == 1 & seq_along(st) == 1
    ep_rows <- D(TIME = st[keepd], EVID = 1L, AMT = amt[keepd], RATE = rate[keepd], kind = "dose")
    if (e > 1) { pure <- pure_prev
      if (pure) { ep_rows <- rbind(D(TIME = ts, EVID = 3L, AMT = 0, RATE = 0, kind = "reset"), ep_rows[ep_rows$TIME > ts, , drop = FALSE]) }
      else { ep_rows$EVID[1] <- 4L; ep_rows$kind[1] <- "reset" } }
    # covariate changes
    nch <- sample(1:4, 1); ct <- sort(runif(nch, ts + 0.05, te - 0.05))
    # cluster some changes shortly after the previous one to exercise short intervals
    if (nch > 1 && runif(1) < 0.7) ct[2] <- ct[1] + 10^runif(1, -3, -0.5)
    ct <- sort(ct[ct < te - 1e-6]); ct <- ct[!(round(ct, 9) %in% round(ep_rows$TIME, 9))]
    ctype <- if (covmode == "obs_only") rep("obs", length(ct)) else sample(c("obs", "evid2"), length(ct), TRUE)
    cs <- D(TIME = c(ts, ct), CRCA = runif(1 + length(ct), 20, 150), AJBW = runif(1 + length(ct), 40, 120))
    sched[[e]] <- cs
    chg_rows <- if (length(ct)) D(TIME = ct, EVID = ifelse(ctype == "evid2", 2L, 0L), AMT = 0, RATE = 0, kind = "chg") else NULL
    otimes <- c(runif(sample(5:9, 1), ts + 0.01, te - 1e-6)); otimes <- otimes[!(round(otimes, 9) %in% round(c(ep_rows$TIME, ct), 9))]
    obs_rows <- D(TIME = otimes, EVID = 0L, AMT = 0, RATE = 0, kind = "obs")
    carry_obs <- if (length(ct)) D(TIME = ct, EVID = 0L, AMT = 0, RATE = 0, kind = "obs") else NULL    # an obs row at every change time
    truth$obs_only_changes <- truth$obs_only_changes + sum(ctype == "obs")
    rows[[e]] <- rbind(ep_rows, chg_rows, obs_rows, carry_obs); rows[[e]]$ep <- e
    ts <- if (!last) tr else NA; pure_prev <- pure_next
  }
  ev <- do.call(rbind, rows); sc <- do.call(rbind, sched)
  idx <- vapply(ev$TIME, function(t) max(which(sc$TIME <= t)), 1L)          # LOCF covariate schedule (episode starts are schedule entries)
  ev$CRCA <- sc$CRCA[idx]; ev$AJBW <- sc$AJBW[idx]
  ev <- ev[order(ev$TIME, runif(nrow(ev))), ]                                # shuffle row order within same-time groups
  data.frame(ID = id, TIME = ev$TIME, EVID = ev$EVID, AMT = ev$AMT, RATE = ev$RATE, CMT = 1L, CRCA = ev$CRCA, AJBW = ev$AJBW,
             TAD = NA_real_, SRC_ROW = NA_integer_, row.names = NULL) -> out
  attr(out, "truth") <- truth; out
}
make_battery <- function(n, id0, allow_reset, covmode) { l <- lapply(seq_len(n), function(i) gen_patient(id0 + i, allow_reset, covmode))
  ev <- do.call(rbind, l); ev$SRC_ROW <- seq_len(nrow(ev)); tr <- do.call(rbind, lapply(l, function(x) as.data.frame(attr(x, "truth")))); tr$ID <- id0 + seq_len(n); list(ev = ev, truth = tr) }
A <- make_battery(150, 0L, FALSE, "obs_only")        # R-COV battery: obs-only changes, short intervals, no resets
B <- make_battery(250, 1000L, TRUE, "mixed")         # R-RESET battery: multi-episode, resets inside/at end of infusions
for (b in list(A, B)) stopifnot(all(b$ev$TIME >= 0))
eta_tab <- function(ev) { ids <- unique(ev$ID); data.frame(PATIENT_ID = ids, eCL = pmax(pmin(rnorm(length(ids), 0, .3), .6), -.6), eV = pmax(pmin(rnorm(length(ids), 0, .3), .6), -.6)) }
etaA <- eta_tab(A$ev); etaB <- eta_tab(B$ev)

# ---- reference on the ORIGINAL raw patient table (reference semantics only) -------------------------------
ref_all <- function(ev, eta, s) do.call(rbind, lapply(split(ev, ev$ID), function(e) { e <- e[order(e$TIME, seq_len(nrow(e))), ]
  et <- eta[eta$PATIENT_ID == e$ID[1], ]; x <- ref_matexp(e[, PKPD02_EV_COLS], s, PHYS[[s]], c(et$eCL, et$eV)); x$PATIENT_ID <- e$ID[1]; x }))

# ---- pipeline under test: split -> canonicalise -> validate -> solve ------------------------------------------
pipeline <- function(ev, eta, s, split = TRUE) {
  tl <- if (split) pkpd02_episode_timelines(ev, clock = "absolute")$timelines else { x <- ev; x$PATIENT_ID <- x$ID; x$EPISODE <- 1L; x$SOLVE_ID <- x$ID; x }
  cz <- pkpd02_canonicalize_events(tl); stopifnot(nrow(cz$rejected) == 0)
  e <- cz$events[, PKPD02_EV_COLS]; ic <- pkpd02_icov_by_patient(tl, eta); pkpd02_assert_shared_eta(ic)
  list(events = e, icov = ic, canon = cz, tl = tl)
}
run_cmp <- function(ev, eta, s, split) {
  p <- pipeline(ev, eta, s, split); stopifnot(length(pkpd02_validate_events(p$events)) == 0)
  o <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], p$events, p$icov[, c("ID", "eCL", "eV")], atol = 1e-12, rtol = 1e-12)
  o$PATIENT_ID <- if (split) o$ID %/% 100L else o$ID
  rf <- ref_all(ev, eta, s); m <- match(keyf(o$PATIENT_ID, o$time), keyf(rf$PATIENT_ID, rf$TIME)); stopifnot(!anyNA(m))
  d <- data.frame(model = s, PATIENT_ID = o$PATIENT_ID, time = o$time, evid = o$evid, C_ode = o$cp, C_ref = rf$cp[m], Ac_ode = o$Ac, Ac_ref = rf$Ac[m])
  d$abs_err <- abs(d$C_ode - d$C_ref); d$rel_err <- ifelse(d$C_ref == 0, ifelse(d$abs_err == 0, 0, Inf), d$abs_err / abs(d$C_ref)); d$tol <- pkpd02_tol(d$C_ref); d$pass <- d$abs_err <= d$tol
  list(cmp = d, p = p) }
unguarded <- function(ev, eta, s) { e <- ev[, PKPD02_EV_COLS]; ic <- data.frame(ID = eta$PATIENT_ID, eCL = eta$eCL, eV = eta$eV)
  o <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], e, ic, validate = FALSE); rf <- ref_all(ev, eta, s)
  m <- match(keyf(o$ID, o$time), keyf(rf$PATIENT_ID, rf$TIME)); er <- abs(o$cp - rf$cp[m]); data.frame(n = length(er), n_exceed = sum(er > pkpd02_tol(rf$cp[m])), max_abs = max(er), min_Ac = min(o$Ac)) }

res <- list(); gate <- list(); info <- list()
for (s in c("1c", "2c")) {
  ra <- run_cmp(A$ev, etaA, s, split = TRUE); rb <- run_cmp(B$ev, etaB, s, split = TRUE)
  res[[paste("A", s)]] <- cbind(battery = "A_obs_only_covariate_changes", ra$cmp); res[[paste("B", s)]] <- cbind(battery = "B_multi_episode_resets", rb$cmp)
  info[[paste("A_unguarded", s)]] <- cbind(battery = "A", model = s, unguarded(A$ev, etaA, s))
  # B unguarded single-timeline route: wrapper must reject those patients with reset inside/at end of infusion
  vB <- lapply(split(B$ev, B$ev$ID), function(e) length(pkpd02_validate_events(pkpd02_canonicalize_events(e)$events[, PKPD02_EV_COLS])) > 0)
  rejected <- as.integer(names(vB)[unlist(vB)]); expected <- B$truth$ID[B$truth$reset_inside | B$truth$reset_at_end]
  gate[[paste0("B_patient_level_rejections_match_independent_expectation_", s)]] <- setequal(rejected, expected)
  info[[paste("B_unguarded_on_rejected_patients", s)]] <- cbind(battery = "B", model = s, unguarded(B$ev[B$ev$ID %in% expected, ], etaB[etaB$PATIENT_ID %in% expected, ], s))
}
cmpall <- do.call(rbind, res); infotab <- do.call(rbind, info)
gate$A_all_rows_pass <- all(cmpall$pass[grepl("^A", cmpall$battery)]); gate$B_all_rows_pass <- all(cmpall$pass[grepl("^B", cmpall$battery)])
gate$A_has_obs_only_changes <- sum(A$truth$obs_only_changes) > 0
gate$B_contains_reset_inside_and_at_end <- any(B$truth$reset_inside) && any(B$truth$reset_at_end)

# ---- canonicalisation does not invent values/times; ordering; same-time permutations -------------------------------
cz <- pkpd02_canonicalize_events(A$ev); ins <- cz$inserted
src_has <- function(id, t, cr, aj) any(A$ev$ID == id & A$ev$TIME == t & A$ev$CRCA == cr & A$ev$AJBW == aj)
gate$inserted_rows_copy_existing_time_and_values <- nrow(ins) > 0 && all(vapply(seq_len(nrow(ins)), function(i) src_has(ins$ID[i], ins$TIME[i], ins$CRCA_to[i], ins$AJBW_to[i]), NA))
indep_count <- sum(vapply(split(A$ev, A$ev$ID), function(e) { e <- e[order(e$TIME), ]; g <- unique(e$TIME); prev <- NULL; n <- 0L
  for (t in g) { r <- e[e$TIME == t, ]; cur <- c(r$CRCA[1], r$AJBW[1]); if (!is.null(prev) && any(cur != prev) && all(r$EVID == 0L)) n <- n + 1L; prev <- cur }; n }, 0L))
gate$inserted_count_equals_independent_count <- nrow(ins) == indep_count
gate$covariate_row_precedes_samples_at_same_time <- all(vapply(split(cz$events, cz$events$ID), function(e) all(vapply(unique(e$TIME), function(t) { r <- e[e$TIME == t, ]; !any(r$EVID == 0L) || all(which(r$EVID == 2L) < min(which(r$EVID == 0L))) }, NA)), NA))
perm <- A$ev[A$ev$ID <= 20, ]; set.seed(SEED + 1L); perm <- perm[order(perm$TIME, runif(nrow(perm))), ]; perm$SRC_ROW <- seq_len(nrow(perm))
e1 <- pkpd02_canonicalize_events(A$ev[A$ev$ID <= 20, ])$events[, PKPD02_EV_COLS]; e2 <- pkpd02_canonicalize_events(perm)$events[, PKPD02_EV_COLS]
key <- function(e) paste(e$ID, sprintf("%.12f", e$TIME), e$EVID, e$AMT, e$CRCA, e$AJBW)
gate$canonical_form_invariant_to_same_time_row_order <- setequal(key(e1), key(e2))
conf <- A$ev[A$ev$ID == 1, ][1:2, ]; conf$TIME <- 7; conf$CRCA <- c(50, 60); conf <- rbind(A$ev[A$ev$ID == 1, ][1, , drop = FALSE], conf); conf$SRC_ROW <- seq_len(nrow(conf))
gate$conflicting_tied_covariates_rejected_not_resolved <- nrow(pkpd02_canonicalize_events(conf)$rejected) == 1

# ---- patient identity ------------------------------------------------------------------------------------------------------
pb <- pipeline(B$ev, etaB, "1c")
gate$episodes_of_a_patient_share_patient_eta <- all(tapply(seq_len(nrow(pb$icov)), pb$icov$PATIENT_ID, function(i) nrow(unique(pb$icov[i, c("eCL", "eV")])) == 1))
gate$multi_episode_patients_present <- any(table(pb$icov$PATIENT_ID) > 1)
bad <- pb$icov; mp <- names(which(table(bad$PATIENT_ID) > 1))[1]; bad$eCL[bad$PATIENT_ID == mp][1] <- bad$eCL[bad$PATIENT_ID == mp][1] + 0.1
gate$per_episode_etas_refused_without_explicit_decision <- inherits(try(pkpd02_assert_shared_eta(bad), silent = TRUE), "try-error") && isTRUE(pkpd02_assert_shared_eta(bad, allow_episode_eta = TRUE))
gate$duplicate_patient_eta_rows_refused <- inherits(try(pkpd02_icov_by_patient(pb$tl, rbind(etaB, etaB[1, ])), silent = TRUE), "try-error")

# ---- (4) reset-tied sample classification: hand-built cases, expected values from independent reference -------------------------------
mk <- function(rows) { x <- do.call(rbind, rows); x$ID <- 1L; x$CMT <- 1L; x$SRC_ROW <- seq_len(nrow(x)); x[, c("ID", "TIME", "EVID", "AMT", "RATE", "CMT", "CRCA", "AJBW", "TAD", "SRC_ROW")] }
rw <- function(t, evid = 0L, amt = 0, rate = 0, crca = 90, aj = 75, tad = NA_real_) D(TIME = t, EVID = evid, AMT = amt, RATE = rate, CRCA = crca, AJBW = aj, TAD = tad)
ep1 <- list(rw(0, 1L, 600, 300, tad = 0), rw(1, tad = 1), rw(3, tad = 3))
tie_cases <- list(
  pre_before_reset_TAD_pos      = list(rows = c(ep1, list(rw(5, tad = 5), rw(5, 4L, 400, 400, 60, 90, 0), rw(6, crca = 60, aj = 90, tad = 1))), reset_evid = 4L, expect = "pre_reset"),
  post_after_reset_TAD_zero     = list(rows = c(ep1, list(rw(5, 4L, 400, 400, 60, 90, 0), rw(5, crca = 60, aj = 90, tad = 0), rw(6, crca = 60, aj = 90, tad = 1))), reset_evid = 4L, expect = "post_reset"),
  indet_before_reset_TAD_zero   = list(rows = c(ep1, list(rw(5, tad = 0), rw(5, 4L, 400, 400, 60, 90, 0), rw(6, crca = 60, aj = 90, tad = 1))), reset_evid = 4L, expect = "indeterminate"),
  indet_after_reset_TAD_pos     = list(rows = c(ep1, list(rw(5, 4L, 400, 400, 60, 90, 0), rw(5, crca = 60, aj = 90, tad = 3), rw(6, crca = 60, aj = 90, tad = 1))), reset_evid = 4L, expect = "indeterminate"),
  indet_TAD_missing             = list(rows = c(ep1, list(rw(5, tad = NA), rw(5, 4L, 400, 400, 60, 90, 0), rw(6, crca = 60, aj = 90, tad = 1))), reset_evid = 4L, expect = "indeterminate"),
  pre_before_EVID3_TAD_pos      = list(rows = c(ep1, list(rw(5, tad = 5), rw(5, 3L, 0, 0, 60, 90, NA), rw(5.5, 1L, 400, 400, 60, 90, 0), rw(6, crca = 60, aj = 90, tad = .5))), reset_evid = 3L, expect = "pre_reset"),
  post_after_EVID3_TAD_zero     = list(rows = c(ep1, list(rw(5, 3L, 0, 0, 60, 90, NA), rw(5, crca = 60, aj = 90, tad = 0), rw(5.5, 1L, 400, 400, 60, 90, 0), rw(6, crca = 60, aj = 90, tad = .5))), reset_evid = 3L, expect = "post_reset"))
trows <- list(); okc <- logical()
for (nm in names(tie_cases)) { tc <- tie_cases[[nm]]; ev <- mk(tc$rows); ties <- pkpd02_classify_reset_ties(ev, "absolute"); tl <- pkpd02_episode_timelines(ev, ties, "absolute")
  cls_ok <- nrow(ties) == 1 && ties$classification == tc$expect
  tied_row <- ev$SRC_ROW[ev$EVID == 0L & ev$TIME == 5]; kept <- tied_row %in% tl$timelines$SRC_ROW; quar <- tied_row %in% tl$quarantined$SRC_ROW
  ep_of <- if (kept) tl$timelines$EPISODE[tl$timelines$SRC_ROW == tied_row] else NA
  place_ok <- switch(tc$expect, pre_reset = kept && ep_of == 1L, post_reset = kept && ep_of == 2L, indeterminate = quar && !kept)
  val_ok <- TRUE; err <- NA_real_
  if (kept) { # expected value from an independent reference built from the relevant episode alone
    refrows <- if (tc$expect == "pre_reset") ev[ev$TIME < 5 | (ev$TIME == 5 & ev$EVID == 0L), ] else { e2 <- ev[ev$TIME >= 5 & !(ev$TIME == 5 & ev$EVID == 0L & ev$SRC_ROW != tied_row), ]; e2 }
    if (tc$expect == "post_reset") { e2 <- ev[ev$TIME >= 5, ]; e2$EVID[e2$EVID == 3L] <- 3L; refrows <- e2 }
    rf <- ref_matexp(refrows[, PKPD02_EV_COLS], "1c", PHYS[["1c"]], c(0, 0)); expv <- rf$cp[rf$TIME == 5]
    cz <- pkpd02_canonicalize_events(tl$timelines); stopifnot(nrow(cz$rejected) == 0)
    sid <- tl$timelines$SOLVE_ID[tl$timelines$SRC_ROW == tied_row]
    o <- pkpd02_solve_ode(plain[["1c"]], "1c", PHYS[["1c"]], cz$events[, PKPD02_EV_COLS], data.frame(ID = unique(cz$events$ID), eCL = 0, eV = 0), atol = 1e-12, rtol = 1e-12)
    got <- o$cp[o$ID == sid & abs(o$time - 5) < 1e-12 & o$evid == 0][1]; err <- abs(got - expv); val_ok <- err <= pkpd02_tol(expv) }
  trows[[nm]] <- data.frame(case = nm, expected = tc$expect, classified = if (nrow(ties)) ties$classification else NA, classification_ok = cls_ok, placement_ok = place_ok, value_ok = val_ok, abs_err = err, quarantined = quar)
  okc <- c(okc, cls_ok && place_ok && val_ok) }
ties_tab <- do.call(rbind, trows); gate$reset_tie_rules_and_values <- all(okc)
code_lines <- sub("#.*$", "", readLines("R/pkpd02_events.R"))
gate$no_DV_reference_in_event_library_code <- !any(grepl("\\bDV\\b", code_lines))
# replay of the preserved S-1r R-COV hazard tables (H1: 0.001 h, H2: 0.1 h obs-only change that reverts): unguarded vs canonicalised
hzev <- read.csv("evidence/s1r_hazard_events.csv"); hzev <- hzev[hzev$ID %in% c(2001, 2002), ]; hzev$CMT <- as.integer(hzev$CMT); hzev$EVID <- as.integer(hzev$EVID)
hzev$TAD <- NA_real_; hzev$SRC_ROW <- seq_len(nrow(hzev)); hz_rows <- list()
for (s in c("1c", "2c")) { et <- data.frame(PATIENT_ID = c(2001, 2002), eCL = 0.10, eV = -0.05)
  rf <- ref_all(transform(hzev, ID = ID), et, s)
  raw <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], hzev[, PKPD02_EV_COLS], data.frame(ID = et$PATIENT_ID, eCL = 0.10, eV = -0.05), validate = FALSE)
  cz <- pkpd02_canonicalize_events(hzev); ce <- cz$events[, PKPD02_EV_COLS]; stopifnot(length(pkpd02_validate_events(ce)) == 0)
  can <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], ce, data.frame(ID = et$PATIENT_ID, eCL = 0.10, eV = -0.05), atol = 1e-12, rtol = 1e-12)
  for (nm in c("raw", "can")) { o <- get(nm); m <- match(keyf(o$ID, o$time), keyf(rf$PATIENT_ID, rf$TIME)); er <- abs(o$cp - rf$cp[m])
    hz_rows[[length(hz_rows) + 1L]] <- data.frame(model = s, route = if (nm == "raw") "unguarded_raw" else "canonicalised", n = length(er), n_exceed = sum(er > pkpd02_tol(rf$cp[m])), max_abs = max(er)) } }
hz_tab <- do.call(rbind, hz_rows)
gate$S1r_hazard_replay_unguarded_disagrees_in_1c <- any(hz_tab$n_exceed[hz_tab$route == "unguarded_raw" & hz_tab$model == "1c"] > 0)
gate$S1r_hazard_replay_canonicalised_all_pass <- all(hz_tab$n_exceed[hz_tab$route == "canonicalised"] == 0)

# ---- battery C: EPISODE-RELATIVE clock (as in the real data: every reset row has TIME 0) ----------------------------------------
gen_rel <- function(id) { n_ep <- sample(1:3, 1); parts <- lapply(seq_len(n_ep), function(e) { x <- gen_patient(1L, FALSE, "mixed"); x$ID <- id; x$ep <- e; x$EVID[x$EVID == 1L & x$TIME == 0] <- if (e > 1) 4L else 1L; x }); do.call(rbind, parts) }
C_ev <- do.call(rbind, lapply(1:200, function(i) gen_rel(3000L + i))); C_ev$TAD <- NA_real_; C_ev$SRC_ROW <- seq_len(nrow(C_ev)); etaC <- eta_tab(C_ev)
ref_episode <- function(ev, eta, s) do.call(rbind, lapply(split(ev, list(ev$ID, ev$ep), drop = TRUE), function(e) { et <- eta[eta$PATIENT_ID == e$ID[1], ]
  x <- ref_matexp(e[, PKPD02_EV_COLS], s, PHYS[[s]], c(et$eCL, et$eV)); x$PATIENT_ID <- e$ID[1]; x$EPISODE <- e$ep[1]; x }))
C_ev$ep <- NULL; C_src <- C_ev; C_src$ep <- pkpd02_episode_index(C_src)
crow <- list(); gap_inv <- list()
for (s in c("1c", "2c")) {
  rf <- ref_episode(C_src, etaC, s); kf <- function(p, e, t) paste(p, e, sprintf("%.9f", t))
  # (i) timelines
  tlC <- pkpd02_episode_timelines(C_src[, c(PKPD02_EV_COLS, "TAD", "SRC_ROW")], clock = "episode_relative"); stopifnot(nrow(tlC$quarantined) == 0)
  cz <- pkpd02_canonicalize_events(tlC$timelines); stopifnot(nrow(cz$rejected) == 0, length(pkpd02_validate_events(cz$events[, PKPD02_EV_COLS])) == 0)
  ic <- pkpd02_icov_by_patient(tlC$timelines, etaC); pkpd02_assert_shared_eta(ic)
  o <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], cz$events[, PKPD02_EV_COLS], ic[, c("ID", "eCL", "eV")], atol = 1e-12, rtol = 1e-12)
  m <- match(kf(o$ID %/% 100L, o$ID %% 100L, o$time), kf(rf$PATIENT_ID, rf$EPISODE, rf$TIME)); stopifnot(!anyNA(m))
  d <- data.frame(model = s, route = "episode_timelines", C_ode = o$cp, C_ref = rf$cp[m]); crow[[length(crow) + 1L]] <- d
  # (ii) patient-level offset table, several gaps
  base <- NULL
  for (gp in c(1e-3, 1, 100)) {
    pl <- pkpd02_patient_level_table(tlC$timelines, gap = gp); czp <- pkpd02_canonicalize_events(pl); stopifnot(nrow(czp$rejected) == 0)
    v <- pkpd02_validate_events(czp$events[, PKPD02_EV_COLS]); stopifnot(length(v) == 0)
    off <- attr(pl, "offsets"); icp <- data.frame(ID = etaC$PATIENT_ID, eCL = etaC$eCL, eV = etaC$eV)
    op <- pkpd02_solve_ode(plain[[s]], s, PHYS[[s]], czp$events[, PKPD02_EV_COLS], icp, atol = 1e-12, rtol = 1e-12)
    epi <- vapply(seq_len(nrow(op)), function(i) { of <- off[off$PATIENT_ID == op$ID[i], ]; max(of$EPISODE[of$offset <= op$time[i] + 1e-12]) }, 1L)
    rel <- op$time - vapply(seq_len(nrow(op)), function(i) off$offset[off$PATIENT_ID == op$ID[i] & off$EPISODE == epi[i]], 0)
    m2 <- match(kf(op$ID, epi, rel), kf(rf$PATIENT_ID, rf$EPISODE, rf$TIME)); stopifnot(!anyNA(m2))
    crow[[length(crow) + 1L]] <- data.frame(model = s, route = paste0("patient_level_offset_gap_", gp), C_ode = op$cp, C_ref = rf$cp[m2])
    if (is.null(base)) base <- op$cp else gap_inv[[length(gap_inv) + 1L]] <- data.frame(model = s, gap = gp, max_abs_diff_vs_smallest_gap = max(abs(op$cp - base[seq_along(op$cp)]))) } }
Cc <- do.call(rbind, crow); Cc$abs_err <- abs(Cc$C_ode - Cc$C_ref); Cc$tol <- pkpd02_tol(Cc$C_ref); Cc$pass <- Cc$abs_err <= Cc$tol
gap_tab <- do.call(rbind, gap_inv)
gate$C_episode_relative_timelines_pass <- all(Cc$pass[Cc$route == "episode_timelines"])
gate$C_patient_level_offset_table_pass_all_gaps <- all(Cc$pass[grepl("patient_level", Cc$route)])
gate$C_results_invariant_to_arbitrary_offset_gap <- all(gap_tab$max_abs_diff_vs_smallest_gap <= 1e-9)
gate$C_multi_episode_patients_present <- length(unique(paste(C_src$ID[C_src$ep > 1]))) > 20
# relative-clock reset ties: obs at TIME 0 AFTER the reset row; TAD==0 -> post_reset (kept), TAD>0 or NA -> indeterminate (quarantined)
relcase <- function(tad) { ev <- mk(c(ep1, list(rw(0, 4L, 400, 400, 60, 90, 0), rw(0, crca = 60, aj = 90, tad = tad), rw(1, crca = 60, aj = 90, tad = 1)))); ev$TIME[4:6] <- c(0, 0, 1)   # episode 2 restarts at TIME 0
  ties <- pkpd02_classify_reset_ties(ev, "episode_relative"); tl <- pkpd02_episode_timelines(ev, ties, "episode_relative"); list(ties = ties, kept = 5L %in% tl$timelines$SRC_ROW, quar = 5L %in% tl$quarantined$SRC_ROW) }
r0 <- relcase(0); r1 <- relcase(1.2); rN <- relcase(NA)
gate$relative_clock_tie_rules <- identical(r0$ties$classification, "post_reset") && r0$kept && identical(r1$ties$classification, "indeterminate") && r1$quar && !r1$kept && identical(rN$ties$classification, "indeterminate") && rN$quar
gate$relative_clock_no_false_ties_across_episodes <- { ev <- mk(c(list(rw(0, 1L, 600, 300, tad = 0), rw(0, tad = 24)), list(rw(5, tad = 5), rw(0, 4L, 400, 400, 60, 90, 0), rw(1, crca = 60, aj = 90, tad = 1)))); ev$TIME[3:5] <- c(5, 0, 1); nrow(pkpd02_classify_reset_ties(ev, "episode_relative")) == 0 }

# ---- verdict & evidence ----------------------------------------------------------------------------------------
VERDICT <- if (all(unlist(gate))) "PASS" else "FAIL"
dir.create("evidence", FALSE); w <- function(x, f) write.csv(x, file.path("evidence", f), row.names = FALSE)
w(cmpall, "s2prep_canonical_comparison_points.csv"); w(cmpall[!cmpall$pass, ], "s2prep_canonical_failing_points.csv"); w(ties_tab, "s2prep_reset_tie_cases.csv"); w(infotab, "s2prep_unguarded_hazard_demonstration.csv")
w(ins, "s2prep_canonicalisation_inserted_rows_battery_A.csv"); w(Cc, "s2prep_battery_C_episode_relative_points.csv"); w(gap_tab, "s2prep_battery_C_gap_invariance.csv"); w(hz_tab, "s2prep_s1r_hazard_replay.csv"); w(B$truth, "s2prep_battery_B_truth.csv")
summ <- do.call(rbind, lapply(split(cmpall, list(cmpall$battery, cmpall$model)), function(d) data.frame(battery = d$battery[1], model = d$model[1], n_rows = nrow(d), n_patients = length(unique(d$PATIENT_ID)), n_fail = sum(!d$pass), max_abs_err = max(d$abs_err), max_rel_err = max(d$rel_err[is.finite(d$rel_err)]), min_Ac_ode = min(d$Ac_ode))))
rownames(summ) <- NULL; w(summ, "s2prep_canonical_summary.csv")
out <- list(test = "s2prep_canonical_equivalence", verdict = VERDICT, synthetic_only = TRUE, seed = SEED, tolerance = pkpd02_tol(1), ode = list(atol = 1e-12, rtol = 1e-12, method = "liblsoda"),
            gates = gate, summary = summ, battery_A = list(n_patients = length(unique(A$ev$ID)), obs_only_changes = sum(A$truth$obs_only_changes), inserted_rows = nrow(ins)),
            battery_B = list(n_patients = length(unique(B$ev$ID)), reset_inside = sum(B$truth$reset_inside), reset_at_end = sum(B$truth$reset_at_end), multi_episode_patients = sum(table(pb$icov$PATIENT_ID) > 1)),
            unguarded_demonstration = infotab, reset_tie_cases = ties_tab)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "string"), "evidence/s2prep_canonical_equivalence.json")
print(summ, row.names = FALSE); cat("\ngates:\n"); print(unlist(gate)); cat("\nunguarded hazard demonstration:\n"); print(infotab, row.names = FALSE); print(ties_tab, row.names = FALSE)
cat("\nCANONICAL EQUIVALENCE VERDICT:", VERDICT, "\n"); quit(status = if (VERDICT == "PASS") 0L else 1L)
