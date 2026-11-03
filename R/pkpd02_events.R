# PKPD-02 v1.2 (DRAFT amendment rev2): canonical event representation.
# METADATA ONLY: no function here reads or needs a measured concentration (DV). Observation rows are
# EVID=0 rows; whether a DV is numeric/positive is not consulted (see audit limitations).
#
# Input ("source-like" table): ID, TIME, EVID, AMT, RATE, CMT, CRCA, AJBW, optional TAD, SRC_ROW
# (original file order). Rules implemented (all by source metadata only):
#  R-COV   canonicalisation: a covariate change carried only by observation rows is RE-CARRIED by an
#          EVID=2 row at the SAME time with the SAME values, placed first in the time group. No value
#          or time is invented; tied rows with conflicting covariates are rejected, never resolved.
#  R-RESET episode-specific timelines: each episode is a separate timeline starting from zero state,
#          so no earlier infusion-stop event can leak across a reset. Patient identity is kept in
#          PATIENT_ID/EPISODE; SOLVE_ID is only the solver key.
#  Tie rule (reset-tied samples): pre-reset iff (file order before the reset row AND TAD>0); post-reset
#          iff (after the reset row AND TAD==0); otherwise INDETERMINATE -> quarantined.

PKPD02_EV_COLS <- c("ID", "TIME", "EVID", "AMT", "RATE", "CMT", "CRCA", "AJBW")

# Effective infusion duration (v1.1 T7/T8). Returns data.frame(dur, rule) per row (NA for non-dose rows).
pkpd02_effective_duration <- function(amt, rate, evid) {
  dose <- evid %in% c(1L, 4L) & !is.na(amt) & amt > 0
  dur <- rep(NA_real_, length(amt)); rule <- rep(NA_character_, length(amt))
  rec <- ifelse(dose & !is.na(rate) & rate > 0, amt / rate, NA_real_)
  zero <- dose & (is.na(rate) | rate == 0)
  bad <- dose & !zero & (rate < 0 | !is.finite(rec) | rec > 24)
  dur[dose] <- rec[dose]; rule[dose] <- "as_recorded"
  dur[zero] <- 0.5; rule[zero] <- "T7_RATE0_half_hour"
  dur[bad] <- 0.5; rule[bad] <- "T8_fallback_half_hour"
  data.frame(dur = dur, rule = rule)
}

pkpd02_episode_index <- function(ev) {            # per-ID episode number; EVID 3/4 start a new episode
  epi <- integer(nrow(ev))
  for (id in unique(ev$ID)) { i <- which(ev$ID == id); st <- ev$EVID[i] %in% c(3L, 4L); st[1] <- TRUE; epi[i] <- cumsum(st) }
  epi
}

# ---- R-COV canonicalisation -------------------------------------------------------------------
pkpd02_canonicalize_events <- function(ev) {
  stopifnot(all(PKPD02_EV_COLS %in% names(ev)))
  if (!"SRC_ROW" %in% names(ev)) ev$SRC_ROW <- seq_len(nrow(ev))
  out <- list(); log <- list(); rej <- list()
  for (id in unique(ev$ID)) {
    x <- ev[ev$ID == id, , drop = FALSE]; prev <- NULL; keep <- vector("list", length(unique(x$TIME))); k <- 0L
    for (tt in unique(x$TIME)) {
      g <- x[x$TIME == tt, , drop = FALSE]; k <- k + 1L
      cv <- g[, c("CRCA", "AJBW")]
      if (anyNA(cv)) { rej[[length(rej) + 1L]] <- data.frame(ID = id, TIME = tt, reason = "NA covariate"); keep[[k]] <- g; next }
      if (nrow(unique(cv)) > 1) { rej[[length(rej) + 1L]] <- data.frame(ID = id, TIME = tt, reason = "tied rows with conflicting covariates (effective value undefined)"); keep[[k]] <- g; next }
      cur <- unlist(cv[1, ])
      if (!is.null(prev) && any(cur != prev) && all(g$EVID == 0L)) {
        add <- g[1, , drop = FALSE]; add$EVID <- 2L; add$AMT <- 0; add$RATE <- 0; add$SRC_ROW <- NA_integer_
        if ("TAD" %in% names(add)) add$TAD <- NA_real_
        log[[length(log) + 1L]] <- data.frame(ID = id, TIME = tt, CRCA_from = prev[1], CRCA_to = cur[1], AJBW_from = prev[2], AJBW_to = cur[2], carried_by_src_row = g$SRC_ROW[1])
        g <- rbind(add, g)                          # covariate row FIRST within the time group
      }
      keep[[k]] <- g; prev <- cur
    }
    out[[length(out) + 1L]] <- do.call(rbind, keep)
  }
  list(events = do.call(rbind, out), inserted = if (length(log)) do.call(rbind, log) else data.frame(),
       rejected = if (length(rej)) do.call(rbind, rej) else data.frame())
}

# ---- reset-tied sample classification (metadata only) ----------------------------------------
# clock = "episode_relative": TIME restarts at each reset (as in the Strathclyde source: every reset row has
#   TIME 0). Episodes are assigned by FILE ORDER, so a tie is an observation of the NEW episode whose TIME equals
#   the reset row's TIME; it follows the reset row in file order and is "post_reset" iff TAD == 0, else
#   indeterminate (TAD NA also indeterminate). Rows before the reset row belong to the earlier episode on a
#   different clock and are never tied.
# clock = "absolute": one continuous clock; an observation at the reset time is pre-reset iff (file order before
#   the reset row AND TAD>0), post-reset iff (after it AND TAD==0), otherwise indeterminate.
pkpd02_classify_reset_ties <- function(ev, clock = c("episode_relative", "absolute")) {
  clock <- match.arg(clock); stopifnot(all(c("SRC_ROW", "TAD") %in% names(ev)))
  epi <- pkpd02_episode_index(ev); res <- list()
  for (id in unique(ev$ID)) { i_id <- which(ev$ID == id); x <- ev[i_id, , drop = FALSE]; xe <- epi[i_id]
    for (j in which(x$EVID %in% c(3L, 4L))) {
      tied <- if (clock == "absolute") which(x$EVID == 0L & x$TIME == x$TIME[j]) else which(x$EVID == 0L & x$TIME == x$TIME[j] & xe == xe[j])
      for (i in tied) {
        before <- x$SRC_ROW[i] < x$SRC_ROW[j]; tad <- x$TAD[i]
        cls <- if (is.na(tad)) "indeterminate" else if (before && tad > 0 && clock == "absolute") "pre_reset" else if (!before && tad == 0) "post_reset" else "indeterminate"
        res[[length(res) + 1L]] <- data.frame(ID = id, TIME = x$TIME[j], obs_src_row = x$SRC_ROW[i], reset_src_row = x$SRC_ROW[j], file_order = if (before) "before_reset_row" else "after_reset_row",
                                              TAD = tad, classification = cls) } } }
  if (length(res)) do.call(rbind, res) else data.frame()
}

# ---- episode-specific timelines ---------------------------------------------------------------
# Returns list(timelines, quarantined, episode_info). timelines carries PATIENT_ID, EPISODE, SOLVE_ID.
pkpd02_episode_timelines <- function(ev, ties = NULL, clock = c("episode_relative", "absolute")) {
  clock <- match.arg(clock)
  if (!"SRC_ROW" %in% names(ev)) ev$SRC_ROW <- seq_len(nrow(ev))
  if (is.null(ties)) ties <- if ("TAD" %in% names(ev)) pkpd02_classify_reset_ties(ev, clock) else data.frame()
  ev$PATIENT_ID <- ev$ID; ev$EPISODE <- pkpd02_episode_index(ev)
  stopifnot(all(ev$EPISODE < 100))
  # File order already places a pre-reset sample (before the reset row) in the earlier episode and a
  # post-reset sample (after it) in the new one; the classification only decides WHETHER the tie is
  # consistent enough to keep. Indeterminate ties are quarantined (excluded from fitting/forecasting).
  drop <- rep(FALSE, nrow(ev))
  if (nrow(ties)) for (r in seq_len(nrow(ties))) { i <- which(ev$ID == ties$ID[r] & ev$SRC_ROW == ties$obs_src_row[r])
    if (ties$classification[r] == "indeterminate") drop[i] <- TRUE }
  quarantined <- ev[drop, , drop = FALSE]; ev <- ev[!drop, , drop = FALSE]
  ev$SOLVE_ID <- ev$PATIENT_ID * 100L + ev$EPISODE
  tl <- list(); info <- list()
  for (sid in unique(ev$SOLVE_ID)) {
    x <- ev[ev$SOLVE_ID == sid, , drop = FALSE]
    x <- x[order(x$TIME, x$EVID == 0L, x$SRC_ROW), ]                    # same-time order: covariate/dose rows, then samples
    # reset rows start the episode: EVID=4 -> dose (EVID=1) ; EVID=3 -> covariate row (EVID=2)
    x$EVID[x$EVID == 4L] <- 1L; x$EVID[x$EVID == 3L] <- 2L
    x$AMT[x$EVID != 1L] <- 0; x$RATE[x$EVID != 1L] <- 0
    tl[[length(tl) + 1L]] <- x
  }
  tl <- do.call(rbind, tl); tl$ID_PATIENT_ORIG <- tl$PATIENT_ID; tl$ID <- tl$SOLVE_ID
  list(timelines = tl, quarantined = quarantined)
}

# ---- explicit patient-level eta mapping (no implicit per-episode etas) --------------------------
# eta_tab: data.frame(PATIENT_ID, eCL, eV) one row per patient. Returns icov keyed by SOLVE_ID.
pkpd02_icov_by_patient <- function(timelines, eta_tab, allow_episode_eta = FALSE) {
  u <- unique(timelines[, c("SOLVE_ID", "PATIENT_ID")])
  if (anyDuplicated(eta_tab$PATIENT_ID)) stop("eta table must have one row per PATIENT_ID (no per-episode etas without an explicit model decision)")
  if (!allow_episode_eta && !all(u$PATIENT_ID %in% eta_tab$PATIENT_ID)) stop("missing patient eta")
  m <- match(u$PATIENT_ID, eta_tab$PATIENT_ID)
  data.frame(ID = u$SOLVE_ID, PATIENT_ID = u$PATIENT_ID, eCL = eta_tab$eCL[m], eV = eta_tab$eV[m])
}
pkpd02_assert_shared_eta <- function(icov, allow_episode_eta = FALSE) {
  if (allow_episode_eta) return(invisible(TRUE))
  bad <- tapply(seq_len(nrow(icov)), icov$PATIENT_ID, function(i) nrow(unique(icov[i, c("eCL", "eV")])) > 1)
  if (any(bad)) stop("episodes of one patient carry different etas: independent per-episode etas need an explicit model decision (allow_episode_eta=TRUE)")
  invisible(TRUE)
}

# ---- patient-level single-ID table (needed so one patient keeps ONE set of etas in nlmixr2) ---------------------
# Input: episode timelines (episode-relative clocks, doses with RATE>0, resets already converted to dose/EVID=2 rows).
# Each episode k>1 is shifted by offset_k = (latest of: last record, last nominal infusion end in the previous
# episode, on the shifted clock) + gap. The offset has NO clinical meaning and the wall-clock gap between real
# episodes is unknown; it only makes TIME monotone and keeps every earlier infusion finished before the reset
# (so the R-RESET hazard cannot occur). Because states are zeroed at a reset, results must be invariant to `gap`
# (tested). The reset row is re-instated as EVID=4 (dose) or EVID=3 (pure reset). Columns EPISODE/PATIENT_ID are kept.
pkpd02_patient_level_table <- function(timelines, gap = 1, reset_evid = NULL) {
  stopifnot(gap > 0); out <- list(); offs <- list()
  for (pid in unique(timelines$PATIENT_ID)) { x <- timelines[timelines$PATIENT_ID == pid, , drop = FALSE]; prev_end <- NA_real_; off <- 0
    for (e in sort(unique(x$EPISODE))) { y <- x[x$EPISODE == e, , drop = FALSE]
      if (e > 1) off <- prev_end + gap
      y$TIME <- y$TIME + off
      if (e > 1) { first <- which(y$TIME == min(y$TIME) & y$EVID %in% c(1L, 2L))[1]
        y$EVID[first] <- if (y$EVID[first] == 1L) 4L else 3L }               # restore the reset marker on the episode's first row
      ends <- ifelse(y$EVID %in% c(1L, 4L) & y$AMT > 0, y$TIME + y$AMT / pmax(y$RATE, 1e-300), -Inf)
      prev_end <- max(y$TIME, ends); out[[length(out) + 1L]] <- y; offs[[length(offs) + 1L]] <- data.frame(PATIENT_ID = pid, EPISODE = e, offset = off) } }
  r <- do.call(rbind, out); r$ID <- r$PATIENT_ID; attr(r, "offsets") <- do.call(rbind, offs); r
}
