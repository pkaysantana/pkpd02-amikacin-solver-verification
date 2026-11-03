# PKPD-02 v1.2 (DRAFT amendment): guarded solver wrapper for the explicit ODE pathway.
# Every call (a) validates the event table, (b) forces LOCF + useLinCmt=FALSE (which also
# overrides options(rxode2.useLinCmt=TRUE)), (c) asserts AFTER solving that the executed
# model contains no linCmt(), has the expected amount states, and ran with the requested
# controls. Any violation is an error (fail closed).

PKPD02_COV <- c("CRCA", "AJBW")
PKPD02_STATES <- list("1c" = "Ac", "2c" = c("Ac", "Ap"))
PKPD02_TOL <- list(abs = 1e-5, rel = 1e-6)   # mg/L; acceptance rule (v1.1 s3)
pkpd02_tol <- function(c_ref) pmax(PKPD02_TOL$abs, PKPD02_TOL$rel * abs(c_ref))

# Event-table contract. Covariate semantics: the value on a row applies from that row's
# TIME forward (LOCF); therefore every covariate change MUST be carried by a row at its
# timestamp, and rows sharing a TIME must carry identical covariates.
pkpd02_validate_events <- function(ev) {
  e <- character()
  need <- c("ID", "TIME", "EVID", "AMT", "RATE", "CMT", PKPD02_COV)
  if (length(miss <- setdiff(need, names(ev)))) return(paste("missing columns:", paste(miss, collapse = ",")))
  if (anyNA(ev[, need])) e <- c(e, "NA in required column")
  if (!all(ev$EVID %in% c(0L, 1L, 2L, 3L, 4L))) e <- c(e, "EVID outside 0-4")
  if (any(ev$CMT != 1)) e <- c(e, "CMT must be 1 (central) for every row")
  dose <- ev$EVID %in% c(1L, 4L) & ev$AMT > 0
  if (any(dose & !(ev$RATE > 0))) e <- c(e, "dose with AMT>0 needs RATE>0 (zero-order infusion; RATE=0 handled by T7 preprocessing, not here)")
  if (any(ev$AMT < 0 | ev$RATE < 0)) e <- c(e, "negative AMT/RATE")
  if (any(!is.finite(ev$AMT / pmax(ev$RATE, 1e-300))[dose])) e <- c(e, "non-finite infusion duration")
  if (any(!ev$EVID %in% c(1L, 4L) & ev$AMT != 0)) e <- c(e, "AMT on a non-dose row")
  for (id in unique(ev$ID)) {
    x <- ev[ev$ID == id, ]
    if (is.unsorted(x$TIME)) e <- c(e, sprintf("ID %s: TIME not non-decreasing", id))
    for (tt in unique(x$TIME[duplicated(x$TIME)])) {
      y <- x[x$TIME == tt, PKPD02_COV]
      if (nrow(unique(y)) > 1) e <- c(e, sprintf("ID %s TIME %s: tied rows carry different covariates", id, format(tt)))
    }
    # R-COV: a covariate change must be carried by a solver-restart row (EVID 1-4), not only by
    # observation rows. Observed in S-1r: a change carried only on an EVID=0 row followed by a short
    # interval was applied with the END-of-interval value (LSODA step overshoot); EVID=2 rows were exact.
    prev <- NULL
    for (tt in unique(x$TIME)) { g <- x[x$TIME == tt, ]; cv <- unlist(g[1, PKPD02_COV])
      if (!is.null(prev) && any(cv != prev) && all(g$EVID == 0L))
        e <- c(e, sprintf("ID %s TIME %s: covariate change carried only by observation rows (R-COV)", id, format(tt)))
      prev <- cv }
    # R-RESET: no reset row (EVID 3/4) while an infusion is active or ends exactly at the reset time.
    # Observed in S-1r: rxode2 still applies the cancelled infusion's stop event after the reset
    # (phantom negative input, negative amounts). Protocol semantics for this case are NOT decided here.
    resets <- which(x$EVID %in% c(3L, 4L))
    for (j in resets) { p0 <- max(c(0L, resets[resets < j]))
      cand <- intersect(seq(max(p0, 1L), j - 1L), which(x$EVID %in% c(1L, 4L) & x$AMT > 0))
      cand <- cand[x$TIME[cand] < x$TIME[j]]
      if (length(cand) && any(x$TIME[cand] + x$AMT[cand] / pmax(x$RATE[cand], 1e-300) >= x$TIME[j] - 1e-12))
        e <- c(e, sprintf("ID %s TIME %s: reset while an infusion is active or ending at the reset time (R-RESET)", id, format(x$TIME[j]))) }
  }
  e
}

pkpd02_atol_rtol_production_candidate <- c(atol = 1e-10, rtol = 1e-10)  # PROPOSAL ONLY; see report s1r

# model: rxUi (function-style) or rxode2 object. params: named numeric; icov: per-ID data.frame
pkpd02_solve_ode <- function(model, structure = c("1c", "2c"), params, events, icov = NULL,
                             atol = 1e-12, rtol = 1e-12, method = "liblsoda",
                             covsInterpolation = "locf", validate = TRUE, ...) {
  structure <- match.arg(structure)
  if (is.null(model)) stop("model is NULL (e.g. ui$foceiModel requires library(nlmixr2est) to be attached); refusing to solve")
  stopifnot(identical(covsInterpolation, "locf"))   # LOCF is the only permitted setting
  if (validate && length(err <- pkpd02_validate_events(events)))
    stop("event table rejected: ", paste(err, collapse = "; "))
  # useLinCmt is an rxSolve.rxUi argument; for plain rxode2 models there is no conversion path
  extra <- if (inherits(model, "rxUi")) list(useLinCmt = FALSE) else list()
  res <- do.call(rxode2::rxSolve, c(list(model, params, events, iCov = icov, covsInterpolation = "locf",
                                         atol = atol, rtol = rtol, method = method, addDosing = TRUE,
                                         cores = 1), extra, list(...)))
  mv <- rxode2::rxModelVars(res); a <- res$env$.args
  chk <- c(
    no_linCmt_flag = isTRUE(unname(mv$flags["linCmt"]) == -100),
    states_are_explicit_amounts = identical(as.character(mv$state), PKPD02_STATES[[structure]]),
    useLinCmt_recorded_FALSE_if_ui = !inherits(model, "rxUi") || identical(a$useLinCmt, FALSE),
    locf_recorded = identical(as.integer(a$covsInterpolation), 1L),
    atol_recorded = isTRUE(all.equal(a$atol, atol)), rtol_recorded = isTRUE(all.equal(a$rtol, rtol)),
    method_recorded = isTRUE(unname(a$method) == unname(rxode2::rxControl(method = method)$method)))
  if (!all(chk)) stop("pathway assertion failed: ", paste(names(chk)[!chk], collapse = ", "))
  out <- as.data.frame(res); names(out) <- sub("^id$", "ID", names(out))
  attr(out, "pathway_checks") <- chk
  out
}
