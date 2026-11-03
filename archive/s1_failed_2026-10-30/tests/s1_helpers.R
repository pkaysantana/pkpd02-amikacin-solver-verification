# PKPD-02 Gate S-1 helpers: synthetic event construction and an INDEPENDENT
# base-R reference solver (piecewise-exact matrix exponential). Nothing here
# uses rxode2/linCmt(); it exists to adjudicate which side is wrong if
# linCmt() and the rxode2 ODE disagree.

S1_TOL_ABS <- 1e-5   # mg/L   (addendum v1.1 s3)
S1_TOL_REL <- 1e-6   # of |C_ODE|

s1_tol <- function(c_ode) pmax(S1_TOL_ABS, S1_TOL_REL * abs(c_ode))

# ---- synthetic population parameters (ARBITRARY; not fitted, not Strathclyde-derived)
S1_PARS <- list(
  "1c" = c(tCL = 4.0, tV = 18.0),
  "2c" = c(tCL = 4.0, tV1 = 14.0, tQ = 3.0, tV2 = 9.0)
)

# ---- event-table construction -------------------------------------------------
# rows : data.frame(time, evid, amt, rate) in desired tie order
# sched: data.frame(time, CRCA, AJBW); covariate value for a row = last sched
#        entry with time <= row time (LOCF; rows at the same time share values).
# EVID=2 rows are auto-added at sched times with no dose/reset row.
s1_events <- function(id, rows, sched) {
  sched <- sched[order(sched$time), ]
  have <- rows$time[rows$evid %in% c(1, 4)]
  need <- setdiff(sched$time, have)
  if (length(need)) {
    rows <- rbind(rows, data.frame(time = need, evid = 2L, amt = 0, rate = 0))
  }
  rows <- rows[order(rows$time, seq_len(nrow(rows))), ]
  idx <- vapply(rows$time, function(t) max(which(sched$time <= t)), 1L)
  data.frame(ID = id, TIME = rows$time, EVID = as.integer(rows$evid),
             AMT = rows$amt, RATE = rows$rate, CMT = 1L,
             CRCA = sched$CRCA[idx], AJBW = sched$AJBW[idx],
             row.names = NULL)
}

# ---- independent reference ---------------------------------------------------
expm_taylor <- function(A) {
  nrm <- max(rowSums(abs(A)))
  s <- max(0L, ceiling(log2(max(nrm, 1e-300))) + 1L)
  A <- A / 2^s
  E <- diag(nrow(A)); Tm <- diag(nrow(A))
  for (k in 1:26) { Tm <- Tm %*% A / k; E <- E + Tm }
  for (i in seq_len(s)) E <- E %*% E
  E
}

ref_coef <- function(model, P, eCL, eV, CRCA, AJBW) {
  CL <- P[["tCL"]] * CRCA / 100 * exp(eCL)
  if (model == "1c") {
    list(CL = CL, V1 = P[["tV"]] * AJBW / 70 * exp(eV), Q = 0, V2 = 1)
  } else {
    list(CL = CL, V1 = P[["tV1"]] * AJBW / 70 * exp(eV), Q = P[["tQ"]],
         V2 = P[["tV2"]] * AJBW / 70)
  }
}

# Advance state x=(Ac,Ap) over dt with constant coefficients and input rate r.
ref_advance <- function(x, cf, r, dt) {
  if (dt <= 0) return(x)
  k10 <- cf$CL / cf$V1; k12 <- cf$Q / cf$V1; k21 <- cf$Q / cf$V2
  # augmented system on (Ac, Ap, 1): dx/dt = [[-(k10+k12), k21],[k12,-k21]] x + (r,0)
  M <- rbind(c(-(k10 + k12), k21, r), c(k12, -k21, 0), c(0, 0, 0))
  as.numeric((expm_taylor(M * dt) %*% c(x, 1))[1:2])
}

# Returns one row per unique event time (all rows at that time share the value).
# Ac_pre/Ap_pre are the left-limit amounts BEFORE row changes at that time.
# interval_cov = "locf": an interval uses coefficients from the record that STARTS it (gating semantics)
#              = "nocb": an interval uses coefficients from the record that ENDS it (diagnostic only)
ref_solve <- function(ev, model, P, eCL, eV, interval_cov = "locf") {
  ev <- ev[order(ev$TIME, seq_len(nrow(ev))), ]
  times <- sort(unique(ev$TIME))
  x <- c(0, 0); tcur <- times[1]
  inf_end <- numeric(0); inf_rate <- numeric(0)
  cf <- NULL
  out <- vector("list", length(times))
  for (g in seq_along(times)) {
    t <- times[g]
    rows <- ev[ev$TIME == t, ]
    cf_next <- ref_coef(model, P, eCL, eV, rows$CRCA[1], rows$AJBW[1])
    if (g > 1) {
      cf_int <- if (interval_cov == "nocb") cf_next else cf
      # advance tcur -> t, splitting at infusion ends
      repeat {
        act <- inf_end > tcur
        nxt <- if (any(act)) min(min(inf_end[act]), t) else t
        r <- sum(inf_rate[act])  # all active infusions run over [tcur, nxt)
        x <- ref_advance(x, cf_int, r, nxt - tcur)
        tcur <- nxt
        if (tcur >= t) break
      }
    }
    pre <- x
    if (any(rows$EVID == 4)) { x <- c(0, 0); inf_end <- numeric(0); inf_rate <- numeric(0) }
    for (k in which(rows$EVID %in% c(1, 4) & rows$AMT > 0)) {
      stopifnot(rows$RATE[k] > 0)  # S-1 uses zero-order infusions only
      inf_end <- c(inf_end, t + rows$AMT[k] / rows$RATE[k])
      inf_rate <- c(inf_rate, rows$RATE[k])
    }
    cf <- cf_next
    out[[g]] <- data.frame(time = t, Ac_pre = pre[1], Ap_pre = pre[2],
                           Ac = x[1], Ap = x[2], V1 = cf$V1,
                           cp = x[1] / cf$V1)
    tcur <- t
  }
  do.call(rbind, out)
}
