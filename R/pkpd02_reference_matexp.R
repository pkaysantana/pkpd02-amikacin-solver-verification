# PKPD-02 v1.2: INDEPENDENT base-R per-interval matrix-exponential reference.
# Independence: written from the addendum equations; no rxode2/nlmixr2 code, no helper from the
# archived S-1 test, own event handling, and a Pade(6,6) scaling-and-squaring exponential
# (the archived reference used a Taylor series). Cross-checked in tests/s1r_reference_selftest.R
# against closed-form solutions and expm::expm.
#
# Semantics implemented (each is a stated convention, tested against the ODE):
#  * Breakpoints = every TIME, every infusion start, every infusion end.
#  * On [a,b) coefficients come from the LAST row with TIME <= a (LOCF); the infusion input is
#    the sum of rates of infusions with start <= a < end.
#  * States are continuous (amounts) across covariate changes, infusion starts and stops.
#  * EVID 3 and EVID 4 reset both amounts to zero AND cancel any active infusion (episode
#    boundary, protocol T1); EVID 4 then starts its own infusion at that time.
#  * Row output at TIME t: state after all events at t (zero-order input => no state jump at a
#    dose start), concentration = Ac / V1 using the covariates carried at t. "Left" states are the
#    limits just before t.

pade_expm <- function(A) {
  n <- nrow(A); I <- diag(n)
  s <- max(0L, ceiling(log2(max(sum(abs(A)), 1e-300))) + 1L)
  A <- A / 2^s
  q <- 6L; cc <- 0.5; X <- A; N <- I + cc * A; D <- I - cc * A; pos <- TRUE
  for (k in 2:q) {
    cc <- cc * (q - k + 1) / (k * (2 * q - k + 1))
    X <- A %*% X
    N <- N + cc * X
    D <- if (pos) D + cc * X else D - cc * X
    pos <- !pos
  }
  E <- solve(D, N)
  for (i in seq_len(s)) E <- E %*% E
  E
}

# par: named vector (1c: tCL,tV ; 2c: tCL,tV1,tQ,tV2) on the PHYSICAL scale; eta: c(CL, V|V1)
ref_matexp <- function(ev, structure, par, eta) {
  stopifnot(all(diff(ev$TIME) >= 0))
  CL_f <- function(crca) par[["tCL"]] * crca / 100 * exp(eta[1])
  V1_f <- function(w) (if (structure == "1c") par[["tV"]] else par[["tV1"]]) * w / 70 * exp(eta[2])
  Q_f  <- if (structure == "2c") par[["tQ"]] else 0
  V2_f <- function(w) if (structure == "2c") par[["tV2"]] * w / 70 else 1

  times <- sort(unique(ev$TIME))
  # active infusions: start, end, rate (reset removes all)
  inf <- data.frame(s = numeric(0), e = numeric(0), r = numeric(0))
  x <- c(0, 0); cov_crca <- NA; cov_w <- NA
  out <- vector("list", length(times))
  prev <- NA
  for (g in seq_along(times)) {
    t <- times[g]
    if (g > 1) {                                   # propagate prev -> t over sub-intervals
      bp <- sort(unique(c(prev, t, inf$e[inf$e > prev & inf$e < t])))
      for (j in seq_len(length(bp) - 1L)) {
        a <- bp[j]; b <- bp[j + 1L]
        rate <- sum(inf$r[inf$s <= a & inf$e > a])
        V1 <- V1_f(cov_w); V2 <- V2_f(cov_w); k10 <- CL_f(cov_crca) / V1
        k12 <- Q_f / V1; k21 <- Q_f / V2
        M <- rbind(c(-(k10 + k12), k21, rate), c(k12, -k21, 0), c(0, 0, 0))
        x <- as.numeric((pade_expm(M * (b - a)) %*% c(x, 1))[1:2])
      }
      inf <- inf[inf$e > t, , drop = FALSE]        # drop finished infusions
    }
    left <- x
    rows <- ev[ev$TIME == t, , drop = FALSE]
    if (any(rows$EVID %in% c(3L, 4L))) { x <- c(0, 0); inf <- inf[0, ] }
    d <- rows[rows$EVID %in% c(1L, 4L) & rows$AMT > 0, , drop = FALSE]
    if (nrow(d)) inf <- rbind(inf, data.frame(s = t, e = t + d$AMT / d$RATE, r = d$RATE))
    cov_crca <- rows$CRCA[1]; cov_w <- rows$AJBW[1]
    out[[g]] <- data.frame(TIME = t, Ac_left = left[1], Ap_left = left[2], Ac = x[1], Ap = x[2],
                           V1 = V1_f(cov_w), cp = x[1] / V1_f(cov_w))
    prev <- t
  }
  do.call(rbind, out)
}
