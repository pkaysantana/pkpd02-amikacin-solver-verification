#!/usr/bin/env Rscript
# Self-test of the independent base-R reference (R/pkpd02_reference_matexp.R).
# Uses closed-form solutions, expm::expm (self-test only) and invariants. No rxode2.
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(normalizePath(file.path(dirname(f), ".."))); source("R/pkpd02_reference_matexp.R")
res <- list(); ok <- function(name, err, tol) { res[[length(res)+1L]] <<- data.frame(check = name, max_abs_err = err, tol = tol, pass = is.finite(err) && err <= tol) }
set.seed(20261102)
# 1. Pade vs expm on random matrices incl. stiff ones
e <- 0; for (i in 1:200) { A <- matrix(rnorm(9, 0, 10^runif(1, -2, 2)), 3); A[3, ] <- 0; e <- max(e, max(abs(pade_expm(A) - expm::expm(A)) / pmax(1, abs(expm::expm(A))))) }
ok("pade_expm vs expm::expm (200 random 3x3, rel)", e, 1e-12)
# 2. 1C closed form (single zero-order infusion, constant coefficients)
mk <- function(times, evid, amt, rate, crca = 100, w = 70) data.frame(ID = 1, TIME = times, EVID = evid, AMT = amt, RATE = rate, CMT = 1L, CRCA = crca, AJBW = w)
P1 <- c(tCL = 4, tV = 18); k <- 4 / 18; R <- 50; Td <- 2
tt <- c(0.5, 1, 2, 3, 7.5, 20); ev <- mk(c(0, tt), c(1L, rep(0L, length(tt))), c(100, rep(0, length(tt))), c(R, rep(0, length(tt))))
r <- ref_matexp(ev, "1c", P1, c(0, 0))
cf <- ifelse(tt <= Td, R / 4 * (1 - exp(-k * tt)), R / 4 * (1 - exp(-k * Td)) * exp(-k * (tt - Td)))
ok("1C single infusion vs closed form (mg/L)", max(abs(r$cp[-1] - cf)), 1e-12)
# 3. 2C closed form via bi-exponential: Ac(t) for infusion rate R, zero initial states
P2 <- c(tCL = 4, tV1 = 14, tQ = 3, tV2 = 9); k10 <- 4 / 14; k12 <- 3 / 14; k21 <- 3 / 9
s <- k10 + k12 + k21; disc <- sqrt(s^2 - 4 * k10 * k21); al <- (s + disc) / 2; be <- (s - disc) / 2
Ac_inf <- function(t, R) { # integral of the bolus impulse response (per unit dose) x R, zero initial states
  R * ((k21 - al) / (be - al) * (1 - exp(-al * t)) / al + (k21 - be) / (al - be) * (1 - exp(-be * t)) / be) }
Ac_f <- function(t) ifelse(t <= Td, Ac_inf(t, R), Ac_inf(t, R) - Ac_inf(pmax(t - Td, 0), R))
ev <- mk(c(0, tt), c(1L, rep(0L, length(tt))), c(100, rep(0, length(tt))), c(R, rep(0, length(tt))))
r <- ref_matexp(ev, "2c", P2, c(0, 0))
ok("2C single infusion vs closed form (mg/L)", max(abs(r$cp[-1] - Ac_f(tt) / 14)), 1e-11)
# 4. mass balance with CL=0: Ac+Ap == delivered amount, overlapping infusions, covariate changes
P0 <- c(tCL = 0, tV1 = 14, tQ = 3, tV2 = 9)
ev <- rbind(mk(c(0, 1), c(1L, 1L), c(300, 200), c(100, 50)), mk(c(0.7, 2.2, 5, 9), c(2L, 0L, 2L, 0L), 0, 0, crca = 100))
ev$AJBW <- c(70, 70, 40, 120, 90, 55); ev <- ev[order(ev$TIME), ]
r <- ref_matexp(ev, "2c", P0, c(0, 0))
deliv <- function(t) pmin(100 * t, 300) + pmin(50 * pmax(t - 1, 0), 200)
ok("2C CL=0 mass balance Ac+Ap == delivered (mg)", max(abs(r$Ac + r$Ap - deliv(r$TIME))), 1e-9)
# 5. splitting an interval with a no-op row must not change anything
ev1 <- mk(c(0, 10), c(1L, 0L), c(500, 0), c(100, 0)); ev2 <- mk(c(0, 3, 4.5, 10), c(1L, 2L, 2L, 0L), c(500, 0, 0, 0), c(100, 0, 0, 0))
ok("no-op covariate rows do not change result (2C)", abs(ref_matexp(ev1, "2c", P2, c(.1, -.1))$cp[2] - tail(ref_matexp(ev2, "2c", P2, c(.1, -.1))$cp, 1)), 1e-12)
# 6. reset zeroes states and cancels an active infusion
ev <- mk(c(0, 1, 1.5, 3), c(1L, 3L, 0L, 0L), c(500, 0, 0, 0), c(100, 0, 0, 0))   # infusion 0-5h cancelled at 1 by EVID=3
r <- ref_matexp(ev, "1c", P1, c(0, 0))
ok("EVID=3 reset cancels active infusion (cp after reset == 0)", max(abs(r$cp[r$TIME >= 1])), 0)
out <- do.call(rbind, res); print(out, row.names = FALSE)
dir.create("evidence", FALSE); write.csv(out, "evidence/s1r_reference_selftest.csv", row.names = FALSE)
quit(status = if (all(out$pass)) 0L else 1L)
