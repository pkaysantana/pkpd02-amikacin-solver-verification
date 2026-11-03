# PKPD-02 v1.2 (DRAFT amendment): explicit amount-based ODE implementations.
# Equations are the v1.1 addendum s2 equations, unchanged:
#   CL(t) = thetaCL * CRCA(t)/100 * exp(eta_CL)
#   V(t)  = thetaV  * AJBW(t)/70  * exp(eta_V)          (1C)
#   V1(t) = thetaV1 * AJBW(t)/70  * exp(eta_V1), Q = thetaQ, V2(t) = thetaV2*AJBW(t)/70   (2C)
#   zero-order IV input R(t) into the central amount; f = Ac / V (central)
#   residual: combined2, v(f) = add.sd^2 + prop.sd^2 f^2 (both arms)
# Typical values are written on the log scale (tcl = log(thetaCL)); this is a
# mathematically identical re-expression used here ONLY for pathway inspection. The
# estimation-scale/bounds decision belongs to the later fit-schema stage (not made here).
# Amounts, not concentrations, are the states, so they stay continuous when CL/V change.
# NO linCmt() appears anywhere in these models.

pkpd02_1c_ode <- function() {
  ini({
    tcl <- log(4.0); tv <- log(18.0)
    eta.cl ~ 0.09; eta.v ~ 0.09
    add.sd <- 0.5; prop.sd <- 0.1
  })
  model({
    cl <- exp(tcl + eta.cl) * CRCA / 100
    v  <- exp(tv + eta.v) * AJBW / 70
    d/dt(Ac) <- -cl / v * Ac
    cp <- Ac / v
    cp ~ add(add.sd) + prop(prop.sd) + combined2()
  })
}

pkpd02_2c_ode <- function() {
  ini({
    tcl <- log(4.0); tv1 <- log(14.0); tq <- log(3.0); tv2 <- log(9.0)
    eta.cl ~ 0.09; eta.v1 ~ 0.09
    add.sd <- 0.5; prop.sd <- 0.1
  })
  model({
    cl <- exp(tcl + eta.cl) * CRCA / 100
    v1 <- exp(tv1 + eta.v1) * AJBW / 70
    q  <- exp(tq)
    v2 <- exp(tv2) * AJBW / 70
    d/dt(Ac) <- -cl / v1 * Ac - q * (Ac / v1 - Ap / v2)
    d/dt(Ap) <-  q * (Ac / v1 - Ap / v2)
    cp <- Ac / v1
    cp ~ add(add.sd) + prop(prop.sd) + combined2()
  })
}

# Plain rxode2 forms of the same equations (used for the bulk equivalence runs; the
# UI forms above are separately shown to compile to the same ODE system).
pkpd02_rx_ode <- list(
  "1c" = function() rxode2::rxode2({
    CL <- tCL * CRCA / 100 * exp(eCL)
    V  <- tV * AJBW / 70 * exp(eV)
    d/dt(Ac) <- -CL / V * Ac
    cp <- Ac / V
  }),
  "2c" = function() rxode2::rxode2({
    CL <- tCL * CRCA / 100 * exp(eCL)
    V1 <- tV1 * AJBW / 70 * exp(eV)
    Q  <- tQ
    V2 <- tV2 * AJBW / 70
    d/dt(Ac) <- -CL / V1 * Ac - Q * (Ac / V1 - Ap / V2)
    d/dt(Ap) <-  Q * (Ac / V1 - Ap / V2)
    cp <- Ac / V1
  }))
