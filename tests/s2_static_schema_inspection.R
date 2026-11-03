#!/usr/bin/env Rscript
# S-2 static fit-schema inspection (NO fitting of any kind, synthetic or real).
# Inspects: explicit-ODE compilation, FOCEI control interface, sensitivity tolerances, covariance options,
# the parameter coordinates declared by the installed packages, and which v1.1 criteria (C1, C4-C9, U2) can
# be evaluated from what. Anything that needs an actual fit object is marked UNVERIFIED.
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)]); setwd(normalizePath(file.path(dirname(f), "..")))
suppressMessages({ library(rxode2); library(nlmixr2est); library(jsonlite) })
source("R/pkpd02_ode_models.R"); source("R/pkpd02_reference_matexp.R")
tc <- function(expr) tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
out <- list(inspection = "S-2 static fit-schema inspection", timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), fitting_performed = FALSE,
            versions = lapply(c(rxode2 = "rxode2", nlmixr2est = "nlmixr2est", nlmixr2 = "nlmixr2", lotri = "lotri", rxode2lincmt = "rxode2lincmt", n1qn1 = "n1qn1", lbfgsb3c = "lbfgsb3c"), function(p) as.character(packageVersion(p))))

# ---- 1. FOCEI control interface -------------------------------------------------------------------------------
fm <- formals(foceiControl); fc <- foceiControl()
ch <- function(n) if (n %in% names(fm) && is.call(fm[[n]])) as.character(eval(fm[[n]])) else if (n %in% names(fm)) as.character(fm[[n]]) else NULL
out$control_interface <- list(
  n_formals = length(fm), covMethod_choices = ch("covMethod"), covMethod_default = "r,s", sensMethod_choices = ch("sensMethod"), derivMethod_choices = ch("derivMethod"), covDerivMethod_choices = ch("covDerivMethod"),
  diagXform_choices = ch("diagXform"), diagXform_default = fc$diagXform, scaleType_choices = ch("scaleType"), boundTol_default = fc$boundTol, interaction_default = fc$interaction, addProp_choices = ch("addProp"), addProp_default = fc$addProp,
  covFallback_default = deparse(fm$covFallback), covFallback_list_accepted = !is.character(tc(foceiControl(covFallback = list()))), covFallback_stored_when_empty_list = deparse(foceiControl(covFallback = list())$covFallback),
  covTryHarder_default = fc$covTryHarder, covFull_default = fc$covFull, covShortcut_default = fc$covShortcut, sigdig_default = fc$sigdig, outerOpt_default_internal_code = fc$outerOpt, innerOpt_default_internal_code = fc$innerOpt,
  outerOpt_choices = ch("outerOpt"), innerOpt_choices = ch("innerOpt"), literalFix_default = fc$literalFix, maxOdeRecalc_default = fc$maxOdeRecalc, odeRecalcFactor_default = fc$odeRecalcFactor,
  tolerance_routes = list(top_level_atol_argument_accepted = !is.character(tc(foceiControl(atol = 1e-10))), via_rxControl_accepted = !is.character(tc(foceiControl(rxControl = rxControl(atol = 1e-10, rtol = 1e-10, atolSens = 1e-10, rtolSens = 1e-10, covsInterpolation = "locf")))),
    stored_when_set_via_rxControl = unlist(foceiControl(rxControl = rxControl(atol = 1e-10, rtol = 1e-10, atolSens = 1e-9, rtolSens = 1e-9, covsInterpolation = "locf"))$rxControl[c("atol", "rtol", "atolSens", "rtolSens", "covsInterpolation")])),
  nlmixr2_default_rxControl = fc$rxControl[c("atol", "rtol", "atolSens", "rtolSens", "maxsteps", "method", "covsInterpolation")],
  covsInterpolation_default_code_1_is_locf = identical(as.integer(fc$rxControl$covsInterpolation), as.integer(rxControl(covsInterpolation = "locf")$covsInterpolation)))

# ---- 2. parameter coordinates declared by the specification (with v1.1 bounds proposed) -----------------------------------------
b1 <- function() { ini({ tcl <- c(log(0.5), log(4), log(50)); tv <- c(log(2), log(18), log(200)); eta.cl ~ 0.09; eta.v ~ 0.09; add.sd <- c(0.001, 0.5, 20); prop.sd <- c(0.001, 0.1, 2) })
  model({ cl <- exp(tcl + eta.cl) * CRCA / 100; v <- exp(tv + eta.v) * AJBW / 70; d/dt(Ac) <- -cl / v * Ac; cp <- Ac / v; cp ~ add(add.sd) + prop(prop.sd) + combined2() }) }
b2 <- function() { ini({ tcl <- c(log(0.5), log(4), log(50)); tv1 <- c(log(2), log(14), log(200)); tq <- c(log(0.01), log(3), log(100)); tv2 <- c(log(0.5), log(9), log(500)); eta.cl ~ 0.09; eta.v1 ~ 0.09; add.sd <- c(0.001, 0.5, 20); prop.sd <- c(0.001, 0.1, 2) })
  model({ cl <- exp(tcl + eta.cl) * CRCA / 100; v1 <- exp(tv1 + eta.v1) * AJBW / 70; q <- exp(tq); v2 <- exp(tv2) * AJBW / 70; d/dt(Ac) <- -cl / v1 * Ac - q * (Ac / v1 - Ap / v2); d/dt(Ap) <- q * (Ac / v1 - Ap / v2); cp <- Ac / v1; cp ~ add(add.sd) + prop(prop.sd) + combined2() }) }
u1 <- suppressMessages(rxode(b1)); u2 <- suppressMessages(rxode(b2))
coord <- function(u) { d <- u$iniDf; d <- d[, c("ntheta", "neta1", "name", "est", "lower", "upper", "fix", "err")]; d$THETA_index_in_compiled_models <- ifelse(is.na(d$ntheta), NA, paste0("THETA[", d$ntheta, "]")); d$ETA_index_in_compiled_models <- ifelse(is.na(d$neta1), NA, paste0("ETA[", d$neta1, "]")); d }
out$parameter_coordinates <- list("1c" = coord(u1), "2c" = coord(u2),
  notes = c("Fixed effects are typical values on the LOG scale (tcl, tv, tq, tv2); v1.1 physical bounds are encoded as log-scale lower/upper (accepted by ini()).",
            "Residual SDs add.sd/prop.sd are on the natural scale with positive bounds; Omega is declared as variances (eta ~ 0.09).",
            "foceiControl diagXform = 'sqrt' means Omega diagonals are optimised on an SD-like (sqrt) coordinate; whether the RETURNED covariance is in that coordinate or on the variance scale is UNVERIFIED until a fit object exists.",
            "boundTol default 0.05 is an ESTIMATION-scale distance; v1.1 C8 requires >0.05 on the estimation log scale: map must be confirmed from the fit object (UNVERIFIED).",
            "Rendering of 'additive SD >= 0.05 mg/L' and 'proportional SD >= 0.02' (C8) needs the natural-scale residual SDs, available from the declared coordinates (verified) but their returned-covariance coordinate is UNVERIFIED."))
cm <- function(u) { m <- u$foceiModel; lapply(m[c("inner", "predOnly", "predNoLhs")], function(x) { mv <- rxModelVars(x); list(states = as.character(mv$state), n_states = length(mv$state), n_sensitivity_states = sum(grepl("^rx__sens_", mv$state)), params = as.character(mv$params), linCmt_flag = unname(mv$flags["linCmt"]),
  contains_linCmt_text = any(grepl("linCmt", rxNorm(x)))) }) }
out$compiled_models <- list("1c" = cm(u1), "2c" = cm(u2))

# ---- 3. fit-object accessors present in the installed nlmixr2est (names only; behaviour needs a fit) ------------------------------------
ns <- asNamespace("nlmixr2est"); acc <- sub("^nmObjGet\\.", "", grep("^nmObjGet\\.", ls(ns), value = TRUE)); out$fit_object_accessors_installed <- sort(acc)
need <- c("atol", "rtol", "covsInterpolation", "control", "cov", "covMethod", "omega", "omegaR", "theta", "eta", "etaMat", "objf", "objective", "parFixed", "parFixedDf", "sigma", "message", "scaleInfo", "dataSav", "ui", "iniUi", "covLvl", "fixef", "seed", "time", "condition", "cor", "coefficients")
out$accessor_presence <- as.list(setNames(need %in% acc, need))
out$accessor_presence_note <- "TRUE = a dedicated nmObjGet.<name> method exists. FALSE does NOT mean the field is absent from a fit object (e.g. fit$cov, fit$objf are data elements, not accessor methods); their names/coordinates are UNVERIFIED until a synthetic fit object is inspected."
out$generics_present <- list(vcov = exists("vcov"), getVarCov = exists("getVarCov"), coef = exists("coef"))

# ---- 4. extraction map: v1.1 criterion -> candidate source -> status ----------------------------------------------------------------------
em <- data.frame(criterion = c("L9/S-1r: solver pathway, LOCF, tolerances in effect", "C1: optimiser success code, boundary warnings", "C4: named, positive-definite covariance with verified r,s method and NO fallback", "C5/C6: structural log-parameter correlation block",
                               "C7: log-scale SEs (CL, V; Q, V2) and residual-SD RSE", "C8: distance from bounds (estimation log scale), SD floors", "C9: IIV variances >= 0.001", "U2: Sigma_z -> xi Jacobian, named block", "Reproducibility: seed / control used"),
  candidate_source = c("fit$atol / fit$rtol / fit$covsInterpolation accessors (present) + rxModelVars of the stored model for linCmt flag", "fit$message / optimiser return code in control+env (field names UNVERIFIED)", "fit$cov + fit$covMethod (names UNVERIFIED); control covFallback=list() to forbid fallback", "named submatrix of fit$cov (UNVERIFIED coordinates)",
                       "fit$cov diagonal / parFixedDf (UNVERIFIED)", "fit$ui$iniDf bounds + boundTol; returned estimates in log scale (UNVERIFIED mapping)", "fit$omega (variance scale expected; UNVERIFIED)", "needs coordinate documentation of fit$cov (sqrt/log xform of Omega, SD of sigma): UNVERIFIED", "fit$control / seed accessors (present)"),
  status = c("PARTIALLY VERIFIED (accessors exist; values unverified)", rep("UNVERIFIED until a synthetic fit object exists", 8)), stringsAsFactors = FALSE)
out$extraction_map <- em

# ---- 5. PILOT (not a fit): inner-model sensitivity accuracy vs finite differences of the independent reference; solve cost ---------------------------------
inner <- u1$foceiModel$inner; arch <- "archive/s1_failed_2026-10-30/evidence"
ev <- read.csv(file.path(arch, "s1_events.csv")); et <- read.csv(file.path(arch, "s1_individual_etas.csv")); ev$CMT <- as.integer(ev$CMT); ev$EVID <- as.integer(ev$EVID)
ids <- c(1, 2, 101:130); sub <- ev[ev$ID %in% ids, ]; etas <- et[et$ID %in% ids, ]
pr <- data.frame(ID = etas$ID, check.names = FALSE); pr[["THETA[1]"]] <- log(4); pr[["THETA[2]"]] <- log(18); pr[["THETA[3]"]] <- 0.5; pr[["THETA[4]"]] <- 0.1; pr[["ETA[1]"]] <- etas$eCL; pr[["ETA[2]"]] <- etas$eV
P <- c(tCL = 4, tV = 18); kf <- function(i, t) paste(i, sprintf("%.9f", t)); cols <- c("ID", "TIME", "EVID", "AMT", "RATE", "CMT", "CRCA", "AJBW")
fd <- function(k, h = 1e-5) do.call(rbind, lapply(split(sub, sub$ID), function(e) { e1 <- etas[etas$ID == e$ID[1], ]; d <- c(0, 0); d[k] <- h
  up <- ref_matexp(e[, cols], "1c", P, c(e1$eCL, e1$eV) + d); dn <- ref_matexp(e[, cols], "1c", P, c(e1$eCL, e1$eV) - d); data.frame(ID = e$ID[1], TIME = up$TIME, dAc = (up$Ac - dn$Ac) / (2 * h)) }))
f1 <- fd(1); f2 <- fd(2)
cfgs <- list(list("nlmixr2_defaults", 1e-6, 1e-3, 1e-6, 1e-3), list("state_1e-10_sens_defaults", 1e-10, 1e-10, 1e-6, 1e-3), list("all_1e-8", 1e-8, 1e-8, 1e-8, 1e-8), list("all_1e-10", 1e-10, 1e-10, 1e-10, 1e-10), list("all_1e-12", 1e-12, 1e-12, 1e-12, 1e-12))
pilot <- do.call(rbind, lapply(cfgs, function(cf) { slv <- function() as.data.frame(rxSolve(inner, pr, sub, covsInterpolation = "locf", atol = cf[[2]], rtol = cf[[3]], atolSens = cf[[4]], rtolSens = cf[[5]], addDosing = TRUE, cores = 1)); o <- slv()
  tm <- median(replicate(5, system.time(slv())[["elapsed"]]))
  names(o) <- sub("^id$", "ID", names(o)); m1 <- match(kf(o$ID, o$time), kf(f1$ID, f1$TIME)); ok <- !is.na(m1); e1 <- abs(o$rx__sens_Ac_BY_ETA_1___[ok] - f1$dAc[m1[ok]]) / pmax(1, abs(f1$dAc[m1[ok]])); e2 <- abs(o$rx__sens_Ac_BY_ETA_2___[ok] - f2$dAc[m1[ok]]) / pmax(1, abs(f2$dAc[m1[ok]]))
  data.frame(config = cf[[1]], atol = cf[[2]], rtol = cf[[3]], atolSens = cf[[4]], rtolSens = cf[[5]], rows = sum(ok), patients = nrow(pr), median_solve_seconds = tm, max_scaled_err_dAc_dEta1 = max(e1), max_scaled_err_dAc_dEta2 = max(e2)) }))
out$pilot_sensitivity_accuracy_nonfit <- list(table = pilot, reference = "central finite difference (h=1e-5) of the independent matrix-exponential reference; its own noise floor is ~2e-7 (scaled), so errors at/below that cannot be resolved",
  note = "Solve of the compiled inner model on 32 archived synthetic patients (915 rows). Not an estimation run; solve cost is per function evaluation only.")
out$blockers_found_statically <- c("covariance coordinates for Omega/sigma/bounded parameters in the returned fit object are undocumented here (need a synthetic fit to read them)",
  "automatic covariance fallback chain r,s -> r -> s is the DEFAULT; C4 forbids undocumented fallback (covFallback = list() accepted by the control, effect on a fit UNVERIFIED)",
  "nlmixr2 default ODE/sensitivity tolerances (1e-6/1e-3) fail the S-1 rule and give ~3e-3 sensitivity errors; must be overridden via rxControl=",
  "tolerances cannot be passed as top-level foceiControl arguments (rejected); only via rxControl=",
  "no way to verify estimation convergence behaviour, runtime of full FOCEI, or multi-start stability without a (not yet approved) synthetic fit")
dir.create("evidence", FALSE); writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null", dataframe = "rows"), "evidence/s2_static_schema_inspection.json")
write.csv(pilot, "evidence/s2_pilot_sensitivity_accuracy_nonfit.csv", row.names = FALSE); write.csv(em, "evidence/s2_extraction_map.csv", row.names = FALSE)
print(pilot, row.names = FALSE); cat("\naccessors present among needed:\n"); print(unlist(out$accessor_presence)); cat("\nSTATIC INSPECTION DONE (no fit performed)\n")
