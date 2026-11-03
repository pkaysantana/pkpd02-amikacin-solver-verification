#!/usr/bin/env Rscript
# S-1r item 5: verify the ACTUAL execution pathway of the explicit-ODE models; show that a
# silent ODE -> linCmt() conversion is possible in this install and that the guarded wrapper
# prevents it. No fitting is performed (the FOCEI models are only compiled, not estimated).
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(normalizePath(file.path(dirname(f), "..")))
suppressMessages({ library(rxode2); library(nlmixr2est); library(jsonlite) })
source("R/pkpd02_ode_models.R"); source("R/pkpd02_solve.R"); source("R/pkpd02_reference_matexp.R")

ARCH <- "archive/s1_failed_2026-10-30/evidence"
ev <- read.csv(file.path(ARCH, "s1_events.csv")); et <- read.csv(file.path(ARCH, "s1_individual_etas.csv"))
ev1 <- ev[ev$ID == 1, ]; ic <- data.frame(ID = 1L, eta.cl = et$eCL[et$ID == 1], eta.v = et$eV[et$ID == 1])
th <- c(tcl = log(4), tv = log(18)); par <- c(tCL = 4, tV = 18)
ui <- suppressMessages(rxode(pkpd02_1c_ode)); ui2 <- suppressMessages(rxode(pkpd02_2c_ode))
tol_ok <- function(x, ref) all(abs(x - ref) <= pkpd02_tol(ref))
ref <- ref_matexp(ev1, "1c", par, c(et$eCL[et$ID == 1], et$eV[et$ID == 1]))
refcp <- function(r) ref$cp[match(sprintf("%.12f", r$time), sprintf("%.12f", ref$TIME))]
slv <- function(...) suppressMessages(rxSolve(ui, th, ev1, iCov = ic, covsInterpolation = "locf", atol = 1e-12, rtol = 1e-12,
                                              addDosing = TRUE, cores = 1, ...))
desc <- function(r) { mv <- rxModelVars(r); d <- as.data.frame(r)
  list(linCmt_flag = unname(mv$flags["linCmt"]), states = as.character(mv$state),
       gate_pass_vs_independent_reference = tol_ok(d$cp, refcp(d)), max_abs_vs_reference = max(abs(d$cp - refcp(d)))) }

out <- list(gate = "S-1r pathway inspection", versions = list(
  rxode2 = as.character(packageVersion("rxode2")), nlmixr2est = as.character(packageVersion("nlmixr2est")),
  rxode2lincmt = as.character(packageVersion("rxode2lincmt"))))

# --- 1. useLinCmt support/defaults in the installed versions ---
out$useLinCmt <- list(
  rxSolve_rxUi_formals_default = deparse(formals(rxode2:::rxSolve.rxUi)$useLinCmt),
  global_option_rxode2.useLinCmt = deparse(getOption("rxode2.useLinCmt")),
  explicit_FALSE_accepted = !inherits(try(slv(useLinCmt = FALSE), silent = TRUE), "try-error"),
  nlmixr2est_functions_referencing_useLinCmt = {
    ns <- asNamespace("nlmixr2est"); h <- Filter(function(n) { fn <- get(n, ns); is.function(fn) &&
      any(grepl("useLinCmt", deparse(fn))) }, ls(ns, all.names = TRUE)); as.list(h) },
  nlmixr2_foceiControl_default_rxControl = foceiControl()$rxControl[c("atol", "rtol", "atolSens", "rtolSens", "method", "covsInterpolation")])

# --- 2. compiled FOCEI / simulation models (compiled, NOT fitted) contain explicit ODE states ---
cm <- function(u) { fm <- u$foceiModel
  l <- list(simulationModel = u$simulationModel, foceiInner = fm$inner, foceiPredOnly = fm$predOnly, foceiPredNoLhs = fm$predNoLhs)
  lapply(l, function(m) { mv <- rxModelVars(m); list(linCmt_flag = unname(mv$flags["linCmt"]),
    states = as.character(mv$state), norm_text_contains_linCmt = any(grepl("linCmt", rxNorm(m)))) }) }
out$compiled_models <- list("1c" = cm(ui), "2c" = cm(ui2))
out$reference_flag_for_genuine_linCmt_model <- unname(rxModelVars(rxode2({ CL <- 4; V <- 18; cp <- linCmt() }))$flags["linCmt"])

# --- 3. solving behaviour under each request (1 patient, canonical events) ---
res <- list(default = desc(slv()), explicit_FALSE = desc(slv(useLinCmt = FALSE)), explicit_TRUE = desc(slv(useLinCmt = TRUE)))
options(rxode2.useLinCmt = TRUE)
res$option_TRUE_default_arg <- desc(slv())
res$option_TRUE_explicit_FALSE <- desc(slv(useLinCmt = FALSE))
guarded <- pkpd02_solve_ode(ui, "1c", th, ev1, ic)
res$option_TRUE_via_guarded_wrapper <- list(pathway_checks = as.list(attr(guarded, "pathway_checks")),
  gate_pass_vs_independent_reference = tol_ok(guarded$cp, refcp(guarded)), max_abs_vs_reference = max(abs(guarded$cp - refcp(guarded))))
options(rxode2.useLinCmt = NULL)
out$solve_behaviour <- res

# --- 4. the actual rxSolve args recorded on a solved object (guarded wrapper) ---
g0 <- pkpd02_solve_ode(ui, "1c", th, ev1, ic)
out$guarded_wrapper_pathway_checks <- as.list(attr(g0, "pathway_checks"))

# --- 5. interpretation flags (pass = pathway behaves as required) ---
chk <- c(
  default_is_explicit_ODE = isTRUE(res$default$linCmt_flag == -100) && identical(res$default$states, "Ac"),
  explicit_FALSE_is_explicit_ODE = isTRUE(res$explicit_FALSE$linCmt_flag == -100) && res$explicit_FALSE$gate_pass_vs_independent_reference,
  explicit_FALSE_overrides_global_option = isTRUE(res$option_TRUE_explicit_FALSE$linCmt_flag == -100) && res$option_TRUE_explicit_FALSE$gate_pass_vs_independent_reference,
  guarded_wrapper_immune_to_global_option = all(unlist(res$option_TRUE_via_guarded_wrapper$pathway_checks)) && res$option_TRUE_via_guarded_wrapper$gate_pass_vs_independent_reference,
  all_compiled_models_free_of_linCmt = all(vapply(unlist(out$compiled_models, recursive = FALSE), function(m) m$linCmt_flag == -100 && !m$norm_text_contains_linCmt, NA)),
  nlmixr2est_never_references_useLinCmt = length(out$useLinCmt$nlmixr2est_functions_referencing_useLinCmt) == 0,
  silent_conversion_hazard_demonstrated = isTRUE(res$explicit_TRUE$linCmt_flag == 1) && !res$explicit_TRUE$gate_pass_vs_independent_reference)
out$checks <- as.list(chk); out$pass <- all(chk)
out$notes <- c("silent_conversion_hazard_demonstrated is a hazard demonstration (useLinCmt=TRUE converts the ODE and reproduces the S-1 error); it must be TRUE for the guard to be meaningful.",
  "nlmixr2est fitting was NOT run; the compiled FOCEI models are the ones fit routines use. The synthetic fit-schema stage must re-assert the pathway from a fit object.")
dir.create("evidence", FALSE)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), "evidence/s1r_pathway_inspection.json")
print(unlist(chk)); cat("PATHWAY INSPECTION:", if (out$pass) "PASS" else "FAIL", "\n")
quit(status = if (out$pass) 0L else 1L)
