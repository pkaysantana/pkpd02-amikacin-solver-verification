#!/usr/bin/env Rscript
# Non-destructive runner for the synthetic-data verification gates.
#
# Why this exists. The original gate scripts write their evidence files in place. Running the S-1 script rewrites
# evidence/s1_*.json (with a new timestamp), and the S-1r runner then correctly reports that the in-place S-1 files no longer match
# their preserved archive, so S-1r fails only because S-1 ran first. The S-1 script and its evidence cannot be edited without
# breaking the hash-manifested failure record, so this runner leaves them alone: each gate runs in its OWN temporary copy of
# the repository. Committed evidence and archive/ are never touched, and gates cannot affect one another.
#
# Usage (repository root):  Rscript tools/verify_public.R [gate ...] [--library DIR] [--keep]
#   gates: s1  s1r  s2prep_canonical  s2_static        (default: all four, about 2.5 minutes)
#   --library DIR   renv library to link into the temporary copies (default: renv/library of this repository)
#   --keep          keep the temporary copies and logs for inspection
# Exit status 0 only if every requested gate behaved as its recorded outcome (note: s1 is EXPECTED to exit 1).
# Needs only base R plus the packages restored by renv. No source data are read by any of these gates.

args <- commandArgs(trailingOnly = TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(self), ".."), mustWork = TRUE)
keep <- "--keep" %in% args
li <- match("--library", args); lib <- if (!is.na(li)) args[li + 1L] else file.path(root, "renv", "library")
want <- setdiff(args, c("--keep", "--library", if (!is.na(li)) args[li + 1L]))

GATES <- list(
  s1               = list(script = "tests/s1_synthetic_solver_equivalence.R", exit = 1L, text = "S-1 VERDICT: FAIL",
                          what = "original solver gate: recorded, intended FAIL"),
  s1r              = list(script = "tests/s1r_run_all.R", exit = 0L, text = "PASS_CONDITIONAL_ON_INPUT_CONTRACT",
                          what = "remediation gate: conditional pass"),
  s2prep_canonical = list(script = "tests/s2prep_canonical_equivalence.R", exit = 0L, text = "CANONICAL EQUIVALENCE VERDICT: PASS",
                          what = "event-representation gate, synthetic batteries"),
  s2_static        = list(script = "tests/s2_static_schema_inspection.R", exit = 0L, text = "STATIC INSPECTION DONE (no fit performed)",
                          what = "static fit-interface inspection, no fit")
)
if (!length(want)) want <- names(GATES)
if (any(!want %in% names(GATES))) stop("unknown gate(s): ", paste(setdiff(want, names(GATES)), collapse = ", "), "; valid: ", paste(names(GATES), collapse = ", "))
if (!dir.exists(lib)) message("note: no renv library at ", lib, " (run renv::restore() first, or pass --library DIR)")

fingerprint <- function(dir) { f <- list.files(dir, recursive = TRUE, all.files = TRUE); f <- f[!grepl("^(\\.git|renv/library|renv/staging)(/|$)", f)]
  setNames(unname(tools::md5sum(file.path(dir, f))), f) }
copy_repo <- function(from, to) { dir.create(to, recursive = TRUE)
  for (e in setdiff(list.files(from, all.files = TRUE, no.. = TRUE), c(".git"))) {
    if (e == "renv") { dir.create(file.path(to, "renv")); for (r in setdiff(list.files(file.path(from, "renv"), all.files = TRUE, no.. = TRUE), c("library", "staging")))
        file.copy(file.path(from, "renv", r), file.path(to, "renv"), recursive = TRUE, copy.date = TRUE) }
    else file.copy(file.path(from, e), to, recursive = TRUE, copy.date = TRUE) }
  if (dir.exists(lib) && !isTRUE(file.symlink(normalizePath(lib), file.path(to, "renv", "library"))))
    file.copy(lib, file.path(to, "renv"), recursive = TRUE) }          # fallback where symlinks are not permitted (unverified path)

# the work area must outlive this R session (tempdir() is deleted on exit), or --keep would keep nothing
before <- fingerprint(root); rows <- list(); work <- file.path(dirname(tempdir()), paste0("pkpd02_verify_", Sys.getpid())); dir.create(work)
for (g in want) { G <- GATES[[g]]; tmp <- file.path(work, g); copy_repo(root, tmp); log <- file.path(work, paste0(g, ".log")); t0 <- Sys.time() - 1
  owd <- setwd(tmp)                                                        # cwd must be the project root so .Rprofile activates renv
  # renv exports RENV_PROJECT into this session; a child would inherit it and silently use THIS repository's project and library, so point it at the temporary copy
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(tmp, G$script)), stdout = log, stderr = log, env = paste0("RENV_PROJECT=", tmp)); setwd(owd)
  txt <- readLines(log, warn = FALSE); ok_text <- any(grepl(G$text, txt, fixed = TRUE))
  ev <- list.files(file.path(tmp, "evidence"), recursive = TRUE, full.names = TRUE); mt <- file.info(ev)$mtime > t0
  newm <- tools::md5sum(ev[mt]); oldm <- tools::md5sum(file.path(root, "evidence", sub("^.*/evidence/", "", ev[mt])))
  differ <- sub("^.*/evidence/", "", names(newm)[is.na(oldm) | newm != oldm])
  rows[[g]] <- data.frame(gate = g, expected_exit = G$exit, exit = status, expected_text = ok_text, evidence_rewritten = sum(mt), differs_from_committed = length(differ),
                          result = if (identical(as.integer(status), G$exit) && ok_text) "AS EXPECTED" else "UNEXPECTED", stringsAsFactors = FALSE)
  cat(sprintf("%-17s %-52s exit %d (expected %d)  %s\n", g, G$what, status, G$exit, rows[[g]]$result))
  if (length(differ)) cat(sprintf("                  evidence that differs from the committed copy (timestamps, timings or environment listings are expected): %s\n", paste(differ, collapse = ", ")))
  if (!identical(rows[[g]]$result, "AS EXPECTED")) { cat("  --- last log lines:\n"); cat(paste0("  ", tail(txt, 6)), sep = "\n") }
}
after <- fingerprint(root); intact <- identical(before, after)
cat(sprintf("\nrepository files unchanged by this run: %s\n", intact))
if (keep) cat("temporary copies and logs kept in:", work, "\n") else unlink(work, recursive = TRUE)
res <- do.call(rbind, rows); ok <- all(res$result == "AS EXPECTED") && intact
cat(if (ok) "ALL REQUESTED GATES BEHAVED AS RECORDED\n" else "AT LEAST ONE GATE DID NOT BEHAVE AS RECORDED\n"); quit(status = if (ok) 0L else 1L)
