#!/usr/bin/env Rscript
# Runs the S-2 preparation scripts and the S-1r regression, verifies ALL preserved archives, and records
# sessionInfo, package versions, seeds, tolerances and file hashes. Exit 0 only if everything passes.
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)]); setwd(normalizePath(file.path(dirname(f), ".."))); suppressMessages(library(jsonlite))
sha <- function(p) digest::digest(p, "sha256", file = TRUE)
arch <- list.dirs("archive", recursive = FALSE); arch_ok <- vapply(arch, function(a) { m <- read.table(file.path(a, "MANIFEST.sha256"), col.names = c("sha", "path"), stringsAsFactors = FALSE); all(vapply(seq_len(nrow(m)), function(i) sha(file.path(a, m$path[i])) == m$sha[i], NA)) }, NA)
run <- function(s) system2("Rscript", file.path("tests", s), stdout = file.path(tempdir(), paste0(s, ".log")), stderr = file.path(tempdir(), paste0(s, ".log")))
codes <- c(canonical_equivalence = run("s2prep_canonical_equivalence.R"), metadata_audit = run("s2prep_metadata_audit.R"), real_structure_domain_check = run("s2prep_real_structure_domain_check.R"),
           static_schema_inspection = run("s2_static_schema_inspection.R"), s1r_regression = run("s1r_run_all.R"))
eq <- fromJSON("evidence/s2prep_canonical_equivalence.json"); dm <- fromJSON("evidence/s2prep_real_structure_domain_check.json"); au <- fromJSON("evidence/s2prep_metadata_audit.json")
pk <- c("rxode2", "nlmixr2est", "nlmixr2", "lotri", "dparser", "rxode2lincmt", "symengine", "expm", "readxl", "jsonlite", "digest"); pkv <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")
overall <- all(codes == 0L) && all(arch_ok)
files <- c(list.files("R", full.names = TRUE), list.files("tests", pattern = "^s2", full.names = TRUE), "protocol/PKPD02_v12_amendment_DRAFT.md", "reports/S2_synthetic_test_design.md", "reports/S2prep_report.md", "renv.lock")
out <- list(stage = "S-2 preparation (no fitting)", overall = if (overall) "PASS" else "FAIL", exit_codes = as.list(codes), archives_verified = as.list(setNames(arch_ok, basename(arch))),
  canonical_equivalence = eq$verdict, real_structure_domain_check = dm$verdict, metadata_audit_DV_loaded = au$audit_trail$DV_column_loaded,
  seeds = list(canonical_equivalence = eq$seed, real_structure_domain_check = dm$seed, failed_S1_generator = 20261101, reference_selftest = 20261102, note = "S-2 static inspection uses archived synthetic scenarios; no RNG"),
  tolerances = list(acceptance = "max(1e-5 mg/L, 1e-6*|C_ref|)", ode_reference = list(atol = 1e-12, rtol = 1e-12), production_candidate = list(atol = 1e-10, rtol = 1e-10, status = "proposal only")),
  fitting_performed = FALSE, strathclyde_concentrations_read = FALSE, forecast_targets_constructed = FALSE, prefit_manifest_created = FALSE, git_push = FALSE,
  package_versions = as.list(pkv), R_version = R.version.string, platform = R.version$platform, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), file_sha256 = setNames(as.list(vapply(files[file.exists(files)], sha, "")), files[file.exists(files)]))
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), "evidence/s2prep_results.json")
writeLines(capture.output(print(sessionInfo()), cat("\nBLAS/LAPACK:\n"), print(extSoftVersion()), cat("\nLa_library:", La_library(), "\nLa_version:", La_version(), "\n")), "evidence/s2prep_sessionInfo.txt")
write.csv(data.frame(package = names(pkv), installed_version = pkv), "evidence/s2prep_package_versions.csv", row.names = FALSE)
for (s in names(codes)) cat(sprintf("%-30s exit %d\n", s, codes[[s]])); cat("archives verified:", paste(basename(arch), arch_ok, collapse = "; "), "\nS-2 PREP OVERALL:", out$overall, "\n"); quit(status = if (overall) 0L else 1L)
