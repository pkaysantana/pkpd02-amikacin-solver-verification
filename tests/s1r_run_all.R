#!/usr/bin/env Rscript
# Runs every S-1r script, verifies the preserved failed-S-1 artefacts are unchanged, and writes the
# combined machine-readable result (evidence/s1r_results.json), sessionInfo and package versions.
# Exit 0 = PASS (conditional on the input contract recorded in the JSON); 1 = FAIL.
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(normalizePath(file.path(dirname(f), ".."))); suppressMessages(library(jsonlite))
ARCH <- "archive/s1_failed_2026-10-30"
man <- read.table(file.path(ARCH, "MANIFEST.sha256"), col.names = c("sha", "path"), stringsAsFactors = FALSE)
sha <- function(p) digest::digest(p, "sha256", file = TRUE)
arch_ok <- all(vapply(seq_len(nrow(man)), function(i) sha(file.path(ARCH, man$path[i])) == man$sha[i], NA))
# original (in-place) failed-S-1 files must still equal the archived copies
orig <- grep("^\\./(tests/s1_|reports/S1_|evidence/s1_)", man$path, value = TRUE)
orig_ok <- all(vapply(orig, function(p) sha(sub("^\\./", "", p)) == man$sha[man$path == p], NA))
run <- function(s) { st <- system2("Rscript", file.path("tests", s), stdout = file.path(tempdir(), paste0(s, ".log")), stderr = file.path(tempdir(), paste0(s, ".log"))); st }
codes <- c(reference_selftest = run("s1r_reference_selftest.R"), pathway_inspection = run("s1r_pathway_inspection.R"), ode_equivalence = run("s1r_ode_equivalence.R"))
eq <- fromJSON("evidence/s1r_equivalence_results.json")
pk <- c("rxode2", "nlmixr2est", "nlmixr2", "lotri", "dparser", "rxode2lincmt", "symengine", "PreciseSums", "n1qn1", "lbfgsb3c", "expm", "jsonlite", "digest")
pkv <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")
rl <- fromJSON("renv.lock", simplifyVector = FALSE)
src <- lapply(setNames(pk, pk), function(p) { e <- rl$Packages[[p]]; if (is.null(e)) NULL else e[intersect(names(e), c("Version", "Source", "Repository", "Hash", "RemoteSha"))] })
overall <- all(codes == 0L) && arch_ok && orig_ok
files <- c(list.files("R", full.names = TRUE), list.files("tests", pattern = "^s1r_", full.names = TRUE), "protocol/PKPD02_v12_amendment_DRAFT.md", "reports/S1r_verification_report.md", "renv.lock",
           "protocol/PKPD02_PreSpecified_Exploratory_Protocol_v0_9.md", "protocol/PKPD02_v11_addendum.md")
out <- list(gate = "S-1r", overall_verdict = if (overall) "PASS_CONDITIONAL_ON_INPUT_CONTRACT" else "FAIL", exit_codes = as.list(codes),
  failed_S1_artefacts_unmodified = list(archive_manifest_verified = arch_ok, in_place_originals_match_archive = orig_ok, manifest_sha256 = sha(file.path(ARCH, "MANIFEST.sha256"))),
  equivalence_verdict = eq$verdict, conditions = eq$conditions, input_contract = eq$input_contract,
  tolerances = list(acceptance = eq$tolerance_rule, primary_ode = eq$primary_ode_settings, production_assessment = eq$production_tolerance_assessment),
  seeds = list(failed_S1_scenario_generator = 20261101, selftest = 20261102, note = "S-1r uses the archived event tables (SHA-verified); no random generation in S-1r scripts"),
  data_scope = "synthetic only; no Strathclyde data, no fitting, no forecast targets",
  package_versions = as.list(pkv), renv_lock_entries = src, R_version = R.version.string, platform = R.version$platform,
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), file_sha256 = setNames(as.list(vapply(files[file.exists(files)], sha, "")), files[file.exists(files)]))
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), "evidence/s1r_results.json")
writeLines(capture.output(print(sessionInfo()), cat("\nBLAS/LAPACK:\n"), print(extSoftVersion()), cat("\nLa_library:", La_library(), "\nLa_version:", La_version(), "\n")), "evidence/s1r_sessionInfo.txt")
write.csv(data.frame(package = names(pkv), installed_version = pkv), "evidence/s1r_package_versions.csv", row.names = FALSE)
for (s in names(codes)) cat(sprintf("%-20s exit %d\n", s, codes[[s]]))
cat("failed-S-1 archive verified:", arch_ok, " originals unmodified:", orig_ok, "\nS-1r OVERALL:", out$overall_verdict, "\n")
quit(status = if (overall) 0L else 1L)
