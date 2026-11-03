#!/usr/bin/env Rscript
# S-2 preparation, decision 5: METADATA-ONLY feasibility audit of the Strathclyde event records for
# R-COV, R-RESET, coincident events, active infusions at reset, and eligibility effects.
# TARGET-FREE: the DV (concentration) column and all other non-needed columns are skipped AT READ TIME;
# no model is fitted; no concentration value is loaded, printed or used. Observation validity (numeric,
# positive DV) therefore cannot be assessed here (documented limitation; counts are upper bounds).
args <- commandArgs(FALSE); f <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(normalizePath(file.path(dirname(f), "..")))
suppressMessages({ library(jsonlite) })
source("R/pkpd02_solve.R"); source("R/pkpd02_events.R"); source("R/pkpd02_source_reader.R")
ev <- pkpd02_read_source_metadata(); audit_trail <- attr(ev, "trail")

# ---- 1. structure ---------------------------------------------------------------------------------------
ev$EPISODE <- pkpd02_episode_index(ev)
dur <- pkpd02_effective_duration(ev$AMT, ev$RATE, ev$EVID); ev$DUR <- dur$dur; ev$DUR_RULE <- dur$rule
ev$IS_DOSE <- ev$EVID %in% c(1L, 4L) & !is.na(ev$AMT) & ev$AMT > 0
struct <- list(rows = nrow(ev), patients = length(unique(ev$ID)), episodes = nrow(unique(ev[, c("ID", "EPISODE")])), evid_counts = as.list(table(ev$EVID)),
  patient_rows_contiguous = all(vapply(split(seq_len(nrow(ev)), ev$ID), function(i) all(diff(i) == 1), NA)),
  time_decreases_within_episode = sum(unlist(lapply(split(ev$TIME, paste(ev$ID, ev$EPISODE)), function(t) sum(diff(t) < 0)))), time_decreases_across_rows_of_a_patient = sum(unlist(lapply(split(ev$TIME, ev$ID), function(t) sum(diff(t) < 0)))),
  reset_rows_with_TIME_zero = sum(ev$EVID %in% c(3L, 4L) & ev$TIME == 0), first_row_TIME_zero_patients = sum(vapply(split(ev$TIME, ev$ID), function(t) t[1] == 0, NA)), TIME_is_episode_relative = all(ev$TIME[ev$EVID %in% c(3L, 4L)] == 0), na_counts = list(TIME = sum(is.na(ev$TIME)), TAD = sum(is.na(ev$TAD)), CRCA = sum(is.na(ev$CRCA)), AJBW = sum(is.na(ev$AJBW))),
  first_row_of_patient_is_reset = sum(vapply(split(ev$EVID, ev$ID), function(e) e[1] %in% c(3L, 4L), NA)),
  dose_rows_with_nonpositive_or_missing_amt = sum(ev$EVID %in% c(1L, 4L) & !ev$IS_DOSE), EVID3_rows = sum(ev$EVID == 3L),
  duration_rules = as.list(table(ev$DUR_RULE, useNA = "no")), recorded_duration_outside_0.45_0.55h = sum(ev$DUR_RULE == "as_recorded" & abs(ev$DUR - 0.5) > 0.05, na.rm = TRUE),
  RATE0_positive_dose_rows = sum(ev$IS_DOSE & (is.na(ev$RATE) | ev$RATE == 0)), RATE0_patients = sort(unique(ev$ID[ev$IS_DOSE & (is.na(ev$RATE) | ev$RATE == 0)])))

# ---- 2. R-COV: covariate changes carried only by observation rows ---------------------------------------------------
grp <- list()
for (key in unique(paste(ev$ID, ev$EPISODE))) { x <- ev[paste(ev$ID, ev$EPISODE) == key, ]; id <- x$ID[1]; prev <- NULL
  for (tt in unique(x$TIME)) { g <- x[x$TIME == tt, ]; cur <- c(g$CRCA[1], g$AJBW[1]); conflict <- nrow(unique(g[, c("CRCA", "AJBW")])) > 1
    chg <- !is.null(prev) && !anyNA(cur) && !anyNA(prev) && any(cur != prev)
    if (chg || conflict || anyNA(cur)) grp[[length(grp) + 1L]] <- data.frame(ID = id, EPISODE = max(g$EPISODE), TIME = tt, n_rows = nrow(g), evids = paste(sort(g$EVID), collapse = "+"), covariate_change = chg,
        carried_by_restart_row = any(g$EVID != 0L), obs_only_change = chg && all(g$EVID == 0L), conflicting_tied_covariates = conflict, na_covariate = anyNA(cur), SRC_ROW_first = min(g$SRC_ROW))
    if (!anyNA(cur)) prev <- cur } }
grp <- do.call(rbind, grp)
rcov <- list(covariate_change_groups = sum(grp$covariate_change), carried_by_restart_row = sum(grp$covariate_change & grp$carried_by_restart_row), obs_only_changes = sum(grp$obs_only_change),
  patients_with_obs_only_change = length(unique(grp$ID[grp$obs_only_change])), episodes_with_obs_only_change = nrow(unique(grp[grp$obs_only_change, c("ID", "EPISODE")])),
  conflicting_tied_covariates = sum(grp$conflicting_tied_covariates), patients_conflicting_ties = length(unique(grp$ID[grp$conflicting_tied_covariates])), na_covariate_groups = sum(grp$na_covariate))

# ---- 3. coincident events ----------------------------------------------------------------------------------------------------
co <- list()
for (key in unique(paste(ev$ID, ev$EPISODE))) { x <- ev[paste(ev$ID, ev$EPISODE) == key, ]; id <- x$ID[1]; for (tt in unique(x$TIME[duplicated(x$TIME)])) { g <- x[x$TIME == tt, ]
  co[[length(co) + 1L]] <- data.frame(ID = id, EPISODE = max(g$EPISODE), TIME = tt, n_rows = nrow(g), evids = paste(g$EVID, collapse = "+"),
    sample_and_dose = any(g$EVID == 0L) && any(g$IS_DOSE), dose_and_dose = sum(g$IS_DOSE) > 1, sample_and_sample = sum(g$EVID == 0L) > 1, reset_and_sample = any(g$EVID %in% c(3L, 4L)) && any(g$EVID == 0L),
    reset_and_other_dose = any(g$EVID %in% c(3L, 4L)) && sum(g$IS_DOSE & g$EVID != 4L) > 0, exact_duplicate_rows = nrow(g) > nrow(unique(g[, c("EVID", "AMT", "RATE", "TAD", "CRCA", "AJBW")])), SRC_ROW_first = min(g$SRC_ROW)) } }
co <- if (length(co)) do.call(rbind, co) else data.frame()
coincident <- list(time_groups_with_multiple_rows = nrow(co), sample_and_dose = sum(co$sample_and_dose), dose_and_dose = sum(co$dose_and_dose), sample_and_sample = sum(co$sample_and_sample),
  reset_and_sample = sum(co$reset_and_sample), reset_and_other_dose = sum(co$reset_and_other_dose), exact_duplicate_row_groups = sum(co$exact_duplicate_rows),
  patients_any = length(unique(co$ID)), sample_dose_ties_TAD_gt0 = NA_integer_)
sd <- co[co$sample_and_dose, ]; coincident$sample_dose_ties_detail <- lapply(seq_len(nrow(sd)), function(i) { g <- ev[ev$ID == sd$ID[i] & ev$EPISODE == sd$EPISODE[i] & ev$TIME == sd$TIME[i], ]; list(ID = sd$ID[i], EPISODE = sd$EPISODE[i], TIME = sd$TIME[i], evids_in_file_order = g$EVID, TAD_of_sample = g$TAD[g$EVID == 0L], file_order_sample_first = which(g$EVID == 0L)[1] < which(g$IS_DOSE)[1], TAD_gt0_but_dose_listed_first = g$TAD[g$EVID == 0L][1] > 0 && which(g$EVID == 0L)[1] > which(g$IS_DOSE)[1]) })

# ---- 4. R-RESET. TIME restarts at every reset (episode-relative clock), so "an infusion active at the reset" CANNOT be
# assessed from recorded times: the wall-clock gap between episodes is unknown. What CAN be determined from metadata:
#   determinately-not-active : every infusion of the previous episode (v1.1 T7/T8 durations) ended by that episode's LAST RECORD
#                              (a reset cannot precede the last record), so no infusion can be interrupted;
#   indeterminate            : some infusion extends beyond the last record; interruption is possible but not evidenced.
rr <- list(); rs <- which(ev$EVID %in% c(3L, 4L))
for (j in rs) { id <- ev$ID[j]; pe <- ev$EPISODE[j] - 1L; x <- ev[ev$ID == id & ev$EPISODE == pe, ]
  d <- x[x$IS_DOSE, ]; last_rec <- max(x$TIME); end_eff <- if (nrow(d)) d$TIME + d$DUR else -Inf
  end_rec <- if (nrow(d)) ifelse(d$DUR_RULE == "as_recorded", d$TIME + d$DUR, d$TIME) else -Inf
  rr[[length(rr) + 1L]] <- data.frame(ID = id, EPISODE_new = ev$EPISODE[j], SRC_ROW = ev$SRC_ROW[j], EVID = ev$EVID[j], reset_TIME = ev$TIME[j], prev_episode_last_record_time = last_rec,
    prev_episode_last_record_is_dose = x$IS_DOSE[nrow(x)], prev_episode_n_doses = nrow(d), max_infusion_end_effective = max(end_eff), slack_h_last_record_minus_max_infusion_end = last_rec - max(end_eff),
    determinately_not_active = all(end_eff <= last_rec + 1e-12), indeterminate_possible_interruption = any(end_eff > last_rec + 1e-12),
    indeterminate_under_recorded_only_durations = any(end_rec > last_rec + 1e-12)) }
rr <- do.call(rbind, rr)
ties <- pkpd02_classify_reset_ties(ev, "episode_relative")
reset <- list(reset_rows = nrow(rr), EVID3 = sum(rr$EVID == 3L), EVID4 = sum(rr$EVID == 4L), R_RESET_assessable_from_recorded_times = FALSE,
  reason_not_assessable = "TIME is episode-relative (all resets at TIME 0); wall-clock gap between episodes is not recorded",
  determinately_not_active_at_reset = sum(rr$determinately_not_active), indeterminate_possible_interruption = sum(rr$indeterminate_possible_interruption),
  indeterminate_under_recorded_only_durations = sum(rr$indeterminate_under_recorded_only_durations), patients_with_indeterminate = length(unique(rr$ID[rr$indeterminate_possible_interruption])),
  previous_episode_ends_with_a_dose_record = sum(rr$prev_episode_last_record_is_dose),
  reset_tied_samples_episode_relative = nrow(ties), reset_tied_by_class = as.list(table(factor(ties$classification, levels = c("pre_reset", "post_reset", "indeterminate")))), patients_with_reset_tied_samples = length(unique(ties$ID)),
  note = "An episode-relative clock means an earlier infusion can never be tied to the next episode's clock; the representation (not the data) decides whether R-RESET can arise.")

# ---- 5. the sanctioned pipeline applied to the real event STRUCTURE (metadata only) ---------------------------------------------------------
tab <- ev; tab$AMT[!tab$IS_DOSE] <- 0; tab$RATE <- ifelse(tab$IS_DOSE, ifelse(tab$DUR_RULE == "as_recorded", ev$RATE, tab$AMT / tab$DUR), 0)
chk <- function(e) { v <- pkpd02_validate_events(e[, PKPD02_EV_COLS]); c(R_COV = sum(grepl("R-COV", v)), R_RESET = sum(grepl("R-RESET", v)), tied_rows_conflicting_covariates = sum(grepl("tied rows carry different", v)),
                                                                              time_not_monotone = sum(grepl("not non-decreasing", v)), other = sum(!grepl("R-COV|R-RESET|tied rows carry different|not non-decreasing", v)), total = length(v)) }
tl <- pkpd02_episode_timelines(tab, ties, "episode_relative"); cz_tl <- pkpd02_canonicalize_events(tl$timelines)
pl <- pkpd02_patient_level_table(tl$timelines, gap = 1); cz_pl <- pkpd02_canonicalize_events(pl)
pipeline <- list(
  naive_patient_level_table_as_recorded = as.list(chk(tab)),
  episode_timelines = list(n = length(unique(tl$timelines$SOLVE_ID)), quarantined_observations = nrow(tl$quarantined), canonicalisation_inserted_rows = nrow(cz_tl$inserted),
    canonicalisation_rejected_groups = nrow(cz_tl$rejected), validator = as.list(chk(cz_tl$events)), patients_with_multiple_episodes = sum(table(unique(tl$timelines[, c("PATIENT_ID", "EPISODE")])$PATIENT_ID) > 1)),
  patient_level_offset_table_gap1h = list(rows = nrow(cz_pl$events), canonicalisation_inserted_rows = nrow(cz_pl$inserted), canonicalisation_rejected_groups = nrow(cz_pl$rejected), validator = as.list(chk(cz_pl$events)),
    note = "episode offsets are a representation device with no clinical meaning; invariance to the gap is verified in tests/s2prep_canonical_equivalence.R"))
tab_inserted <- cz_tl$inserted

# ---- 6. eligibility (target-free; DV validity not assessable -> upper bounds) -------------------------------------------------------------------
eligible <- function(d, drop_rows = integer(0), mode = "one") {
  d <- d[!(d$SRC_ROW %in% drop_rows), ]; out <- list()
  for (key in unique(paste(d$ID, d$EPISODE))) { x <- d[paste(d$ID, d$EPISODE) == key, ]; ds <- x$TIME[x$IS_DOSE]; ob <- x[x$EVID == 0L, ]
    if (!nrow(ob) || !length(ds)) next
    usable <- ob$TIME > min(ds) & vapply(ob$TIME, function(t) any(ds < t), NA)
    if (mode == "one") { if (!any(usable)) next; i0 <- which(usable)[1]; t0 <- ob$TIME[i0]
      tg <- which(seq_len(nrow(ob)) > i0 & ob$TIME > t0 & vapply(ob$TIME, function(t) any(ds > t0 & ds < t), NA)); if (length(tg)) out[[length(out) + 1L]] <- data.frame(ID = x$ID[1], EPISODE = x$EPISODE[1], origin_row = ob$SRC_ROW[i0], target_row = ob$SRC_ROW[tg[1]]) }
    else { ok <- NULL; for (k in seq_len(nrow(ob))) { prior <- which(usable & ob$TIME < ob$TIME[k] & seq_len(nrow(ob)) < k); if (length(prior) >= 2) { b <- prior[length(prior)]
        if (any(ds > ob$TIME[b] & ds < ob$TIME[k])) { ok <- k; break } } }
      if (!is.null(ok)) out[[length(out) + 1L]] <- data.frame(ID = x$ID[1], EPISODE = x$EPISODE[1], origin_row = NA, target_row = ob$SRC_ROW[ok]) } }
  if (length(out)) do.call(rbind, out) else data.frame(ID = integer(), EPISODE = integer(), origin_row = integer(), target_row = integer()) }
cnt <- function(e) list(episodes = nrow(e), patients = length(unique(e$ID)))
e1 <- eligible(ev, mode = "one"); e2 <- eligible(ev, mode = "two")
pre <- fromJSON(Sys.getenv("PKPD02_PREFLIGHT_JSON", "PKPD02_v11_metadata_preflight.json"))
repro <- list(one_prior = cnt(e1), two_prior = cnt(e2), v11_preflight = list(one_prior = list(episodes = pre$one_prior_eligible_episodes, patients = pre$one_prior_eligible_patients),
  two_prior = list(episodes = pre$two_prior_eligible_episodes, patients = pre$two_prior_eligible_patients)),
  matches_v11_one_prior = identical(c(nrow(e1), length(unique(e1$ID))), c(pre$one_prior_eligible_episodes, pre$one_prior_eligible_patients)),
  matches_v11_two_prior = identical(c(nrow(e2), length(unique(e2$ID))), c(pre$two_prior_eligible_episodes, pre$two_prior_eligible_patients)))
q_rows <- tl$quarantined$SRC_ROW
scen <- list(S0_as_recorded = list(one = cnt(e1), two = cnt(e2)),
  S1_quarantine_indeterminate_reset_tied_samples = list(one = cnt(eligible(ev, q_rows, "one")), two = cnt(eligible(ev, q_rows, "two")), quarantined_rows = length(q_rows)),
  S2_exclude_episodes_with_indeterminate_possible_interruption_if_resets_were_treated_as_unsafe = { bad <- rr[rr$indeterminate_possible_interruption, ]; keyb <- paste(bad$ID, bad$EPISODE_new - 1L)
     dd <- ev[!(paste(ev$ID, ev$EPISODE) %in% keyb), ]; list(one = cnt(eligible(dd, q_rows, "one")), two = cnt(eligible(dd, q_rows, "two")), excluded_episodes = nrow(bad)) },
  S3_exclude_episodes_with_obs_only_covariate_change_if_not_canonicalised = { keyc <- unique(paste(grp$ID[grp$obs_only_change], grp$EPISODE[grp$obs_only_change])); dd <- ev[!(paste(ev$ID, ev$EPISODE) %in% keyc), ]
     list(one = cnt(eligible(dd, q_rows, "one")), two = cnt(eligible(dd, q_rows, "two")), excluded_episodes = length(keyc)) })
used_rows <- c(e1$origin_row, e1$target_row, e2$target_row); affected_primary <- list(
  primary_one_prior_forecasts_whose_origin_or_target_is_a_quarantined_sample = sum(e1$origin_row %in% q_rows | e1$target_row %in% q_rows),
  primary_one_prior_episodes_with_obs_only_covariate_change = sum(paste(e1$ID, e1$EPISODE) %in% paste(grp$ID[grp$obs_only_change], grp$EPISODE[grp$obs_only_change])),
  primary_one_prior_episodes_followed_by_a_reset_with_indeterminate_interruption = sum(paste(e1$ID, e1$EPISODE) %in% paste(rr$ID[rr$indeterminate_possible_interruption], rr$EPISODE_new[rr$indeterminate_possible_interruption] - 1L)))
first_per_patient <- function(e) e[!duplicated(e$ID), ]; scen$first_eligible_episode_per_patient_S0_one_prior <- cnt(first_per_patient(e1))
# one-prior example rows from the v1.1 preflight (metadata only: patient, origin row, target row) as an algorithm cross-check
ex <- pre$one_prior_examples; ex_ok <- vapply(seq_len(nrow(ex)), function(i) any(e1$ID == ex[i, 1] & e1$origin_row + 2L == ex[i, 2] & e1$target_row + 2L == ex[i, 3]), NA)   # preflight rows are Excel row numbers (2 header rows)
repro$preflight_one_prior_example_rows_reproduced <- sum(ex_ok); repro$preflight_one_prior_example_rows_total <- length(ex_ok)
repro$preflight_examples_vs_this_audit_origin_target_rows <- NULL

# ---- evidence ---------------------------------------------------------------------------------------------------------------------------------
dir.create("evidence", FALSE); w <- function(x, f) write.csv(x, file.path("evidence", f), row.names = FALSE)
w(grp, "s2prep_audit_covariate_change_groups.csv"); w(co, "s2prep_audit_coincident_event_groups.csv"); w(rr, "s2prep_audit_reset_rows.csv"); w(ties, "s2prep_audit_reset_tied_samples.csv")
w(if (nrow(tl$quarantined)) tl$quarantined[, c("ID", "TIME", "TAD", "EVID", "SRC_ROW")] else data.frame(), "s2prep_audit_quarantined_observations.csv")
w(e1, "s2prep_audit_eligible_one_prior_rows.csv"); w(e2, "s2prep_audit_eligible_two_prior_rows.csv")
ep <- do.call(rbind, lapply(split(ev, paste(ev$ID, ev$EPISODE)), function(x) data.frame(ID = x$ID[1], EPISODE = x$EPISODE[1], n_rows = nrow(x), n_doses = sum(x$IS_DOSE), n_obs = sum(x$EVID == 0L), start_time = min(x$TIME), end_time = max(x$TIME))))
w(ep, "s2prep_audit_episode_summary.csv")
res <- list(audit = "S-2 prep metadata-only feasibility audit", timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), audit_trail = audit_trail, structure = struct, R_COV = rcov, coincident_events = coincident, R_RESET = reset,
  sanctioned_pipeline_on_real_event_structure = pipeline, eligibility_reproduction = repro, eligibility_scenarios = scen, forecasts_touching_flagged_events = affected_primary,
  script_sha256 = digest::digest("tests/s2prep_metadata_audit.R", "sha256", file = TRUE), library_sha256 = list(events = digest::digest("R/pkpd02_events.R", "sha256", file = TRUE), solve = digest::digest("R/pkpd02_solve.R", "sha256", file = TRUE)),
  limitations = c("Observation validity (numeric positive DV) not assessed: DV never read; eligibility counts are upper bounds on that dimension.", "Effective infusion durations apply v1.1 T7/T8; recorded-only sensitivity reported for R-RESET.",
                  "SEX/AGE/WT/IBW/HT/SCREAT/CRCL not read; covariate-change analysis covers only the model covariates CRCA and AJBW."))
writeLines(toJSON(res, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"), "evidence/s2prep_metadata_audit.json")
cat(toJSON(res[c("structure", "R_COV", "coincident_events", "R_RESET", "sanctioned_pipeline_on_real_event_structure", "eligibility_reproduction", "eligibility_scenarios", "forecasts_touching_flagged_events")], auto_unbox = TRUE, pretty = TRUE, na = "null"))
