# Reads ONLY event/covariate metadata from the Strathclyde workbook. The DV (concentration) column and all
# other columns are skipped at read time; the function stops if DV were ever loaded. SHA-256 is verified.
pkpd02_read_source_metadata <- function(xlsx = Sys.getenv("PKPD02_SOURCE_XLSX", "Amikacin_Data_File.xlsx"),
                                        expect_sha = "521daef9e7d753419006259ed15e3ac4f26a1345c290c79c1cacf42e11c767cb") {
  suppressMessages(library(readxl)); sha <- digest::digest(xlsx, "sha256", file = TRUE); stopifnot(identical(sha, expect_sha))
  hdr <- as.character(unlist(read_excel(xlsx, sheet = 1, range = cell_rows(2), col_names = FALSE, col_types = "text", .name_repair = "minimal")))
  keep <- c("#ID", "TIME", "TAD", "EVID", "AMT", "RATE", "AJBW", "CRCA")
  raw <- suppressMessages(read_excel(xlsx, sheet = 1, skip = 2, col_names = hdr[hdr %in% keep], col_types = ifelse(hdr %in% keep, "text", "skip")))
  stopifnot(!"DV" %in% names(raw), ncol(raw) == length(keep))
  num <- function(v) suppressWarnings(as.numeric(v))
  ev <- data.frame(ID = as.integer(num(raw[["#ID"]])), TIME = num(raw$TIME), TAD = num(raw$TAD), EVID = as.integer(num(raw$EVID)), AMT = num(raw$AMT), RATE = num(raw$RATE), CRCA = num(raw$CRCA), AJBW = num(raw$AJBW))
  ev$SRC_ROW <- seq_len(nrow(ev)); ev$CMT <- 1L
  attr(ev, "trail") <- list(source_file = basename(xlsx), sha256 = sha, sha256_matches_protocol = TRUE, columns_loaded = names(raw), DV_column_loaded = FALSE, other_columns_skipped_at_read = setdiff(hdr, keep),
    n_rows_loaded = nrow(ev), nonnumeric_token_counts = as.list(vapply(raw, function(v) sum(!is.na(v) & is.na(num(v))), 1L)),
    statement = "No concentration (DV) value was read, printed, summarised or used; no model was fitted; no forecast target was constructed.")
  ev
}
