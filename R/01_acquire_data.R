# Acquire and validate the Prostate Cancer Study data.
#
# Run from the project root:
#   source("R/01_acquire_data.R")
#
# Inputs:
#   lbreg::PCS (public data distributed with the lbreg package)
#
# Outputs (ignored by Git):
#   data/processed/prostate_modern.rds
#   data/processed/prostate_2018.rds
#   data/processed/acquisition_metadata.txt

check_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

project_markers <- c(
  "prostatic-capsular-penetration.Rproj",
  "renv.lock"
)

check_true(
  all(file.exists(project_markers)),
  paste(
    "Run this script from the project root.",
    "Expected to find:",
    paste(project_markers, collapse = ", ")
  )
)

check_true(
  requireNamespace("lbreg", quietly = TRUE),
  paste(
    "Package 'lbreg' is required but not installed.",
    "Restore the project environment with renv::restore()."
  )
)

data("PCS", package = "lbreg", envir = environment())

check_true(exists("PCS", inherits = FALSE), "Could not load lbreg::PCS.")
check_true(is.data.frame(PCS), "lbreg::PCS must be a data frame.")

expected_names <- c(
  "id", "tumor", "age", "race", "dpros",
  "dcaps", "psa", "vol", "gleason"
)

check_true(nrow(PCS) == 380L, "Expected 380 observations in lbreg::PCS.")
check_true(ncol(PCS) == 9L, "Expected 9 variables in lbreg::PCS.")
check_true(
  identical(names(PCS), expected_names),
  paste(
    "Unexpected variable names or order in lbreg::PCS.",
    "Expected:",
    paste(expected_names, collapse = ", ")
  )
)

check_true(!anyDuplicated(PCS$id), "Record identifiers must be unique.")
check_true(
  identical(PCS$id, seq_len(380L)),
  "Record identifiers must be the integers 1 through 380 in order."
)

expected_outcome_counts <- c(`0` = 227L, `1` = 153L)
observed_outcome_counts <- table(PCS$tumor, useNA = "ifany")

check_true(
  identical(names(observed_outcome_counts), names(expected_outcome_counts)) &&
    identical(
      as.integer(observed_outcome_counts),
      unname(expected_outcome_counts)
    ),
  "Unexpected tumor-penetration outcome counts."
)

expected_missing_counts <- c(
  id = 0L,
  tumor = 0L,
  age = 0L,
  race = 3L,
  dpros = 0L,
  dcaps = 0L,
  psa = 0L,
  vol = 0L,
  gleason = 0L
)
observed_missing_counts <- colSums(is.na(PCS))

check_true(
  identical(names(observed_missing_counts), names(expected_missing_counts)) &&
    identical(
      as.integer(observed_missing_counts),
      unname(expected_missing_counts)
    ),
  "Unexpected missing-value counts in lbreg::PCS."
)

expected_missing_race_rows <- c(22L, 46L, 252L)
observed_missing_race_rows <- which(is.na(PCS$race))

check_true(
  identical(observed_missing_race_rows, expected_missing_race_rows),
  "Unexpected locations for missing race values."
)

check_true(
  identical(sort(unique(PCS$tumor)), 0:1),
  "The outcome must contain only 0 and 1."
)
check_true(
  identical(sort(unique(stats::na.omit(PCS$race))), 1:2),
  "Race must contain only codes 1 and 2, apart from missing values."
)
check_true(
  identical(sort(unique(PCS$dpros)), 1:4),
  "Digital rectal examination results must contain codes 1 through 4."
)
check_true(
  identical(sort(unique(PCS$dcaps)), 1:2),
  "Capsular-involvement examination results must contain codes 1 and 2."
)
check_true(
  all(PCS$gleason >= 0L & PCS$gleason <= 10L),
  "Gleason scores must be between 0 and 10."
)

prostate_modern <- PCS
names(prostate_modern)[names(prostate_modern) == "tumor"] <- "capsule"

historical_names <- c(
  "id", "capsule", "age", "race",
  "dpros", "dcaps", "psa", "gleason"
)
prostate_2018 <- prostate_modern[historical_names]

check_true(
  identical(prostate_modern$capsule, PCS$tumor),
  "Outcome values changed unexpectedly during renaming."
)
check_true(
  identical(prostate_2018$id, PCS$id),
  "Observation order changed while deriving the 2018 dataset."
)
check_true(
  identical(names(prostate_2018), historical_names),
  "The derived 2018 dataset has unexpected variables."
)
check_true(
  nrow(prostate_modern) == 380L && nrow(prostate_2018) == 380L,
  "Acquisition must retain all 380 observations."
)

output_dir <- file.path("data", "processed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

modern_path <- file.path(output_dir, "prostate_modern.rds")
historical_path <- file.path(output_dir, "prostate_2018.rds")
metadata_path <- file.path(output_dir, "acquisition_metadata.txt")

saveRDS(prostate_modern, modern_path, version = 3)
saveRDS(prostate_2018, historical_path, version = 3)

output_md5 <- tools::md5sum(c(modern_path, historical_path))

metadata <- c(
  "Prostate Cancer Study data acquisition",
  paste("Acquisition date:", Sys.Date()),
  "Source dataset: PCS",
  "Source package: lbreg",
  paste("Source package version:", as.character(utils::packageVersion("lbreg"))),
  paste("R version:", R.version.string),
  "Canonical dimensions: 380 rows x 9 columns",
  "Historical dimensions: 380 rows x 8 columns",
  "Outcome counts: 0 = 227; 1 = 153",
  "Missing race rows: 22, 46, 252",
  "Modern transformation: tumor renamed to capsule; vol retained",
  "Historical transformation: tumor renamed to capsule; vol excluded",
  paste("MD5", names(output_md5), unname(output_md5), sep = ": ")
)

writeLines(metadata, metadata_path, useBytes = TRUE)

message("Data acquisition and validation completed successfully.")
message("Created: ", modern_path)
message("Created: ", historical_path)
message("Created: ", metadata_path)
