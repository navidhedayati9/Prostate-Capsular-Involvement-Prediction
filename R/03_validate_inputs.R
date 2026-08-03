# Validate inputs for the prespecified modern analysis.
#
# Run from the project root after acquiring the data:
#   source("R/01_acquire_data.R")
#   source("R/03_validate_inputs.R")
#
# Input (generated and ignored by Git):
#   data/processed/prostate_modern.rds
#   data/processed/acquisition_metadata.txt
#
# Output:
#   results/modern/input_validation.csv
#
# This script validates data only. It does not fit or evaluate models.

check_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

project_markers <- c(
  "prostatic-capsular-penetration.Rproj",
  "renv.lock",
  "report/modern-analysis-plan.md",
  "report/sample-size-assessment.md"
)

check_true(
  all(file.exists(project_markers)),
  paste(
    "Run this script from the project root.",
    "Expected to find:",
    paste(project_markers, collapse = ", ")
  )
)

input_path <- file.path("data", "processed", "prostate_modern.rds")
metadata_path <- file.path(
  "data",
  "processed",
  "acquisition_metadata.txt"
)
output_dir <- file.path("results", "modern")
output_path <- file.path(output_dir, "input_validation.csv")

check_true(
  file.exists(input_path),
  paste(
    "Modern dataset not found:",
    input_path,
    "Run source(\"R/01_acquire_data.R\") first."
  )
)
check_true(
  file.exists(metadata_path),
  paste(
    "Acquisition metadata not found:",
    metadata_path,
    "Run source(\"R/01_acquire_data.R\") first."
  )
)

modern_data <- tryCatch(
  readRDS(input_path),
  error = function(error) {
    stop(
      paste("Could not read", input_path, "-", conditionMessage(error)),
      call. = FALSE
    )
  }
)

metadata <- readLines(metadata_path, warn = FALSE)
expected_md5_line_prefix <- paste0("MD5: ", input_path, ": ")
md5_line <- metadata[startsWith(metadata, expected_md5_line_prefix)]

expected_md5 <- if (length(md5_line) == 1L) {
  sub(expected_md5_line_prefix, "", md5_line, fixed = TRUE)
} else {
  NA_character_
}
observed_md5 <- unname(tools::md5sum(input_path))

validation_rows <- list()

format_value <- function(value) {
  if (length(value) == 0L) {
    return("<none>")
  }
  if (length(value) > 12L) {
    value <- c(value[seq_len(12L)], "...")
  }
  paste(value, collapse = ", ")
}

add_check <- function(
  check_id,
  category,
  description,
  condition,
  observed,
  expected
) {
  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    check_id = check_id,
    category = category,
    description = description,
    status = if (isTRUE(condition)) "PASS" else "FAIL",
    observed = format_value(observed),
    expected = format_value(expected),
    stringsAsFactors = FALSE
  )
}

expected_names <- c(
  "id", "capsule", "age", "race", "dpros",
  "dcaps", "psa", "vol", "gleason"
)
expected_classes <- c(
  id = "integer",
  capsule = "integer",
  age = "integer",
  race = "integer",
  dpros = "integer",
  dcaps = "integer",
  psa = "numeric",
  vol = "numeric",
  gleason = "integer"
)

add_check(
  "file_md5",
  "integrity",
  "Dataset MD5 agrees with acquisition metadata.",
  !is.na(expected_md5) && identical(observed_md5, expected_md5),
  observed_md5,
  expected_md5
)
add_check(
  "data_frame",
  "structure",
  "Input is a data frame.",
  is.data.frame(modern_data),
  class(modern_data)[1],
  "data.frame"
)
add_check(
  "row_count",
  "structure",
  "Dataset contains the expected number of observations.",
  nrow(modern_data) == 380L,
  nrow(modern_data),
  380L
)
add_check(
  "column_count",
  "structure",
  "Dataset contains the expected number of variables.",
  ncol(modern_data) == 9L,
  ncol(modern_data),
  9L
)
add_check(
  "column_names",
  "structure",
  "Variable names and order match the modern data specification.",
  identical(names(modern_data), expected_names),
  names(modern_data),
  expected_names
)

structure_is_valid <- is.data.frame(modern_data) &&
  identical(names(modern_data), expected_names)

if (structure_is_valid) {
  observed_classes <- vapply(modern_data, class, character(1))

  add_check(
    "column_classes",
    "structure",
    "Variable storage classes match the acquisition specification.",
    identical(observed_classes, expected_classes),
    paste(names(observed_classes), observed_classes, sep = "="),
    paste(names(expected_classes), expected_classes, sep = "=")
  )
  add_check(
    "finite_numeric_values",
    "values",
    "All observed numeric values are finite.",
    all(vapply(
      modern_data,
      function(x) all(is.finite(x[!is.na(x)])),
      logical(1)
    )),
    "checked all variables",
    "all finite apart from documented NA values"
  )
  add_check(
    "unique_ids",
    "identifiers",
    "Record identifiers are unique.",
    !anyDuplicated(modern_data$id),
    length(unique(modern_data$id)),
    380L
  )
  add_check(
    "sequential_ids",
    "identifiers",
    "Record identifiers are the integers 1 through 380 in order.",
    identical(modern_data$id, seq_len(380L)),
    range(modern_data$id),
    c(1L, 380L)
  )

  outcome_counts <- table(modern_data$capsule, useNA = "always")
  add_check(
    "outcome_codes",
    "outcome",
    "Capsule penetration contains only binary codes 0 and 1.",
    identical(sort(unique(modern_data$capsule)), 0:1),
    sort(unique(modern_data$capsule)),
    0:1
  )
  add_check(
    "outcome_counts",
    "outcome",
    "Capsule-penetration counts match the verified source.",
    identical(
      as.integer(table(modern_data$capsule)),
      c(227L, 153L)
    ),
    paste(names(outcome_counts), outcome_counts, sep = "="),
    c("0=227", "1=153", "<NA>=0")
  )

  observed_missing <- colSums(is.na(modern_data))
  expected_missing <- c(
    id = 0L,
    capsule = 0L,
    age = 0L,
    race = 3L,
    dpros = 0L,
    dcaps = 0L,
    psa = 0L,
    vol = 0L,
    gleason = 0L
  )
  add_check(
    "missing_counts",
    "missingness",
    "Missing-value counts match the verified source.",
    identical(
      as.integer(observed_missing),
      unname(expected_missing)
    ),
    paste(names(observed_missing), observed_missing, sep = "="),
    paste(names(expected_missing), expected_missing, sep = "=")
  )
  add_check(
    "missing_race_ids",
    "missingness",
    "Race is missing only for the three documented records.",
    identical(
      modern_data$id[is.na(modern_data$race)],
      c(22L, 46L, 252L)
    ),
    modern_data$id[is.na(modern_data$race)],
    c(22L, 46L, 252L)
  )
  add_check(
    "primary_predictors_complete",
    "missingness",
    "All primary-model predictors are complete.",
    !anyNA(
      modern_data[c("age", "dpros", "psa", "gleason")]
    ),
    colSums(is.na(
      modern_data[c("age", "dpros", "psa", "gleason")]
    )),
    c(age = 0L, dpros = 0L, psa = 0L, gleason = 0L)
  )

  add_check(
    "race_codes",
    "codes",
    "Observed race values use only documented codes 1 and 2.",
    identical(sort(unique(stats::na.omit(modern_data$race))), 1:2),
    sort(unique(stats::na.omit(modern_data$race))),
    1:2
  )
  add_check(
    "dpros_codes",
    "codes",
    "Digital rectal examination uses all four documented codes.",
    identical(sort(unique(modern_data$dpros)), 1:4),
    sort(unique(modern_data$dpros)),
    1:4
  )
  add_check(
    "dcaps_codes",
    "codes",
    "Capsular-involvement examination uses documented codes 1 and 2.",
    identical(sort(unique(modern_data$dcaps)), 1:2),
    sort(unique(modern_data$dcaps)),
    1:2
  )

  add_check(
    "age_range",
    "ranges",
    "Age values match the verified historical range.",
    identical(range(modern_data$age), c(43L, 79L)),
    range(modern_data$age),
    c(43L, 79L)
  )
  add_check(
    "psa_positive",
    "transformations",
    "All PSA values are positive, permitting log2 transformation.",
    all(modern_data$psa > 0),
    range(modern_data$psa),
    "minimum > 0"
  )
  log2_psa <- log2(modern_data$psa)
  add_check(
    "log2_psa_finite",
    "transformations",
    "The planned log2(PSA) transformation produces finite values.",
    all(is.finite(log2_psa)),
    range(log2_psa),
    "all finite"
  )

  add_check(
    "gleason_range",
    "ranges",
    "Gleason scores lie within the documented range 0 through 10.",
    all(modern_data$gleason >= 0L & modern_data$gleason <= 10L),
    range(modern_data$gleason),
    c(0L, 10L)
  )
  add_check(
    "gleason_zero_records",
    "edge_cases",
    "The two prespecified zero Gleason scores occur at expected records.",
    identical(modern_data$id[modern_data$gleason == 0L], c(282L, 357L)),
    modern_data$id[modern_data$gleason == 0L],
    c(282L, 357L)
  )

  add_check(
    "volume_nonnegative",
    "ranges",
    "Tumor-volume values are nonnegative.",
    all(modern_data$vol >= 0),
    range(modern_data$vol),
    "minimum >= 0"
  )
  add_check(
    "volume_zero_count",
    "edge_cases",
    "The number of recorded zero tumor volumes matches the plan.",
    sum(modern_data$vol == 0) == 167L,
    sum(modern_data$vol == 0),
    167L
  )

  volume_positive <- as.integer(modern_data$vol > 0)
  volume_log_positive <- ifelse(
    modern_data$vol > 0,
    log2(modern_data$vol),
    0
  )
  add_check(
    "volume_transform_finite",
    "transformations",
    "The planned two-part tumor-volume representation is finite.",
    all(is.finite(volume_positive)) &&
      all(is.finite(volume_log_positive)),
    c(
      paste0("indicator=", paste(range(volume_positive), collapse = ":")),
      paste0(
        "log-positive=",
        paste(round(range(volume_log_positive), 4), collapse = ":")
      )
    ),
    "all finite"
  )

  spline_basis <- splines::ns(log2_psa, df = 3)
  add_check(
    "psa_spline_basis",
    "transformations",
    "The planned natural-spline PSA basis has three finite columns.",
    identical(dim(spline_basis), c(380L, 3L)) &&
      all(is.finite(spline_basis)),
    dim(spline_basis),
    c(380L, 3L)
  )

  primary_factor_counts <- c(
    table(factor(modern_data$dpros, levels = 1:4)),
    table(modern_data$capsule)
  )
  add_check(
    "primary_category_support",
    "model_readiness",
    "Every dpros level and outcome category has observed support.",
    all(primary_factor_counts > 0L),
    primary_factor_counts,
    "all counts > 0"
  )
}

validation_report <- do.call(rbind, validation_rows)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  validation_report,
  output_path,
  row.names = FALSE,
  na = ""
)

failed_checks <- validation_report$check_id[
  validation_report$status == "FAIL"
]

if (length(failed_checks) > 0L) {
  stop(
    paste(
      "Modern input validation failed.",
      "Review:",
      output_path,
      "Failed checks:",
      paste(failed_checks, collapse = ", ")
    ),
    call. = FALSE
  )
}

message(
  "Modern input validation passed: ",
  nrow(validation_report),
  " checks."
)
message("Created: ", output_path)
