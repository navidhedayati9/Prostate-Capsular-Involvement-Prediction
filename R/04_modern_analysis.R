# Run the prespecified primary modern analysis.
#
# Run from the project root after acquiring and validating the data:
#   source("R/01_acquire_data.R")
#   source("R/03_validate_inputs.R")
#   source("R/04_modern_analysis.R")
#
# This script compares:
#   1. the historical benchmark,
#   2. the prespecified primary modern model, and
#   3. the prespecified nonlinear PSA model.
#
# The dcaps, tumor-volume, and zero-Gleason sensitivity analyses are reserved
# for R/05_sensitivity_analyses.R.

check_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

project_markers <- c(
  "prostatic-capsular-penetration.Rproj",
  "renv.lock",
  "report/modern-analysis-plan.md",
  "report/sample-size-assessment.md",
  "R/03_validate_inputs.R"
)

check_true(
  all(file.exists(project_markers)),
  paste(
    "Run this script from the project root.",
    "Expected to find:",
    paste(project_markers, collapse = ", ")
  )
)

# Re-run the approved validation checks immediately before analysis.
source("R/03_validate_inputs.R", local = TRUE)

validation_path <- file.path(
  "results",
  "modern",
  "input_validation.csv"
)
validation_report <- utils::read.csv(
  validation_path,
  stringsAsFactors = FALSE
)
check_true(
  nrow(validation_report) == 27L &&
    !anyDuplicated(validation_report$check_id) &&
    all(validation_report$status == "PASS"),
  "The modern input-validation report must contain 27 unique passing checks."
)

input_path <- file.path("data", "processed", "prostate_modern.rds")
modern_data <- readRDS(input_path)

results_dir <- file.path("results", "modern")
figures_dir <- file.path("figures", "modern")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

set.seed_value <- 20260803L
number_folds <- 10L
number_repeats <- 20L

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
set.seed(set.seed_value)

model_formulas <- list(
  historical_benchmark =
    capsule ~ factor(dpros) + psa + gleason,
  primary_modern =
    capsule ~ age + factor(dpros) + log2(psa) + gleason,
  nonlinear_psa =
    capsule ~ age + factor(dpros) +
      splines::ns(log2(psa), df = 3) + gleason
)

model_labels <- c(
  historical_benchmark = "Historical benchmark",
  primary_modern = "Primary modern",
  nonlinear_psa = "Nonlinear PSA"
)

make_stratified_folds <- function(outcome, folds) {
  assignment <- integer(length(outcome))

  for (outcome_level in sort(unique(outcome))) {
    row_indices <- which(outcome == outcome_level)
    shuffled_indices <- sample(row_indices, length(row_indices))
    assignment[shuffled_indices] <- rep(
      seq_len(folds),
      length.out = length(shuffled_indices)
    )
  }

  assignment
}

resampling_assignments <- do.call(
  rbind,
  lapply(seq_len(number_repeats), function(repeat_number) {
    data.frame(
      id = modern_data$id,
      capsule = modern_data$capsule,
      repeat_id = repeat_number,
      fold = make_stratified_folds(
        modern_data$capsule,
        number_folds
      )
    )
  })
)

resampling_assignments <- resampling_assignments[
  order(
    resampling_assignments$repeat_id,
    resampling_assignments$fold,
    resampling_assignments$id
  ),
]
row.names(resampling_assignments) <- NULL

assignment_counts <- stats::aggregate(
  id ~ repeat_id + fold + capsule,
  data = resampling_assignments,
  FUN = length
)
names(assignment_counts)[names(assignment_counts) == "id"] <- "observations"

check_true(
  nrow(resampling_assignments) ==
    nrow(modern_data) * number_repeats,
  "Unexpected number of resampling-assignment rows."
)
check_true(
  all(table(resampling_assignments$id) == number_repeats),
  "Every observation must appear once in every repeat."
)
check_true(
  all(assignment_counts$observations > 0L),
  "Every assessment fold must contain both outcome categories."
)

write_csv(
  resampling_assignments,
  file.path(results_dir, "resampling_assignments.csv")
)
write_csv(
  assignment_counts,
  file.path(results_dir, "resampling_fold_counts.csv")
)

fit_logistic_model <- function(formula, analysis_data) {
  fit <- stats::glm(
    formula,
    family = stats::binomial(link = "logit"),
    data = analysis_data,
    control = stats::glm.control(maxit = 100L)
  )

  check_true(fit$converged, "A cross-validation model did not converge.")
  check_true(
    !anyNA(stats::coef(fit)),
    "A cross-validation model produced unidentified coefficients."
  )

  fit
}

prediction_rows <- vector(
  mode = "list",
  length = length(model_formulas) * number_repeats * number_folds
)
diagnostic_rows <- vector(
  mode = "list",
  length = length(prediction_rows)
)
output_index <- 0L

for (repeat_number in seq_len(number_repeats)) {
  repeat_assignment <- resampling_assignments[
    resampling_assignments$repeat_id == repeat_number,
  ]

  for (fold_number in seq_len(number_folds)) {
    assessment_ids <- repeat_assignment$id[
      repeat_assignment$fold == fold_number
    ]
    analysis_data <- modern_data[
      !modern_data$id %in% assessment_ids,
      ,
      drop = FALSE
    ]
    assessment_data <- modern_data[
      modern_data$id %in% assessment_ids,
      ,
      drop = FALSE
    ]

    check_true(
      nrow(analysis_data) + nrow(assessment_data) == nrow(modern_data),
      "Analysis and assessment portions do not reconstruct the dataset."
    )
    check_true(
      length(intersect(analysis_data$id, assessment_data$id)) == 0L,
      "Analysis and assessment portions overlap."
    )

    for (model_name in names(model_formulas)) {
      output_index <- output_index + 1L
      fit <- fit_logistic_model(
        model_formulas[[model_name]],
        analysis_data
      )
      predicted_probability <- stats::predict(
        fit,
        newdata = assessment_data,
        type = "response"
      )

      check_true(
        length(predicted_probability) == nrow(assessment_data) &&
          all(is.finite(predicted_probability)) &&
          all(
            predicted_probability > 0 &
              predicted_probability < 1
          ),
        "A cross-validation model produced invalid probabilities."
      )

      prediction_rows[[output_index]] <- data.frame(
        id = assessment_data$id,
        capsule = assessment_data$capsule,
        repeat_id = repeat_number,
        fold = fold_number,
        model = model_name,
        predicted_probability = unname(predicted_probability)
      )

      diagnostic_rows[[output_index]] <- data.frame(
        repeat_id = repeat_number,
        fold = fold_number,
        model = model_name,
        analysis_observations = nrow(analysis_data),
        analysis_events = sum(analysis_data$capsule == 1L),
        assessment_observations = nrow(assessment_data),
        assessment_events = sum(assessment_data$capsule == 1L),
        converged = fit$converged,
        iterations = fit$iter,
        minimum_prediction = min(predicted_probability),
        maximum_prediction = max(predicted_probability)
      )
    }
  }
}

out_of_fold_predictions <- do.call(rbind, prediction_rows)
resampling_diagnostics <- do.call(rbind, diagnostic_rows)

out_of_fold_predictions <- out_of_fold_predictions[
  order(
    out_of_fold_predictions$model,
    out_of_fold_predictions$repeat_id,
    out_of_fold_predictions$id
  ),
]
row.names(out_of_fold_predictions) <- NULL

expected_prediction_rows <-
  nrow(modern_data) * number_repeats * length(model_formulas)

check_true(
  nrow(out_of_fold_predictions) == expected_prediction_rows,
  "Unexpected number of out-of-fold predictions."
)

prediction_key <- interaction(
  out_of_fold_predictions$model,
  out_of_fold_predictions$repeat_id,
  out_of_fold_predictions$id,
  drop = TRUE
)
check_true(
  !anyDuplicated(prediction_key),
  "Each model-repeat-observation combination must be unique."
)
check_true(
  all(
    table(
      out_of_fold_predictions$model,
      out_of_fold_predictions$repeat_id
    ) == nrow(modern_data)
  ),
  "Every model must produce 380 out-of-fold predictions per repeat."
)
check_true(
  all(resampling_diagnostics$converged),
  "All cross-validation models must converge."
)

write_csv(
  out_of_fold_predictions,
  file.path(results_dir, "out_of_fold_predictions.csv")
)
write_csv(
  resampling_diagnostics,
  file.path(results_dir, "resampling_diagnostics.csv")
)

calculate_auc <- function(outcome, probability) {
  positive <- outcome == 1L
  number_positive <- sum(positive)
  number_negative <- sum(!positive)
  probability_ranks <- rank(probability, ties.method = "average")

  (
    sum(probability_ranks[positive]) -
      number_positive * (number_positive + 1) / 2
  ) / (number_positive * number_negative)
}

calculate_calibration <- function(outcome, probability) {
  clipped_probability <- pmin(
    pmax(probability, .Machine$double.eps),
    1 - .Machine$double.eps
  )
  linear_predictor <- stats::qlogis(clipped_probability)

  intercept_fit <- stats::glm(
    outcome ~ 1,
    offset = linear_predictor,
    family = stats::binomial(link = "logit")
  )
  slope_fit <- stats::glm(
    outcome ~ linear_predictor,
    family = stats::binomial(link = "logit")
  )

  c(
    calibration_intercept = unname(stats::coef(intercept_fit)[1]),
    calibration_slope = unname(stats::coef(slope_fit)[2])
  )
}

calculate_metrics <- function(outcome, probability) {
  clipped_probability <- pmin(
    pmax(probability, .Machine$double.eps),
    1 - .Machine$double.eps
  )
  calibration <- calculate_calibration(outcome, clipped_probability)

  c(
    auc = calculate_auc(outcome, clipped_probability),
    brier = mean((outcome - clipped_probability)^2),
    log_loss = -mean(
      outcome * log(clipped_probability) +
        (1 - outcome) * log(1 - clipped_probability)
    ),
    calibration
  )
}

metric_rows <- list()
metric_index <- 0L

for (model_name in names(model_formulas)) {
  for (repeat_number in seq_len(number_repeats)) {
    repeat_predictions <- out_of_fold_predictions[
      out_of_fold_predictions$model == model_name &
        out_of_fold_predictions$repeat_id == repeat_number,
    ]
    metrics <- calculate_metrics(
      repeat_predictions$capsule,
      repeat_predictions$predicted_probability
    )

    for (metric_name in names(metrics)) {
      metric_index <- metric_index + 1L
      metric_rows[[metric_index]] <- data.frame(
        model = model_name,
        repeat_id = repeat_number,
        metric = metric_name,
        estimate = unname(metrics[[metric_name]])
      )
    }
  }
}

resampling_performance <- do.call(rbind, metric_rows)
check_true(
  nrow(resampling_performance) ==
    length(model_formulas) * number_repeats * 5L,
  "Unexpected number of repeat-level performance estimates."
)
check_true(
  all(is.finite(resampling_performance$estimate)),
  "Performance estimates must be finite."
)

summarize_values <- function(values) {
  c(
    mean = mean(values),
    sd = stats::sd(values),
    median = stats::median(values),
    resampling_p025 = unname(stats::quantile(values, 0.025)),
    resampling_p975 = unname(stats::quantile(values, 0.975)),
    minimum = min(values),
    maximum = max(values)
  )
}

performance_summary <- do.call(
  rbind,
  lapply(
    split(
      resampling_performance,
      list(
        resampling_performance$model,
        resampling_performance$metric
      ),
      drop = TRUE
    ),
    function(group) {
      summaries <- summarize_values(group$estimate)
      data.frame(
        model = group$model[1],
        metric = group$metric[1],
        as.list(summaries),
        repeats = nrow(group),
        interval_note = paste(
          "2.5th and 97.5th percentiles describe resampling variability;",
          "they are not formal confidence limits."
        ),
        row.names = NULL
      )
    }
  )
)
row.names(performance_summary) <- NULL
performance_summary <- performance_summary[
  order(performance_summary$metric, performance_summary$model),
]

comparison_pairs <- list(
  c("primary_modern", "historical_benchmark"),
  c("nonlinear_psa", "primary_modern")
)

paired_rows <- list()
paired_index <- 0L

for (pair in comparison_pairs) {
  for (metric_name in unique(resampling_performance$metric)) {
    model_a <- resampling_performance[
      resampling_performance$model == pair[1] &
        resampling_performance$metric == metric_name,
      c("repeat_id", "estimate")
    ]
    model_b <- resampling_performance[
      resampling_performance$model == pair[2] &
        resampling_performance$metric == metric_name,
      c("repeat_id", "estimate")
    ]
    paired <- merge(
      model_a,
      model_b,
      by = "repeat_id",
      suffixes = c("_model_a", "_model_b"),
      sort = TRUE
    )

    paired_index <- paired_index + 1L
    paired_rows[[paired_index]] <- data.frame(
      model_a = pair[1],
      model_b = pair[2],
      metric = metric_name,
      repeat_id = paired$repeat_id,
      model_a_estimate = paired$estimate_model_a,
      model_b_estimate = paired$estimate_model_b,
      raw_difference_a_minus_b =
        paired$estimate_model_a - paired$estimate_model_b,
      interpretation = paste(
        "Raw paired difference; preferred direction depends on the metric."
      )
    )
  }
}

paired_metric_differences <- do.call(rbind, paired_rows)

write_csv(
  resampling_performance,
  file.path(results_dir, "resampling_performance.csv")
)
write_csv(
  performance_summary,
  file.path(results_dir, "performance_summary.csv")
)
write_csv(
  paired_metric_differences,
  file.path(results_dir, "paired_metric_differences.csv")
)

mean_out_of_fold_predictions <- stats::aggregate(
  predicted_probability ~ id + capsule + model,
  data = out_of_fold_predictions,
  FUN = mean
)
names(mean_out_of_fold_predictions)[
  names(mean_out_of_fold_predictions) == "predicted_probability"
] <- "mean_out_of_fold_probability"

prediction_variability <- stats::aggregate(
  predicted_probability ~ id + capsule + model,
  data = out_of_fold_predictions,
  FUN = stats::sd
)
names(prediction_variability)[
  names(prediction_variability) == "predicted_probability"
] <- "sd_out_of_fold_probability"

mean_out_of_fold_predictions <- merge(
  mean_out_of_fold_predictions,
  prediction_variability,
  by = c("id", "capsule", "model"),
  sort = TRUE
)

write_csv(
  mean_out_of_fold_predictions,
  file.path(results_dir, "mean_out_of_fold_predictions.csv")
)

# Fit models to the full dataset only after cross-validation is complete.
full_model_fits <- lapply(
  model_formulas,
  function(formula) fit_logistic_model(formula, modern_data)
)

full_coefficient_rows <- lapply(
  names(full_model_fits),
  function(model_name) {
    fit <- full_model_fits[[model_name]]
    coefficient_matrix <- summary(fit)$coefficients
    confidence_multiplier <- stats::qnorm(0.975)

    data.frame(
      model = model_name,
      term = rownames(coefficient_matrix),
      estimate = coefficient_matrix[, "Estimate"],
      std_error = coefficient_matrix[, "Std. Error"],
      z_value = coefficient_matrix[, "z value"],
      p_value = coefficient_matrix[, "Pr(>|z|)"],
      wald_ci_lower =
        coefficient_matrix[, "Estimate"] -
          confidence_multiplier * coefficient_matrix[, "Std. Error"],
      wald_ci_upper =
        coefficient_matrix[, "Estimate"] +
          confidence_multiplier * coefficient_matrix[, "Std. Error"],
      row.names = NULL
    )
  }
)
full_model_coefficients <- do.call(rbind, full_coefficient_rows)

full_model_summary <- do.call(
  rbind,
  lapply(names(full_model_fits), function(model_name) {
    fit <- full_model_fits[[model_name]]
    full_probability <- stats::fitted(fit)
    metrics <- calculate_metrics(modern_data$capsule, full_probability)

    data.frame(
      model = model_name,
      formula = paste(
        deparse(stats::formula(fit)),
        collapse = " "
      ),
      observations = stats::nobs(fit),
      predictor_parameters = length(stats::coef(fit)) - 1L,
      aic = stats::AIC(fit),
      converged = fit$converged,
      apparent_auc = unname(metrics[["auc"]]),
      apparent_brier = unname(metrics[["brier"]]),
      apparent_log_loss = unname(metrics[["log_loss"]]),
      apparent_calibration_intercept =
        unname(metrics[["calibration_intercept"]]),
      apparent_calibration_slope =
        unname(metrics[["calibration_slope"]]),
      performance_note = paste(
        "Apparent full-data performance;",
        "cross-validated performance is primary."
      )
    )
  })
)

write_csv(
  full_model_coefficients,
  file.path(results_dir, "full_model_coefficients.csv")
)
write_csv(
  full_model_summary,
  file.path(results_dir, "full_model_summary.csv")
)

calculate_roc_points <- function(outcome, probability) {
  thresholds <- c(
    Inf,
    sort(unique(probability), decreasing = TRUE),
    -Inf
  )
  data.frame(
    threshold = thresholds,
    false_positive_rate = vapply(
      thresholds,
      function(threshold) {
        mean(probability[outcome == 0L] >= threshold)
      },
      numeric(1)
    ),
    true_positive_rate = vapply(
      thresholds,
      function(threshold) {
        mean(probability[outcome == 1L] >= threshold)
      },
      numeric(1)
    )
  )
}

model_colors <- c(
  historical_benchmark = "#666666",
  primary_modern = "#2166AC",
  nonlinear_psa = "#B2182B"
)

grDevices::png(
  filename = file.path(figures_dir, "cross_validated_roc.png"),
  width = 1200,
  height = 1000,
  res = 150
)
graphics::plot(
  c(0, 1),
  c(0, 1),
  type = "n",
  xlim = c(0, 1),
  ylim = c(0, 1),
  xaxs = "i",
  yaxs = "i",
  xlab = "False-positive rate",
  ylab = "True-positive rate",
  main = "Cross-validated ROC curves"
)
graphics::abline(0, 1, lty = 2, col = "grey70")

roc_legend <- character(length(model_formulas))
roc_index <- 0L
for (model_name in names(model_formulas)) {
  roc_index <- roc_index + 1L
  model_predictions <- mean_out_of_fold_predictions[
    mean_out_of_fold_predictions$model == model_name,
  ]
  roc_points <- calculate_roc_points(
    model_predictions$capsule,
    model_predictions$mean_out_of_fold_probability
  )
  graphics::lines(
    roc_points$false_positive_rate,
    roc_points$true_positive_rate,
    col = model_colors[[model_name]],
    lwd = 3
  )
  roc_auc <- calculate_auc(
    model_predictions$capsule,
    model_predictions$mean_out_of_fold_probability
  )
  roc_legend[roc_index] <- paste0(
    model_labels[[model_name]],
    " (AUC ",
    format(round(roc_auc, 3), nsmall = 3),
    ")"
  )
}
graphics::legend(
  "bottomright",
  legend = roc_legend,
  col = model_colors[names(model_formulas)],
  lwd = 3,
  bty = "n"
)
grDevices::dev.off()

assign_calibration_groups <- function(probability, groups = 10L) {
  ordered_rows <- order(probability)
  group <- integer(length(probability))
  group[ordered_rows] <- ceiling(
    seq_along(probability) * groups / length(probability)
  )
  pmin(group, groups)
}

grDevices::png(
  filename = file.path(figures_dir, "cross_validated_calibration.png"),
  width = 1800,
  height = 650,
  res = 150
)
graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

for (model_name in names(model_formulas)) {
  model_predictions <- mean_out_of_fold_predictions[
    mean_out_of_fold_predictions$model == model_name,
  ]
  calibration_group <- assign_calibration_groups(
    model_predictions$mean_out_of_fold_probability
  )
  grouped_calibration <- stats::aggregate(
    cbind(
      predicted = model_predictions$mean_out_of_fold_probability,
      observed = model_predictions$capsule
    ),
    by = list(group = calibration_group),
    FUN = mean
  )
  smooth_calibration <- stats::loess(
    capsule ~ mean_out_of_fold_probability,
    data = model_predictions,
    span = 0.75,
    degree = 1
  )
  probability_grid <- seq(
    min(model_predictions$mean_out_of_fold_probability),
    max(model_predictions$mean_out_of_fold_probability),
    length.out = 200L
  )
  smooth_observed <- stats::predict(
    smooth_calibration,
    newdata = data.frame(
      mean_out_of_fold_probability = probability_grid
    )
  )

  graphics::plot(
    c(0, 1),
    c(0, 1),
    type = "n",
    xlab = "Predicted probability",
    ylab = "Observed proportion",
    main = model_labels[[model_name]],
    asp = 1
  )
  graphics::abline(0, 1, lty = 2, col = "grey70")
  graphics::lines(
    probability_grid,
    pmin(pmax(smooth_observed, 0), 1),
    col = model_colors[[model_name]],
    lwd = 3
  )
  graphics::points(
    grouped_calibration$predicted,
    grouped_calibration$observed,
    pch = 19,
    col = model_colors[[model_name]]
  )
  graphics::rug(
    model_predictions$mean_out_of_fold_probability,
    side = 1,
    col = grDevices::adjustcolor("grey30", alpha.f = 0.25)
  )
}
grDevices::dev.off()

grDevices::png(
  filename = file.path(
    figures_dir,
    "cross_validated_probability_distributions.png"
  ),
  width = 1800,
  height = 650,
  res = 150
)
graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

for (model_name in names(model_formulas)) {
  model_predictions <- mean_out_of_fold_predictions[
    mean_out_of_fold_predictions$model == model_name,
  ]
  graphics::boxplot(
    mean_out_of_fold_probability ~ factor(capsule),
    data = model_predictions,
    col = c("#D9EAF7", "#F4CCCC"),
    xlab = "Capsule penetration",
    ylab = "Mean out-of-fold probability",
    names = c("No", "Yes"),
    main = model_labels[[model_name]],
    ylim = c(0, 1)
  )
  graphics::stripchart(
    mean_out_of_fold_probability ~ factor(capsule),
    data = model_predictions,
    method = "jitter",
    vertical = TRUE,
    add = TRUE,
    pch = 19,
    cex = 0.45,
    col = grDevices::adjustcolor("black", alpha.f = 0.30)
  )
}
grDevices::dev.off()

metric_order <- c(
  "auc",
  "brier",
  "log_loss",
  "calibration_intercept",
  "calibration_slope"
)
metric_labels <- c(
  auc = "ROC AUC",
  brier = "Brier score",
  log_loss = "Log loss",
  calibration_intercept = "Calibration intercept",
  calibration_slope = "Calibration slope"
)
performance_axis_labels <- c(
  historical_benchmark = "Historical",
  primary_modern = "Primary",
  nonlinear_psa = "Nonlinear PSA"
)

grDevices::png(
  filename = file.path(
    figures_dir,
    "paired_resampling_performance.png"
  ),
  width = 2100,
  height = 1000,
  res = 150
)
graphics::par(mfrow = c(2, 3), mar = c(7, 4.5, 3, 1))

for (metric_name in metric_order) {
  metric_data <- resampling_performance[
    resampling_performance$metric == metric_name,
  ]
  metric_data$model <- factor(
    metric_data$model,
    levels = names(model_formulas)
  )
  graphics::boxplot(
    estimate ~ model,
    data = metric_data,
    names = unname(performance_axis_labels[names(model_formulas)]),
    las = 2,
    col = unname(model_colors[names(model_formulas)]),
    xlab = "",
    ylab = metric_labels[[metric_name]],
    main = metric_labels[[metric_name]]
  )

  for (repeat_number in seq_len(number_repeats)) {
    repeat_data <- metric_data[
      metric_data$repeat_id == repeat_number,
    ]
    repeat_data <- repeat_data[
      order(repeat_data$model),
    ]
    graphics::lines(
      seq_along(model_formulas),
      repeat_data$estimate,
      col = grDevices::adjustcolor("black", alpha.f = 0.15)
    )
  }

  if (metric_name == "calibration_intercept") {
    graphics::abline(h = 0, lty = 2, col = "grey30")
  }
  if (metric_name == "calibration_slope") {
    graphics::abline(h = 1, lty = 2, col = "grey30")
  }
}
graphics::plot.new()
graphics::legend(
  "center",
  legend = c(
    "Each line connects models evaluated",
    "on the same resampling repeat."
  ),
  bty = "n"
)
grDevices::dev.off()

analysis_metadata <- c(
  "Primary modern analysis",
  paste("Run date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste("Input:", input_path),
  paste("Input MD5:", unname(tools::md5sum(input_path))),
  paste("Random seed:", set.seed_value),
  paste("Folds:", number_folds),
  paste("Repeats:", number_repeats),
  "Resampling: outcome-stratified repeated 10-fold cross-validation",
  paste(
    "Historical benchmark:",
    paste(deparse(model_formulas$historical_benchmark), collapse = " ")
  ),
  paste(
    "Primary modern:",
    paste(deparse(model_formulas$primary_modern), collapse = " ")
  ),
  paste(
    "Nonlinear PSA:",
    paste(deparse(model_formulas$nonlinear_psa), collapse = " ")
  ),
  "Performance: ROC AUC, Brier score, log loss, calibration intercept, calibration slope",
  "Primary estimates: repeat-level out-of-fold performance",
  "Percentile ranges describe resampling variability, not confidence intervals",
  "Full-data performance is apparent and secondary",
  "Sensitivity and exploratory analyses are excluded from this script"
)
writeLines(
  analysis_metadata,
  file.path(results_dir, "modern_analysis_metadata.txt"),
  useBytes = TRUE
)

message("Primary modern analysis completed successfully.")
message(
  "Created ",
  nrow(out_of_fold_predictions),
  " out-of-fold predictions."
)
message("Results: ", results_dir)
message("Figures: ", figures_dir)
