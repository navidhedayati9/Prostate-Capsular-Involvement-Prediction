# Run the prespecified sensitivity and exploratory analyses.
#
# Run from the project root after the primary modern analysis:
#   source("R/01_acquire_data.R")
#   source("R/03_validate_inputs.R")
#   source("R/04_modern_analysis.R")
#   source("R/05_sensitivity_analyses.R")
#
# This script evaluates:
#   1. the primary model plus dcaps,
#   2. the primary model after excluding two zero Gleason scores, and
#   3. the primary model plus a two-part tumor-volume representation.
#
# The saved primary-analysis resampling assignments are reused unchanged.

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
  "R/03_validate_inputs.R",
  "R/04_modern_analysis.R"
)

check_true(
  all(file.exists(project_markers)),
  paste(
    "Run this script from the project root.",
    "Expected to find:",
    paste(project_markers, collapse = ", ")
  )
)

source("R/03_validate_inputs.R", local = TRUE)

results_dir <- file.path("results", "modern")
figures_dir <- file.path("figures", "modern")
input_path <- file.path("data", "processed", "prostate_modern.rds")
assignment_path <- file.path(results_dir, "resampling_assignments.csv")
primary_prediction_path <- file.path(
  results_dir,
  "out_of_fold_predictions.csv"
)

check_true(
  file.exists(assignment_path),
  paste(
    "Saved resampling assignments not found:",
    assignment_path,
    "Run source(\"R/04_modern_analysis.R\") first."
  )
)
check_true(
  file.exists(primary_prediction_path),
  paste(
    "Primary out-of-fold predictions not found:",
    primary_prediction_path,
    "Run source(\"R/04_modern_analysis.R\") first."
  )
)

modern_data <- readRDS(input_path)
resampling_assignments <- utils::read.csv(
  assignment_path,
  stringsAsFactors = FALSE
)
primary_predictions_all <- utils::read.csv(
  primary_prediction_path,
  stringsAsFactors = FALSE
)

expected_assignment_names <- c(
  "id", "capsule", "repeat_id", "fold"
)
check_true(
  identical(names(resampling_assignments), expected_assignment_names),
  "Saved resampling assignments have unexpected columns or order."
)
check_true(
  nrow(resampling_assignments) == 7600L,
  "Expected 7,600 saved resampling-assignment rows."
)
check_true(
  !anyDuplicated(
    interaction(
      resampling_assignments$id,
      resampling_assignments$repeat_id,
      drop = TRUE
    )
  ),
  "Each observation must occur exactly once per saved repeat."
)
check_true(
  all(table(resampling_assignments$id) == 20L),
  "Every observation must occur in all 20 saved repeats."
)
check_true(
  identical(sort(unique(resampling_assignments$fold)), 1:10),
  "Saved fold identifiers must be 1 through 10."
)
check_true(
  identical(sort(unique(resampling_assignments$repeat_id)), 1:20),
  "Saved repeat identifiers must be 1 through 20."
)
check_true(
  all(
    resampling_assignments$capsule ==
      modern_data$capsule[
        match(resampling_assignments$id, modern_data$id)
      ]
  ),
  "Saved assignment outcomes do not match the validated dataset."
)

primary_predictions <- primary_predictions_all[
  primary_predictions_all$model == "primary_modern",
]
check_true(
  nrow(primary_predictions) == 7600L,
  "Expected 7,600 primary-model out-of-fold predictions."
)
check_true(
  !anyDuplicated(
    interaction(
      primary_predictions$id,
      primary_predictions$repeat_id,
      drop = TRUE
    )
  ),
  "Primary predictions must be unique by observation and repeat."
)

# Deterministic two-part representation prespecified for tumor volume.
modern_data$vol_positive <- as.integer(modern_data$vol > 0)
modern_data$log2_vol_positive <- ifelse(
  modern_data$vol > 0,
  log2(modern_data$vol),
  0
)

check_true(
  all(is.finite(modern_data$vol_positive)) &&
    all(is.finite(modern_data$log2_vol_positive)),
  "The two-part tumor-volume representation must be finite."
)
check_true(
  sum(modern_data$vol_positive == 0L) == 167L,
  "Expected 167 zero tumor-volume records."
)

zero_gleason_ids <- c(282L, 357L)
check_true(
  identical(
    modern_data$id[modern_data$gleason == 0L],
    zero_gleason_ids
  ),
  "Unexpected records have zero Gleason score."
)

sensitivity_definitions <- list(
  dcaps_sensitivity = list(
    label = "Primary + dcaps",
    formula =
      capsule ~ age + factor(dpros) + log2(psa) +
        gleason + factor(dcaps),
    eligible_ids = modern_data$id,
    classification = "Sensitivity"
  ),
  zero_gleason_exclusion = list(
    label = "Exclude zero Gleason",
    formula =
      capsule ~ age + factor(dpros) + log2(psa) + gleason,
    eligible_ids = setdiff(modern_data$id, zero_gleason_ids),
    classification = "Data-quality sensitivity"
  ),
  volume_exploratory = list(
    label = "Primary + tumor volume",
    formula =
      capsule ~ age + factor(dpros) + log2(psa) +
        gleason + vol_positive + log2_vol_positive,
    eligible_ids = modern_data$id,
    classification = "Exploratory"
  )
)

fit_logistic_model <- function(formula, analysis_data) {
  fit <- stats::glm(
    formula,
    family = stats::binomial(link = "logit"),
    data = analysis_data,
    control = stats::glm.control(maxit = 100L)
  )

  check_true(fit$converged, "A sensitivity model did not converge.")
  check_true(
    !anyNA(stats::coef(fit)),
    "A sensitivity model produced unidentified coefficients."
  )

  fit
}

prediction_rows <- vector(
  mode = "list",
  length = length(sensitivity_definitions) * 20L * 10L
)
diagnostic_rows <- vector(
  mode = "list",
  length = length(prediction_rows)
)
output_index <- 0L

for (analysis_name in names(sensitivity_definitions)) {
  definition <- sensitivity_definitions[[analysis_name]]
  eligible_data <- modern_data[
    modern_data$id %in% definition$eligible_ids,
    ,
    drop = FALSE
  ]
  eligible_assignments <- resampling_assignments[
    resampling_assignments$id %in% definition$eligible_ids,
  ]

  check_true(
    all(table(eligible_assignments$id) == 20L),
    paste("Incomplete saved assignments for", analysis_name)
  )

  eligible_fold_counts <- table(
    eligible_assignments$repeat_id,
    eligible_assignments$fold,
    eligible_assignments$capsule
  )
  check_true(
    all(eligible_fold_counts > 0L),
    paste("An assessment fold lacks an outcome category for", analysis_name)
  )

  for (repeat_number in 1:20) {
    repeat_assignments <- eligible_assignments[
      eligible_assignments$repeat_id == repeat_number,
    ]

    for (fold_number in 1:10) {
      assessment_ids <- repeat_assignments$id[
        repeat_assignments$fold == fold_number
      ]
      analysis_data <- eligible_data[
        !eligible_data$id %in% assessment_ids,
        ,
        drop = FALSE
      ]
      assessment_data <- eligible_data[
        eligible_data$id %in% assessment_ids,
        ,
        drop = FALSE
      ]

      check_true(
        nrow(analysis_data) + nrow(assessment_data) ==
          nrow(eligible_data),
        paste("Fold portions do not reconstruct", analysis_name)
      )
      check_true(
        length(intersect(analysis_data$id, assessment_data$id)) == 0L,
        paste("Fold portions overlap for", analysis_name)
      )

      fit <- fit_logistic_model(definition$formula, analysis_data)
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
        paste("Invalid sensitivity probabilities for", analysis_name)
      )

      output_index <- output_index + 1L
      prediction_rows[[output_index]] <- data.frame(
        id = assessment_data$id,
        capsule = assessment_data$capsule,
        repeat_id = repeat_number,
        fold = fold_number,
        analysis = analysis_name,
        classification = definition$classification,
        predicted_probability = unname(predicted_probability)
      )
      diagnostic_rows[[output_index]] <- data.frame(
        repeat_id = repeat_number,
        fold = fold_number,
        analysis = analysis_name,
        classification = definition$classification,
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

sensitivity_predictions <- do.call(rbind, prediction_rows)
sensitivity_diagnostics <- do.call(rbind, diagnostic_rows)

sensitivity_predictions <- sensitivity_predictions[
  order(
    sensitivity_predictions$analysis,
    sensitivity_predictions$repeat_id,
    sensitivity_predictions$id
  ),
]
row.names(sensitivity_predictions) <- NULL

expected_prediction_count <- 7600L + 7560L + 7600L
check_true(
  nrow(sensitivity_predictions) == expected_prediction_count,
  "Unexpected number of sensitivity out-of-fold predictions."
)
check_true(
  !anyDuplicated(
    interaction(
      sensitivity_predictions$analysis,
      sensitivity_predictions$repeat_id,
      sensitivity_predictions$id,
      drop = TRUE
    )
  ),
  "Sensitivity predictions must be unique by analysis, repeat, and ID."
)
check_true(
  nrow(sensitivity_diagnostics) == 600L &&
    all(sensitivity_diagnostics$converged),
  "Expected 600 converged sensitivity-model fits."
)

write_csv(
  sensitivity_predictions,
  file.path(results_dir, "sensitivity_out_of_fold_predictions.csv")
)
write_csv(
  sensitivity_diagnostics,
  file.path(results_dir, "sensitivity_resampling_diagnostics.csv")
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

calculate_metrics <- function(outcome, probability) {
  clipped_probability <- pmin(
    pmax(probability, .Machine$double.eps),
    1 - .Machine$double.eps
  )
  linear_predictor <- stats::qlogis(clipped_probability)
  calibration_intercept_fit <- stats::glm(
    outcome ~ 1,
    offset = linear_predictor,
    family = stats::binomial(link = "logit")
  )
  calibration_slope_fit <- stats::glm(
    outcome ~ linear_predictor,
    family = stats::binomial(link = "logit")
  )

  c(
    auc = calculate_auc(outcome, clipped_probability),
    brier = mean((outcome - clipped_probability)^2),
    log_loss = -mean(
      outcome * log(clipped_probability) +
        (1 - outcome) * log(1 - clipped_probability)
    ),
    calibration_intercept =
      unname(stats::coef(calibration_intercept_fit)[1]),
    calibration_slope =
      unname(stats::coef(calibration_slope_fit)[2])
  )
}

calculate_repeat_performance <- function(
  predictions,
  analysis_column,
  analysis_names
) {
  performance_rows <- list()
  performance_index <- 0L

  for (analysis_name in analysis_names) {
    for (repeat_number in 1:20) {
      keep <- predictions[[analysis_column]] == analysis_name &
        predictions$repeat_id == repeat_number
      repeat_predictions <- predictions[keep, ]
      metrics <- calculate_metrics(
        repeat_predictions$capsule,
        repeat_predictions$predicted_probability
      )

      for (metric_name in names(metrics)) {
        performance_index <- performance_index + 1L
        performance_rows[[performance_index]] <- data.frame(
          analysis = analysis_name,
          repeat_id = repeat_number,
          metric = metric_name,
          estimate = unname(metrics[[metric_name]])
        )
      }
    }
  }

  do.call(rbind, performance_rows)
}

sensitivity_performance <- calculate_repeat_performance(
  sensitivity_predictions,
  "analysis",
  names(sensitivity_definitions)
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

sensitivity_performance_summary <- do.call(
  rbind,
  lapply(
    split(
      sensitivity_performance,
      list(
        sensitivity_performance$analysis,
        sensitivity_performance$metric
      ),
      drop = TRUE
    ),
    function(group) {
      summaries <- summarize_values(group$estimate)
      data.frame(
        analysis = group$analysis[1],
        classification =
          sensitivity_definitions[[group$analysis[1]]]$classification,
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
row.names(sensitivity_performance_summary) <- NULL
sensitivity_performance_summary <- sensitivity_performance_summary[
  order(
    sensitivity_performance_summary$metric,
    sensitivity_performance_summary$analysis
  ),
]

check_true(
  nrow(sensitivity_performance) == 300L &&
    all(is.finite(sensitivity_performance$estimate)),
  "Expected 300 finite repeat-level sensitivity performance estimates."
)

write_csv(
  sensitivity_performance,
  file.path(results_dir, "sensitivity_performance.csv")
)
write_csv(
  sensitivity_performance_summary,
  file.path(results_dir, "sensitivity_performance_summary.csv")
)

# Recalculate primary-reference metrics on the applicable population.
reference_definitions <- list(
  dcaps_sensitivity = modern_data$id,
  zero_gleason_exclusion =
    setdiff(modern_data$id, zero_gleason_ids),
  volume_exploratory = modern_data$id
)

reference_rows <- list()
reference_index <- 0L
for (comparison_name in names(reference_definitions)) {
  eligible_ids <- reference_definitions[[comparison_name]]

  for (repeat_number in 1:20) {
    reference_predictions <- primary_predictions[
      primary_predictions$id %in% eligible_ids &
        primary_predictions$repeat_id == repeat_number,
    ]
    metrics <- calculate_metrics(
      reference_predictions$capsule,
      reference_predictions$predicted_probability
    )

    for (metric_name in names(metrics)) {
      reference_index <- reference_index + 1L
      reference_rows[[reference_index]] <- data.frame(
        analysis = comparison_name,
        repeat_id = repeat_number,
        metric = metric_name,
        reference_estimate = unname(metrics[[metric_name]])
      )
    }
  }
}
reference_performance <- do.call(rbind, reference_rows)

sensitivity_paired_differences <- merge(
  sensitivity_performance,
  reference_performance,
  by = c("analysis", "repeat_id", "metric"),
  sort = TRUE
)
names(sensitivity_paired_differences)[
  names(sensitivity_paired_differences) == "estimate"
] <- "sensitivity_estimate"

sensitivity_paired_differences$raw_difference_sensitivity_minus_primary <-
  sensitivity_paired_differences$sensitivity_estimate -
  sensitivity_paired_differences$reference_estimate

calculate_preferred_improvement <- function(
  metric,
  sensitivity_estimate,
  reference_estimate
) {
  if (metric == "auc") {
    return(sensitivity_estimate - reference_estimate)
  }
  if (metric %in% c("brier", "log_loss")) {
    return(reference_estimate - sensitivity_estimate)
  }
  if (metric == "calibration_intercept") {
    return(
      abs(reference_estimate) - abs(sensitivity_estimate)
    )
  }
  if (metric == "calibration_slope") {
    return(
      abs(reference_estimate - 1) -
        abs(sensitivity_estimate - 1)
    )
  }
  stop("Unknown metric.", call. = FALSE)
}

sensitivity_paired_differences$preferred_improvement <- mapply(
  calculate_preferred_improvement,
  sensitivity_paired_differences$metric,
  sensitivity_paired_differences$sensitivity_estimate,
  sensitivity_paired_differences$reference_estimate
)
sensitivity_paired_differences$interpretation <- paste(
  "Positive preferred_improvement favors the sensitivity or exploratory",
  "analysis; negative values favor the applicable primary-model reference."
)

write_csv(
  sensitivity_paired_differences,
  file.path(results_dir, "sensitivity_paired_differences.csv")
)

# Fit full-data models only after cross-validation is complete.
full_sensitivity_fits <- lapply(
  names(sensitivity_definitions),
  function(analysis_name) {
    definition <- sensitivity_definitions[[analysis_name]]
    eligible_data <- modern_data[
      modern_data$id %in% definition$eligible_ids,
      ,
      drop = FALSE
    ]
    fit_logistic_model(definition$formula, eligible_data)
  }
)
names(full_sensitivity_fits) <- names(sensitivity_definitions)

full_sensitivity_coefficients <- do.call(
  rbind,
  lapply(names(full_sensitivity_fits), function(analysis_name) {
    fit <- full_sensitivity_fits[[analysis_name]]
    coefficient_matrix <- summary(fit)$coefficients
    confidence_multiplier <- stats::qnorm(0.975)

    data.frame(
      analysis = analysis_name,
      classification =
        sensitivity_definitions[[analysis_name]]$classification,
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
  })
)

full_sensitivity_summary <- do.call(
  rbind,
  lapply(names(full_sensitivity_fits), function(analysis_name) {
    fit <- full_sensitivity_fits[[analysis_name]]
    eligible_ids <- sensitivity_definitions[[analysis_name]]$eligible_ids
    eligible_data <- modern_data[
      modern_data$id %in% eligible_ids,
      ,
      drop = FALSE
    ]
    probability <- stats::fitted(fit)
    metrics <- calculate_metrics(eligible_data$capsule, probability)

    data.frame(
      analysis = analysis_name,
      classification =
        sensitivity_definitions[[analysis_name]]$classification,
      formula = paste(deparse(stats::formula(fit)), collapse = " "),
      observations = stats::nobs(fit),
      events = sum(eligible_data$capsule == 1L),
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
  full_sensitivity_coefficients,
  file.path(results_dir, "sensitivity_full_model_coefficients.csv")
)
write_csv(
  full_sensitivity_summary,
  file.path(results_dir, "sensitivity_full_model_summary.csv")
)

mean_sensitivity_predictions <- stats::aggregate(
  predicted_probability ~ id + capsule + analysis,
  data = sensitivity_predictions,
  FUN = mean
)
names(mean_sensitivity_predictions)[
  names(mean_sensitivity_predictions) == "predicted_probability"
] <- "mean_out_of_fold_probability"

analysis_colors <- c(
  dcaps_sensitivity = "#5E3C99",
  zero_gleason_exclusion = "#E66101",
  volume_exploratory = "#1B7837"
)
analysis_axis_labels <- c(
  dcaps_sensitivity = "+ dcaps",
  zero_gleason_exclusion = "Exclude zeros",
  volume_exploratory = "+ volume"
)

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

grDevices::png(
  filename = file.path(
    figures_dir,
    "sensitivity_performance_comparison.png"
  ),
  width = 2100,
  height = 1000,
  res = 150
)
graphics::par(mfrow = c(2, 3), mar = c(6.5, 4.5, 3, 1))

for (metric_name in metric_order) {
  plot_data <- sensitivity_paired_differences[
    sensitivity_paired_differences$metric == metric_name,
  ]
  plot_data$analysis <- factor(
    plot_data$analysis,
    levels = names(sensitivity_definitions)
  )

  graphics::boxplot(
    preferred_improvement ~ analysis,
    data = plot_data,
    names = unname(
      analysis_axis_labels[names(sensitivity_definitions)]
    ),
    las = 2,
    col = unname(
      analysis_colors[names(sensitivity_definitions)]
    ),
    xlab = "",
    ylab = "Improvement",
    main = metric_labels[[metric_name]]
  )
  graphics::abline(h = 0, lty = 2, col = "grey30")
}
graphics::plot.new()
graphics::legend(
  "center",
  legend = c(
    "Positive values favor the sensitivity",
    "or exploratory analysis.",
    "Negative values favor the primary reference."
  ),
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
  filename = file.path(figures_dir, "sensitivity_calibration.png"),
  width = 1800,
  height = 650,
  res = 150
)
graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

for (analysis_name in names(sensitivity_definitions)) {
  plot_data <- mean_sensitivity_predictions[
    mean_sensitivity_predictions$analysis == analysis_name,
  ]
  calibration_group <- assign_calibration_groups(
    plot_data$mean_out_of_fold_probability
  )
  grouped_calibration <- stats::aggregate(
    cbind(
      predicted = plot_data$mean_out_of_fold_probability,
      observed = plot_data$capsule
    ),
    by = list(group = calibration_group),
    FUN = mean
  )
  smooth_calibration <- stats::loess(
    capsule ~ mean_out_of_fold_probability,
    data = plot_data,
    span = 0.75,
    degree = 1
  )
  probability_grid <- seq(
    min(plot_data$mean_out_of_fold_probability),
    max(plot_data$mean_out_of_fold_probability),
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
    main = sensitivity_definitions[[analysis_name]]$label
  )
  graphics::abline(0, 1, lty = 2, col = "grey70")
  graphics::lines(
    probability_grid,
    pmin(pmax(smooth_observed, 0), 1),
    col = analysis_colors[[analysis_name]],
    lwd = 3
  )
  graphics::points(
    grouped_calibration$predicted,
    grouped_calibration$observed,
    pch = 19,
    col = analysis_colors[[analysis_name]]
  )
  graphics::rug(
    plot_data$mean_out_of_fold_probability,
    side = 1,
    col = grDevices::adjustcolor("grey30", alpha.f = 0.25)
  )
}
grDevices::dev.off()

# Exploratory full-data volume relationship at prespecified reference values.
positive_volumes <- modern_data$vol[modern_data$vol > 0]
volume_grid <- exp(
  seq(
    log(min(positive_volumes)),
    log(max(positive_volumes)),
    length.out = 200L
  )
)
reference_age <- stats::median(modern_data$age)
reference_dpros <- as.integer(
  names(which.max(table(modern_data$dpros)))
)
reference_psa <- stats::median(modern_data$psa)
reference_gleason <- stats::median(modern_data$gleason)

volume_prediction_data <- data.frame(
  vol = c(0, volume_grid),
  age = reference_age,
  dpros = reference_dpros,
  psa = reference_psa,
  gleason = reference_gleason
)
volume_prediction_data$vol_positive <- as.integer(
  volume_prediction_data$vol > 0
)
volume_prediction_data$log2_vol_positive <- ifelse(
  volume_prediction_data$vol > 0,
  log2(volume_prediction_data$vol),
  0
)
volume_prediction_data$predicted_probability <- stats::predict(
  full_sensitivity_fits$volume_exploratory,
  newdata = volume_prediction_data,
  type = "response"
)
volume_prediction_data$interpretation <- paste(
  "Exploratory apparent prediction at age",
  reference_age,
  "dpros",
  reference_dpros,
  "PSA",
  reference_psa,
  "and Gleason",
  reference_gleason
)

write_csv(
  volume_prediction_data,
  file.path(results_dir, "exploratory_volume_prediction_curve.csv")
)

grDevices::png(
  filename = file.path(
    figures_dir,
    "exploratory_volume_predictions.png"
  ),
  width = 1200,
  height = 900,
  res = 150
)
graphics::plot(
  volume_grid,
  volume_prediction_data$predicted_probability[-1],
  type = "l",
  log = "x",
  lwd = 3,
  col = analysis_colors[["volume_exploratory"]],
  xlab = "Positive recorded tumor volume (cm³; log scale)",
  ylab = "Apparent predicted probability",
  ylim = c(0, 1),
  main = "Exploratory tumor-volume relationship"
)
graphics::points(
  min(positive_volumes),
  volume_prediction_data$predicted_probability[1],
  pch = 18,
  cex = 1.5,
  col = "#B2182B"
)
graphics::legend(
  "topright",
  legend = c(
    "Positive recorded volume",
    "Zero recorded volume shown at left boundary"
  ),
  col = c(analysis_colors[["volume_exploratory"]], "#B2182B"),
  lwd = c(3, NA),
  pch = c(NA, 18),
  bty = "n"
)
graphics::mtext(
  paste(
    "Reference: age",
    reference_age,
    "| dpros",
    reference_dpros,
    "| PSA",
    reference_psa,
    "| Gleason",
    reference_gleason
  ),
  side = 3,
  line = 0.3,
  cex = 0.8
)
grDevices::dev.off()

sensitivity_metadata <- c(
  "Sensitivity and exploratory analyses",
  paste("Run date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste("Input:", input_path),
  paste("Input MD5:", unname(tools::md5sum(input_path))),
  paste("Saved assignments:", assignment_path),
  "Resampling: reused primary outcome-stratified 10-fold assignments across 20 repeats",
  paste(
    "dcaps sensitivity:",
    paste(
      deparse(sensitivity_definitions$dcaps_sensitivity$formula),
      collapse = " "
    )
  ),
  paste(
    "Zero-Gleason sensitivity:",
    paste(
      deparse(sensitivity_definitions$zero_gleason_exclusion$formula),
      collapse = " "
    )
  ),
  "Zero-Gleason exclusions: IDs 282 and 357",
  paste(
    "Volume exploratory:",
    paste(
      deparse(sensitivity_definitions$volume_exploratory$formula),
      collapse = " "
    )
  ),
  "Volume representation: indicator(vol > 0) plus log2(vol) among positive values",
  "Performance: ROC AUC, Brier score, log loss, calibration intercept, calibration slope",
  "Primary estimates: repeat-level out-of-fold performance",
  "Positive preferred_improvement favors the sensitivity or exploratory analysis",
  "Percentile ranges describe resampling variability, not confidence intervals",
  "Full-data performance and the volume prediction curve are apparent and secondary",
  "The primary modern model remains the prespecified main analysis"
)
writeLines(
  sensitivity_metadata,
  file.path(results_dir, "sensitivity_analysis_metadata.txt"),
  useBytes = TRUE
)

message("Sensitivity and exploratory analyses completed successfully.")
message(
  "Created ",
  nrow(sensitivity_predictions),
  " sensitivity out-of-fold predictions."
)
message("Results: ", results_dir)
message("Figures: ", figures_dir)
