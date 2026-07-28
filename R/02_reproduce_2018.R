# Reproduce the principal analyses reported in the 2018 SDSU project.
#
# Run from the project root after acquiring the data:
#   source("R/01_acquire_data.R")
#   source("R/02_reproduce_2018.R")
#
# Input (generated and ignored by Git):
#   data/processed/prostate_2018.rds
#
# Outputs:
#   results/historical/
#   figures/historical/

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
  "report/historical-reproduction-plan.md"
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
  requireNamespace("MASS", quietly = TRUE),
  "Package 'MASS' is required. Restore the project with renv::restore()."
)

input_path <- file.path("data", "processed", "prostate_2018.rds")

check_true(
  file.exists(input_path),
  paste(
    "Historical dataset not found:",
    input_path,
    "Run source(\"R/01_acquire_data.R\") first."
  )
)

prostate_all <- readRDS(input_path)

expected_names <- c(
  "id", "capsule", "age", "race",
  "dpros", "dcaps", "psa", "gleason"
)

check_true(is.data.frame(prostate_all), "Historical input must be a data frame.")
check_true(nrow(prostate_all) == 380L, "Expected 380 historical observations.")
check_true(
  identical(names(prostate_all), expected_names),
  "Historical input has unexpected variables or variable order."
)
check_true(
  identical(prostate_all$id, seq_len(380L)),
  "Historical record identifiers must be 1 through 380 in order."
)
check_true(
  identical(which(is.na(prostate_all$race)), c(22L, 46L, 252L)),
  "Historical race values must be missing only in rows 22, 46, and 252."
)
check_true(
  identical(as.integer(table(prostate_all$capsule)), c(227L, 153L)),
  "Historical capsule outcome counts must be 227 and 153."
)

complete_rows <- stats::complete.cases(prostate_all)
prostate_complete <- prostate_all[complete_rows, , drop = FALSE]

check_true(nrow(prostate_complete) == 377L, "Expected 377 complete cases.")
check_true(
  identical(prostate_complete$id, setdiff(seq_len(380L), c(22L, 46L, 252L))),
  "Unexpected observations were removed by complete-case deletion."
)

results_dir <- file.path("results", "historical")
figures_dir <- file.path("figures", "historical")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

options(contrasts = c("contr.treatment", "contr.poly"))

# Dataset description --------------------------------------------------------

dataset_counts <- data.frame(
  measure = c(
    "Original observations",
    "Capsule penetration: no",
    "Capsule penetration: yes",
    "Missing race values",
    "Complete cases"
  ),
  value = c(
    nrow(prostate_all),
    unname(table(prostate_all$capsule)[["0"]]),
    unname(table(prostate_all$capsule)[["1"]]),
    sum(is.na(prostate_all$race)),
    nrow(prostate_complete)
  )
)

continuous_variables <- c("age", "psa", "gleason")

descriptive_statistics <- do.call(
  rbind,
  lapply(continuous_variables, function(variable) {
    x <- prostate_all[[variable]]
    data.frame(
      variable = variable,
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      minimum = min(x, na.rm = TRUE),
      p25 = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
      p75 = unname(stats::quantile(x, 0.75, na.rm = TRUE)),
      maximum = max(x, na.rm = TRUE)
    )
  })
)

write_csv(dataset_counts, file.path(results_dir, "dataset_counts.csv"))
write_csv(
  descriptive_statistics,
  file.path(results_dir, "descriptive_statistics.csv")
)

contingency_to_data <- function(variable) {
  tab <- table(
    predictor = prostate_complete[[variable]],
    capsule = prostate_complete$capsule
  )
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(out)[1] <- variable
  out
}

contingency_tables <- do.call(
  rbind,
  lapply(c("race", "dpros", "dcaps"), function(variable) {
    out <- contingency_to_data(variable)
    out$variable <- variable
    names(out)[1] <- "level"
    out[c("variable", "level", "capsule", "Freq")]
  })
)

names(contingency_tables)[names(contingency_tables) == "Freq"] <- "count"
write_csv(
  contingency_tables,
  file.path(results_dir, "contingency_tables.csv")
)

# Exploratory figures --------------------------------------------------------

grDevices::png(
  filename = file.path(figures_dir, "continuous_boxplots.png"),
  width = 1800,
  height = 650,
  res = 150
)
graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
graphics::boxplot(
  psa ~ factor(capsule),
  data = prostate_complete,
  varwidth = TRUE,
  col = "cyan",
  xlab = "Penetrated the capsule",
  ylab = "PSA",
  main = "PSA"
)
graphics::boxplot(
  age ~ factor(capsule),
  data = prostate_complete,
  varwidth = TRUE,
  col = "cyan",
  xlab = "Penetrated the capsule",
  ylab = "Age (years)",
  main = "Age"
)
graphics::boxplot(
  gleason ~ factor(capsule),
  data = prostate_complete,
  varwidth = TRUE,
  col = "cyan",
  xlab = "Penetrated the capsule",
  ylab = "Gleason score",
  main = "Gleason score"
)
grDevices::dev.off()

historical_correlations <- stats::cor(prostate_complete)

grDevices::png(
  filename = file.path(figures_dir, "historical_correlation_plot.png"),
  width = 1050,
  height = 900,
  res = 150
)
graphics::par(mar = c(7, 7, 3, 2))
cor_palette <- grDevices::colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)
graphics::image(
  x = seq_len(ncol(historical_correlations)),
  y = seq_len(nrow(historical_correlations)),
  z = t(historical_correlations[nrow(historical_correlations):1, ]),
  col = cor_palette,
  zlim = c(-1, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Historical numeric-code correlation plot"
)
graphics::axis(
  1,
  at = seq_len(ncol(historical_correlations)),
  labels = colnames(historical_correlations),
  las = 2
)
graphics::axis(
  2,
  at = seq_len(nrow(historical_correlations)),
  labels = rev(rownames(historical_correlations)),
  las = 2
)
graphics::box()
grDevices::dev.off()

write_csv(
  data.frame(
    variable = rownames(historical_correlations),
    historical_correlations,
    row.names = NULL,
    check.names = FALSE
  ),
  file.path(results_dir, "historical_correlation_matrix.csv")
)

# Univariate analyses --------------------------------------------------------

t_test_results <- do.call(
  rbind,
  lapply(c("psa", "age", "gleason"), function(variable) {
    test <- stats::t.test(
      prostate_complete[[variable]][prostate_complete$capsule == 0L],
      prostate_complete[[variable]][prostate_complete$capsule == 1L]
    )
    data.frame(
      variable = variable,
      mean_capsule_0 = unname(test$estimate[[1]]),
      mean_capsule_1 = unname(test$estimate[[2]]),
      mean_difference_0_minus_1 = unname(diff(rev(test$estimate))),
      statistic = unname(test$statistic),
      degrees_freedom = unname(test$parameter),
      p_value = test$p.value,
      confidence_lower = test$conf.int[[1]],
      confidence_upper = test$conf.int[[2]],
      method = test$method
    )
  })
)

write_csv(t_test_results, file.path(results_dir, "two_sample_t_tests.csv"))

chi_square_results <- do.call(
  rbind,
  lapply(c("race", "dpros", "dcaps"), function(variable) {
    test <- suppressWarnings(
      stats::chisq.test(
        x = prostate_complete[[variable]],
        y = prostate_complete$capsule
      )
    )
    data.frame(
      variable = variable,
      statistic = unname(test$statistic),
      degrees_freedom = unname(test$parameter),
      p_value = test$p.value,
      method = test$method
    )
  })
)

write_csv(chi_square_results, file.path(results_dir, "chi_square_tests.csv"))

univariate_glm_results <- do.call(
  rbind,
  lapply(c("psa", "age", "race", "dpros", "dcaps", "gleason"), function(variable) {
    formula <- stats::reformulate(variable, response = "capsule")
    fit <- stats::glm(
      formula,
      family = stats::binomial(link = "logit"),
      data = prostate_complete
    )
    coefficients <- summary(fit)$coefficients
    data.frame(
      predictor = variable,
      term = rownames(coefficients),
      estimate = coefficients[, "Estimate"],
      std_error = coefficients[, "Std. Error"],
      z_value = coefficients[, "z value"],
      p_value = coefficients[, "Pr(>|z|)"],
      row.names = NULL
    )
  })
)

write_csv(
  univariate_glm_results,
  file.path(results_dir, "univariate_logistic_models.csv")
)

# Historical model selection ------------------------------------------------

candidate_formula <- capsule ~
  dcaps + dpros + psa + gleason +
  dpros:dcaps + dpros:psa + dpros:gleason +
  dcaps:psa + dcaps:gleason + psa:gleason

candidate_model <- stats::glm(
  candidate_formula,
  family = stats::binomial(link = "logit"),
  data = prostate_complete
)

step_model <- MASS::stepAIC(
  candidate_model,
  direction = "both",
  trace = FALSE
)

reported_final_formula <- capsule ~ factor(dpros) + psa + gleason
reported_final_model <- stats::glm(
  reported_final_formula,
  family = stats::binomial(link = "logit"),
  data = prostate_complete
)

formula_text <- function(model) {
  paste(deparse(stats::formula(model)), collapse = " ")
}

model_selection_summary <- data.frame(
  model = c("candidate", "stepAIC_selected", "reported_final"),
  formula = c(
    formula_text(candidate_model),
    formula_text(step_model),
    formula_text(reported_final_model)
  ),
  degrees_freedom = c(
    attr(stats::logLik(candidate_model), "df"),
    attr(stats::logLik(step_model), "df"),
    attr(stats::logLik(reported_final_model), "df")
  ),
  aic = c(
    stats::AIC(candidate_model),
    stats::AIC(step_model),
    stats::AIC(reported_final_model)
  ),
  observations = c(
    stats::nobs(candidate_model),
    stats::nobs(step_model),
    stats::nobs(reported_final_model)
  )
)

step_matches_reported <- identical(
  attr(stats::terms(step_model), "term.labels"),
  attr(stats::terms(reported_final_model), "term.labels")
)

write_csv(
  model_selection_summary,
  file.path(results_dir, "model_selection_summary.csv")
)

# Reported final model -------------------------------------------------------

coefficient_matrix <- summary(reported_final_model)$coefficients
wald_multiplier <- stats::qnorm(0.975)
wald_confidence <- cbind(
  coefficient_matrix[, "Estimate"] -
    wald_multiplier * coefficient_matrix[, "Std. Error"],
  coefficient_matrix[, "Estimate"] +
    wald_multiplier * coefficient_matrix[, "Std. Error"]
)
profile_confidence <- suppressMessages(
  stats::confint(reported_final_model, level = 0.95)
)

final_model_results <- data.frame(
  term = rownames(coefficient_matrix),
  estimate = coefficient_matrix[, "Estimate"],
  std_error = coefficient_matrix[, "Std. Error"],
  z_value = coefficient_matrix[, "z value"],
  p_value = coefficient_matrix[, "Pr(>|z|)"],
  wald_coefficient_ci_lower = wald_confidence[, 1],
  wald_coefficient_ci_upper = wald_confidence[, 2],
  profile_coefficient_ci_lower = profile_confidence[, 1],
  profile_coefficient_ci_upper = profile_confidence[, 2],
  odds_ratio = exp(coefficient_matrix[, "Estimate"]),
  odds_ratio_ci_lower = exp(profile_confidence[, 1]),
  odds_ratio_ci_upper = exp(profile_confidence[, 2]),
  row.names = NULL
)

write_csv(
  final_model_results,
  file.path(results_dir, "final_logistic_model.csv")
)

reported_benchmarks <- data.frame(
  term = c(
    "(Intercept)",
    "factor(dpros)2",
    "factor(dpros)3",
    "factor(dpros)4",
    "psa",
    "gleason"
  ),
  reported_estimate = c(-8.14, 0.77, 1.55, 1.43, 0.03, 1.00),
  reported_coefficient_ci_lower = c(-10.22, 0.07, 0.83, 0.55, 0.01, 0.68),
  reported_coefficient_ci_upper = c(-6.07, 1.47, 2.28, 2.31, 0.05, 1.31),
  reported_odds_ratio = c(0.00, 2.17, 4.73, 4.18, 1.03, 2.71),
  reported_or_ci_lower = c(0.00, 1.09, 2.32, 1.75, 1.01, 2.00),
  reported_or_ci_upper = c(0.00, 4.43, 10.00, 10.25, 1.05, 3.76),
  reported_p_value = c(0.00, 0.31, 0.00, 0.00, 0.00, 0.00)
)

final_comparison <- merge(
  reported_benchmarks,
  final_model_results,
  by = "term",
  all.x = TRUE,
  sort = FALSE
)

final_comparison$estimate_status <- ifelse(
  round(final_comparison$estimate, 2) == final_comparison$reported_estimate,
  "Reproduced",
  "Different"
)
final_comparison$coefficient_ci_status <- ifelse(
  round(final_comparison$wald_coefficient_ci_lower, 2) ==
    final_comparison$reported_coefficient_ci_lower &
    round(final_comparison$wald_coefficient_ci_upper, 2) ==
      final_comparison$reported_coefficient_ci_upper,
  "Reproduced",
  "Different"
)
final_comparison$odds_ratio_status <- ifelse(
  round(final_comparison$odds_ratio, 2) ==
    final_comparison$reported_odds_ratio,
  "Reproduced",
  "Different"
)
final_comparison$ci_status <- ifelse(
  round(final_comparison$odds_ratio_ci_lower, 2) ==
    final_comparison$reported_or_ci_lower &
    round(final_comparison$odds_ratio_ci_upper, 2) ==
      final_comparison$reported_or_ci_upper,
  "Reproduced",
  "Different"
)
final_comparison$p_value_note <- ifelse(
  final_comparison$term == "factor(dpros)2",
  paste(
    "The report lists 0.31; the model-derived p-value is",
    format(final_comparison$p_value, digits = 6)
  ),
  "The report rounded the p-value to two decimals."
)

write_csv(
  final_comparison,
  file.path(results_dir, "reported_vs_reproduced_final_model.csv")
)

model_level_comparison <- data.frame(
  measure = c("Observations", "AIC"),
  reported_value = c(377, 393.22),
  reproduced_value = c(
    stats::nobs(reported_final_model),
    stats::AIC(reported_final_model)
  ),
  status = c(
    ifelse(stats::nobs(reported_final_model) == 377L, "Reproduced", "Different"),
    ifelse(
      round(stats::AIC(reported_final_model), 2) == 393.22,
      "Reproduced",
      "Different"
    )
  )
)

write_csv(
  model_level_comparison,
  file.path(results_dir, "reported_vs_reproduced_model_level.csv")
)

# Diagnostics and conservative sensitivity -------------------------------------

fitted_probabilities <- stats::fitted(reported_final_model)
standardized_pearson <- stats::rstandard(
  reported_final_model,
  type = "pearson"
)
cook_distance <- stats::cooks.distance(reported_final_model)
leverage <- stats::hatvalues(reported_final_model)
number_parameters <- length(stats::coef(reported_final_model))
number_observations <- stats::nobs(reported_final_model)

diagnostic_data <- data.frame(
  row = seq_len(nrow(prostate_complete)),
  id = prostate_complete$id,
  fitted_probability = fitted_probabilities,
  standardized_pearson = standardized_pearson,
  cook_distance = cook_distance,
  leverage = leverage,
  residual_flag = abs(standardized_pearson) > 2,
  cook_flag = cook_distance > 4 / number_observations,
  leverage_flag = leverage > 2 * number_parameters / number_observations
)
diagnostic_data$diagnostic_flag <- with(
  diagnostic_data,
  residual_flag | cook_flag | leverage_flag
)

write_csv(
  diagnostic_data,
  file.path(results_dir, "model_diagnostics.csv")
)

flagged_data <- diagnostic_data[diagnostic_data$diagnostic_flag, , drop = FALSE]
write_csv(
  flagged_data,
  file.path(results_dir, "diagnostically_flagged_observations.csv")
)

grDevices::png(
  filename = file.path(figures_dir, "final_model_diagnostics.png"),
  width = 1500,
  height = 700,
  res = 150
)
graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
graphics::plot(
  fitted_probabilities,
  standardized_pearson,
  xlab = "Estimated probability",
  ylab = "Standardized Pearson residual",
  main = "Residuals vs. estimated probability",
  pch = 19,
  col = ifelse(diagnostic_data$diagnostic_flag, "#B2182B", "#2166AC")
)
graphics::abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = "grey40")
graphics::plot(
  leverage,
  cook_distance,
  xlab = "Leverage",
  ylab = "Cook's distance",
  main = "Cook's distance vs. leverage",
  pch = 19,
  col = ifelse(diagnostic_data$diagnostic_flag, "#B2182B", "#2166AC")
)
graphics::abline(
  h = 4 / number_observations,
  v = 2 * number_parameters / number_observations,
  lty = 2,
  col = "grey40"
)
grDevices::dev.off()

# This intentionally conservative stress test removes the union of all three
# diagnostic screening sets. A flag is not, by itself, a reason to exclude an
# observation from the primary historical reproduction.
retained_rows <- !diagnostic_data$diagnostic_flag
combined_flag_sensitivity_model <- stats::update(
  reported_final_model,
  data = prostate_complete[retained_rows, , drop = FALSE]
)

full_coefficients <- stats::coef(reported_final_model)
sensitivity_coefficients <- stats::coef(combined_flag_sensitivity_model)

combined_flag_sensitivity <- data.frame(
  term = names(full_coefficients),
  full_estimate = unname(full_coefficients),
  sensitivity_estimate = unname(
    sensitivity_coefficients[names(full_coefficients)]
  ),
  absolute_difference = abs(
    unname(full_coefficients) -
      unname(sensitivity_coefficients[names(full_coefficients)])
  ),
  relative_difference = abs(
    (
      unname(sensitivity_coefficients[names(full_coefficients)]) -
        unname(full_coefficients)
    ) / unname(full_coefficients)
  ),
  full_observations = stats::nobs(reported_final_model),
  sensitivity_observations = stats::nobs(combined_flag_sensitivity_model),
  interpretation = paste(
    "Conservative stress test removing the union of all diagnostic flags;",
    "not a recommendation to exclude these observations."
  )
)

write_csv(
  combined_flag_sensitivity,
  file.path(results_dir, "conservative_combined_flag_sensitivity.csv")
)

# Reproduction summary -------------------------------------------------------

reproduced_coefficient_count <- sum(
  final_comparison$estimate_status == "Reproduced"
)
reproduced_coefficient_ci_count <- sum(
  final_comparison$coefficient_ci_status == "Reproduced"
)
reproduced_or_count <- sum(
  final_comparison$odds_ratio_status == "Reproduced"
)
reproduced_ci_count <- sum(final_comparison$ci_status == "Reproduced")

summary_lines <- c(
  "Historical reproduction summary",
  paste("Run date:", Sys.Date()),
  paste("R version:", R.version.string),
  "Input: data/processed/prostate_2018.rds",
  paste("Original observations:", nrow(prostate_all)),
  paste("Complete-case observations:", nrow(prostate_complete)),
  paste("Candidate model AIC:", format(stats::AIC(candidate_model), digits = 10)),
  paste("stepAIC-selected formula:", formula_text(step_model)),
  paste("stepAIC-selected AIC:", format(stats::AIC(step_model), digits = 10)),
  paste("Reported final formula:", formula_text(reported_final_model)),
  paste(
    "Stepwise model matches reported final model:",
    ifelse(step_matches_reported, "yes", "no")
  ),
  paste(
    "Reported final model AIC:",
    format(stats::AIC(reported_final_model), digits = 10)
  ),
  paste(
    "Reported final AIC reproduces at two decimals:",
    ifelse(
      round(stats::AIC(reported_final_model), 2) == 393.22,
      "yes",
      "no"
    )
  ),
  paste(
    "Coefficient estimates reproduced at two decimals:",
    paste0(reproduced_coefficient_count, "/", nrow(final_comparison))
  ),
  paste(
    "Wald coefficient confidence intervals reproduced at two decimals:",
    paste0(reproduced_coefficient_ci_count, "/", nrow(final_comparison))
  ),
  paste(
    "Odds ratios reproduced at two decimals:",
    paste0(reproduced_or_count, "/", nrow(final_comparison))
  ),
  paste(
    "Odds-ratio confidence intervals reproduced at two decimals:",
    paste0(reproduced_ci_count, "/", nrow(final_comparison))
  ),
  paste(
    "Observations with |standardized Pearson residual| > 2:",
    sum(diagnostic_data$residual_flag)
  ),
  paste(
    "Observations with Cook's distance > 4/n:",
    sum(diagnostic_data$cook_flag)
  ),
  paste(
    "Observations with leverage > 2p/n:",
    sum(diagnostic_data$leverage_flag)
  ),
  paste(
    "Observations flagged by at least one diagnostic screen:",
    nrow(flagged_data)
  ),
  paste(
    paste(
      "Largest absolute coefficient change in conservative",
      "combined-flag removal stress test:"
    ),
    format(max(combined_flag_sensitivity$absolute_difference), digits = 8)
  ),
  paste(
    "Diagnostic interpretation:",
    "screening flags identify observations for review and do not, by themselves,"
  ),
  "justify excluding observations from the primary historical reproduction.",
  paste(
    "Important discrepancy: reported p-value for factor(dpros)2 = 0.31;",
    "the model-derived p-value is",
    format(
      final_model_results$p_value[
        final_model_results$term == "factor(dpros)2"
      ],
      digits = 8
    )
  ),
  paste(
    "Interpretation boundary:",
    "results concern capsule penetration among patients with prostate cancer,"
  ),
  "not initial prostate cancer detection."
)

writeLines(
  summary_lines,
  file.path(results_dir, "reproduction_summary.txt"),
  useBytes = TRUE
)

message("Historical reproduction completed successfully.")
message("Results: ", results_dir)
message("Figures: ", figures_dir)
