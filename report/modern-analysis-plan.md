# Modern Analysis Plan

## Status and purpose

**Approved:** August 3, 2026

This document prespecifies the modern extension of the 2018 SDSU analysis. It is written before modern model fitting or performance comparison begins.

The purpose is to reduce avoidable analytical flexibility and to make later decisions, deviations, and results transparent. Any material change made after model results are examined will be documented as a deviation from this plan.

## Research question

The primary research question is:

> How accurately can baseline clinical measurements predict prostatic capsule penetration among patients diagnosed with prostate cancer?

The analysis concerns capsule penetration among patients who already have prostate cancer. It does not address cancer screening, initial diagnosis, prognosis, treatment selection, or clinical deployment.

## Objectives

The modern analysis will:

1. evaluate the out-of-sample predictive performance of the historical model;
2. evaluate a prespecified modern logistic regression model;
3. examine whether limited nonlinear modeling of PSA improves prediction;
4. evaluate discrimination, calibration, and overall predictive error;
5. assess the effect of including `dcaps` in a prespecified sensitivity analysis;
6. explore whether ultrasound-derived tumor volume (`vol`) adds predictive information; and
7. describe model stability and the uncertainty created by the small historical dataset.

## Data source

The analysis will use:

```text
data/processed/prostate_modern.rds
```

This file is generated reproducibly from the public `PCS` dataset distributed by the R package `lbreg`.

The dataset contains 380 observations. Capsule penetration is present for 153 observations and absent for 227 observations.

The available variables are:

- `id`: record identifier;
- `capsule`: capsule-penetration outcome;
- `age`: age in years;
- `race`: recorded race category;
- `dpros`: four-level digital rectal examination result;
- `dcaps`: baseline detection of capsular involvement during rectal examination;
- `psa`: prostate-specific antigen;
- `vol`: ultrasound-derived tumor volume; and
- `gleason`: total Gleason score.

The published values were modified to protect patient confidentiality. The analysis therefore cannot be treated as a contemporary clinical validation study.

## Outcome

The binary outcome is `capsule`:

- `0`: no prostatic capsule penetration;
- `1`: prostatic capsule penetration.

The predicted quantity is the probability that `capsule = 1`.

Outcome prevalence will be reported with a confidence interval, but the observed prevalence will not be altered through over-sampling, under-sampling, synthetic sampling, or class weighting.

## Analysis population

All 380 observations will be eligible for the primary modern analysis.

The primary predictors—age, `dpros`, PSA, and Gleason score—contain no missing values. Consequently, the primary analysis requires neither complete-case deletion nor imputation.

The three missing race values will not cause otherwise complete observations to be discarded because race is not included in the primary predictive models.

## Predictor policy

Predictors will be chosen according to the prespecified scientific roles below. They will not be screened using univariate p-values, stepwise selection, or examination of the modern validation results.

### Record identifier

`id` will never be used as a predictor.

### Age

Age will be modeled as a continuous linear term in the primary model. It will not be categorized.

### Digital rectal examination result

`dpros` will be modeled as a categorical factor with four levels. Level 1 will be the reference category.

Treating `dpros` as categorical avoids assuming equal spacing or a constant linear effect across its four clinical categories.

### PSA

PSA is positive but strongly right-skewed, ranging from 0.3 to 139.7. The primary model will use:

```r
log2(psa)
```

This permits interpretation in terms of a doubling of PSA and reduces the influence of the highest recorded values. No offset is required because all PSA values are greater than zero.

The historical benchmark will retain untransformed PSA so that it represents the 2018 model faithfully.

### Gleason score

Gleason score will be modeled as a continuous linear term, preserving the parsimonious representation used historically.

The dataset documentation permits values from 0 through 10, and the recorded data contain two zero values. These values will be retained in the primary analysis because they are part of the documented public dataset. A data-quality sensitivity analysis will refit the primary model after excluding the two zero-score observations.

### Race

Race will not be used as a predictor in the modern models.

The historical dataset records race using only two categories, contains only 36 observations in category 2, has three missing values, and was modified for confidentiality. These limitations prevent a responsible interpretation of race as a portable clinical predictor or an adequate assessment of model fairness.

Race will remain in the descriptive data documentation and historical reproduction.

### Baseline detection of capsular involvement

`dcaps` will be excluded from the primary model comparison because it directly records whether capsular involvement was detected during the baseline rectal examination.

It will be added in a prespecified sensitivity analysis to quantify how much this direct assessment changes predictive performance.

### Ultrasound-derived tumor volume

`vol` will be excluded from the primary analysis because it was unavailable in the 2018 SDSU dataset and is designated as a secondary exploratory predictor.

There are 167 recorded zero values. To avoid treating the difference between zero and a positive value as equivalent to proportional changes among positive volumes, the exploratory volume analysis will use a two-part representation:

1. an indicator for `vol > 0`; and
2. `log2(vol)` among positive values, coded as zero when `vol = 0`.

This representation uses two predictor parameters and will be clearly labeled exploratory.

## Prespecified models

No stepwise selection or automated predictor elimination will be used.

### Model 1: historical benchmark

```r
capsule ~ factor(dpros) + psa + gleason
```

This model reproduces the final model reported in 2018 and provides the historical benchmark.

### Model 2: primary modern model

```r
capsule ~ age + factor(dpros) + log2(psa) + gleason
```

This is the prespecified primary modern model. It updates the PSA representation, retains the clinically interpretable predictors, and includes age without relying on historical univariate significance testing.

### Model 3: limited nonlinear PSA model

```r
capsule ~ age + factor(dpros) + splines::ns(log2(psa), df = 3) + gleason
```

The spline will be a natural cubic spline with three degrees of freedom. Knot placement will be determined from each analysis portion within resampling, not from the held-out assessment portion.

This model tests a single prespecified nonlinear extension. It will not be expanded through a search over spline degrees of freedom, transformations, or interactions.

### `dcaps` sensitivity model

The primary modern model will be extended by adding categorical `dcaps`:

```r
capsule ~ age + factor(dpros) + log2(psa) + gleason + factor(dcaps)
```

Its results will be reported separately from the primary comparison.

### Exploratory tumor-volume model

The primary modern model will be extended with the prespecified two-part representation of `vol`.

This model will be labeled exploratory and will not replace the primary model based solely on an observed performance difference.

## Model estimation

All models will use logistic regression.

Predictor terms will be retained regardless of their p-values. The analysis will focus on predicted probabilities and model performance rather than statistical significance as a variable-selection mechanism.

The primary analysis will use maximum-likelihood logistic regression. Model convergence, extreme fitted probabilities, and evidence of separation will be checked. If instability or separation is detected, a ridge-penalized sensitivity analysis will be added and explicitly documented; it will not silently replace the primary estimator.

## Internal validation

The dataset is too small to support an efficient one-time division into training and test sets. All observations will therefore contribute to model development while out-of-sample performance is estimated using resampling.

The primary internal-validation procedure will be:

- stratified 10-fold cross-validation;
- 20 complete repeats;
- identical folds for every compared model; and
- a fixed, recorded random seed.

Within each repeat, every observation will receive one prediction from a model that was fitted without that observation. Metrics will be calculated from the complete set of out-of-fold predictions for that repeat.

All data-dependent preprocessing, including spline construction, will be estimated using only the analysis portion of each fold.

The distribution of each metric across the 20 repeats will be reported as resampling variability. Percentile ranges across repeats will not be described as formal confidence intervals because the repeated estimates are correlated.

After validation is completed, each prespecified model may be fitted to the full dataset to report its final coefficients and prediction equation. Apparent full-data performance will be clearly distinguished from cross-validated performance.

## Performance measures

### Primary performance measures

The principal evaluation will include:

- area under the receiver operating characteristic curve (ROC AUC) for discrimination;
- Brier score for overall probabilistic error;
- logarithmic loss for the quality of predicted probabilities;
- calibration intercept; and
- calibration slope.

A calibration intercept near 0 and calibration slope near 1 indicate better calibration. These quantities will be interpreted together with their resampling variability rather than as pass/fail tests.

### Graphical evaluation

The analysis will produce:

- an ROC curve based on out-of-fold predictions;
- a calibration plot comparing observed and predicted risk;
- a distribution plot of predicted probabilities by outcome; and
- paired comparisons of model metrics across resampling repeats.

Calibration graphics will display the observed event proportion, the 45-degree reference line, and a smooth calibration estimate. Sparse regions of predicted risk will be interpreted cautiously.

### Classification measures

Accuracy, sensitivity, specificity, positive predictive value, and negative predictive value depend on a probability threshold.

Because no clinically justified decision threshold is available for this historical educational dataset, these measures will not be used to choose the primary model. If shown for illustration, results will be reported across several explicitly labeled thresholds and will remain secondary.

Decision-curve analysis will not be presented because the project does not define a clinical decision, intervention, or defensible threshold range.

## Model comparison

The models will be compared using paired results from identical resampling splits.

The historical benchmark, primary modern model, and nonlinear PSA model will all be reported. The primary modern model is designated in advance and will not be replaced simply because another model achieves the numerically highest AUC.

The nonlinear PSA model will be considered meaningfully helpful only if it shows a consistent improvement in probabilistic performance or calibration without an important loss of discrimination. Small, unstable differences will be described as inconclusive.

The `dcaps` and tumor-volume extensions will be reported separately and will not participate in selecting the primary model.

## Model stability

Model stability will be assessed using the fitted coefficients and predicted risks across resamples.

The analysis will examine:

- variability of coefficient estimates;
- variability of model performance;
- the distribution of individual out-of-fold predicted probabilities; and
- whether a small number of observations produce unusually large changes.

Diagnostic observations will be investigated but will not be automatically deleted. Any exclusion analysis will be labeled as a sensitivity analysis.

## Sample-size and complexity assessment

The analysis contains 380 observations and 153 capsule-penetration outcomes. Model complexity will be expressed as the total number of estimated predictor parameters rather than by counting variable names.

A formal development sample-size assessment will be recorded before final model fitting, using the binary-outcome framework of Riley and colleagues. Because anticipated model performance is uncertain, the assessment will include a conservative range rather than a single optimistic assumption.

This assessment will inform the interpretation of model stability and precision. It will not be used retrospectively to search for a model that passes a chosen rule.

## Missing-data policy

No primary predictor is missing, so the primary analysis will not impute data.

If a later approved analysis introduces a predictor with missing values, imputation must occur separately within each resampling analysis fold. Imputation performed once using the complete dataset would constitute information leakage and is prohibited.

The three missing race values will not be imputed because race is excluded from modern predictive modeling.

## Reproducibility and leakage controls

The modern workflow will:

- use a fixed random seed;
- save the resampling assignments;
- apply identical resamples to all models;
- perform preprocessing within resampling folds;
- preserve outcome prevalence rather than altering class balance;
- avoid stepwise selection and univariate screening;
- avoid using the validation results to create additional unreported candidate models;
- save machine-readable performance estimates and out-of-fold predictions;
- save figures through scripts rather than interactively; and
- record package versions through `renv`.

Model selection, transformation choices, and performance measures will not be changed after results are examined unless the change is documented as exploratory or as a protocol deviation.

## Planned outputs

The modern analysis will be implemented in separate scripts, anticipated to include:

```text
R/03_validate_inputs.R
R/04_modern_analysis.R
R/05_sensitivity_analyses.R
```

Machine-readable outputs will be stored under:

```text
results/modern/
```

Figures will be stored under:

```text
figures/modern/
```

A later findings document will distinguish prespecified primary results, sensitivity results, exploratory results, and any deviations from this protocol.

## Interpretation limits

This is an internal validation study using a small historical dataset from a single clinical setting. It does not provide external validation, contemporary transportability, or evidence that the models would improve patient outcomes.

The data were modified for confidentiality, provide limited demographic information, and do not support a meaningful fairness assessment. Any apparent predictive accuracy may be specific to this dataset and clinical context.

The project is intended for statistical education and methodological demonstration. It is not a medical device and must not be used for clinical decisions.

## Reporting guidance

The modern analysis and final report will be organized with reference to:

1. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*. 2024;385:e078378. <https://doi.org/10.1136/bmj-2023-078378>
2. Moons KGM, Damen JAA, Kaul T, et al. PROBAST+AI: an updated quality, risk of bias, and applicability assessment tool for prediction models using regression or artificial intelligence methods. *BMJ*. 2025;388:e082505. <https://doi.org/10.1136/bmj-2024-082505>
3. Riley RD, Snell KIE, Ensor J, et al. Minimum sample size for developing a multivariable prediction model: Part II—binary and time-to-event outcomes. *Statistics in Medicine*. 2019;38:1276–1296. <https://doi.org/10.1002/sim.7992>

These references guide transparent reporting, risk-of-bias awareness, validation, and model-complexity decisions. They do not imply that this historical educational dataset is suitable for clinical deployment.
