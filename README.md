# Predicting Prostatic Capsule Penetration

## A Modern Extension of a 2018 Analysis

## Project overview

This project reproduces and extends a logistic regression analysis originally completed in 2018 for STAT 696: Statistical Communication in Data Science at San Diego State University.

The study examines whether baseline clinical measurements can predict prostatic capsule penetration among patients diagnosed with prostate cancer.

The project has two principal components:

1. a faithful reproduction of the 2018 analysis; and
2. a prespecified modern extension emphasizing out-of-sample prediction, calibration, model stability, and reproducibility.

The modern extension follows a documented analysis plan created before model fitting. It uses repeated cross-validation and compares a limited set of interpretable logistic regression models.

## Final report

The completed integrated report is available in three formats:

- [Self-contained HTML report](report/final-report.html)
- [PDF report](report/final-report.pdf)
- [Quarto source](report/final-report.qmd)

The HTML file contains its figures and supporting resources in a single file.
When viewing the repository on GitHub, download the HTML file and open it in a
web browser. GitHub can preview the PDF directly. The Quarto source is the
authoritative, reproducible report source.

## Research question

> How accurately can baseline clinical measurements predict prostatic capsule penetration among patients diagnosed with prostate cancer?

### Secondary questions

1. How does including baseline detection of capsular involvement (`dcaps`) affect predictive performance?
2. Does ultrasound-derived tumor volume (`vol`) provide additional predictive information?

The `dcaps` analysis is treated as a sensitivity analysis because the measurement directly assesses possible capsular involvement during the baseline rectal examination.

Tumor volume is treated as a secondary exploratory analysis because it was present in the canonical public dataset but not in the dataset used for the 2018 SDSU project.

## Project objectives

1. Reproduce and verify the 2018 logistic regression analysis.
2. Develop a fully documented and reproducible R workflow.
3. Estimate predictive performance using repeated cross-validation.
4. Evaluate discrimination, calibration, and overall predictive error.
5. Compare the historical model with a limited number of interpretable modern alternatives.
6. Assess the effect of including `dcaps`.
7. Examine the stability of results after excluding two unusual zero Gleason scores.
8. Explore the additional predictive value of ultrasound-derived tumor volume.

## Data source

The project uses the Prostate Cancer Study dataset described by Hosmer and Lemeshow in *Applied Logistic Regression*.

The data were collected by Dr. Donn Young at The Ohio State University Comprehensive Cancer Center. Published values were modified to protect patient confidentiality.

The canonical dataset is publicly distributed as `PCS` in the R package `lbreg`. It contains 380 observations and nine variables:

| Variable | Meaning |
|---|---|
| `id` | Record identification code |
| `tumor` | Tumor penetration of the prostatic capsule |
| `age` | Age in years |
| `race` | Recorded race category |
| `dpros` | Digital rectal examination result |
| `dcaps` | Detection of capsular involvement during rectal examination |
| `psa` | Prostate-specific antigen |
| `vol` | Tumor volume obtained by ultrasound |
| `gleason` | Total Gleason score |

During reproducible acquisition, `tumor` is renamed to `capsule`.

The eight-variable SDSU dataset is an exact subset of the canonical dataset after:

1. renaming `tumor` to `capsule`; and
2. excluding `vol`.

All remaining values, observation order, and missing-value locations match the public dataset.

Generated data files are not committed to Git. They are recreated from the documented public source by `R/01_acquire_data.R`.

## Historical reproduction

The historical reproduction used the same eight-variable structure available in 2018 and preserved the original complete-case population of 377 observations.

The central numerical results were highly reproducible:

- dataset counts and missing-value locations matched;
- descriptive statistics matched;
- principal univariate conclusions matched;
- all six final-model coefficients reproduced;
- all coefficient confidence intervals reproduced;
- all odds ratios and odds-ratio confidence intervals reproduced; and
- the reported final-model AIC of 393.22 reproduced.

Two important qualifications were identified.

### Corrected p-value

The 2018 report listed a p-value of `0.31` for `dpros = 2`. The reproduced model gives:

```text
p = 0.02997
```

Because the corresponding coefficient, standard error, confidence interval, and odds ratio all reproduce, the reported value is almost certainly a typographical error.

### Model-selection ambiguity

The reproduced bidirectional stepwise procedure selected:

```r
capsule ~ dpros + psa + gleason + dpros:psa
```

The final model reported in 2018 was:

```r
capsule ~ factor(dpros) + psa + gleason
```

The available historical materials do not document why the interaction was removed or why `dpros` was changed from a numeric term to a categorical factor. The final model is reproducible, but it is not the direct result of the documented stepwise procedure.

Detailed historical findings are available in:

- [Historical reproduction plan](report/historical-reproduction-plan.md)
- [Historical reproduction findings](report/historical-reproduction-findings.md)

## Modern analysis

The modern analysis used all 380 observations and compared three prespecified models.

### Historical benchmark

```r
capsule ~ factor(dpros) + psa + gleason
```

### Primary modern model

```r
capsule ~ age + factor(dpros) + log2(psa) + gleason
```

### Nonlinear PSA model

```r
capsule ~ age +
  factor(dpros) +
  splines::ns(log2(psa), df = 3) +
  gleason
```

Predictive performance was estimated using outcome-stratified 10-fold cross-validation repeated 20 times. Every model used the same saved resampling assignments.

### Cross-validated performance

| Model | ROC AUC | Brier score | Log loss | Calibration intercept | Calibration slope |
|---|---:|---:|---:|---:|---:|
| Historical benchmark | 0.809 | 0.17315 | 0.51878 | −0.002 | 0.920 |
| Primary modern | 0.810 | 0.17310 | 0.51925 | −0.003 | 0.914 |
| Nonlinear PSA | 0.805 | 0.17461 | 0.52266 | −0.003 | 0.884 |

The historical benchmark and primary modern model performed almost identically. The modern model did not materially outperform the simpler historical benchmark.

The nonlinear PSA model increased complexity and prediction instability without improving out-of-sample performance.

Calibration intercepts were close to zero, while calibration slopes below one indicated modest overfitting. Overfitting was greatest for the nonlinear model.

## Sensitivity and exploratory findings

### Adding `dcaps`

Adding `dcaps` did not provide meaningful or stable incremental predictive value. Small improvements in some metrics were offset by slightly worse log loss and calibration slope.

### Excluding zero Gleason scores

The dataset contains two records with Gleason score equal to zero. Excluding these records produced essentially no change in discrimination, prediction error, calibration, or fitted coefficients.

The primary conclusions are not sensitive to these two observations.

### Adding tumor volume

The exploratory volume model used a two-part representation:

1. an indicator distinguishing zero from positive recorded volume; and
2. `log2(volume)` among positive values.

The model produced small average improvements in AUC, Brier score, and log loss but worse calibration slope. The fitted relationship was nearly flat among positive volumes, and the meaning of the 167 recorded zero values is uncertain.

The tumor-volume findings are therefore inconclusive and remain exploratory.

Detailed modern results are available in:

- [Modern analysis plan](report/modern-analysis-plan.md)
- [Sample-size and model-complexity assessment](report/sample-size-assessment.md)
- [Modern analysis findings](report/modern-analysis-findings.md)

## Main conclusion

Digital rectal examination findings, PSA, and Gleason score contain meaningful predictive information about prostatic capsule penetration in this historical dataset.

However, the analysis does not show that a more complicated model predicts better. The historical benchmark and prespecified primary modern model performed essentially equally, while nonlinear PSA modeling performed somewhat worse.

The results favor a parsimonious and interpretable logistic regression model.

The reported performance represents internal validation only. The models have not been externally validated and should not be considered contemporary clinical prediction tools.

## Reproducing the analysis

### Requirements

- R
- Quarto for rendering the integrated report
- the package versions recorded in `renv.lock`
- internet access during initial package restoration and data acquisition

Clone the repository and open the R project:

```text
prostatic-capsular-penetration.Rproj
```

Restore the recorded package environment:

```r
renv::restore()
```

Run the scripts from the project root in numerical order:

```r
source("R/01_acquire_data.R")
source("R/02_reproduce_2018.R")
source("R/03_validate_inputs.R")
source("R/04_modern_analysis.R")
source("R/05_sensitivity_analyses.R")
```

Alternatively, the scripts can be run non-interactively:

```sh
Rscript --vanilla R/01_acquire_data.R
Rscript --vanilla R/02_reproduce_2018.R
Rscript --vanilla R/03_validate_inputs.R
Rscript --vanilla R/04_modern_analysis.R
Rscript --vanilla R/05_sensitivity_analyses.R
```

The workflow will stop if expected files, dataset structure, values, transformations, resampling assignments, or model outputs fail validation.

Render both reader-ready report formats after completing the analysis:

```sh
quarto render report/final-report.qmd
```

The rendered outputs are written to:

```text
report/final-report.html
report/final-report.pdf
```

The repository includes approved rendered copies for convenient reading.
The `.qmd` file remains the authoritative source.

## Repository structure

```text
.
├── R/
│   ├── 01_acquire_data.R
│   ├── 02_reproduce_2018.R
│   ├── 03_validate_inputs.R
│   ├── 04_modern_analysis.R
│   └── 05_sensitivity_analyses.R
├── data/
│   ├── raw/
│   │   └── README.md
│   └── processed/
│       └── README.md
├── figures/
│   ├── historical/
│   └── modern/
├── report/
│   ├── final-report.qmd
│   ├── final-report.html
│   ├── final-report.pdf
│   ├── historical-reproduction-plan.md
│   ├── historical-reproduction-findings.md
│   ├── modern-analysis-plan.md
│   ├── sample-size-assessment.md
│   └── modern-analysis-findings.md
├── results/
│   ├── historical/
│   └── modern/
├── renv/
├── renv.lock
└── prostatic-capsular-penetration.Rproj
```

## Reproducibility design

The project separates the workflow into documented stages:

1. public data acquisition;
2. historical reproduction;
3. modern input validation;
4. primary modern analysis; and
5. sensitivity and exploratory analyses.

Additional reproducibility features include:

- generated data excluded from version control;
- source and transformation checks;
- stored dataset checksums;
- explicit input-validation results;
- fixed and recorded random seed;
- saved resampling assignments;
- identical folds for model comparisons;
- preprocessing performed within resampling where required;
- machine-readable predictions and performance estimates;
- scripted, non-interactive figures;
- deterministic repeated runs; and
- package versions recorded through `renv`.

## Scope and limitations

This study evaluates prediction of prostatic capsule penetration among patients already diagnosed with prostate cancer.

It does not evaluate:

- prostate cancer screening;
- initial prostate cancer diagnosis;
- treatment selection;
- prognosis;
- cancer-specific mortality;
- fairness across contemporary populations;
- external validity;
- clinical utility; or
- deployment of a prediction tool.

The dataset is historical, relatively small, and derived from a single clinical setting. Published values were modified for confidentiality. Race is recorded using only two categories, and the meaning of zero tumor-volume values is uncertain.

These limitations restrict generalizability and prevent responsible clinical deployment.

## Clinical-use limitation

This repository is intended for statistical education and methodological demonstration.

It is not a medical device and must not be used to make clinical decisions.

## Licensing

The repository uses separate licenses for software and original project
content:

- R scripts and other software code are licensed under the MIT License.
- Documentation, report text, original figures, and project-produced
  analytical results are licensed under the Creative Commons Attribution 4.0
  International License.

The source dataset and other third-party materials remain governed by their
original licenses and terms. See [LICENSE.md](LICENSE.md) for details.

## Historical context and suggested citation

The original report was completed by Navid Hedayati in 2018 as part of STAT 696 at San Diego State University.

The historical work may be cited as:

> Hedayati, N. (2018). *Detection of Prostate Cancer: Data Analysis Report 2*. Unpublished course report, STAT 696, San Diego State University.

Citation metadata for the complete modern repository is provided in
[`CITATION.cff`](CITATION.cff).

## References

Hosmer, D. W., & Lemeshow, S. (2000). *Applied Logistic Regression* (2nd ed.). John Wiley & Sons.

Collins, G. S., Moons, K. G. M., Dhiman, P., et al. (2024). TRIPOD+AI statement: Updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*, 385, e078378. <https://doi.org/10.1136/bmj-2023-078378>

Moons, K. G. M., Damen, J. A. A., Kaul, T., et al. (2025). PROBAST+AI: An updated quality, risk of bias, and applicability assessment tool for prediction models using regression or artificial intelligence methods. *BMJ*, 388, e082505. <https://doi.org/10.1136/bmj-2024-082505>

Riley, R. D., Snell, K. I. E., Ensor, J., et al. (2019). Minimum sample size for developing a multivariable prediction model: Part II—binary and time-to-event outcomes. *Statistics in Medicine*, 38, 1276–1296. <https://doi.org/10.1002/sim.7992>

## Project status

The project is complete and includes:

- verified public data provenance and reproducible acquisition;
- historical reproduction of the 2018 analysis;
- documented corrections and model-selection ambiguity;
- a prespecified modern analysis plan;
- formal sample-size and model-complexity assessment;
- modern input validation;
- repeated cross-validation of the primary model comparison;
- sensitivity and exploratory analyses;
- machine-readable results and programmatically generated figures;
- a complete integrated Quarto report;
- an approved self-contained HTML report;
- a recorded and restorable R package environment; and
- a final repository audit for public release.

The reported model performance represents internal validation of a historical
teaching dataset. The models are not intended for clinical use.
