# Modern Analysis Findings

## Purpose

The modern analysis evaluated how accurately baseline clinical measurements predicted prostatic capsule penetration among patients diagnosed with prostate cancer.

Unlike the 2018 analysis, which focused primarily on in-sample associations, the modern analysis emphasized out-of-sample predicted probabilities, discrimination, calibration, overall prediction error, and model stability.

The analysis followed the prespecified modern analysis plan and formal sample-size assessment. No stepwise selection, univariate predictor screening, class balancing, or post hoc search across additional models was performed.

## Analysis population

The modern dataset contained 380 observations:

- 153 patients with capsule penetration;
- 227 patients without capsule penetration; and
- an observed capsule-penetration prevalence of 40.3%.

The primary predictors—age, digital rectal examination result, PSA, and Gleason score—had no missing values. All 380 observations were therefore included in the primary modern analysis.

Race was excluded from modern predictive modeling for the reasons specified in the analysis plan. The record identifier was never used as a predictor.

## Internal-validation procedure

Predictive performance was evaluated using outcome-stratified 10-fold cross-validation repeated 20 times.

The same saved folds were used for every model comparison. Within each repeat, every observation received one prediction from a model that had been fitted without that observation.

The analysis generated:

- 7,600 saved resampling assignments;
- 22,800 primary-analysis out-of-fold predictions;
- 600 fitted primary-analysis cross-validation models; and
- 300 repeat-level primary performance estimates.

All models converged, all predicted probabilities were finite and strictly between zero and one, and every observation received the expected predictions.

The percentile ranges across repeats describe resampling variability. They are not formal confidence intervals because estimates from repeated cross-validation are correlated.

## Models compared

### Historical benchmark

```r
capsule ~ factor(dpros) + psa + gleason
```

This reproduced the structure of the final model reported in 2018.

### Primary modern model

```r
capsule ~ age + factor(dpros) + log2(psa) + gleason
```

This model added age and represented PSA on a base-2 logarithmic scale. The PSA coefficient therefore describes the association with a doubling of PSA.

### Nonlinear PSA model

```r
capsule ~ age +
  factor(dpros) +
  splines::ns(log2(psa), df = 3) +
  gleason
```

This model allowed a modestly curved PSA relationship using two additional parameters relative to the primary model.

## Cross-validated performance

Mean performance across the 20 cross-validation repeats was:

| Model | ROC AUC | Brier score | Log loss | Calibration intercept | Calibration slope |
|---|---:|---:|---:|---:|---:|
| Historical benchmark | 0.809 | 0.17315 | 0.51878 | −0.002 | 0.920 |
| Primary modern | 0.810 | 0.17310 | 0.51925 | −0.003 | 0.914 |
| Nonlinear PSA | 0.805 | 0.17461 | 0.52266 | −0.003 | 0.884 |

Higher AUC is preferable, while lower Brier score and log loss are preferable. A calibration intercept of zero and calibration slope of one are ideal.

## Historical benchmark versus primary modern model

The primary modern model and historical benchmark produced almost identical out-of-sample performance.

On average, the primary model differed from the historical benchmark by:

- `+0.00045` in AUC;
- `−0.00005` in Brier score;
- `+0.00047` in log loss; and
- `−0.0063` in calibration slope.

The primary model therefore had a negligibly higher AUC and slightly lower Brier score, while the historical benchmark had slightly better log loss and a calibration slope closer to one.

These differences are too small to support a claim that the primary modern model materially outperformed the historical benchmark. Adding age and replacing raw PSA with `log2(PSA)` changed the model’s interpretation but did not produce a meaningful improvement in out-of-sample prediction.

The primary modern model remains the designated primary analysis because it was prespecified before performance was examined. The historical benchmark remains a similarly performing, somewhat simpler alternative.

## Nonlinear PSA model

The nonlinear PSA model did not improve prediction.

Compared with the primary modern model, it had:

- AUC lower by approximately `0.0045`;
- Brier score higher by approximately `0.0015`;
- log loss higher by approximately `0.0034`; and
- calibration slope lower by approximately `0.030`.

Its mean calibration slope was 0.884, compared with 0.914 for the primary model. Its individual out-of-fold predictions also varied more across resampling repeats.

The mean standard deviation of an observation’s predicted probability across repeats was:

| Model | Mean prediction SD |
|---|---:|
| Historical benchmark | 0.0160 |
| Primary modern | 0.0173 |
| Nonlinear PSA | 0.0197 |

The additional flexibility therefore increased instability without improving out-of-sample performance. The results do not support replacing the linear `log2(PSA)` term with the three-degree-of-freedom spline in this dataset.

## Calibration

All three models had mean calibration intercepts close to zero. This indicates little average tendency to systematically overpredict or underpredict the overall event rate.

The calibration slopes were below one:

- 0.920 for the historical benchmark;
- 0.914 for the primary modern model; and
- 0.884 for the nonlinear PSA model.

Slopes below one indicate that predictions were somewhat too extreme when applied to held-out observations, consistent with modest overfitting. The lower slope for the nonlinear model indicates greater overfitting from its additional flexibility.

The calibration plots showed generally reasonable agreement between predicted and observed probabilities, with local deviations in parts of the risk range. These deviations should be interpreted cautiously because some regions contained relatively few observations.

## Full-data primary-model coefficients

After cross-validation was completed, the primary modern model was fitted to all 380 observations.

| Predictor | Coefficient | Odds ratio | 95% Wald odds-ratio interval |
|---|---:|---:|---:|
| Age, per year | −0.012 | 0.99 | 0.95 to 1.03 |
| `dpros = 2` versus 1 | 0.768 | 2.15 | 1.07 to 4.33 |
| `dpros = 3` versus 1 | 1.527 | 4.60 | 2.22 to 9.56 |
| `dpros = 4` versus 1 | 1.484 | 4.41 | 1.85 to 10.52 |
| PSA, per doubling | 0.349 | 1.42 | 1.15 to 1.75 |
| Gleason score, per point | 0.962 | 2.62 | 1.89 to 3.62 |

Holding the other predictors constant:

- higher digital rectal examination categories were associated with greater odds of capsule penetration;
- each doubling of PSA was associated with approximately 42% higher odds of capsule penetration; and
- each one-point increase in Gleason score was associated with approximately 2.62 times the odds of capsule penetration.

Age contributed little additional information after accounting for the other predictors.

These full-data coefficients describe associations within the historical dataset. They do not establish causality, and their apparent full-data performance is not the primary estimate of predictive accuracy.

## `dcaps` sensitivity analysis

The prespecified `dcaps` sensitivity model added baseline detection of capsular involvement during the rectal examination to the primary model.

Compared with the primary model, adding `dcaps` produced average changes of approximately:

- `+0.00017` in AUC;
- `−0.00018` in Brier score;
- `+0.00059` in log loss; and
- `−0.0106` in preferred calibration-slope performance.

The changes were small and did not consistently favor the extended model. The full-data `dcaps` coefficient was also imprecisely estimated.

Adding `dcaps` therefore did not provide meaningful or stable incremental predictive value beyond age, `dpros`, PSA, and Gleason score.

This finding also reduces concern that excluding `dcaps` from the primary model omitted a major source of predictive information.

## Zero-Gleason sensitivity analysis

The dataset contains two observations with Gleason score equal to zero:

- record 282; and
- record 357.

Both observations had no capsule penetration. The primary model was refitted and revalidated after excluding these records while preserving the saved fold assignments for the remaining 378 observations.

For a matched comparison on the same 378 observations:

- AUC rankings were unchanged;
- mean Brier-score difference was less than `0.000001`;
- mean log-loss difference was less than `0.000002`;
- calibration intercept was effectively unchanged; and
- calibration slope was effectively unchanged.

The full-data coefficients were also very similar.

The two zero Gleason records therefore do not materially affect the primary conclusions. Retaining them in the primary analysis, as required by the documented public dataset, is supported by this stability result.

## Exploratory tumor-volume analysis

The exploratory volume model added:

1. an indicator distinguishing zero from positive recorded tumor volume; and
2. `log2(volume)` among observations with positive recorded volume.

Compared with the primary model, the volume model showed average changes of approximately:

- `+0.0032` in AUC;
- `−0.0011` in Brier score;
- `−0.0007` in log loss; and
- `−0.0213` in preferred calibration-slope performance.

The small improvements in discrimination and overall error were accompanied by worse calibration slope. The improvements were also not uniform across resampling repeats.

In the full-data exploratory model:

- the coefficient distinguishing zero from positive recorded volume was imprecise;
- the coefficient describing changes among positive volumes was effectively zero; and
- the fitted relationship was nearly flat across the positive-volume range.

The apparent predictive information may therefore arise from the distinction between zero and positive recorded volume rather than from increasing positive tumor volume.

Because the meaning of the 167 zero values is uncertain, this distinction cannot be given a confident biological interpretation. Zero might indicate no measurable volume, a measurement below detection, an unrecorded measurement encoded as zero, or another historical coding practice.

The tumor-volume findings are therefore classified as **inconclusive and exploratory**. They do not justify adding volume to the primary model.

## Model stability

The resampling analysis showed that the historical and primary models were relatively stable, although individual predictions still varied across data partitions.

The nonlinear and exploratory models required additional parameters and generally displayed greater instability or worse calibration. This pattern is consistent with the formal sample-size assessment, which found that the more complex models were not adequately supported under the most conservative anticipated-performance scenario.

No observations were removed based on diagnostic results, and no additional models were introduced after performance was examined.

## Overall findings

The principal modern findings are:

1. The historical benchmark achieved useful but imperfect discrimination, with a mean cross-validated AUC of approximately 0.81.
2. The primary modern model performed almost identically to the historical benchmark.
3. Transforming PSA to `log2(PSA)` provides a clearer doubling-based interpretation but did not materially improve prediction when combined with age.
4. Limited nonlinear PSA modeling increased complexity and instability without improving out-of-sample performance.
5. Calibration intercepts were close to zero, but calibration slopes below one indicated modest overfitting.
6. Adding `dcaps` did not provide meaningful incremental predictive value.
7. Excluding the two zero Gleason records did not change the conclusions.
8. Tumor volume showed small, inconsistent exploratory improvements accompanied by worse calibration and substantial uncertainty about zero-volume coding.

## Conclusion

The modern extension confirms that digital rectal examination findings, PSA, and Gleason score contain meaningful predictive information about prostatic capsule penetration in this historical dataset.

However, the analysis does not demonstrate that a more complicated model predicts better. The simple historical benchmark and prespecified primary modern model performed essentially equally, while the nonlinear PSA model performed somewhat worse.

The results favor parsimony. In this dataset, adding flexibility or additional predictors did not produce a clear, stable improvement over a straightforward logistic regression model.

The reported performance represents internal validation only. The models have not been externally validated and should not be considered contemporary clinical prediction tools.

## Limitations

Important limitations include:

- a sample of only 380 observations;
- data from a single historical clinical setting;
- values modified to protect confidentiality;
- limited and incomplete demographic information;
- no external validation dataset;
- uncertain meaning of zero tumor-volume values;
- two unusual zero Gleason scores;
- correlated estimates from repeated cross-validation;
- possible instability in individual predicted probabilities; and
- no evaluation of clinical consequences or decision thresholds.

The analysis cannot establish contemporary transportability, fairness, clinical utility, or improvement in patient outcomes.

## Reproducibility record

Input validation is implemented in:

```text
R/03_validate_inputs.R
```

The primary modern analysis is implemented in:

```text
R/04_modern_analysis.R
```

The sensitivity and exploratory analyses are implemented in:

```text
R/05_sensitivity_analyses.R
```

Machine-readable results are stored in:

```text
results/modern/
```

Figures are stored in:

```text
figures/modern/
```

The modern analysis followed the approved protocol. No material protocol deviations were introduced during model fitting or evaluation.
