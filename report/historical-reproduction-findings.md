# Historical Reproduction Findings

## Purpose

This historical reproduction evaluated whether the principal numerical results and conclusions from the 2018 SDSU project could be reproduced using the verified historical dataset and a documented R workflow.

The reproduction was designed to assess the fidelity of the earlier analysis, not whether its statistical methods represented current best practice. Methodological improvements will be addressed separately in the modern extension.

## Historical analysis population

The verified historical dataset contains 380 observations and eight variables. Prostatic capsule penetration is recorded for 153 patients, while 227 patients did not experience capsule penetration.

Race is missing for records 22, 46, and 252. Following the 2018 workflow, complete-case deletion was applied across the dataset, producing an analysis population of 377 observations. This deletion was preserved for historical comparability, although it is unnecessary for models that do not include race.

The outcome represents penetration of the prostatic capsule among patients already diagnosed with prostate cancer. It does not represent initial prostate cancer detection or screening.

## Descriptive findings

The descriptive statistics from the historical dataset were reproduced successfully.

| Variable | Mean | Standard deviation | Minimum | 25th percentile | 75th percentile | Maximum |
|---|---:|---:|---:|---:|---:|---:|
| Age | 66.04 | 6.53 | 43.0 | 62.0 | 71.0 | 79.0 |
| PSA | 15.41 | 20.00 | 0.3 | 5.0 | 17.13 | 139.7 |
| Gleason score | 6.38 | 1.09 | 0.0 | 6.0 | 7.0 | 9.0 |

These values agree with the descriptive results presented in the 2018 report at the reported precision.

## Univariate findings

The historical two-sample comparisons showed strong evidence that PSA and Gleason score differed according to capsule-penetration status.

Mean PSA was approximately 10.00 among patients without capsule penetration and 23.12 among patients with capsule penetration. The estimated mean difference was −13.12 when calculated as the mean for patients without penetration minus the mean for patients with penetration (`p = 4.50 × 10⁻⁸`).

Mean Gleason score was approximately 5.98 among patients without capsule penetration and 6.98 among patients with capsule penetration. The corresponding mean difference was −1.00 (`p = 2.30 × 10⁻²⁰`).

Age did not show evidence of a difference between the two outcome groups. Mean age was approximately 66.28 among patients without capsule penetration and 65.74 among those with capsule penetration (`p = 0.426`).

The categorical comparisons produced the following results:

| Predictor | Chi-square statistic | Degrees of freedom | P-value | Interpretation |
|---|---:|---:|---:|---|
| Race | Approximately 0.00 | 1 | Approximately 1.000 | No evidence of association |
| Digital rectal examination (`dpros`) | 38.74 | 3 | `1.97 × 10⁻⁸` | Evidence of association |
| Baseline detection of capsular involvement (`dcaps`) | 21.16 | 1 | `4.22 × 10⁻⁶` | Evidence of association |

These univariate findings support the 2018 conclusions that PSA, Gleason score, digital rectal examination findings, and baseline detection of capsular involvement were associated with capsule penetration, while age and race were not.

These comparisons describe unadjusted associations. They should not be interpreted as measures of predictive performance or as evidence of causal effects.

## Reproduction of the model-selection process

The historical candidate logistic regression model included `dpros`, `dcaps`, PSA, Gleason score, and the documented two-way interactions. Its AIC was 398.73.

Applying bidirectional stepwise AIC selection to this candidate model produced:

```r
capsule ~ dpros + psa + gleason + dpros:psa
```

The selected model had an AIC of 391.94.

The final model presented in the 2018 report was instead:

```r
capsule ~ factor(dpros) + psa + gleason
```

This model had an AIC of 393.22.

The exact stepwise-selected model therefore does not match the model ultimately reported in 2018. The historical code appears to have manually fitted the reported model after the stepwise procedure, but the available materials do not document why the interaction was removed or why `dpros` was changed from a numeric term to a categorical factor.

This workflow discrepancy is classified as **Ambiguous**. The reported final model itself remains reproducible, but it cannot be described as the direct output of the documented stepwise procedure.

## Reproduction of the reported final model

The reported final logistic regression model was fitted to the 377 complete observations. Its AIC of 393.22 was reproduced.

All six reported coefficient estimates, coefficient confidence intervals, odds ratios, and odds-ratio confidence intervals were reproduced at the precision presented in the 2018 report.

| Predictor | Coefficient | P-value | Odds ratio | 95% odds-ratio CI |
|---|---:|---:|---:|---:|
| Intercept | −8.145 | `< 0.001` | `< 0.001` | `< 0.001` to `0.002` |
| `dpros = 2` versus `dpros = 1` | 0.773 | 0.030 | 2.17 | 1.09 to 4.43 |
| `dpros = 3` versus `dpros = 1` | 1.553 | `< 0.001` | 4.73 | 2.32 to 10.00 |
| `dpros = 4` versus `dpros = 1` | 1.429 | 0.001 | 4.18 | 1.75 to 10.25 |
| PSA | 0.027 | 0.004 | 1.03 | 1.01 to 1.05 |
| Gleason score | 0.995 | `< 0.001` | 2.71 | 2.00 to 3.76 |

Holding the other variables in the model constant, higher digital rectal examination categories were associated with greater odds of capsule penetration relative to category 1.

Each one-unit increase in PSA was associated with approximately 2.8% higher odds of capsule penetration. Because PSA has a wide and strongly skewed distribution, this per-unit interpretation should be treated cautiously and reconsidered in the modern analysis.

Each one-point increase in Gleason score was associated with approximately 2.71 times the odds of capsule penetration, conditional on PSA and `dpros`.

These estimates describe associations within the historical dataset. They do not establish causality, and the historical analysis did not measure out-of-sample predictive performance.

## Corrected p-value

The 2018 regression table reports a p-value of `0.31` for `dpros = 2`. The reproduced model gives:

```text
p = 0.02997
```

The corresponding coefficient, standard error, confidence interval, and odds ratio all reproduce successfully. In addition, the reproduced confidence interval excludes zero on the coefficient scale and excludes one on the odds-ratio scale.

The reported value of `0.31` is therefore almost certainly a typographical error. This result is classified as **Corrected**. The historical report remains unchanged, while the reproducible output records the model-derived p-value.

## Diagnostic findings

The historical final model was evaluated using standardized Pearson residuals, Cook’s distance, and leverage.

The diagnostic screens identified:

- 13 observations with an absolute standardized Pearson residual greater than 2;
- 25 observations with Cook’s distance greater than `4/n`;
- 24 observations with leverage greater than `2p/n`; and
- 43 observations meeting at least one of these screening criteria.

These are diagnostically flagged observations, not observations that should automatically be excluded. A diagnostic flag identifies a record for further examination and does not, by itself, demonstrate an error or justify deletion.

As a deliberately conservative stress test, the model was refitted after simultaneously removing the union of all 43 flagged observations. This reduced the analysis population from 377 to 334 and produced substantial coefficient changes. The largest absolute change was approximately 4.28 for the intercept; meaningful changes also occurred for the predictor coefficients.

Because this stress test removes every observation meeting any of three broad screening rules, it should not be interpreted as a preferred alternative analysis. Nevertheless, it does not support an unqualified claim that removing all diagnostically flagged observations would leave the fitted model materially unchanged.

The primary historical reproduction therefore retains all 377 complete observations.

## Classification of the historical results

| Component | Classification | Finding |
|---|---|---|
| Dataset counts and missing-value locations | **Reproduced** | Exact agreement |
| Descriptive statistics | **Reproduced** | Agreement at the reported precision |
| Principal univariate findings | **Reproduced** | Same substantive conclusions |
| Final-model coefficients | **Reproduced** | All six estimates reproduced |
| Coefficient confidence intervals | **Reproduced** | All six intervals reproduced |
| Odds ratios and confidence intervals | **Reproduced** | All six reproduced |
| Final-model AIC | **Reproduced** | AIC = 393.22 |
| P-value for `dpros = 2` | **Corrected** | Reported as 0.31; reproduced as 0.02997 |
| Transition from stepwise selection to the reported model | **Ambiguous** | The available workflow does not explain the change |
| Claim concerning removal of flagged observations | **Not supported without qualification** | Conservative simultaneous removal materially changed the estimates |

## Overall conclusion

The central numerical results of the 2018 SDSU analysis are highly reproducible. The dataset description, principal univariate findings, final logistic regression coefficients, confidence intervals, odds ratios, and AIC were all recovered from the verified historical data.

Two qualifications emerged. First, one p-value in the historical report is almost certainly a typographical error. Second, the model selected by the documented stepwise procedure differs from the model ultimately reported, and the reason for this change is not documented.

The reproduction also shows why the historical model should not yet be treated as a validated prediction model. The 2018 analysis estimated associations using the same observations employed to fit the model, relied on stepwise selection, did not report out-of-sample discrimination or calibration, and offered limited evidence about model stability.

The historical conclusions are therefore reproducible as in-sample statistical associations. The next stage of the project will evaluate whether these predictors provide stable and well-calibrated out-of-sample predictions using a prespecified modern modeling workflow.

## Reproducibility record

The historical reproduction is implemented in:

```text
R/02_reproduce_2018.R
```

Machine-readable results are stored in:

```text
results/historical/
```

Reproduced figures are stored in:

```text
figures/historical/
```

The original 2018 report remains unchanged as a historical artifact.
