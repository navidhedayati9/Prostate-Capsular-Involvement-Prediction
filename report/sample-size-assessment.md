# Formal Sample-Size and Model-Complexity Assessment

## Purpose

This assessment evaluates whether the available dataset is reasonably adequate for the models prespecified in the modern analysis plan.

It was completed before fitting or comparing the modern models. Its purpose is not to prove that the sample size is sufficient for clinical model development, but to identify which planned analyses are reasonably supported and where overfitting or instability remains a concern.

## Dataset information

The modern dataset contains:

- 380 observations;
- 153 capsule-penetration outcomes;
- 227 observations without capsule penetration; and
- an observed outcome prevalence of 40.3%.

Because the primary predictors contain no missing values, all 380 observations are available for the modern analysis.

## Predictor-parameter counts

Model complexity was measured by the number of estimated predictor parameters, excluding the intercept.

| Model | Predictor parameters | Explanation |
|---|---:|---|
| Historical benchmark | 5 | Three `dpros` indicators, PSA, and Gleason score |
| Primary modern model | 6 | Age, three `dpros` indicators, `log2(PSA)`, and Gleason score |
| Nonlinear PSA model | 8 | Primary terms, with three spline parameters replacing one linear PSA parameter |
| `dcaps` sensitivity model | 7 | Primary model plus one `dcaps` parameter |
| Tumor-volume exploratory model | 8 | Primary model plus two tumor-volume parameters |

The intercept is estimated in every model but is not included in the predictor-parameter count used by the sample-size framework.

## Events per predictor parameter

For descriptive context, the dataset provides:

| Model | Predictor parameters | Events per predictor parameter |
|---|---:|---:|
| Historical benchmark | 5 | 30.6 |
| Primary modern model | 6 | 25.5 |
| `dcaps` sensitivity model | 7 | 21.9 |
| Nonlinear PSA model | 8 | 19.1 |
| Tumor-volume exploratory model | 8 | 19.1 |

These ratios are not treated as definitive adequacy rules. Modern sample-size guidance recommends considering anticipated model performance, overfitting, and precision rather than relying only on a fixed events-per-parameter threshold.

## Anticipated model performance

The reproduced historical model has an apparent Cox–Snell \(R^2\) of approximately 0.285 in the 377-observation historical analysis population.

Because apparent performance is measured in the same data used to fit the model, using 0.285 as the sole planning assumption would be optimistic. Three scenarios were therefore evaluated:

| Scenario | Anticipated Cox–Snell \(R^2\) | Rationale |
|---|---:|---|
| Conservative | 0.142 | 50% of the historical apparent value |
| Intermediate | 0.214 | 75% of the historical apparent value |
| Optimistic boundary | 0.285 | 100% of the historical apparent value |

The third scenario is called an optimistic boundary because it assumes that the apparent historical explanatory strength is fully retained. It is not the preferred planning assumption.

## Assessment method

Minimum sample sizes were calculated using the binary-outcome framework of Riley and colleagues, as implemented in `pmsampsize` version 1.1.3.

The calculations used:

- outcome prevalence of 0.403;
- target global shrinkage of at least 0.90;
- maximum acceptable difference of 0.05 between apparent and adjusted \(R^2\); and
- maximum margin of error of 0.05 for estimation of the overall outcome risk.

The final required sample size is the largest requirement produced by the framework’s three criteria.

## Required sample sizes

| Model | Parameters | Conservative \(R^2=0.142\) | Intermediate \(R^2=0.214\) | Optimistic boundary \(R^2=0.285\) |
|---|---:|---:|---:|---:|
| Historical benchmark | 5 | 370 | 370 | 370 |
| Primary modern model | 6 | 370 | 370 | 370 |
| `dcaps` sensitivity model | 7 | 407 | 370 | 370 |
| Nonlinear PSA model | 8 | 465 | 370 | 370 |
| Tumor-volume exploratory model | 8 | 465 | 370 | 370 |

The available sample size is 380.

## Interpretation by model

### Historical benchmark

The five-parameter historical benchmark is supported under all three scenarios. It remains subject to internal-validation uncertainty and should not be interpreted as externally validated.

### Primary modern model

The six-parameter primary model is also supported under all three scenarios. Under the conservative scenario, the required sample size is 370, only 10 observations fewer than the available 380.

The primary model is therefore reasonably supported, but not generously supported. Its calibration, coefficients, and individual predictions may still be unstable. This reinforces the decision to avoid interactions, automated variable selection, and additional unplanned predictors.

### Nonlinear PSA model

The eight-parameter nonlinear PSA model is supported under the intermediate and optimistic scenarios but not under the conservative scenario, which requires approximately 465 observations.

This model should remain a limited secondary comparison. Any apparent improvement over the primary model must be consistent across resamples and large enough to justify the two additional PSA parameters. An unstable or very small improvement should be interpreted as inconclusive.

### `dcaps` sensitivity model

The seven-parameter `dcaps` model requires approximately 407 observations under the conservative scenario, exceeding the available sample by 27 observations. It is supported under the intermediate and optimistic scenarios.

This result does not require abandoning the prespecified sensitivity analysis. It does require labeling its estimates as potentially more susceptible to overfitting under weak-signal conditions.

### Tumor-volume exploratory model

The eight-parameter tumor-volume model requires approximately 465 observations under the conservative scenario. As with the nonlinear PSA model, it is supported only under the intermediate and optimistic scenarios.

The volume model should remain exploratory. Its added complexity is especially important because 167 observations have recorded tumor volume equal to zero, requiring a two-part representation.

## Implications for internal validation

A one-time training/test split would reduce the effective model-development sample substantially and would produce an imprecise test-set evaluation. The approved repeated cross-validation strategy remains appropriate for using the available observations efficiently.

Nevertheless, repeated cross-validation does not create additional independent patients. It estimates internal performance while exposing variability caused by different data partitions, but it cannot eliminate uncertainty caused by the limited original sample.

Particular attention should be paid to:

- variability in calibration slope and calibration intercept;
- variability in individual predicted probabilities;
- instability of spline and exploratory-model coefficients;
- extreme predictions or convergence problems; and
- performance differences that are small relative to resampling variability.

## Decision

The formal assessment supports proceeding with the prespecified analysis subject to the following interpretation rules:

1. The six-parameter primary modern model is adequately supported for an educational internal-validation study, although its sample-size margin is narrow under conservative assumptions.
2. The historical benchmark can be evaluated as planned.
3. The nonlinear PSA model will remain a secondary comparison and will not automatically replace the primary model.
4. The `dcaps` analysis will remain a separately labeled sensitivity analysis.
5. The tumor-volume analysis will remain exploratory.
6. No additional interactions, spline terms, predictors, or data-driven selection procedures will be introduced without documenting a protocol deviation.
7. Results will not be described as externally validated or ready for clinical use.

## Overall conclusion

The dataset is sufficient to proceed with the prespecified primary logistic regression model and its internal validation. It is not large enough to support unrestricted modeling flexibility.

The sample-size assessment therefore reinforces the project’s intentionally limited model set: a historical benchmark, one primary modern model, one modest nonlinear comparison, and clearly separated sensitivity and exploratory extensions.

Greater uncertainty is expected for the nonlinear PSA, `dcaps`, and tumor-volume models under conservative assumptions. Their results should be interpreted as supporting evidence rather than as independent clinical model-development claims.

## Reference

Riley RD, Snell KIE, Ensor J, et al. Minimum sample size for developing a multivariable prediction model: Part II—binary and time-to-event outcomes. *Statistics in Medicine*. 2019;38:1276–1296. <https://doi.org/10.1002/sim.7992>
