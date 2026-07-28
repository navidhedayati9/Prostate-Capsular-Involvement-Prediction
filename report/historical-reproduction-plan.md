# Historical reproduction plan

## Purpose

The historical reproduction will answer:

> Can the principal numerical results and conclusions reported in the 2018 SDSU project be reproduced from the verified historical dataset and documented R workflow?

This is different from asking whether the original approach was optimal.

## Analysis population

The reproduction will use `prostate_2018.rds`:

- 380 observations
- 8 variables
- `vol` excluded
- `tumor` renamed to `capsule`

To reproduce the original workflow, complete-case deletion will be applied across the dataset:

- Records 22, 46, and 252 will be removed because `race` is missing.
- The historical modeling population will contain 377 observations.

This deletion is unnecessary for models that do not contain `race`, but it will be preserved for historical comparability.

## Results to reproduce

### 1. Dataset description

Verify:

- 380 original observations
- 153 capsule-penetration outcomes
- three missing race values
- 377 complete cases

Reproduce the descriptive statistics for:

- age
- PSA
- Gleason score

The calculated values will be compared with Table 1 of the original report.

### 2. Exploratory figures

Reproduce:

- PSA by capsule-penetration status
- age by capsule-penetration status
- Gleason score by capsule-penetration status
- the original correlation plot

The figures should preserve the original analytical intent but do not need to be pixel-for-pixel copies. Scientifically questionable choices, such as correlations involving category codes, will be identified separately.

### 3. Univariate analyses

Reproduce the original comparisons:

- PSA: two-sample t-test
- age: two-sample t-test
- Gleason score: two-sample t-test
- race: chi-square test
- digital rectal examination (`dpros`): chi-square test
- baseline detection of capsular involvement (`dcaps`): chi-square test

Full p-values will be reported rather than replacing small values with `0.00`.

### 4. Historical model-selection process

Reproduce the original candidate logistic model containing:

- `dpros`
- `dcaps`
- PSA
- Gleason score
- documented two-way interactions

Then reproduce the bidirectional stepwise AIC procedure.

A potential inconsistency must be investigated: the report says interactions were considered, but the final model appears to have been manually refitted after the stepwise procedure.

The reproduction will record:

- the exact model selected by `stepAIC`,
- the manually specified final model in the original script, and
- whether those models are identical.

### 5. Reported final logistic model

Fit the historical final model:

```r
capsule ~ factor(dpros) + psa + gleason
```

Reproduce:

- regression coefficients
- standard errors
- test statistics
- p-values
- 95% confidence intervals
- odds ratios
- odds-ratio confidence intervals
- number of observations
- AIC

This model is the central result of the 2018 report.

### 6. Model diagnostics

Reproduce the main historical diagnostics:

- standardized residuals versus fitted probabilities
- Cook's distance and leverage
- identification of potentially influential observations
- sensitivity of coefficients to influential observations, where the original report claims their removal would not materially change the result

Diagnostic generation will be non-interactive while preserving the historical method.

### 7. Conclusions

Compare the reproduced results with the original conclusions about:

- digital rectal examination findings
- PSA
- Gleason score
- age
- race
- `dcaps`

The reproduction will distinguish association with capsule penetration from prostate cancer detection.

## Comparison standards

| Result | Reproduction criterion |
|---|---|
| Observation and category counts | Exact match |
| Missing-value locations | Exact match |
| Descriptive statistics | Match before report rounding |
| Regression coefficients | Numerical agreement within `1e-6` |
| Odds ratios and confidence intervals | Match when rounded to two decimals |
| AIC | Match when rounded to two decimals |
| P-values | Compare unrounded calculated values |
| Figures | Same data and analytical meaning; not pixel-identical |
| Written conclusions | Classify as supported, overstated, or contradicted |

## Discrepancy classifications

Each result will receive one status:

- **Reproduced:** agrees with the original report.
- **Rounding difference:** underlying result agrees, presentation differs.
- **Corrected:** original value or label is demonstrably wrong.
- **Not reproducible:** available materials do not generate the reported result.
- **Ambiguous:** original workflow does not provide enough information.

The historical report itself will remain unchanged.

## Proposed implementation

The reproduction script will be:

```text
R/02_reproduce_2018.R
```

It will:

- read only the processed historical dataset,
- reproduce calculations without manual p-values,
- avoid `attach()`, `file.choose()`, and workspace deletion,
- avoid interactive diagnostics,
- save structured tables and figures to designated output folders, and
- stop if the expected historical dataset is unavailable or invalid.
