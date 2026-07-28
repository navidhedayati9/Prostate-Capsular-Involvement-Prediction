# Predicting Prostatic Capsule Penetration

## A Modern Extension of a 2018 Analysis

## Project overview

This project reproduces and extends a logistic regression analysis originally completed in 2018 for STAT 696: Statistical Communication in Data Science at San Diego State University.

The original analysis examined whether baseline clinical measurements could predict whether a prostate tumor had penetrated the prostatic capsule. It identified digital rectal examination findings, prostate-specific antigen (PSA), and Gleason score as important predictors.

The present project will first reproduce the 2018 analysis and verify its reported results. It will then apply a modern, reproducible modeling workflow that evaluates out-of-sample predictive performance, calibration, and model stability.

## Central research question

> How accurately can baseline clinical measurements predict prostatic capsule penetration among patients diagnosed with prostate cancer?

## Secondary research questions

1. How does including the baseline assessment of capsular involvement (`dcaps`) affect predictive performance?
2. Does including ultrasound-derived tumor volume (`vol`) provide additional predictive information?

The `dcaps` analysis is prespecified as a sensitivity analysis because this measurement directly assesses possible capsular involvement during the baseline rectal examination.

Tumor volume is treated as a secondary exploratory analysis because it is present in the canonical dataset but was not included in the dataset used for the 2018 SDSU project.

## Objectives

1. Reproduce and verify the 2018 logistic regression analysis.
2. Develop a fully documented and reproducible R workflow.
3. Estimate predictive performance using appropriate resampling methods.
4. Evaluate both model discrimination and calibration.
5. Compare the historical model with a limited number of interpretable modern alternatives.
6. Assess the effect of including `dcaps` in a sensitivity analysis.
7. Explore the additional predictive value of `vol`.

## Historical context

The original report was completed by Navid Hedayati in 2018 as part of STAT 696 at San Diego State University.

The historical analysis used an eight-variable version of the Prostate Cancer Study dataset and selected a logistic regression model containing digital rectal examination findings (`dpros`), PSA, and Gleason score.

The original report will be retained as a historical artifact. The reproduction will distinguish between:

- results reported in 2018,
- results obtained by faithfully reproducing the original workflow,
- corrections or discrepancies identified during reproduction, and
- results from the modern extension.

The original work may be referenced as:

> Hedayati, N. (2018). *Detection of Prostate Cancer: Data Analysis Report 2*. Unpublished course report, STAT 696, San Diego State University.

## Data source

The project uses the Prostate Cancer Study dataset described by Hosmer and Lemeshow in *Applied Logistic Regression*.

The data were collected by Dr. Donn Young at The Ohio State University Comprehensive Cancer Center. The published values were modified to protect patient confidentiality.

The canonical dataset contains 380 observations and nine variables:

- `id`: record identification code
- `tumor`: tumor penetration of the prostatic capsule
- `age`: age in years
- `race`: recorded race category
- `dpros`: digital rectal examination result
- `dcaps`: detection of capsular involvement during the rectal examination
- `psa`: prostate-specific antigen value
- `vol`: tumor volume obtained by ultrasound
- `gleason`: total Gleason score

The canonical data are publicly distributed as the `PCS` dataset in the R package `lbreg`.

The eight-variable SDSU dataset is an exact subset of the canonical data after:

1. renaming `tumor` to `capsule`, and
2. excluding `vol`.

All remaining values, observation order, and missing-value locations match exactly.

### Reference

Hosmer, D. W., & Lemeshow, S. (2000). *Applied Logistic Regression* (2nd ed.). John Wiley & Sons.

## Analysis design

The project separates the analysis into three related components.

### 1. Historical reproduction

The historical reproduction will use the same eight variables available in the 2018 SDSU project. It will preserve the original complete-case approach where necessary so that the reproduced estimates can be compared fairly with the reported results.

### 2. Modern primary analysis

The primary modern analysis will:

- use `capsule` as the binary outcome,
- exclude `id` from all predictive models,
- use predictors available in the 2018 analysis,
- exclude `dcaps` from the primary model comparison,
- use model-appropriate missing-data handling,
- estimate performance using resampling, and
- evaluate both discrimination and calibration.

Decisions about transformations of PSA and representation of Gleason score will be documented before model fitting.

### 3. Sensitivity and exploratory analyses

Two additional analyses are planned:

- a prespecified sensitivity analysis that includes `dcaps`, and
- a secondary exploratory analysis that includes ultrasound tumor volume (`vol`).

These results will be reported separately from the primary model comparison.

## Scope and limitations

This study evaluates prediction of prostatic capsule penetration among patients already diagnosed with prostate cancer.

It does not evaluate:

- prostate cancer screening,
- initial prostate cancer diagnosis,
- treatment selection,
- patient prognosis,
- cancer-specific mortality, or
- clinical deployment of a prediction tool.

The dataset is historical, contains a relatively small sample from a single clinical setting, and records race using only two categories. These limitations restrict the generalizability of any results.

This repository is intended for statistical education and methodological demonstration. It is not a medical device and must not be used to make clinical decisions.

## Reproducibility

The analysis will be developed as a reproducible R project. Data preparation, historical reproduction, modern modeling, evaluation, and reporting will be separated into documented stages.

Package versions and computational requirements will be recorded as the project develops.

## Project status

This project is being developed collaboratively and incrementally.

Current decisions completed:

- project title
- central and secondary research questions
- project objectives
- report structure
- dataset provenance investigation
- comparison of the public and SDSU datasets
- initial data and analysis specification

Analysis code and results have not yet been produced.
