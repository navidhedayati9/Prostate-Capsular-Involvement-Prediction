# Raw data

This directory is reserved for unmodified source data used by the project.

Raw data files are intentionally excluded from Git version control. Only this README should be committed from the `data/raw/` directory.

## Canonical dataset

The project uses the Prostate Cancer Study dataset described by Hosmer and Lemeshow in *Applied Logistic Regression*.

The data were collected by Dr. Donn Young at The Ohio State University Comprehensive Cancer Center. Published descriptions state that the observed values were modified to protect patient confidentiality.

The canonical dataset is publicly distributed as `PCS` in the R package `lbreg`.

### Reference

Hosmer, D. W., & Lemeshow, S. (2000). *Applied Logistic Regression* (2nd ed.). John Wiley & Sons.

## Canonical structure

The public `PCS` dataset contains 380 observations and nine variables:

- `id`: record identification code
- `tumor`: tumor penetration of the prostatic capsule
- `age`: age in years
- `race`: recorded race category
- `dpros`: digital rectal examination result
- `dcaps`: detection of capsular involvement during the rectal examination
- `psa`: prostate-specific antigen value
- `vol`: tumor volume obtained by ultrasound
- `gleason`: total Gleason score

There are three missing values in `race`, occurring in records 22, 46, and 252. The other canonical variables are complete.

## Relationship to the 2018 SDSU dataset

The dataset used for the 2018 SDSU analysis contains 380 observations and eight variables.

A row-by-row comparison established that it is an exact subset of the public `PCS` data after:

1. renaming `tumor` to `capsule`, and
2. excluding `vol`.

All remaining values, data types, observation order, and missing-value locations match the public version.

The historical reproduction will derive this eight-variable structure from the canonical data rather than treating the course copy as a separate data source.

## Data acquisition

A reproducible data-preparation script will obtain the canonical `PCS` data from the documented public R package source and create the required analysis datasets.

The acquisition and preparation process will record:

- the source package and version,
- the date of acquisition,
- the original dimensions and variable names,
- integrity checks,
- all renaming and variable-selection operations, and
- the locations of missing values.

Until that script is developed, no raw dataset should be added to this directory.

## Data-handling rules

Files placed in this directory must be treated as immutable source data.

Do not:

- edit raw values manually,
- overwrite the source dataset with cleaned data,
- store derived variables here,
- remove observations from the raw file, or
- commit raw data files to Git.

Cleaned and analysis-ready datasets, if retained, should be created programmatically and stored under `data/processed/`.

## Redistribution status

Although the dataset is publicly available through educational and software sources, an explicit data-specific redistribution license has not yet been confirmed.

For that reason, this repository will document how to obtain the canonical public data rather than directly republishing the raw dataset.

## Clinical-use limitation

This historical dataset is used for statistical education and methodological demonstration. It must not be used for clinical decision-making or represented as contemporary clinical evidence.
