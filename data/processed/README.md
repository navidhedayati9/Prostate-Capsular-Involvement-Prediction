# Processed data

This directory is reserved for datasets generated reproducibly from the canonical source data.

Generated data files are intentionally excluded from Git version control. Only this README should be committed from the `data/processed/` directory.

## Purpose

Processed datasets provide stable inputs for the historical reproduction and modern analysis without altering the canonical source data.

Every file in this directory must be created by a tracked script. Manual edits are not permitted.

## Planned datasets

The data-acquisition script creates:

- `prostate_modern.rds`: the canonical nine-variable dataset with `tumor` renamed to `capsule` and `vol` retained.
- `prostate_2018.rds`: the exact eight-variable structure used for the 2018 SDSU analysis, with `tumor` renamed to `capsule` and `vol` excluded.
- `acquisition_metadata.txt`: a human-readable record of the source package, package version, acquisition date, validation results, transformations, and output checksums.

Both analysis datasets retain all 380 observations. Missing values are preserved at this stage.

## Processing rules

The acquisition process must:

1. obtain `PCS` from the documented public `lbreg` package source,
2. verify the expected dimensions and column names,
3. verify unique sequential record identifiers,
4. verify outcome counts and missing-value locations,
5. preserve the original observation order,
6. rename the outcome without changing its values,
7. derive the historical dataset by excluding only `vol`, and
8. stop with an error if any validation check fails.

No observations are removed during acquisition.

Categorical encoding, missing-data handling, transformations, resampling, and model-specific preprocessing will be performed in later analysis stages.

## Regeneration

From the project root, run:

```r
source("R/01_acquire_data.R")
```

The generated files may be deleted and recreated at any time by rerunning the tracked script.

## Version-control policy

Do not commit generated `.rds`, `.csv`, `.txt`, or other data artifacts from this directory.

The top-level `.gitignore` excludes all contents of `data/processed/` except this README.

## Clinical-use limitation

These datasets are derived from historical data for statistical education and methodological demonstration. They must not be used for clinical decision-making.
