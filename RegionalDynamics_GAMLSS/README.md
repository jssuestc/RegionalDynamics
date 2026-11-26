# Global RAD: GAMLSS (SN2)

Install once: `install.packages(c("gamlss", "gamlss.dist", "R.matlab"))`. Open `Main_GAMLSS.R` in RStudio and click Source, or run `Rscript Main_GAMLSS.R`. Edit all input paths, MAT variable names and output folder names in the main script only. The four modules are invoked automatically.

## Input

The default file is `IndividualGlobalRAD/GlobalRAD.mat`, relative to `Main_GAMLSS.R`. It contains the already harmonized whole-cortex RAD measure. `VARS` explicitly maps `BayesDataN` (RAD), `AgeN`, `GenderN`, `SitesN`. Each must contain one aligned value per participant. Use MAT v7 rather than v7.3/HDF5. No ComBat or global averaging is performed. The bundled MAT file contains participant-level RAD values and covariates; no patient files are included.

`INPUT_WITHOUT_GROUP <- "HC"` assumes the file is HC-only when `VARS$group` is NULL. For a combined cohort file, explicitly set `VARS$group` to its diagnosis field and confirm `GROUP_MAP`. Set `VARS$id` if IDs exist; otherwise row-based IDs are used. There is no automatic variable-name guessing. Optional `PATIENT_FILES` contains named SZ/MDD paths; use `PATIENT_VARS` to specify different field names for individual patient files if needed. Fixed-site predictions require every patient site to be represented in HC. No missing data are imputed or observations silently removed.

## Statistical definitions

SN2 is prespecified from the study-wide mean BIC comparison across global and seven network measures. No distribution selection is performed in this global-only example.

| Module | Method |
|---|---|
| 01 Normative modeling | `mu ~ pb(age) + gender + site`, `sigma ~ pb(age)`, `nu ~ 1`; categorical fixed effects. Standardized conditional centiles are weighted averages over empirical joint HC gender/site cells, not mixture quantiles. |
| 02 Normative deviation | HC site-stratified 10-fold CV; patients use the full HC fit. `qnorm(clip(F(y),1e-10,1-1e-10))`, extreme if abs(z) >= 1.96. Calibration KS p-values are descriptive. |
| 03 Age effect | Remove age from mu only, retain sigma age and nu intercept. Partial R² from SSE of fitted mu, signed using mean mu derivative. Effective-df chi-square p-value is approximate. |
| 04 Age derivative | 1000 attempted subject resamples within site x gender; refit same model, original age grid and covariate weights; finite-difference P50 derivative and pointwise percentile 95% CIs. Withhold CIs if fewer than 80% of attempts succeed. No separate P50 CI plot. |

## Output

The default `OUTPUT_ROOT` is `results/GlobalRAD/` within this module. This directory is ignored by Git; change `OUTPUT_ROOT` if results should be stored outside the repository. Results go directly to the four configured folders (no automatic timestamp): centile CSV/MAT/PNG and model summary; deviation CSV/MAT and HC calibration plots; age-effect CSV/MAT and mu curve; slope CSV/MAT/PNG, bootstrap log and draws. Existing filenames are overwritten; choose a new root for each new analysis. `SAVE_MODEL <- FALSE` by default. Treat individual results, scatter plots, fitted models and bootstrap draws as private unless authorized for release. Routine console prints and administrative run-level files were removed; convergence, invalid prediction and bootstrap-success checks remain. All numerical settings and seeds are unchanged. For a smoke test only, reduce `N_BOOTSTRAP` to 20 and restore 1000 for inference.

The updated paths and input dimensions were checked statically; the revised workflow has not been rerun end to end in R/GAMLSS.
