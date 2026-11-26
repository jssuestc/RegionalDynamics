# Regional dynamics: group-map analyses

Run `Main_GroupMaps.m` in MATLAB R2021a or later with Statistics and Machine Learning Toolbox. Edit all input paths and output directory names in that main script. The three previous separate main scripts have been consolidated. To run a single analysis, execute the configuration/setup sections first, then its analysis section.

## Inputs

All cortical maps follow numeric Brainnetome ROI order 1:210.

| Input | MAT fields |
|---|---|
| Pooled RAD/iRTD mean maps | `dataROI`, 210 values (not an average of age-bin means) |
| Neuromaps | `roi_matrix` (210 x 5 or 5 x 210), `file_names`: T1T2ratio, G1_fMRI, G1_dMRI, PET_AG, Cortical_Thickness |
| Neurosynth | `roi_matrix` (124 x 210 or 210 x 124), `file_names`, optional `roi_labels` |
| Individual PLS weights | `PLS1w_Z_individual`, N x 124, same term order as Neurosynth |
| Optional age | `Age` or `AgeN`, N values in exactly the same subject row order as each weight file |
| Four symptom-map files | `Y_all`, 210 x 3; `symNames`, three labels |

The input MAT files and BN annotation/cortex labels are included in this module. The main script uses the bundled `Neuromaps/Neuromaps_BN_210.mat` and `NeursynthMaps/NeuroSynth124_BN_210.mat`. Registered fsaverage spheres (`lh.sphere.reg`, `rh.sphere.reg`, 163842 vertices each) are not included and are required only if a valid 210 x 5000 spin cache is unavailable. By default, `sphereDir` points to `fsaverage/surf/` at the repository root. Set it to your local sphere directory, or set `spinFile` to reuse an existing validated cache.

Each bundled `Individual_pls_RAD/PLS1w_Z.mat` and `Individual_pls_iRTD/PLS1w_Z.mat` contains its own `Age` vector, so `radAgeFile` and `irtAgeFile` default to their respective weight files. Set either to `''` to skip its age fits. A missing configured age file warns and skips age fits. No subject rows are automatically deleted. Verify participant and term alignment yourself.

## Analyses and outputs

- **Dominance:** all-age pooled maps, five predictors, average incremental ordinary R² over all subset sizes. MAT, CSV and PNG outputs; weights sum to full-model R².
- **Individual clustering:** cluster 124 term profiles (not participants), row-wise normalization across participants, squared-Euclidean k-means++ for K=2:8, 100 restarts, max 1000 iterations, seed 2026. Retain the best mean silhouette candidate, as in the uploaded tested archive. Quadratic age fits use the normalized term profiles and centroids. `Cluster_results.mat`, `Cluster_Data_ForPython.mat`, member/silhouette CSVs, word clouds and optional age trajectories are retained.
- **Symptom PLS:** MDD_RAD, MDD_iRTD, SZ_RAD and SZ_iRTD, three symptom columns jointly, standardized across ROIs. 20-component variance curve, PLS1/2 inference, 5000 paired ROI bootstrap resamples and 5000 bilateral spins of all Y columns together, X fixed. Refit PLS for every spin. Upper-tail plus-one variance p-values and two-sided absolute Pearson score-correlation plus-one p-values. The PLS2 spatial null is not conditional on PLS1. Bootstrap Z ratios are not spatially corrected. Retain `PLS_results.mat`, PLS1/2 weights/loadings, `PLS_summary.csv`, `Variance_curve.csv`, `Spin_null_variance.csv`, `Spin_null_scores.csv` and the aggregate summary. Legacy fields are unchanged.

All input paths default to the module directory. The default output root is `results/` within this module. All output directories (including separate RAD/iRTD clustering and the four symptom models) are editable in `Main_GroupMaps.m`. Existing filenames are overwritten; use a new output root for a new run. Only routine prints and duplicate checks were removed; important map/ROI/term-order and resampling checks remain. No statistical algorithms, seeds or resampling counts were changed. Original FreeSurfer reader attributions are retained. Verify redistribution permissions before public release.

The updated paths and input dimensions were checked statically; the revised workflow has not been rerun end to end in MATLAB.
