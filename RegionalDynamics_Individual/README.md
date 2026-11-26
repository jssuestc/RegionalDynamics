# Regional dynamics: single-participant demo

Run `Run_Individual.m` in MATLAB R2021a or later with Statistics and Machine Learning Toolbox and SPM12. Edit the **Input paths** and **Output paths** sections of that script only. No separate configuration file is needed. The main script runs all five steps; to rerun a subset, comment out the other calls after the configuration and setup. Recompute upstream files when inputs change; the simplified functions do not compare saved provenance against the current configuration.

## Inputs

- Preprocessed 4D resting-state fMRI, grey-matter mask and Brainnetome 246 atlas on the same 5-mm isotropic grid. No preprocessing or reslicing is performed.
- Full-resolution fsaverage registered `lh.sphere.reg` and `rh.sphere.reg` (163842 vertices each). BN210 annotations and cortex labels are included in `atlas/`. Use actual registered spheres, not the anatomical coordinates in cortex labels.
- Neurosynth MAT: `roi_matrix` (124 x 210, or 124 x 246 with cortical IDs in columns 1:210), `file_names` (124 terms), optional `roi_labels` in numeric BN order.

The example scan is in `individual_data/`, and the mask, volumetric atlas and parcellated Neurosynth maps are in `atlas_mask/`. The registered spheres are not bundled. By default, both MATLAB modules look for them in `fsaverage/surf/` at the repository root; alternatively, set `cfg.SphereDir` to your local sphere directory. Set `cfg.SPMDir` if SPM12 is not already on the MATLAB path. All input and output directory paths, including intermediate FCD/metric files and the spin cache, are configured in the main script. Output filenames are unchanged from the supplied tested archive.

## Analyses

| Step | Method and principal outputs |
|---|---|
| 01 FCD | Raw counts of positive Pearson r > 0.6 across all mask voxels, 50-volume windows, step 1, self-correlation included for nonconstant signals. No normalization. `FCD_windows.mat`, window summary. |
| 02 RAD/iRTD | BN cortical IDs 1:210. RAD = sample SD of window-wise parcel means (exclude zero FCD values); iRTD = 1 minus mean Pearson spatial correlation across distinct valid window pairs on fixed voxel support. `RAD_iRTD_results.mat`, ROI CSV, NIfTI maps. |
| 03 Coupling | Spearman spatial correlation; 5000 Vasa-style bilateral spins of iRTD, seed 1, two-sided plus-one p-value. Coupling MAT/CSV/plots and spin indices. |
| 04 Balance | `atan2d(tiedrank(RAD),tiedrank(iRTD))`; >45 degrees denotes relative RAD predominance. Balance MAT/CSV/maps/plot. |
| 05 PLS | Separate RAD/iRTD PLS1, standardized X/Y across ROIs, 5000 paired ROI bootstrap resamples with sign alignment, 5000 spatial-null refits per metric. `*_PLS1_results.mat`, weights CSV, score/variance null CSV and summary. |

Constant time series yield zero undefined FC counts. RAD requires a nonempty regional mean in every window; iRTD excludes undefined constant spatial patterns. At least three voxels per parcel and two valid patterns are required for iRTD. All 210 paired regional metrics must be finite for downstream analyses. BN spins are one-to-one and hemisphere-preserving (odd IDs left, even IDs right), using arithmetic mean registered-sphere coordinates. The Step03 spin matrix is reused by Step05. Explained variance uses an upper-tail plus-one spin p-value; absolute fitted X/Y score correlation uses a two-sided plus-one spin p-value. Bootstrap term Z ratios are not spatially corrected.

## Output and scope

The default root is `results/sub01/` relative to `Run_Individual.m`; change it and each step folder in `Run_Individual.m`. Existing filenames are overwritten. Use a new root for a different participant or analysis version. Routine console messages were removed; errors and meaningful QC warnings remain. MAT/CSV/NIfTI/PNG results and legacy fields are retained. The example input scan is included; generated results are not included. Check permissions before publishing data or third-party atlas/reader files; original FreeSurfer attributions are retained.

Paths are relative to the main script by default. The numerical definitions, random seeds and resampling counts are unchanged. The updated paths and input dimensions were checked statically; the revised workflow has not been rerun end to end in MATLAB/SPM.
