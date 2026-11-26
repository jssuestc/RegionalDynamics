# Regional functional dynamics: RAD and iRTD

## System requirements and dependencies

The analysis scripts are implemented in MATLAB and R and can be executed on a standard CPU-based workstation without GPU acceleration.

### MATLAB environment

- MATLAB R2021a or later (recommended)
- Statistics and Machine Learning Toolbox
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) for reading and writing NIfTI images
- [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) fsaverage registered spherical surfaces for spatial permutation tests

Spatial permutation tests require the following files from the full-resolution fsaverage template:

```text
fsaverage/surf/lh.sphere.reg
fsaverage/surf/rh.sphere.reg
```

The registered spheres are not included. Place them in `fsaverage/surf/` at the repository root, or edit `cfg.SphereDir` in `Run_Individual.m` and `sphereDir` in `Main_GroupMaps.m`. A complete FreeSurfer installation is not required if these files are available. The Brainnetome cortical annotation files and spatial permutation functions are included. Set `cfg.SPMDir` in `Run_Individual.m` if SPM12 is not already on the MATLAB path.

### R environment

R 4.3.x or later is recommended, with the following packages:

- `gamlss`
- `gamlss.dist`
- `R.matlab`

The required packages can be installed from CRAN:

```r
install.packages(c("gamlss", "gamlss.dist", "R.matlab"))
```

RStudio is optional but recommended for interactive execution. Input MAT files should be saved in MATLAB v7 format rather than v7.3 format.

### Hardware requirements

A standard workstation with at least 16 GB of RAM is recommended. Additional memory (32 GB or more) may be beneficial for voxel-wise FCD estimation. Computational time depends on data size and the number of bootstrap and spatial permutation iterations.

## Running the code

The repository contains three independent analysis modules:

| Module | Main script | Description |
|---|---|---|
| Individual-level analysis | [`Run_Individual.m`](RegionalDynamics_Individual/RegionalDynamics_Individual/Run_Individual.m) | FCD, RAD, iRTD, spatial coupling, relative balance, and individual-level PLS analysis |
| Group-level analysis | [`Main_GroupMaps.m`](RegionalDynamics_Group/Main_GroupMaps.m) | Dominance analysis, clustering, and symptom-related PLS analysis |
| Normative modelling | [`Main_GAMLSS.R`](RegionalDynamics_GAMLSS/Main_GAMLSS.R) | GAMLSS fitting, normative deviation estimation, age-effect analysis, and first-derivative estimation |

Bundled input files are located relative to the corresponding main script, so the repository can be moved without changing the default data paths. Input paths, external software paths and output directory names can be edited in the main scripts. Results are written to each module's `results/` directory by default; existing files with the same names may be overwritten. Generated results are excluded by the module-level `.gitignore` files.

The MATLAB scripts can be executed directly in MATLAB. The R workflow can be executed in RStudio or from the command line:

```bash
Rscript RegionalDynamics_GAMLSS/Main_GAMLSS.R
```

The normative modelling demonstration uses the prespecified SN2 distribution, as selected in the full study. Distribution selection is not repeated in this demonstration.

See the module-specific README files for input formats and statistical definitions. The software versions and hardware configurations above are recommendations, not verified minimum requirements.
