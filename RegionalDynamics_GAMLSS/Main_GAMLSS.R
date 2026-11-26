# Global RAD GAMLSS (SN2). Run with source() or Rscript.
# install.packages(c("gamlss", "gamlss.dist", "R.matlab"))

# Find the code folder when sourced or run with Rscript.
script_root <- function() {
  paths <- lapply(sys.frames(),function(f) f$ofile)
  paths <- Filter(function(x) is.character(x) && length(x)==1L && file.exists(x),paths)
  if (length(paths)) return(dirname(normalizePath(tail(paths,1L)[[1L]],winslash="/")))
  arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=","",arg[1L]),winslash="/")))
  getwd()
}
CODE_DIR <- script_root()

# Input paths and MAT variable names
INPUT_FILE <- file.path(CODE_DIR,"IndividualGlobalRAD","GlobalRAD.mat")
PATIENT_FILES <- c() # Optional: named SZ/MDD paths, e.g. c(SZ=file.path(CODE_DIR,"SZ.mat"))
INPUT_WITHOUT_GROUP <- "HC"
VARS <- list(y="BayesDataN", age="AgeN", gender="GenderN", site="SitesN", group=NULL, id=NULL)
# Set group/id to explicit MAT field names when present. NULL uses the configured cohort and row IDs.
GROUP_MAP <- c("1"="HC", "2"="SZ", "3"="MDD", HC="HC", SZ="SZ", MDD="MDD", SCZ="SZ", CONTROL="HC")
# Optional per-file variable names (default: same as VARS).
PATIENT_VARS <- list() # e.g. list(SZ=list(y="data",age="age",gender="gender",site="sites"))

# Output paths: choose a new root for a new run; existing filenames are overwritten.
OUTPUT_ROOT <- file.path(CODE_DIR,"results","GlobalRAD")
MODEL_DIR <- file.path(OUTPUT_ROOT,"01_Normative_Modeling")
DEVIATION_DIR <- file.path(OUTPUT_ROOT,"02_Normative_Deviation")
AGE_EFFECT_DIR <- file.path(OUTPUT_ROOT,"03_Age_Effect")

# Analysis settings (unchanged)
FAMILY <- "SN2"
INCLUDE_SITE <- TRUE
N_FOLDS <- 10L
N_AGE_GRID <- 200L
MAX_CYCLES <- 100L
CV_SEED <- 20260714L
SAVE_MODEL <- FALSE # Full model objects may contain individual data.

# Load analysis modules from this script
for (f in c("functions.R","01_Normative_Modeling.R","02_Normative_Deviation.R",
            "03_Age_Effect.R")) source(file.path(CODE_DIR,f),local=TRUE)
require_packages()

main_gamlss <- function() {
  cfg <- list(input_file=INPUT_FILE,patient_files=PATIENT_FILES,patient_vars=PATIENT_VARS,
    input_without_group=INPUT_WITHOUT_GROUP,vars=VARS,group_map=GROUP_MAP,
    include_site=INCLUDE_SITE,family=FAMILY,n_folds=N_FOLDS,
    n_age_grid=N_AGE_GRID,max_cycles=MAX_CYCLES,
    cv_seed=CV_SEED,save_model=SAVE_MODEL)
  stopifnot(N_FOLDS>=2L,N_AGE_GRID>=3L)
  cohorts <- load_global_cohorts(cfg)
  dat <- analytic_data(cohorts$HC)
  for (out in c(MODEL_DIR,DEVIATION_DIR,AGE_EFFECT_DIR))
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
  fit <- run_normative_modeling(dat,FAMILY,cfg,MODEL_DIR)
  run_normative_deviation(cohorts,fit,cfg,DEVIATION_DIR)
  run_age_effect(dat,fit,cfg,AGE_EFFECT_DIR)
  cat("Results: ",OUTPUT_ROOT,"\n",sep="")
  invisible(OUTPUT_ROOT)
}
main_gamlss()
