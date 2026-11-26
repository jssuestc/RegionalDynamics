# Step 1: Full HC GAMLSS model and standardized conditional centiles.
run_normative_modeling <- function(dat,family,cfg,out) {
  fit <- fit_model(dat,family,cfg$include_site,max_cycles=cfg$max_cycles)
  grid <- seq(min(dat$age),max(dat$age),length.out=cfg$n_age_grid)
  curve <- standardized_curve(fit,grid)
  write_csv(curve,file.path(out,"GlobalRAD_centiles.csv"))
  save_numeric_mat(file.path(out,"GlobalRAD_centiles.mat"),curve)
  plot_centiles(dat,curve,file.path(out,"GlobalRAD_centiles.png"),"Global RAD: HC normative trajectories")
  write_csv(data.frame(family=family,n=nrow(dat),BIC=fit$model$sbc,
                       AIC=fit$model$aic,global_deviance=fit$model$G.deviance,
                       df_fit=fit$model$df.fit,converged=fit$model$converged,
                       warning=paste(fit$warnings,collapse=" | ")),
            file.path(out,"model_summary.csv"))
  if (cfg$save_model) saveRDS(fit,file.path(out,"GlobalRAD_HC_model_PRIVATE.rds"))
  fit
}
