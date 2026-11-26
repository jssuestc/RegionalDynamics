# Step 3: Full/reduced model comparison; age removed from mu only.
run_age_effect <- function(dat,full,cfg,out) {
  reduced <- fit_model(dat,full$family,cfg$include_site,reduced=TRUE,max_cycles=cfg$max_cycles)
  stopifnot(identical(as.numeric(full$model$y),as.numeric(reduced$model$y)))
  sse_full <- sum((dat$y-full$model$mu.fv)^2)
  sse_reduced <- sum((dat$y-reduced$model$mu.fv)^2)
  partial_R2 <- if (is.finite(sse_reduced) && sse_reduced>1e-12) (sse_reduced-sse_full)/sse_reduced else NA_real_
  grid <- seq(min(dat$age),max(dat$age),length.out=cfg$n_age_grid)
  cr <- standardized_curve(full,grid,probs=numeric(0))
  mean_derivative <- mean(cr$mu_derivative)
  signed_R2 <- if (mean_derivative<0) -partial_R2 else partial_R2
  delta_dev <- reduced$model$G.deviance-full$model$G.deviance
  delta_df <- full$model$df.fit-reduced$model$df.fit
  p_age <- if (is.finite(delta_df) && delta_df>0 && is.finite(delta_dev)) pchisq(max(0,delta_dev),df=delta_df,lower.tail=FALSE) else NA_real_
  result <- data.frame(family=full$family,n=nrow(dat),SSE_full=sse_full,SSE_reduced=sse_reduced,
    partial_R2_age=partial_R2,mean_mu_derivative=mean_derivative,signed_partial_R2_age=signed_R2,
    BIC_full=full$model$sbc,BIC_reduced=reduced$model$sbc,
    deviance_full=full$model$G.deviance,deviance_reduced=reduced$model$G.deviance,
    df_full=full$model$df.fit,df_reduced=reduced$model$df.fit,
    deviance_difference=delta_dev,df_difference=delta_df,p_age_approx=p_age,
    reduced_warning=paste(reduced$warnings,collapse=" | "))
  write_csv(result,file.path(out,"GlobalRAD_AgeEffect.csv"))
  save_numeric_mat(file.path(out,"GlobalRAD_AgeEffect.mat"),result)
  write_csv(cr,file.path(out,"GlobalRAD_mu_curve.csv"))
  if (cfg$save_model) saveRDS(reduced,file.path(out,"GlobalRAD_reduced_model_PRIVATE.rds"))
  invisible(result)
}
