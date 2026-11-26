# Step 2: HC out-of-fold deviation; patients use the full HC model.
run_normative_deviation <- function(cohorts,fit,cfg,out) {
  hc <- cohorts$HC
  dat <- analytic_data(hc)
  folds <- make_folds(hc$site,cfg$n_folds,cfg$cv_seed)
  cv <- vector("list",cfg$n_folds)
  for (k in seq_len(cfg$n_folds)) {
    test <- which(folds==k)
    f <- fit_model(dat[folds!=k,,drop=FALSE],fit$family,cfg$include_site,max_cycles=cfg$max_cycles)
    d <- calculate_deviation(f,dat[test,,drop=FALSE])
    cv[[k]] <- cbind(row=test,fold=k,d)
  }
  hc_dev <- do.call(rbind,cv)
  hc_dev <- hc_dev[order(hc_dev$row),,drop=FALSE]
  stopifnot(identical(as.integer(hc_dev$row),seq_len(nrow(hc))))
  for (g in names(cohorts)) {
    d <- cohorts[[g]]
    dev <- if (g=="HC") hc_dev[,setdiff(names(hc_dev),"row"),drop=FALSE] else
      cbind(fold=NA_integer_,calculate_deviation(fit,analytic_data(d)))
    meta <- d[,c("subject_id","input_row","input_file","group","age","gender","site"),drop=FALSE]
    result <- cbind(meta,dev)
    write_csv(result,file.path(out,paste0(g,"_GlobalRAD_Deviation.csv")))
    # Keep original input order within each group; CSV holds IDs and site/sex labels.
    mat <- data.frame(input_row=d$input_row,age=d$age,fold=dev$fold,
                      deviation_z=dev$deviation_z,deviationZ=dev$deviation_z,
                      extreme_95=dev$extreme_95,probability=dev$probability,
                      percentile=dev$percentile,conditional_median=dev$conditional_median)
    save_numeric_mat(file.path(out,paste0(g,"_GlobalRAD_Deviation.mat")),mat)
  }
  calibration <- function(v,s) {
    ks <- suppressWarnings(ks.test(v,"pnorm"))
    data.frame(site=s,n=length(v),z_mean=mean(v),z_sd=sd(v),
               proportion_extreme=mean(abs(v)>=1.96),KS_D=as.numeric(ks$statistic),
               KS_p_descriptive=ks$p.value)
  }
  z <- hc_dev$deviation_z
  cal <- do.call(rbind,c(list(calibration(z,"ALL")),lapply(levels(hc$site),function(s) calibration(z[hc$site==s],s))))
  write_csv(cal,file.path(out,"HC_CV_calibration.csv"))
  save_plot(file.path(out,"HC_CV_deviation_histogram.png"),{
    hist(z,breaks=30,probability=TRUE,col="grey80",border="white",main="HC out-of-fold deviation",xlab="Deviation Z")
    curve(dnorm(x),add=TRUE,col="steelblue4",lwd=2)
    abline(v=c(-1.96,1.96),lty=2)
  })
  save_plot(file.path(out,"HC_CV_deviation_QQ.png"),{
    qqnorm(z,main="HC out-of-fold deviation: normal Q-Q plot"); qqline(z,col="steelblue4",lwd=2)
  })
  invisible(cal)
}
