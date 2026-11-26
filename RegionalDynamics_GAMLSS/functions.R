# Shared GAMLSS implementation. Input is already harmonized; no ComBat is fitted here.
require_packages <- function() {
  pkgs <- c("gamlss", "gamlss.dist", "R.matlab")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install packages: ", paste(missing, collapse=", "))
  suppressPackageStartupMessages(library(gamlss))
  suppressPackageStartupMessages(library(gamlss.dist))
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  utils::write.csv(x, path, row.names=FALSE, na="")
}

family_object <- function(name) {
  f <- getExportedValue("gamlss.dist",name)()
  if (!inherits(f,"gamlss.family")) stop("Invalid family: ",name)
  f
}

fit_model <- function(dat, family="SN2", include_site=TRUE, reduced=FALSE, max_cycles=100L, crit=.001) {
  dat <- droplevels(dat)
  if (nrow(dat)<50L || length(unique(dat$y))<5L || length(unique(dat$age))<5L) stop("Insufficient observations or unique outcomes.")
  if (nlevels(dat$gender)<2L || (include_site && nlevels(dat$site)<2L)) stop("Insufficient factor levels in training data.")
  mu_formula <- stats::as.formula(paste("y ~", paste(c(if(!reduced) "pb(age)","gender",if(include_site) "site"),collapse=" + ")))
  fam <- family_object(family)
  args <- list(formula=mu_formula, family=fam, data=dat,
               control=gamlss::gamlss.control(c.crit=crit,n.cyc=max_cycles,trace=FALSE))
  if (fam$nopar>=2L) args$sigma.formula <- ~pb(age)
  if (fam$nopar>=3L) args$nu.formula <- ~1
  if (fam$nopar>=4L) args$tau.formula <- ~1
  warns <- character()
  m <- withCallingHandlers(do.call(gamlss::gamlss,args),warning=function(w) {
    warns <<- c(warns,conditionMessage(w)); invokeRestart("muffleWarning")
  })
  if (!isTRUE(m$converged) || length(m$sbc)!=1L || !is.finite(m$sbc)) stop("GAMLSS did not converge: ",paste(unique(warns),collapse=" | "))
  list(model=m,data=dat,family=family,include_site=include_site,reduced=reduced,warnings=unique(warns))
}

predict_parameters <- function(fit, nd) {
  nd <- as.data.frame(nd[,c("age","gender","site"),drop=FALSE])
  if (any(!as.character(nd$gender) %in% levels(fit$data$gender))) stop("Unseen gender level in prediction.")
  if (fit$include_site && any(!as.character(nd$site) %in% levels(fit$data$site)))
    stop("Patient/held-out site absent from HC training data; cannot estimate its fixed effect.")
  nd$gender <- factor(as.character(nd$gender),levels=levels(fit$data$gender))
  if (fit$include_site) nd$site <- factor(as.character(nd$site),levels=levels(fit$data$site))
  else nd$site <- NULL # Site is not a predictor in the optional no-site model.
  pars <- intersect(c("mu","sigma","nu","tau"),fit$model$parameters)
  ans <- setNames(lapply(pars,function(p) as.numeric(stats::predict(fit$model,what=p,newdata=nd,type="response",data=fit$data))),pars)
  if (any(vapply(ans,function(v) length(v)!=nrow(nd) || any(!is.finite(v)),logical(1)))) stop("Invalid parameter predictions.")
  ans
}

family_quantile <- function(fit, pars, p) {
  as.numeric(do.call(getExportedValue("gamlss.dist",paste0("q",fit$family)),c(list(p=rep(p,length(pars$mu))),pars)))
}

standardization_cells <- function(dat) {
  tab <- as.data.frame(table(gender=dat$gender,site=dat$site),stringsAsFactors=FALSE)
  tab <- tab[tab$Freq>0,,drop=FALSE]
  tab$weight <- tab$Freq/sum(tab$Freq)
  tab
}

age_derivative <- function(age,y) {
  n <- length(age)
  stopifnot(n>=3L,length(y)==n,all(diff(age)>0))
  c((y[2]-y[1])/(age[2]-age[1]),(y[3:n]-y[1:(n-2)])/(age[3:n]-age[1:(n-2)]),(y[n]-y[n-1])/(age[n]-age[n-1]))
}

standardized_curve <- function(fit, grid, cells=standardization_cells(fit$data), probs=c(.025,.05,.25,.5,.75,.95,.975)) {
  # Weighted average of conditional centiles, not a quantile of a mixture distribution.
  k <- nrow(cells)
  nd <- data.frame(age=rep(grid,each=k),gender=rep(cells$gender,length(grid)),site=rep(cells$site,length(grid)))
  pars <- predict_parameters(fit,nd)
  avg <- function(v) {
    if (length(v)!=nrow(nd) || any(!is.finite(v))) stop("Invalid curve prediction.")
    as.numeric(crossprod(cells$weight,matrix(v,nrow=k)))
  }
  out <- data.frame(age=grid,mu=avg(pars$mu))
  for (p in setdiff(names(pars),"mu")) out[[p]] <- avg(pars[[p]])
  for (p in probs) out[[if (abs(100*p-round(100*p))<1e-8) sprintf("p%02d",as.integer(round(100*p))) else paste0("p",format(100*p,trim=TRUE,scientific=FALSE))]] <- avg(family_quantile(fit,pars,p))
  if ("p50" %in% names(out)) out$P50_derivative <- age_derivative(grid,out$p50)
  out$mu_derivative <- age_derivative(grid,out$mu)
  out
}

calculate_deviation <- function(fit, dat, eps=1e-10, threshold=1.96) {
  pars <- predict_parameters(fit,dat)
  p <- as.numeric(do.call(getExportedValue("gamlss.dist",paste0("p",fit$family)),c(list(q=dat$y),pars)))
  if (length(p)!=nrow(dat) || any(!is.finite(p)) || any(p<0 | p>1)) stop("Invalid CDF predictions.")
  z <- stats::qnorm(pmin(pmax(p,eps),1-eps))
  med <- family_quantile(fit,pars,.5)
  data.frame(y=dat$y,conditional_mu=pars$mu,conditional_median=med,raw_from_median=dat$y-med,
             probability=p,percentile=100*p,deviation_z=z,extreme_95=as.integer(abs(z)>=threshold))
}

make_folds <- function(site,k,seed) {
  set.seed(seed)
  counts <- table(site)
  if (any(counts<2L)) stop("Each HC site needs at least two participants for site-stratified CV.")
  f <- integer(length(site))
  for (s in names(counts)) {
    idx <- which(site==s)
    idx <- idx[sample.int(length(idx))]
    labels <- rep(seq_len(k),length.out=length(idx))
    f[idx] <- labels[sample.int(length(labels))]
  }
  if (!setequal(unique(f),seq_len(k))) stop("Some CV folds are empty.")
  f
}

plot_centiles <- function(dat, curve, path, title) {
  grDevices::png(path, width=1800, height=1200, res=200)
  on.exit(grDevices::dev.off(), add=TRUE)
  plot(dat$age, dat$y, pch=16, cex=.45, col=adjustcolor("grey40",alpha.f=.25),
       xlab="Age (years)", ylab="Metric value", main=title,
       ylim=range(c(dat$y,curve$p05,curve$p95)))
  for (nm in c("p05","p25","p75","p95")) lines(curve$age,curve[[nm]],col="steelblue",lty=2)
  lines(curve$age,curve$p50,col="steelblue4",lwd=2)
}

# ---- Global MAT input (one value per participant) ----
# Explicit MAT field names are configured in Main_GAMLSS.R.
mat_vector <- function(v, label, n=NULL, numeric_only=FALSE) {
  # Numeric row/column vectors and MATLAB cell arrays of scalar strings.
  if (is.list(v)) {
    parts <- lapply(v, function(z) {
      z <- unlist(z,recursive=TRUE,use.names=FALSE)
      if (is.character(z) && length(z)>1L && all(nchar(z)==1L)) z <- paste0(z,collapse="")
      if (length(z)!=1L) stop(label,": non-scalar MATLAB cell entry.")
      z
    })
    v <- unlist(parts,use.names=FALSE)
  } else if (is.character(v) && !is.null(n) && length(dim(v))==2L &&
             nrow(v)==n && ncol(v)>1L && all(nchar(v)==1L)) {
    v <- apply(v,1L,paste0,collapse="")
  } else if (!is.null(dim(v)) && sum(dim(v)>1L)>1L) {
    stop(label,": expected a single row/column vector, not a multi-column matrix.")
  }
  v <- as.vector(v)
  if (!is.null(n) && length(v)!=n) stop(label,": expected ",n," entries; found ",length(v))
  if (numeric_only) {
    if (is.complex(v)) stop(label,": complex values are not supported.")
    v <- suppressWarnings(as.numeric(v))
    if (any(!is.finite(v))) stop(label,": non-finite/non-numeric entries; clean input first.")
  } else {
    if (is.numeric(v) && any(!is.finite(v))) stop(label,": non-finite entries.")
    v <- trimws(as.character(v))
    if (anyNA(v) || any(!nzchar(v))) stop(label,": missing/empty entries.")
  }
  v
}

read_global_mat <- function(path, vars, group_map, default_group=NULL) {
  if (!file.exists(path)) stop("MAT file not found: ",path)
  x <- tryCatch(R.matlab::readMat(path,fixNames=FALSE),error=function(e) {
    stop("Cannot read MAT file: ",conditionMessage(e),
         ". For MATLAB v7.3/HDF5 files, resave the required vectors using save(...,'-v7').")
  })
  keys <- vars
  required <- c("y","age","gender","site")
  for (k in c(required,"group","id")) {
    field <- keys[[k]]
    if (is.null(field) && !k %in% required) next
    if (length(field)!=1L || !field %in% names(x))
      stop("Missing MAT field for ",k," in ",path,"; edit VARS in Main_GAMLSS.R.")
  }
  y <- mat_vector(x[[keys$y]],"RAD",numeric_only=TRUE)
  n <- length(y)
  if (!n) stop("Empty RAD vector.")
  age <- mat_vector(x[[keys$age]],"age",n,TRUE)
  gender <- mat_vector(x[[keys$gender]],"gender",n)
  site <- mat_vector(x[[keys$site]],"site",n)
  if (any(age<0)) stop("Negative age: verify age units and values.")
  if (is.null(keys$group)) {
    if (is.null(default_group) || !default_group %in% c("HC","SZ","MDD"))
      stop("No diagnosis field. Explicitly specify the group of this file.")
    group <- rep(default_group,n)
  } else {
    raw <- toupper(mat_vector(x[[keys$group]],"group",n))
    map <- setNames(as.character(group_map),toupper(names(group_map)))
    group <- unname(map[raw])
    if (anyNA(group) || any(!group %in% c("HC","SZ","MDD")))
      stop("Unrecognized diagnosis codes: ",paste(unique(raw[is.na(group) | !group %in% c("HC","SZ","MDD")]),collapse=", "),
           ". Edit GROUP_MAP in Main_GAMLSS.R.")
  }
  id <- if (is.null(keys$id)) paste0(group,"_row_",seq_len(n)) else mat_vector(x[[keys$id]],"subject ID",n)
  if (anyDuplicated(paste(group,id,sep="|"))) stop("Duplicate subject IDs within group.")
  data.frame(y=y,age=age,gender=factor(gender),site=factor(site),group=group,
             subject_id=id,input_row=seq_len(n),input_file=basename(path),stringsAsFactors=FALSE)
}

load_global_cohorts <- function(cfg) {
  all <- read_global_mat(cfg$input_file,cfg$vars,cfg$group_map,cfg$input_without_group)
  for (g in names(cfg$patient_files)) {
    if (!g %in% c("SZ","MDD")) stop("PATIENT_FILES names must be SZ or MDD.")
    if (g %in% all$group) stop(g," already present in main file; do not supply it twice.")
    vars <- cfg$patient_vars[[g]]
    if (is.null(vars)) vars <- cfg$vars
    d <- read_global_mat(cfg$patient_files[[g]],vars,cfg$group_map,g)
    if (any(d$group!=g)) stop("Separate ",g," file contains another diagnosis.")
    all <- rbind(all,d)
  }
  cohorts <- lapply(split(all,all$group),droplevels)
  hc <- cohorts$HC
  if (is.null(hc)) stop("No HC observations. HC are required for normative fitting.")
  if (nrow(hc)<max(56L,cfg$n_folds)) stop("Too few HC for the configured model and CV.")
  for (g in setdiff(names(cohorts),"HC")) {
    d <- cohorts[[g]]
    if (any(!as.character(d$gender) %in% levels(hc$gender))) stop(g,": gender absent from HC.")
    if (cfg$include_site && any(!as.character(d$site) %in% levels(hc$site)))
      stop(g,": sites absent from HC: ",paste(setdiff(as.character(d$site),levels(hc$site)),collapse=", "))
    if (any(d$age<min(hc$age) | d$age>max(hc$age))) warning(g,": age extrapolation beyond HC age range.")
  }
  cohorts
}

analytic_data <- function(d) droplevels(d[,c("y","age","gender","site"),drop=FALSE])

save_plot <- function(path, expr) {
  grDevices::png(path,width=1800,height=1200,res=200)
  on.exit(grDevices::dev.off(),add=TRUE)
  force(expr)
}

save_numeric_mat <- function(path, d) {
  # MATLAB-compatible numeric columns only; original IDs/labels remain in CSV.
  cols <- d[vapply(d,is.numeric,logical(1)) | vapply(d,is.logical,logical(1))]
  do.call(R.matlab::writeMat,c(list(con=path),lapply(cols,function(v) matrix(as.numeric(v),ncol=1))))
}
