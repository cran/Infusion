.calc_init_smooth_fac <- function(object, given=NULL) {
  means <- object$jointdens@parameters@mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  nc <- ncol(means)
  smooth_fac <- numeric(nc)
  colnames(means) <- names(smooth_fac) <- object$jointdens@varNames
  for (stat_st in object$jointdens@statNames ) smooth_fac[stat_st] <- max(dist(means[,stat_st])) # some elements may remain 0
  if (is.null(given)) {
    means <- object$pardens@parameters@mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  } else {
    condpardens <- .conditional_Rmixmod(object$pardens, given=given, expansion=1) 
    means <- condpardens@parameters@mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  }
  for (par_st in colnames(means)) smooth_fac[par_st] <- max(dist(means[,par_st])) # some elements may remain 0
  smooth_fac <- smooth_fac*smooth_fac # we will correct the *cov* matrix
  smooth_fac
}

.smooth_opt_Rmixmod <- function(init, 
                                objfn, # must have a 'solve_t_chol_sigma_lists' argument
                                object, lower, upper, 
                                smooth_fac, # initial vector of smoothing variances 
                                maxit=4L, rate=5, # iterative smoothing parameters
                                ...) {
  nstats <- length(object$jointdens@statNames)
  for (it in 0:maxit) {
    p_smoothing_mat <- diag(x=head(smooth_fac,-nstats)/(rate^it), nrow=length(smooth_fac)-nstats)
    j_smoothing_mat <- diag(x=c(smooth_fac/(rate^it)) , nrow=length(smooth_fac))
    solve_t_chol_sigma_lists <- list(
      pardens= lapply(object$pardens@parameters["variance"], .solve_t_cholfn, smoothing_mat=p_smoothing_mat) ,
      jointdens= lapply(object$jointdens@parameters["variance"], .solve_t_cholfn, smoothing_mat=j_smoothing_mat)
    ) 
    optr <- .safe_opt(init=init, objfn, lower=lower, upper=upper, verbose=FALSE, LowUp=list(), 
                      solve_t_chol_sigma_lists=solve_t_chol_sigma_lists, ...)
    print(optr)
    init <- optr$solution
  }
  solve_t_chol_sigma_lists <- object$solve_t_chol_sigma_lists
  optr <- .safe_opt(init=init, objfn, lower=lower, upper=upper, verbose=FALSE, LowUp=list(), 
                    solve_t_chol_sigma_lists=solve_t_chol_sigma_lists, ...)
  optr
}

.smoothoptim_mlogL_newobs <- function(object, 
                                      givenpars=NULL, # optional profile parameters, but necessary only if their value differ from the MSLE
                                      newobs, # 1-row matrix in projected summstats space 
                                      init=NULL,
                                      lower, upper, # in restricted space if profile
                                      solve_t_chol_sigma_lists=object$solve_t_chol_sigma_lists,
                                      log=TRUE) {
  if (is.null(init)) init <- object$MSL$MSLE[names(lower)]
  template <- object$MSL$MSLE
  if (!is.null(givenpars)) template[names(givenpars)] <- givenpars
  if (inherits(object$jointdens,"Mclust")) { # code not tested (_F I X M E_)
    prof_mlogLfn_Dsim <- function(parmv, solve_t_chol_sigma_lists) { # 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .predict_SLik_j_mclust(object=object, newdata=fullpar, tstat.obs=newobs, log=log, which="",
                               solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)
    }
  } else {
    prof_mlogLfn_Dsim <- function(parmv, solve_t_chol_sigma_lists) { 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .pointpredict.Rmixmod(object=object, X=fullpar, tstat.obs=newobs, log=log, which="", 
                              solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)}
  }
  smooth_fac <- .calc_init_smooth_fac(object, given=givenpars) 
  .smooth_opt_Rmixmod(init=init, objfn=prof_mlogLfn_Dsim, object, lower=lower,upper=upper, smooth_fac)
}


