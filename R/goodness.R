.get_gof_stats <- function(projectors, rawstatempdist, 
                           gdimfn=.Infusion.data$options$gof_nstats_fn, 
                           importance_fn=importance
                           ) {
  stats <- unique(unlist(lapply(projectors, attr, which="stats"))) ## all stats used by all projectors
  rawstatempdist <- rawstatempdist[,stats, drop=FALSE] # some simulated raw stats might not have been used in the projections
  
  ## We need to check for linear redundancies in the raw stats 
  rawstatempdist <- check_raw_stats(rawstatempdist,statNames = stats,remove = TRUE, verbose=0.5) # silent if not problem detected, but still verbose otherwise
  stats <- colnames(rawstatempdist)
  ##
  projnames <- ls(projectors)
  # rel_imps are importances for each projector, standardized by the RMSE of oob prediction for each projector
  rel_imps <- try(lapply(projnames, function(v) importance_fn(projectors[[v]])/projectors[[v]]$prediction.error), silent=TRUE)
  if ( ! inherits(rel_imps,"try-error")) { 
    for (it in seq_along(rel_imps)) rel_imps[[it]] <- rel_imps[[it]][intersect(names(rel_imps[[it]]),stats)]
    rel_imps <- do.call(unlist,list(rel_imps)) # list() important here...
    nr <- nrow(rawstatempdist)
    gdim <- gdimfn(nr=nr, nstats=length(stats))
    # find the minimum importance of each raw stat over projectors, and sort raw stats by increasing importance:
    stats <- names(head(sort(sapply(unique(names(rel_imps)), 
                                    function(v) {min(rel_imps[grep(v, names(rel_imps), fixed=TRUE)] )})),
                        n = gdim))
  }
  stats
}

.check_importance_method <- function(projnames, projectors, object) {
  #allranger <- all(sapply(projnames, function(v) inherits( projectors[[v]],"ranger")))
  imp_methods <- unlist(lapply(projnames, function(v) projectors[[v]]$"importance.mode"))
  if (length(imp_methods)==length(projnames)) {
    if (length(u_imp <- unique(imp_methods))==1L) {
      if (u_imp=="none") {
        message("Recomputing projections with importance='permutation' for the goodness-of-fit test...")
        reproject(object, methodArgs=list(importance="permutation"))
      } else if ( ! u_imp %in% c("impurity_corrected","permutation")) {
        warning("Importance method not 'permutation' (preferred) nor 'impurity_corrected': user is responsible to make sure that GoF test will be meaningful.",
                immediate.=TRUE)
      } # In particular, the importance measure must indeed increase with importance...
    } else stop("Different importance methods used for different projectors, cannot be used together for GoF test.")
  } else stop("Variable importance method not found for all projectors.")
}


goftest <- function(object, nsim=99L, method="", stats=NULL, plot.=TRUE,
                    nb_cores=NULL, Simulate=get_from(object,"Simulate"),
                    control.Simulate=get_from(object,"control.Simulate"),
                    packages=get_from(object,"packages"), env=get_from(object,"env"),
                    verbose=interactive(),
                    cl_seed=.update_seed(object),
                    get_gof_stats=.get_gof_stats
) {
  feasible <- FALSE
  colTypes <- object$colTypes
  stat_obs <- t(get_from(object,"stat.obs"))
  if (method=="mixture") { # Result won't be a valid test of goodness of fit !
    statempdist <- simulate(object, nsim=nsim)
    plotframe <- as.data.frame(rbind(statempdist,stat_obs))
    plotmain <- "Obs. vs. mixture model"
  } else { # actual method: requires resimulation of processus
    if (is.null(Simulate)) {
      stop("'Simulate' function no available. Hint: its default value in goftest() call\n     assumes it was made available through the original add_reftable() call.")
    }
    MSLErep <- t(object$MSL$MSLE)
    MSLErep <- as.data.frame(MSLErep[rep(1,nsim),])
    MSLErep <- cbind(MSLErep,object$colTypes$fixedPars) ## add fixedPars for simulation
    MSLErep <- MSLErep[,colTypes$allPars,drop=FALSE] ## column reordering and remove polluting things (cumul_iter in fixedPars!: _F I X M E_)
    if (inherits(object,"SLik_j")) {
      rawstatempdist <- add_reftable(Simulate=Simulate, parsTable=MSLErep, verbose=verbose,
                                control.Simulate=control.Simulate,
                                nb_cores=nb_cores, packages=packages$add_simulation, env=env,
                                cl_seed=cl_seed)     
    } else {
      rawstatempdist <- add_simulation(Simulate=Simulate, parsTable=MSLErep,verbose=verbose,
                                  control.Simulate=control.Simulate,
                                  nb_cores=nb_cores, packages=packages$add_simulation, env=env,
                                  cl_seed=cl_seed)   
    }
    rawstatempdist <- na.omit(rawstatempdist)
    if ( ! is.null(projectors <- object$projectors)) {
      projnames <- ls(projectors)
      .check_importance_method(projnames, projectors, object)
      plotmain <- "Obs. vs. process: summ.stat. residuals."
      feasible <- TRUE
      statempdist <- .project_reftable_raw(rawstatempdist,projectors=projectors, ext_projdata=object$projdata) # if composite param were fitted, the projectors are predictors of composite params
      statempdist <- as.matrix(statempdist[,object$colTypes$statNames,drop=FALSE])
      
      # A good GoF stats should be defined by consideration of some alternative model, but this concept is not available here
      if (is.null(stats)) stats <- .get_gof_stats(projectors, rawstatempdist=rawstatempdist) # trying to get variable with least importance
      if (verbose) cat(paste("The following raw summary statistics were retained for the GoF test:\n",
                                  paste(stats,collapse=", "),"\n"))
      #
      rawstatempdist <- rawstatempdist[,stats,drop=FALSE]
      statsFitToProj <- qr.solve(statempdist,rawstatempdist) # regression to projected stats (potentially predictors of composite params)
      colnames(statsFitToProj) <- stats
      gofStats <- rawstatempdist - statempdist %*% statsFitToProj # residuals from fit to projected stats
      gofObsStat <- attr(stat_obs,"raw_data")[stats] - stat_obs %*% statsFitToProj # residuals idem
      plotframe <- as.data.frame(rbind(gofStats,gofObsStat))
      colnames(plotframe) <- paste0("Res(",colnames(plotframe),")")
    } else {
      plotmain <- "Obs. vs. process"
      statempdist <- rawstatempdist[,plot.,drop=FALSE]
      plotframe <- as.data.frame(rbind(statempdist,stat_obs))
    }
  }
  # interpret plot. as character strings belonging to colnames(plotframe):
  plotnames <- colnames(plotframe)
  if (is.logical(plot.)) {
    if (plot.) {
      plot. <- plotnames[seq(min(8L,length(plotnames)))]
    } else plot. <- c()
  } else if (is.character(plot.)) {
    plot. <- intersect(plot., plotnames)
  } else if (is.numeric(plot.)) { # simplest way to force more than 8:
    plot. <- plotnames[plot.]
  }
  if (length(plot.)) {
    if (length(plot.)<length(plotnames)) message(paste0("Plotting ",length(plot.)," GoF statistics out of ", length(plotnames),":"))
    chk <- try(plot(plotframe,
                    main=plotmain,
                    col=c(rep("black",nsim),"orange"),
                    pch=c(rep(20,nsim),21),
                    bg="orange"))
  }
  if (feasible) {
    safe_nbCluster <- get_nbCluster_range(projdata=gofStats) # this checks that one cluster (or more) can be fitted (.get_gof_stats() presumably handles that, but it provides only a default value.)
    using <- Infusion.getOption("mixturing")
    gofStats <- as.data.frame(gofStats)
    if (using=="mclust") {
      if ( ! "package:mclust" %in% search()) stop("'mclust' should be loaded first.")
      models <- vector("list",length(safe_nbCluster))
      for (it in seq_along(safe_nbCluster)) {
        models[[it]] <- .do_call_wrap("densityMclust",
                                      list(data=gofStats,modelNames=Infusion.getOption("mclustModel"), 
                                           G=safe_nbCluster[it], verbose=FALSE, plot=FALSE),
                                      pack="mclust")
      }
      gofdens <- .get_best_mclust_by_IC(models) 
      # predict(gofdens...) will call predict.densityMclust() which will ignore the undefined solve_t_chol_sigma
    } else {
      gofdens <- .do_call_wrap(
        "mixmodCluster",
        list(data=gofStats, 
             nbCluster=safe_nbCluster, 
             models=.do_call_wrap("mixmodGaussianModel",
                                  list(listModels=Infusion.getOption("mixmodGaussianModel"))), 
             seed=123, 
             strategy=.Infusion.data$options$get_mixModstrategy(nc=ncol(gofStats))))
      gofdens <- .get_best_mixmod_by_IC(gofdens) 
      solve_t_chol_sigma_list <- lapply(gofdens@parameters["variance"], .solve_t_cholfn)
    }
    emp_logLs <- predict(gofdens,newdata=gofStats,solve_t_chol_sigma_list=solve_t_chol_sigma_list,log=TRUE) ## pas de binFactor!
    obs_logL <- predict(gofdens,newdata =gofObsStat,solve_t_chol_sigma_list=solve_t_chol_sigma_list,log=TRUE)
    pval <- (sum(obs_logL>emp_logLs)+1L)/(nsim+1L)
    resu <- list(pval=pval, plotframe=plotframe)
  } else {
    message("No feasible test")
    resu <- list(pval=NULL)
  }
  class(resu) <- c("goftest", class(resu))
  return(resu)
}

summary.goftest <- function(object, ...) utils::str(object)

print.goftest <- function(x, ...) utils::str(x)