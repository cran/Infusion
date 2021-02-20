goftest <- function(object, nsim=99L, method="", stats=NULL, plot.=TRUE,
                    nb_cores=NULL, Simulate=attr(object$logLs,"Simulate"),
                    packages=attr(object$logLs,"packages"), env=attr(object$logLs,"env"),
                    verbose=interactive()
                    ) {
  feasible <- FALSE
  stat_obs <- t(attr(object$logLs,"stat.obs"))
  if (method=="mixture") { # not a valid test of goodness of fit !
    statdens <- .conditional_Rmixmod(object$jointdens,#fittedPars=object$colTypes$statNames,
                                     given=object$MSL$MSLE, expansion=1) # stat dens|ML parameter estimates
    plotmain <- "Obs. vs. mixture model"
    statempdist <- .simulate.MixmodResults(statdens, nsim=1L, size=nsim, drop=TRUE) # directly in projected space
    plotframe <- as.data.frame(rbind(statempdist,stat_obs))
  } else { # 
    MSLErep <- t(object$MSL$MSLE)
    MSLErep <- as.data.frame(MSLErep[rep(1,nsim),])
    MSLErep <- cbind(MSLErep,object$colTypes$fixedPars) ## add fixedPars for simulation
    MSLErep <- MSLErep[,object$colTypes$allPars,drop=FALSE] ## column reordering and remove polluting things (cumul_iter in fixedPars!: _F I X ME_)
    if (inherits(object,"SLik_j")) {
      statempdist <- add_reftable(Simulate=Simulate, par.grid=MSLErep, verbose=verbose,
                                control.Simulate=attr(object$logLs,"control.Simulate"),
                                nb_cores=nb_cores, packages=packages$add_simulation, env=env)     
    } else {
      statempdist <- add_simulation(Simulate=Simulate, par.grid=MSLErep,verbose=verbose,
                                  control.Simulate=attr(object$logLs,"control.Simulate"),
                                  nb_cores=nb_cores, packages=packages$add_simulation, env=env)   
    }
    if ( ! is.null(projectors <- object$projectors)) {
      plotmain <- "Obs. vs. process: summ.stat. residuals."
      feasible <- TRUE
      rawstatempdist <- statempdist
      statempdist <- project(rawstatempdist,projectors=projectors)
      statempdist <- as.matrix(statempdist[,object$colTypes$statNames,drop=FALSE])
      if (is.null(stats)) stats <- unique(unlist(lapply(projectors, attr, which="stats"))) ## all stats used by all projectors
      rawstatempdist <- rawstatempdist[,stats,drop=FALSE]
      statsFitToProj <- qr.solve(statempdist,rawstatempdist)
      colnames(statsFitToProj) <- stats
      gofStats <- rawstatempdist - statempdist %*% statsFitToProj
      gofObsStat <- attr(stat_obs,"raw_data")[stats] - stat_obs %*% statsFitToProj
      plotframe <- as.data.frame(rbind(gofStats,gofObsStat))
      colnames(plotframe) <- paste0("Res(",colnames(plotframe),")")
    } else {
      plotmain <- "Obs. vs. process"
      statempdist <- statempdist[,plot.,drop=FALSE]
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
    plot(plotframe,
         main=plotmain,
         col=c(rep("black",nsim),"orange"),
         pch=c(rep(20,nsim),21),
         bg="orange")
  }
  if (feasible) {
    nbCluster <- eval(Infusion.getOption("nbCluster"),list2env(list(data=gofStats)))
    using <- Infusion.getOption("mixturing")
    gofStats <- as.data.frame(gofStats)
    if (using=="mclust") {
      if ( ! "package:mclust" %in% search()) stop("'mclust' should be loaded first.")
      models <- vector("list",length(nbCluster))
      for (it in seq_along(nbCluster)) {
        models[[it]] <- .do_call_wrap("densityMclust",
                                      list(data=gofStats,modelNames=Infusion.getOption("mclustModel"), 
                                           G=nbCluster[it], verbose=FALSE),
                                      pack="mclust")
      }
      gofdens <- .get_best_mclust_by_IC(models) 
    } else {
      gofdens <- .do_call_wrap("mixmodCluster",list(data=gofStats, 
                                                    nbCluster=nbCluster, 
                                                    models=.do_call_wrap("mixmodGaussianModel",list(listModels=Infusion.getOption("mixmodGaussianModel"))), 
                                                    seed=123, 
                                                    strategy=eval(Infusion.getOption("strategy"))))
      gofdens <- .get_best_mixmod_by_IC(gofdens) 
    }
    solve_t_chol_sigma <- lapply(gofdens@parameters["variance"], function(mat) solve(t(chol(mat))))
    emp_logLs <- predict(gofdens,newdata=gofStats,solve_t_chol_sigma=solve_t_chol_sigma,log=TRUE) ## pas de binFactor!
    obs_logL <- predict(gofdens,newdata =gofObsStat,solve_t_chol_sigma=solve_t_chol_sigma,log=TRUE)
    pval <- (sum(obs_logL>emp_logLs)+1L)/(nsim+1L)
    return(list(pval=pval))
  } else {
    message("No feasible test")
    return(list(pval=NULL))
  }
}
