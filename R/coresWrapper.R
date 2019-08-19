.check_nb_cores <- local({
  nb_cores_warned <- FALSE
  cores_avail_warned <- FALSE
  function(nb_cores=NULL) {
    if (is.null(nb_cores)) nb_cores <- Infusion.getOption("nb_cores") ## may be NULL
    machine_cores <- parallel::detectCores()
    if (is.null(nb_cores)) {
      nb_cores <- 1L ## default
      if (machine_cores>1L && interactive()) {
        if (! cores_avail_warned) {
          message(paste(machine_cores,"cores are available for parallel computation\n(you may be allowed to fewer of them on a shared cluster).\n"))
          message("Change 'nb_cores' argument to use some of them.\nUse Infusion.options(nb_cores=<n>) to control nb_cores globally.")
          cores_avail_warned <<- TRUE
        }
      }
    } else if (nb_cores>machine_cores) {
      if (! nb_cores_warned) {
        warning("More cores were requested than found by parallel::detectCores(). Check Infusion.getOption(\"nb_cores\") argument.
                Number of availlable cores automatically downsized to the number of cores found by parallel::detectCores(). I continue.")
        nb_cores_warned <<- TRUE
      }
      nb_cores <- machine_cores
    }
    return(nb_cores)
  }
})

.init_cores <- local({
  doSNOW_warned <- FALSE
  function(nb_cores=NULL, ## passing explicit value from user
           ...) {  ## ... are arguments used by functions called by the loc_calc_logL function
    nb_cores <- .check_nb_cores(nb_cores=nb_cores)
    cores_info <- list(nb_cores=nb_cores)
    if (nb_cores > 1L) {
      cores_info$cl <- parallel::makeCluster(nb_cores) 
      dotenv <- list2env(list(...))
      parallel::clusterExport(cl=cores_info$cl, as.list(ls(dotenv)),envir=dotenv) 
      ## foreach is NOT a parallel backend so there is no point using it if doSNOW is not available
      if (cores_info$has_doSNOW <- (isNamespaceLoaded("doSNOW"))) {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        ## allows progressbar but then requires foreach
        assign(".Random.seed", R.seed, envir = .GlobalEnv) # loading (?) the namespace of 'snow' changes the global RNG state!
        fn <- get("registerDoSNOW", asNamespace("doSNOW"))
        do.call(fn,list(cl=cores_info$cl)) 
      } else {
        if ( ! doSNOW_warned) {
          message("If the 'doSNOW' package were attached, better load-balancing might be possible.")
          doSNOW_warned <<- TRUE
        } 
      }
    }
    return(cores_info)
  }
})

.eval_Sobs_densities <- function(method, object, cores_info, packages=NULL, stat.obs,logLname,
                                 verbose # list
                                 ) {
  #prevmsglength <- 0 ## no longer used ?
  if (cores_info$nb_cores > 1L) {
    #blackboxOptions <- blackbox.options() ## FIXME: none of the package options are passed to the child processed
    #InfusionOptions <- Infusion.options()
    packages <- c("Infusion","blackbox",packages)
    parallel::clusterExport(cores_info$cl, list(method),envir=environment()) ## passes useks
    parallel::clusterExport(cores_info$cl, list("packages"),envir=environment()) ## passes the list of packages to load
    abyss <- parallel::clusterEvalQ(cores_info$cl, {sapply(packages,library,character.only=TRUE)}) ## snif
    if (cores_info$has_doSNOW) {
      show_pb <- (verbose$most && ! isTRUE(getOption('knitr.in.progress')))
      if (show_pb) {
        pb <- txtProgressBar(max = length(object), style = 3, char="P")
        progress <- function(n) setTxtProgressBar(pb, n)
        parallel::clusterExport(cl=cores_info$cl, list("progress"),envir=environment()) ## slow! why?
        .options.snow = list(progress = progress)
      } else .options.snow = NULL
      ii <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'ii'
      foreach_args <- list(
        ii = seq_len(length(object)), 
        .packages= packages,
        .options.snow = .options.snow,
        .inorder = TRUE, .errorhandling = "remove"
        #                                 "pass"## "pass" to see error messages
      )
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      Sobs.densities <- foreach::`%dopar%`(foreach_blob,
                                     do.call(method,c(list(object[[ii]]),
                                                      list(stat.obs=stat.obs,logLname=logLname,verbose=verbose) )))
      if (show_pb) close(pb)
    } else {
      pbopt <- pboptions(nout=min(100,2*length(object)),type="timer", char="p")
      Sobs.densities <- pblapply(X=object, FUN = method, cl=cores_info$cl, stat.obs=stat.obs,logLname=logLname,verbose=verbose)
      pboptions(pbopt)
    }
  } else { 
    pbopt <- pboptions(nout=min(100,2*length(object)),type="timer", char="s")
    Sobs.densities <- pblapply(X=object, FUN = method, cl=NULL, stat.obs=stat.obs,logLname=logLname,verbose=verbose)
    pboptions(pbopt)
  }
  if (verbose$final) {
    areValid <- sapply(Sobs.densities,`[`,name="isValid")
    nInvalid <- sum(! areValid)
    if (nInvalid>0L) message(paste(nInvalid,"distributions tagged as 'outlier'(s))"))
  }
  return(Sobs.densities)
}

.close_cores <- function(cores_info) {
  if ( cores_info$nb_cores > 1L) {
    if (cores_info$has_doSNOW) foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
    parallel::stopCluster(cores_info$cl)
  }
}
