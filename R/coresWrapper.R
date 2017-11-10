.check_nb_cores <- function(nb_cores=NULL) {
  if (is.null(nb_cores)) nb_cores <- Infusion.getOption("nb_cores") ## may be NULL
  machine_cores <- parallel::detectCores()
  if (is.null(nb_cores)) {
    nb_cores <- 1L ## default
    if (machine_cores>1L && interactive()) {
      if (! identical(Infusion.getOption("cores_avail_warned"),TRUE)) {
        message(paste(machine_cores,"cores are available for parallel computation\n(you may be allowed to fewer of them on a shared cluster).\n"))
        message("Change 'nb_cores' argument to use some of them.\nUse Infusion.options(nb_cores=<n>) to control nb_cores globally.")
        Infusion.options(cores_avail_warned=TRUE)
      }
    }
  } else if (nb_cores>machine_cores) {
    if (! identical(Infusion.getOption("nb_cores_warned"),TRUE)) {
      warning("More cores were requested than found by parallel::detectCores(). Check Infusion.getOption(\"nb_cores\") argument.
                Number of availlable cores automatically downsized to the number of cores found by parallel::detectCores(). I continue.")
      Infusion.options(nb_cores_warned=TRUE)
    }
    nb_cores <- machine_cores
  }
  return(nb_cores)
  }

.init_cores <- function(nb_cores=NULL, ## passing explicit value from user
                        ...) {  ## ... are arguments used by functions called by the loc_calc_logL function
  nb_cores <- .check_nb_cores(nb_cores=nb_cores)
  if (nb_cores > 1L) {
    cl <- parallel::makeCluster(nb_cores) 
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    if (has_doSNOW <- ("doSNOW" %in% .packages() )) { ## allows progressbar but then requires foreach
      # loading (?) the namespace of 'snow' changes the global RNG state!
      assign(".Random.seed", R.seed, envir = .GlobalEnv)
      eval(as.call(c(quote(registerDoSNOW),list(cl=cl)))) 
    } else if ( ! identical(Infusion.getOption("doSNOW_warned"),TRUE)) {
      message("If the 'doSNOW' package were attached, the progress of the computation could be reported.")
      Infusion.options(doSNOW_warned=TRUE)
    } 
    dotenv <- list2env(list(...))
    parallel::clusterExport(cl=cl, as.list(ls(dotenv)),envir=dotenv) 
    return(list(cl=cl,nb_cores=nb_cores,has_doSNOW=has_doSNOW))
  } else return(list(nb_cores=1L)) ## $has_doSNOW and $cl are NULL
  
}

.run_cores <- function(method, object, cores_info, packages=NULL, stat.obs,logLname,verbose) {
  ii <- 0 ## 'global definition' (!)
  prevmsglength <- 0
  if (cores_info$nb_cores > 1L) {
    InfusionOptions <- Infusion.options()
    packages <- c("Infusion","blackbox",packages)
    if (cores_info$has_doSNOW) {
      pb <- txtProgressBar(max = length(object), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      parallel::clusterExport(cl=cores_info$cl, list("progress",method),envir=environment()) 
      `%foreachdopar%` <- foreach::`%dopar%`
      Sobs.densities <- foreach::foreach(
        ii = 1:length(object),
        #.combine = "rbind", ## leave as a list
        .inorder = TRUE,
        .packages = packages,
        .errorhandling = "remove", ## "pass" to get the error objects in the list: useful for debugging
        .options.snow = opts
      ) %foreachdopar% {
        do.call(method, list(EDF=object[[ii]], stat.obs=stat.obs,logLname=logLname,verbose=verbose))
      }
      close(pb)
    } else {
      parallel::clusterExport(cores_info$cl, list(method),envir=environment()) ## passes useks
      parallel::clusterExport(cores_info$cl, list("packages"),envir=environment()) ## passes the list of packages to load
      abyss <- parallel::clusterEvalQ(cores_info$cl, {sapply(packages,library,character.only=TRUE)}) ## snif
      Sobs.densities <- parallel::parLapply(cores_info$cl,X=object, fun = method, stat.obs=stat.obs,logLname=logLname,verbose=verbose)
    }
    if (verbose$final) {
      areValid <- sapply(Sobs.densities,`[`,name="isValid")
      nInvalid <- sum(! areValid)
      if (nInvalid>0L) message(paste(nInvalid,"distributions tagged as 'outlier'(s))"))
    }
  } else {
    prevmsglength <- 0L
    nInvalid <- 0L
    lit <- 0L
    Sobs.densities <- lapply(object, function(element) {
      par_logL_indic <- do.call(method, list(EDF=element, stat.obs=stat.obs,logLname=logLname,verbose=verbose))
      # par_logL_indic is vector of pars + logL + isValid
      if (! par_logL_indic["isValid"]) nInvalid <<- nInvalid+1L
      lit <<- lit+1L
      if (verbose$most) {
        msg <- paste("Already ", lit, " distributions smoothed",sep="")
        if (nInvalid>0L) msg <- paste(msg," (",nInvalid," tagged as 'outlier'(s))",sep="")
        prevmsglength <<- .overcat(msg, prevmsglength)
      }
      return(par_logL_indic)
    }) 
    if (verbose$final) {
      msg <- paste(lit, " distributions smoothed",sep="")
      if (nInvalid>0L) msg <- paste(msg," (",nInvalid," tagged as 'outlier'(s))",sep="")
      msg <- paste(msg,"            ")
      prevmsglength <<- .overcat(msg, prevmsglength)
    }
    if (any(unlist(verbose))) cat("\n")
  }
  return(Sobs.densities)
}