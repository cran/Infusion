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

.run_cores <- function(method, object, cl, packages=NULL, stat.obs,logLname,verbose) {
  ii <- 0 ## 'global definition' (!)
  prevmsglength <- 0
  if ( ! is.null(cl)) {
    InfusionOptions <- Infusion.options()
    packages <- c("Infusion","blackbox",packages)
    parallel::clusterExport(cl, list(method),envir=environment()) ## passes useks
    parallel::clusterExport(cl, list("packages"),envir=environment()) ## passes the list of packages to load
    abyss <- parallel::clusterEvalQ(cl, {sapply(packages,library,character.only=TRUE)}) ## snif
  } 
  pbopt <- pboptions(nout=min(100,2*length(object)),type="timer")
  Sobs.densities <- pblapply(X=object, FUN = method, cl=cl, stat.obs=stat.obs,logLname=logLname,verbose=verbose)
  pboptions(pbopt)
  if (verbose$final) {
    areValid <- sapply(Sobs.densities,`[`,name="isValid")
    nInvalid <- sum(! areValid)
    if (nInvalid>0L) message(paste(nInvalid,"distributions tagged as 'outlier'(s))"))
  }
  return(Sobs.densities)
}