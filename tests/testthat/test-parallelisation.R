#library("Infusion")
enorm <- function() {Sys.sleep(1); c(r=rnorm(n=1,mean = 0,s=1))}
                                                                         # no more than 2 cores in R checks...
add_reftable(Simulate=enorm, parsTable=data.frame(matrix(ncol=0,nrow=3)), nb_cores=2, cl_seed=123)  # replicable
add_reftable(Simulate=enorm, parsTable=data.frame(matrix(ncol=0,nrow=3)), nb_cores=2, cl_seed=123)  # replicable
if (.Platform$OS.type != "windows") {
  add_reftable(Simulate=enorm, parsTable=data.frame(matrix(ncol=0,nrow=8)), nb_cores=2, cl_seed=123, cluster_args=list(type="FORK"))  # replicable
  add_reftable(Simulate=enorm, parsTable=data.frame(matrix(ncol=0,nrow=8)), nb_cores=2, cl_seed=123, cluster_args=list(type="FORK"))  # replicable
  if (suppressWarnings(do.call("require",list(package="progressr", quietly = TRUE)))) {
    add_reftable(Simulate=enorm, parsTable=data.frame(matrix(ncol=0,nrow=8)), nb_cores=2, cl_seed=123, cluster_args=list(type="FORK"))  # replicable
    add_reftable(Simulate=enorm, parsTable=data.frame(matrix(ncol=0,nrow=8)), nb_cores=2, cl_seed=123, cluster_args=list(type="FORK"))  # replicable
    unloadNamespace("progressr")
  } else {message(paste("'progressr' required but not available.",sep=""))}
}


if (FALSE) {
  # template for devel 
  xs <- 1:8
  with_progress({
    p <- progressor(steps=4)
    yyy <- parallel::mclapply(xs, function(x) {
      res <- enorm()
      p()
      res
    }, mc.cores=2, mc.preschedule = FALSE)
  })
  str(unlist(yyy,use.names = FALSE))
}

