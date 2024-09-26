\name{get_workflow_design}
\alias{get_workflow_design}
\title{
Workflow design
}
\description{

\code{get_workflow_design} provides default control values for the simulation plan, ideally chosen for good performance. The default design is ilustrated in the examples.  
}
\usage{
get_workflow_design(npar, n_proj_stats=npar, n_latent=0L,
                    final_reft_size=NULL,
                    refine_blocksize=NULL, subblock_nbr=NULL,
                    version=packageVersion("Infusion"),
                    cumn_over_maxit = NULL,
                    test_fac=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npar}{
  Number of fitted parameters of the statistical model.
}
\item{n_proj_stats}{number of projected summary statistics.}
\item{n_latent}{Number of latent variables to be predicted.}

\code{npar}, \code{n_proj_stats} and \code{n_latent} are here distinguished for clarity, but only their sum is currently used to define the sampling design.  

  \item{final_reft_size}{
  NULL, or integer specifying a non-default value of the final reference table size.
}
  \item{refine_blocksize}{
  NULL, or integer specifying a non-default value of the number of points added per refine() call (except perhaps the first refine call in a workflow).
}
  \item{subblock_nbr}{
  NULL, or integer specifying a non-default value of the target number of iterations per refine() call (actual number of iterations may differ).
}
  \item{version}{
  A version number for \pkg{Infusion}. This is intended to allow reproducibility of results of some past versions, as return values have changed over versions.
}
\item{cumn_over_maxit}{
  logical; Whether to stop iteration when the target cumulative number of points is added to the reference table, or when the target number of iterations is first reached.
}
\item{test_fac}{
optional numeric value; for testing a workflow, it may be useful to run it with smaller reference table sizes. \code{test_fac} specifies a reduction factor for these sizes, relative to the default design.
}
}
\value{
\code{get_workflow_design} returns a list of control values, with elements \code{final_reft_size}, \code{init_reft_size}, \code{refine_blocksize}, \code{reftable_sizes}, and \code{subblock_nbr}, which can be used as shown in the Examples; as well as \code{cumn_over_maxit}, \code{first_refine_ntot}, and possibly other elements.
}
\examples{
## This shows how the get_workflow_design() may be used, 
## but in most cases one does not need to manipulate it.

\dontrun{
  blurred <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    s <- exp(s/4)
    return(c(mean=mean(s),var=var(s)))
  }
  
  ## Simulated data:
  set.seed(123)
  dSobs <- blurred(mu=4,s2=1,sample.size=40)
  
  workflow_design <- get_workflow_design(npar=2L)

  ## Construct initial reference table:
  
  # Sample its parameters:
  if (IMPLICIT <- TRUE) { # use implicit control by get_workflow_design()
    parsp_j <- init_reftable(lower=c(mu=2.5,s2=0.25,sample.size=40),
                           upper=c(mu=5.2,s2=2.4,sample.size=40))
    # => get_workflow_design() has been called internally with default values 
    #    to provide the dimension of the initial reference table.  
    #    The following syntax provides a more explicit control:
  } else {
    parsp_j <- init_reftable(lower=c(mu=2.5,s2=0.25,sample.size=40),
                             upper=c(mu=5.2,s2=2.4,sample.size=4,
                             nUnique=workflow_design$init_reft_size))
  }
  # and yet another way to the same result could be
  #
  # parsp_j <- data.frame(mu=runif(init_reft_size,min=2.5,max=5.2),
  #                       s2=runif(init_reft_size,min=0.25,max=2.4),
  #                       sample.size=40)
  
  # Generate the initial simulations:
  dsimuls <- add_reftable(, Simulate="blurred", parsTable=parsp_j, verbose=FALSE)
  
  ## Construct projections
  mufit <- project("mu",stats=c("mean","var"),data=dsimuls,verbose=FALSE)
  s2fit <- project("s2",stats=c("mean","var"),data=dsimuls,verbose=FALSE)
  dprojectors <- list(MEAN=mufit,VAR=s2fit)
  
  ## Apply projections on simulated statistics and 'data':
  dprojSimuls <- project(dsimuls,projectors=dprojectors,verbose=FALSE)
  dprojSobs <- project(dSobs,projectors=dprojectors)
  
  ## Summary-likelihood inference:
  # Initial Inference of log-likelihood surface
  slik_j <- infer_SLik_joint(dprojSimuls, stat.obs=dprojSobs, verbose=TRUE)
  # Find maximum, confidence intervals...
  slik_j <- MSL(slik_j, eval_RMSEs=FALSE, CIs=FALSE)
  
  ## Refinements over iterations
  # Here, with only two estimated parameters, workflow_design$final_reft_size
  # suggests a final reference table of 5000 simulations, attained through
  # 6 refine() calls with intermediate sizes given by 
  (workflow_design$reftable_sizes)
  # here 500, 1000, 2000, 3000, 4000, 5000.
  # 
  if (IMPLICIT) { # Again using implicit control by get_workflow_design()
    # Essentially, it suffices to call 
    for (it in seq(6)) slik_j <- 
      refine(slik_j, eval_RMSEs= it==6L, CIs= it==6L)
    # to run the default workflow. Again, the following syntax, 
    # showing how successive table sizes are controlled internally, 
    # provides a more explicit control:
  } else {
    reftable_sizes <- workflow_design$reftable_sizes
    init_reft_size <- workflow_design$init_reft_size
    refine_sizes <- diff(c(init_reft_size, reftable_sizes))
    maxit <- workflow_design$subblock_nbr
    for(it in seq_len(length(refine_sizes)-1L)) {
      add_it <- refine_sizes[it]
      slik_j <- refine(slik_j, ntot=add_it, maxit=maxit, 
                       eval_RMSEs=FALSE, CIs=FALSE)
    }
    add_it <- tail(refine_sizes,1L)
    slik_j <- refine(slik_j, ntot=add_it, maxit=maxit, 
                     CIs=add_it, eval_RMSEs=add_it)
  }
}
  
}
