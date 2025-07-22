\name{infer_SLik_joint}
\alias{infer_SLik_joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer a (summary) likelihood surface from a simulation table
}
\description{
  This infers the likelihood surface from a simulation table where each simulated data set is drawn for a distinct (vector-valued) parameter, as is usual for reference tables in other forms of simulation-based inference such as Approximate Bayesian Computation.
  A parameter density is inferred, as well as a joint density of parameters and summary statistics, and the likelihood surface is inferred from these two densities. 
}
\usage{
infer_SLik_joint(data, stat.obs, logLname = Infusion.getOption("logLname"), 
                Simulate = attr(data, "Simulate"), 
                nbCluster= seq_nbCluster(nr=nrow(data), nc=ncol(data)),
                using = Infusion.getOption("mixturing"), 
                verbose = list(most=interactive(),pedantic=FALSE,final=FALSE),
                marginalize = TRUE,
                constr_crits=NULL,
                projectors=NULL,
                is_trainset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame, whose each row contains a vector of parameters and one realization of the summary statistics for these parameters. Typically this holds the projected reference table (but see \code{projectors} argument for an experimental alternative). 
  }
  \item{stat.obs}{
  Named numeric vector of observed values of summary statistics. Typically this holds the projected values (but see \code{projectors} argument for an experimental alternative).
}
  \item{logLname}{
 The name to be given to the log Likelihood in the return object, or the root of the latter name in case of conflict with other names in this object.
}
  \item{Simulate}{Either NULL or the name of the simulation function if it can be called from the R session (see \code{\link{add_reftable}} for more information on this function). 
}
  \item{nbCluster}{controls the \code{nbCluster} argument of \code{Rmixmod::mixmodCluster}%\code{\link[Rmixmod]{mixmodCluster}}
  ; a vector of integers, or \code{"max"} which is interpreted as the maximum of the default \code{nbCluster} value.   
}
  \item{using}{Either \code{"Rmixmod"}, \code{"EMCluster"} or \code{"mclust"} to select the clustering methods used. 
}
\item{marginalize}{Boolean; whether to derive the clustering of fitted parameters by marginalization of the joint clustering; if not, a distinct call to a clustering function is performed. It is strongly advised not to change the default. This argument might be deprecated in future versions.
}
\item{constr_crits}{NULL, or quoted expression specifying a constraints on parameters, beyond the ones defined by the ranges over each parameter: see \code{\link{constr_crits}} for details.  This will control the parameter space both for maximization of the summary-likelihood, and for generation of new parameter points when \code{refine()} is called on the return object. See Examples section for a nice artificial toy example.
}
  \item{verbose}{
A list as shown by the default, or simply a vector of booleans, indicating respectively
  whether to display (1) some information about progress; (2) more information whose importance is not clear to me; (3) a final summary of the results after all elements of \code{simuls} have been processed.
}
  \item{projectors}{
  if not NULL, this argument may be passed to \code{project} in the case the \code{data} or \code{stat.obs} do not yet contain the projected statistics. This \bold{experimental} feature aims to remove the two user-level calls to \code{project} in the inference workflow.
}
  \item{is_trainset}{
  Passed to \code{project} in the case the \code{projectors} argument is used.
}
}
\value{
  An object of class \code{SLik_j}, which is a list including an \code{Rmixmod::mixmodCluster}%\code{\link[Rmixmod]{mixmodCluster}}, 
  object (or equivalent objects produced by non-default methods), and additional members not documented here. If projection was used, the list includes a  data.frame \code{reftable_raw} of cumulated unprojected simulations.
}
\seealso{ \code{\link{constr_crits}} for using \code{infer_SLik_joint} in a model 
with non-box constraints on parameters.}
\examples{
if (Infusion.getOption("example_maxtime")>2) {
  myrnorm <- function(mu,s2,sample.size=40L) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    return(c(mean=mean(s),var=var(s)))
  } # simulate means and variances of normal samples of size 'sample.size'
  set.seed(123)
  # simulated data with stands for the actual data to be analyzed:  
  Sobs <- myrnorm(mu=4,s2=1) 
  # Uniform sampling in parameter space:
  parsp <- init_reftable(lower=c(mu=2.8,s2=0.4),
                           upper=c(mu=5.2,s2=2.4))
  # Build simulation table:
  simuls <- add_reftable(Simulate="myrnorm", parsTable=parsp)
  # Infer surface:
  densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
  # Usual workflow using inferred surface:
  slik_j <- MSL(densv, eval_RMSEs=FALSE) ## find the maximum of the log-likelihood surface
  slik_j <- refine(slik_j, eval_RMSEs=FALSE)
  plot(slik_j)
  # etc:
  profile(slik_j,c(mu=4)) ## profile summary logL for given parameter value
  confint(slik_j,"mu") ## compute 1D confidence interval for given parameter
  plot1Dprof(slik_j,pars="s2",gridSteps=40) ## 1D profile
  
  # See also help for 'constr_crits'
}
}
