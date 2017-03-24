\name{infer_SLik_joint}
\alias{infer_SLik_joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer as (summary) likelihood surface from a simulation table
}
\description{
  This infers the likelihood surface from a simulation table where each simulated data set is drawn for a distinct (vector-valued) parameter, as is usual for reference tables in ABC.
  A parameter density is inferred, as well as a joint density of parameters and summary statistics, and the likelihood surface is inferred from these two densities.    This is not yet extensively tested, nor the code has been optimized. 
}
\usage{
infer_SLik_joint(data, stat.obs, logLname = Infusion.getOption("logLname"), 
                Simulate = attr(data, "Simulate"), 
                nbCluster= Infusion.getOption("nbCluster"),
                verbose = list(most = interactive(), final = FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame, whose rows contain a vector of parameters and one realization of the summary statistics for these parameters.
  }
  \item{stat.obs}{
  Named numeric vector of observed values of summary statistics.
}
  \item{logLname}{
 The name to be given to the log Likelihood in the return object, or the root of the latter name in case of conflict with other names in this object.
}
  \item{Simulate}{Either NULL or the name of the simulation function if it can be called from the R session. 
}
  \item{nbCluster}{\code{nbCluster} argument of \code{\link[Rmixmod]{mixmodCluster}} 
}
  \item{verbose}{
A list as shown by the default, or simply a vector of booleans, indicating respectively
  whether to display (1) some information about progress; (2) a final summary of the results after all elements of \code{simuls} have been processed.
}
}
\value{
  An object of class \code{SLik_j}, which is a list including two \code{\link[Rmixmod]{mixmodCluster}} objects, and additional members not documented here.
}
\examples{
if (Infusion.getOption("example_maxtime")>50) {
  myrnorm <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    return(c(mean=mean(s),var=var(s)))
  } # simulate means and variances of normal samples of size 'sample.size'
  set.seed(123)
  # pseudo-sample with stands for the actual data to be analyzed:  
  ssize <- 40
  Sobs <- myrnorm(mu=4,s2=1,sample.size=ssize) 
  # Uniform sampling in parameter space:
  npoints <- 600
  parsp <- data.frame(mu=runif(npoints,min=2.8,max=5.2),
                      s2=runif(npoints,min=0.4,max=2.4),sample.size=ssize)
  # Build simulation table:
  simuls <- add_reftable(Simulate="myrnorm",par.grid=parsp)
  # Infer surface:
  densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
  # Usual workflow using inferred suface:
  slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
  slik_j <- refine(slik_j,maxit=5)
  plot(slik_j)
  # etc:
  profile(slik_j,c(mu=4)) ## profile summary logL for given parameter value
  confint(slik_j,"mu") ## compute 1D confidence interval for given parameter
  plot1Dprof(slik_j,pars="s2",gridSteps=40) ## 1D profile
}
}
