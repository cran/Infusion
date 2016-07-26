\name{infer_SLik_joint}
\alias{infer_SLik_joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer as (summary) likelihood surface from simulation table
}
\description{
  This infers the likelihood surface from a simulation table where each realization of the summary statistic is simulated under a distinct parameter vector (as is usual for reference tables in ABC).
  A parameter density is inferred, as well as a joint density of parameters and summary statistics, and the likelihood surface is inferred from these two densities.
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
  A data frame, which rows contain a vector of parameters and one realization of the summary statistics for these parameters.
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
myrnorm <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  return(c(mean=mean(s),var=var(s)))
} # simulate means and variances of normal samples of size 'sample.size'
set.seed(123)
# pseudo-sample with stands for the actual data to be analyzed:  
Sobs <- myrnorm(mu=4,s2=1,sample.size=20) 
# Uniform sampling in parameter space:
npoints <- 3000
parsp <- data.frame(mu=runif(npoints,min=2.8,max=5.2),
                    s2=runif(npoints,min=0.4,max=2.4),sample.size=20)
# Build simulation table:
simuls <- add_reftable(Simulate="myrnorm",par.grid=parsp)
# Infer surface (nbCluster argument for fast demo only):
densv <- infer_SLik_joint(simuls,stat.obs=Sobs,nbCluster=2)
# Usual workflow using inferred suface:
slik <- MSL(densv) ## find the maximum of the log-likelihood surface
plot(slik)
# More workflow and other auxiliary functions (not run):
#slik <- refine(slik) ## refine iteratively
#profile(slik,c(mu=4)) ## profile summary logL for given parameter value
#confint(slik,"mu") ## compute confidence interval for given parameter
#plot1Dprof(slik,pars="s2",gridSteps=40) ## 1D profile
}
