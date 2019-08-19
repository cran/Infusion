\name{example_raw}
\alias{example_raw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Workflow for originally described method}
\description{
Example of the workflow with \code{add_simulation}), implementing the method described in the original publication 
(Rousset et al. 2017 <doi:10.1111/1755-0998.12627>).
}
\examples{
## The following example illustrates the workflow.
## However, most steps run longer than accepted by the CRAN checks,
## So by default they will not run.
##
## (1) The user must provide the function for simulation of summary statistics
myrnorm <- function(mu,s2,sample.size) {
 s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
 return(c(mean=mean(s),var=var(s)))
} # simulate means and variances of normal samples of size 'sample.size'
#
## pseudo-sample:  
set.seed(123)
Sobs <- myrnorm(mu=4,s2=1,sample.size=40) ## stands for the actual data to be analyzed
#
## (2) Generate, and simulate distributions for, 
##        an irregular grid of parameter values, with some replicates
if (Infusion.getOption("example_maxtime")>45) {
  parsp <- init_grid(lower=c(mu=2.8,s2=0.2,sample.size=40),
                     upper=c(mu=5.2,s2=3,sample.size=40))
  simuls <- add_simulation(NULL,Simulate="myrnorm",par.grid=parsp)
  
  ## (3) infer logL(pars,stat.obs) for each simulated 'pars'
  # Relatively slow, hence saved as data 'densv'
  densv <- infer_logLs(simuls,stat.obs=Sobs)
} else {
  data(densv)
  .Random.seed <- saved_seed
}
#
## (4) infer a log-likelihood surface and its maximum;
##       plot and extract various information. 
if (Infusion.getOption("example_maxtime")>17) {
 slik <- infer_surface(densv)
 slik <- MSL(slik) ## find the maximum of the log-likelihood surface
 plot(slik)
 profile(slik,c(mu=4)) ## profile summary logL for given parameter value
 confint(slik,"mu") ## compute confidence interval for given parameter
 plot1Dprof(slik,pars="s2",gridSteps=40) ## 1D profile
}
#
## (5) ## refine iteratively
if (Infusion.getOption("example_maxtime")>68) {
 slik <- refine(slik) 
}
}
