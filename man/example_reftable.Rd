\name{example_reftable}
\alias{example_reftable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Workflow for method with reference table}
\description{
Example of workflow with a reference table produced by \code{add_reftable}, possibly faster in many applications than the originally described method.
}
\examples{
if ((Infusion.getOption("example_maxtime")>116)) {
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
  
  ## trivial projections that should produce an y=x regression:
  mufit <- project("mu",stats=c("mean","var"),data=simuls)
  s2fit <- project("s2",stats=c("mean","var"),data=simuls)
  
  ## additional plots for some projection method
  if (inherits(mufit,"HLfit")) mapMM(mufit,map.asp=1,
        plot.title=title(main="prediction of normal mean",xlab="mean",ylab="var"))
  if (inherits(s2fit,"HLfit")) mapMM(s2fit,map.asp=1,
        plot.title=title(main="prediction of normal var",xlab="mean",ylab="var"))
  
  ## apply projections on simulated statistics
  corrSobs <- project(Sobs,projectors=list("MEAN"=mufit,"VAR"=s2fit))
  corrSimuls <- project(simuls,projectors=list("MEAN"=mufit,"VAR"=s2fit))
  
  
  # Infer surface:
  densv <- infer_SLik_joint(corrSimuls,stat.obs=corrSobs)
  # Usual workflow using inferred surface:
  slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
  slik_j <- refine(slik_j,maxit=5)
  plot(slik_j)
  # etc:
  profile(slik_j,c(mu=4)) ## profile summary logL for given parameter value
  confint(slik_j,"mu") ## compute 1D confidence interval for given parameter
  plot1Dprof(slik_j,pars="s2",gridSteps=40) ## 1D profile
}
}
