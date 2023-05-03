\name{example_reftable}
\alias{example_reftable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Workflow for method with reference table}
\description{
Examples of workflow with a reference table produced by \code{add_reftable}, possibly faster in many applications than the originally described method.
}
\examples{
if (Infusion.getOption("example_maxtime")>46) {

  ## Normal(mu,sd) model, with inefficient raw summary statistics:
  ## To illustrate that case we transform normal random deviates rnorm(,mu,sd)
  ## so that the mean of transformed sample is not sufficient for mu,
  ## and the variance of transformed sample is not sufficient for sd.
  blurred <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    s <- exp(s/4)
    return(c(mean=mean(s),var=var(s)))
  }
  
  ## pseudo-sample which stands for the actual data to be analyzed:  
  set.seed(123)
  dSobs <- blurred(mu=4,s2=1,sample.size=40)
  
  ## Construct reference table:
  parsp_j <- data.frame(mu=runif(600L,min=2.8,max=5.2),
                        s2=runif(600L,min=0.4,max=2.4),sample.size=40)
  dsimuls <- add_reftable(,Simulate="blurred",par.grid=parsp_j,verbose=FALSE)
  
  #- When no 'Simulate' function is provided, 
  #- but only a data.frame 'toydf' of simulations,
  #- a formal reference table can be produced by  
  # dsimuls <- structure(toydf, LOWER=c(mu=2,s2=0,sample.size=40))
  # dsimuls <- add_reftable(dsimuls)
  #- where the 'LOWER' attribute tells 
  #- the parameters apart from the summary statistics.
  
  ## Construct projections
  mufit <- project("mu",stats=c("mean","var"),data=dsimuls,verbose=FALSE)
  s2fit <- project("s2",stats=c("mean","var"),data=dsimuls,verbose=FALSE)
  dprojectors <- list(MEAN=mufit,VAR=s2fit)
  
  ## Apply projections on simulated statistics and 'data':
  dprojSimuls <- project(dsimuls,projectors=dprojectors,verbose=FALSE)
  dprojSobs <- project(dSobs,projectors=dprojectors)
  
  ## Summary-likelihood inference:
  # Infer log-likelihood surface
  slik_j <- infer_SLik_joint(dprojSimuls,stat.obs=dprojSobs,verbose=TRUE)
  # Find maximum, confidence intervals...
  slik_j <- MSL(slik_j)
  
  # Convenience function for plotting projections...
  plot_proj(slik_j, parm="mu", proj="MEAN")
  
  # ... and for computing likelihoods for new parameters and/or data:
  summLik(slik_j, parm=slik_j$MSL$MSLE+0.1)
  
  ## refine estimates iteratively
  slik_j <- refine(slik_j,maxit=5, update_projectors=TRUE)
  
  if (Infusion.getOption("example_maxtime")>99) { # Post-fit procedures,
    #                                   all with distinct documentation: 
  
    plot(slik_j)
    profile(slik_j,c(mu=4)) ## profile summary logL for given parameter value
    confint(slik_j,"mu") ## compute 1D confidence interval for given parameter
    plot1Dprof(slik_j,pars="s2",gridSteps=40) ## 1D profile
    summary(slik_j) # or print()
    logLik(slik_j)
    
    SLRT(slik_j, h0=slik_j$MSL$MSLE+0.1, nsim = 100L) # LRT
    SLRT(slik_j, h0=slik_j$MSL$MSLE[1]+0.1, nsim = 100L) # profile LRT
    
    goftest(slik_j) # goodness of fit test
    
    # Low-level predict() method (rarely directly used, otherwise see its documentation!)
    predict(slik_j, newdata = slik_j$MSL$MSLE)         # the 'data' are here parameters!

    
  }
}
}
