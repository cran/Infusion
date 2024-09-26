\name{example_reftable}
\alias{example_reftable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Workflow for method with reference table}
\description{
Examples of workflow with a reference table produced by \code{add_reftable}, to be preferred the primitive method described in the first publication of \pkg{Infusion}.
}
\examples{
if (Infusion.getOption("example_maxtime")>56) {

  ## Normal(mu,sd) model, with inefficient raw summary statistics:
  ## To illustrate that case we transform normal random deviates rnorm(,mu,sd)
  ## so that the mean of transformed sample is not sufficient for mu,
  ## and the variance of transformed sample is not sufficient for sd.
  blurred <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    s <- exp(s/4)
    return(c(mean=mean(s),var=var(s)))
  }
  
  ## simulated data which stands for the actual data to be analyzed:  
  set.seed(123)
  dSobs <- blurred(mu=4,s2=1,sample.size=40)
  
  ## Construct initial reference table:
  parsp_j <- init_reftable(lower=c(mu=2.5, s2=0.25, sample.size=40), 
                           upper=c(mu=5.2, s2=2.4, sample.size=40))
  dsimuls <- add_reftable(,Simulate="blurred", parsTable=parsp_j,verbose=FALSE)

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
  plot_proj(slik_j)
  plot_importance(slik_j, parm="mu")
  
  # ... and for computing likelihoods for new parameters and/or data:
  summLik(slik_j, parm=slik_j$MSL$MSLE+0.1)
  
  ## refine estimates iteratively
  maxrefines <- 2L
  # See get_workflow_design() for a suggested number of refine() calls, 
  # typically more than the 2 refine calls shown in this small example.
  for(it in 1:maxrefines) slik_j <- 
    refine(slik_j, eval_RMSEs=it==maxrefines, CIs=it==maxrefines)
  ##

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
    
    # Low-level predict() method (rarely directly used, but documented)
    predict(slik_j, newdata = slik_j$MSL$MSLE)   # the 'data' are here parameters!

    # 'ranger' projections can take a lot of memory. One can reduce them by... 
    # deforest_projectors(slik_j)
    # ...before...
    # save(slik_j, file="slik_j")
    # They will be rebuilt on the fly if needed for further iterations.
  }
}
}
