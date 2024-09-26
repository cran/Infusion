cat(cli::col_yellow("test reftable:\n"))

myrnorm <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  return(c(mean=mean(s),var=var(s)))
} # simulate means and variances of normal samples of size 'sample.size'
set.seed(123)
# simulated data with stands for the actual data to be analyzed:  
ssize <- 40L
(Sobs <- myrnorm(mu=4,s2=1,sample.size=ssize) )
## Construct initial reference table:
# Uniform sampling in parameter space:
parsp <- init_reftable(lower=c(mu=2.8, s2=0.4, sample.size=ssize), 
                         upper=c(mu=5.2, s2=2.4, sample.size=ssize))
# Build simulation table:
# set.seed(456) 
simuls <- add_reftable(Simulate="myrnorm", parsTable=parsp)

## trivial projections which should produce an y=x regression:
#mufit <- project("mu",stats=c("mean","var"),data=simuls, method="randomForest")
#s2fit <- project("s2",stats=c("mean","var"),data=simuls, method="randomForest")
mufit <- project("mu",stats=c("mean","var"),data=simuls)
s2fit <- project("s2",stats=c("mean","var"),data=simuls)

## apply projections on simulated statistics
projectors <- list("MEAN"=mufit,"VAR"=s2fit)
projectors <- list2env(projectors)
corrSobs <- project(Sobs,projectors=projectors)
corrSimuls <- project(simuls,projectors=projectors)

# Infer surface:
densv <- infer_SLik_joint(corrSimuls,stat.obs=corrSobs)
# Usual workflow using inferred surface:
slik_1 <- MSL(densv) ## find the maximum of the log-likelihood surface
#slik_j <- refine(slik_j,maxit=5, update_projectors=TRUE)

# Convenience function for plotting projections...
plot_importance(slik_1, parm="mu")


slik_j <- refine(slik_1)
# slik_j <- refine(slik_j)
plot(slik_j)
testthat::test_that(
  'Warning expected for plot_proj(slik_1, parm="mu") as projectors were updated in the input slik_1 object.', 
  {
    warn <- ""
    withCallingHandlers(
      plot_proj(slik_1, parm="mu"), warning=function(w) {
        warn <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
    testthat::expect_true(substr(warn,1,23)=="Projectors were updated")
  }
)
# No warning expected on this one:
plot_proj(slik_j)

# etc:
profile(slik_j,c(mu=4)) ## profile summary logL for given parameter value
confint(slik_j,"mu") ## compute 1D confidence interval for given parameter
plot1Dprof(slik_j,pars="s2",gridSteps=40) ## 1D profile
# goftest(slik_j, nsim = 300L) # goodness of fit test



if (requireNamespace("xLLiM", quietly=TRUE)) { # workflow with xLLiM::gllim
  densvx <- infer_SLik_joint(corrSimuls,stat.obs=corrSobs, using="xLLiM")
  # Usual workflow using inferred surface:
  slik_jx <- MSL(densvx, 
                 eval_RMSEs = FALSE # had to add this bc of hard to fix issue in RMSE computation.
                 # I need to allow browsing from errors in bootstrap replicates to debug this. (___F I X M E___)
                ) ## find the maximum of the log-likelihood surface
  #slik_j <- refine(slik_j,maxit=5, update_projectors=TRUE)
  slik_jx <- refine(slik_jx)
  SLRT(slik_jx, h0=slik_jx$MSL$MSLE+0.1, nsim = 100L) # LRT
  goftest(slik_jx, nsim = 300L) # goodness of fit test
} else warning("package 'xLLiM' not available for testing.")

if (FALSE) { # example of distinct trainsample
  trainsample <- sample(nrow(.get_reft_raw(slik_j)),1000)
  mufit <- project("mu",stats=c("mean","var"),data=.get_reft_raw(slik_j)[trainsample,])
  s2fit <- project("s2",stats=c("mean","var"),data=.get_reft_raw(slik_j)[trainsample,])
  
  ## apply projections on simulated statistics
  projectors <- list("MEAN"=mufit,"VAR"=s2fit)
  projectors <- list2env(projectors)
  corrSobs <- project(Sobs,projectors=projectors)
  corrSimuls <- project(.get_reft_raw(slik_j)[-trainsample,],projectors=projectors)
  
  
  # Infer surface:
  densv <- infer_SLik_joint(corrSimuls,stat.obs=corrSobs)
  # Usual workflow using inferred surface:
  slik_j <- MSL(densv) 
  plot(slik_j)
}


if (FALSE) { # example of reprojecting accumulated simulations
  remufit <- project("mu",stats=c("mean","var"),data=.get_reft_raw(slik_j), method="REML", train_cP_size=400, trainingsize=1000)
  res2fit <- project("s2",stats=c("mean","var"),data=.get_reft_raw(slik_j), method="REML", train_cP_size=400, trainingsize=1000)
  reprojectors <- list("MEAN"=remufit,"VAR"=res2fit)
  reprojectors <- list2env(reprojectors)
  recorrSobs <- project(Sobs,projectors=reprojectors)
  recorrSimuls <- project(.get_reft_raw(slik_j),projectors=reprojectors)
  # Infer surface:
  redensv <- infer_SLik_joint(recorrSimuls,stat.obs=recorrSobs)
  # Usual workflow using inferred suface:
  reslik_j <- MSL(redensv) ## find the maximum of the log-likelihood surface
  plot(reslik_j)
}

if (FALSE) { # 1-parameter example
  parsp1 <- data.frame(mu=4,
                      s2=runif(npoints,min=0.4,max=2.4),sample.size=ssize)
  # Build simulation table:
  simuls1 <- add_reftable(Simulate="myrnorm", parsTable=parsp1)
  s2fit1 <- project("s2",stats=c("mean","var"),data=simuls1)
  ## apply projections on simulated statistics
  projectors1 <- list("VAR"=s2fit1)
  projectors1 <- list2env(projectors1)
  corrSobs1 <- project(Sobs,projectors=projectors1)
  corrSimuls1 <- project(simuls1,projectors=projectors1)
  # Infer surface:
  densv1 <- infer_SLik_joint(corrSimuls1,stat.obs=corrSobs1)
  # Usual workflow using inferred surface:
  slik_j1 <- MSL(densv1) ## find the maximum of the log-likelihood surface
  slik_j1 <- refine(slik_j1,maxit=5, update_projectors=TRUE)
  plot(slik_j1)
  # etc:
  confint(slik_j1,"s2") ## compute 1D confidence interval for given parameter
  plot1Dprof(slik_j1,pars="s2",gridSteps=40) ## 1D profile
}
