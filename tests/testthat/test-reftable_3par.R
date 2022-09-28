cat(crayon::yellow("test 3 parameters, with tricky NA's:\n"))

wout_refine <- 32
if (Infusion.getOption("example_maxtime")>wout_refine) {
  ## (1) function for sampling from 3-parameter gaussian mixture
  myrnorm2 <- function(mu1,mu2,s2,sample.size) {
    sam1 <- rnorm(n=sample.size,mean=mu1,sd=sqrt(s2))
    sam2 <- rnorm(n=sample.size,mean=mu2,sd=sqrt(s2))
    s <- c(sam1,sam2)
    e_mu <- mean(s)
    e_s2 <- var(s)
    if (abs(mu1-mu2)<1 && runif(1L) <0.1) {
      return(c(mean=NA,var=NA,kurt=NA)) 
    } else return(c(mean=e_mu,var=e_s2,kurt=sum((s-e_mu)^4)/e_s2^2))
  } 
  #
  ## pseudo-sample, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- myrnorm2(mu1=4,mu2=2,s2=1,sample.size=40) ## 
  #
  npoints <- 600
  parsp <- data.frame(mu1=runif(npoints,min=2.8,max=5.2),
                      mu2=runif(npoints,min=1,max=3),
                      s2=runif(npoints,min=0.2,max=3),sample.size=40)
  simuls <- add_reftable(Simulate="myrnorm2",par.grid=parsp)
  
  # verif that projection works with missing data
  abyss <- (project("mu1", data=simuls, stats=c("mean","var","kurt")))
  # (but the projection is not used below)
  

  densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
  # Usual workflow using inferred surface:
  slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
  #slik_j <- refine(slik_j,maxit=5, update_projectors=TRUE)
  slik_j <- refine(slik_j,maxit=2,update_projectors=TRUE)
  plot(slik_j)
  # etc:
  profile(slik_j,c(mu1=4)) ## profile summary logL for given parameter value
  confint(slik_j,"mu1") ## compute 1D confidence interval for given parameter
  plot1Dprof(slik_j,pars="s2",gridSteps=40) ## 1D profile
  plot2Dprof(slik_j)
  
}
