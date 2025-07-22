cat(cli::col_yellow("test 3 parameters (primitive workflow !):\n"))

wout_refine <- 85
if (Infusion.getOption("example_maxtime")>wout_refine) { # ___F I X M E___ but in fact it takes more time (incl mixmodCluster failures)
  ## (1) function for sampling from 3-parameter gaussian mixture
  myrnorm2 <- function(mu1,mu2,s2,sample.size) {
    sam1 <- rnorm(n=sample.size,mean=mu1,sd=sqrt(s2))
    sam2 <- rnorm(n=sample.size,mean=mu2,sd=sqrt(s2))
    s <- c(sam1,sam2)
    e_mu <- mean(s)
    e_s2 <- var(s)
    return(c(mean=e_mu,var=e_s2,kurt=sum((s-e_mu)^4)/e_s2^2))
  } 
  #
  ## simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- myrnorm2(mu1=4,mu2=2,s2=1,sample.size=40) ## 
  #
  ## (2) Generate, and simulate distributions for, 
  ##        an irregular grid of parameter values, with some replicates
  
  parsp <- init_grid(lower=c(mu1=2.8,mu2=1,s2=0.2,sample.size=40),
                     upper=c(mu1=5.2,mu2=3,s2=3,sample.size=40))
  simuls <- add_simulation(NULL,Simulate="myrnorm2", parsTable=parsp)
  
  ## (3) infer logL(pars,stat.obs) for each simulated 'pars'
  # Relatively slow, hence saved as data 'densv'
  densv <- infer_logLs(simuls,stat.obs=Sobs)
  
  #
  ## (4) infer a log-likelihood surface and its maximum;
  ##       plot and extract various information. 
  slik <- infer_surface(densv)
  slik <- MSL(slik) ## find the maximum of the log-likelihood surface
  plot(slik)
  profile(slik,c(mu1=4)) ## profile summary logL for given parameter value
  confint(slik,"mu1") ## compute confidence interval for given parameter
  plot1Dprof(slik,pars="s2",gridSteps=40) ## 1D profile 
  #
  ## (5) ## refine iteratively
  if (Infusion.getOption("example_maxtime")>(wout_refine+116)) {
    slik <- refine(slik) 
  }
}
