if (FALSE) {
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
  parsp <- init_grid(lower=c(mu=2.8,s2=0.2,sample.size=40),
                     upper=c(mu=5.2,s2=3,sample.size=40))
  simuls <- add_simulation(NULL,Simulate="myrnorm",par.grid=parsp)
  
  ## (3) infer logL(pars,stat.obs) for each simulated 'pars'
  # Relatively slow, hence saved as data 'densv'
  densv <- infer_logLs(simuls,stat.obs=Sobs)
  saved_seed <- .Random.seed
  save(densv,saved_seed,file="./package/data/densv.RData")
}
  
if (FALSE) {  
  ## Transform normal random deviates rnorm(,mu,sd)
  ## so that the mean of transformed sample is not sufficient for mu,
  ## and that variance of transformed sample is not sufficient for sd,
  blurred <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    s <- exp(s/4)
    return(c(mean=mean(s),var=var(s)))
  }
  
  set.seed(123)
  dSobs <- blurred(mu=4,s2=1,sample.size=20) ## stands for the actual data to be analyzed
  
  ## Sampling design as in canonical example 
  parsp <- init_grid(lower=c(mu=2.8,s2=0.4,sample.size=20),
                     upper=c(mu=5.2,s2=2.4,sample.size=20))
  # simulate distributions
  dsimuls <- add_simulation(,Simulate="blurred", par.grid=parsp) 
  
  ## Use projection to construct better summary statistics for each each parameter 
  mufit <- project("mu",stats=c("mean","var"),data=dsimuls)
  s2fit <- project("s2",stats=c("mean","var"),data=dsimuls)
  
  ## plots
  mapMM(mufit,map.asp=1,
        plot.title=title(main="prediction of normal mean",xlab="exp mean",ylab="exp var"))
  mapMM(s2fit,map.asp=1,
        plot.title=title(main="prediction of normal var",xlab="exp mean",ylab="exp var"))
  
  ## apply projections on simulated statistics
  corrSobs <- project(dSobs,projectors=list("MEAN"=mufit,"VAR"=s2fit))
  corrSimuls <- project(dsimuls,projectors=list("MEAN"=mufit,"VAR"=s2fit))
  
  ## Analyze 'projected' data as any data (cf canonical example)
  densb <- infer_logLs(corrSimuls,stat.obs=corrSobs) 
  saved_seed <- .Random.seed
  save(densb,saved_seed,file="./package/data/densb.RData")
  
}