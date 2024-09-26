cat(cli::col_yellow("test likelihood with narrow peak among local maxima:\n"))

# Example from Wilkinson's ABC tutorial
# Where Infusion.options(mixmodGaussianModel) matters... 
nips <- function(theta) {
  s <- rnorm(n=1,mean=2*theta*(theta+2)*(theta-2),sd=sqrt(0.1+theta^2))
  return(c(D=s))
}

lik <- function(theta) dnorm(x=2,mean=2*theta*(theta+2)*(theta-2),sd=sqrt(0.1+theta^2),log=FALSE)
liks <- sapply(seqx <- seq(-4,4,0.001),lik)


set.seed(123)
Dobs <- c(D=2) ## stands for the actual data to be analyzed
# plot(seqx, liks, xlab="x", ylab="likelihood(D=2)", col="red", type="l")

npoints <- get_workflow_design(npar=1)$init_reft_size
parsp <- data.frame(theta=runif(npoints,min=-4,max=4))
simuls <- add_reftable(Simulate=nips, parsTable=parsp)
#
# library(abcrf)
# rfabc <- regAbcrf(theta~D, simuls)
# densityPlot(rfabc, data.frame(D=2), simuls)
# predict(rfabc,obs=data.frame(D=2), training=simuls)
#
densv <- infer_SLik_joint(simuls,stat.obs=Dobs)
currMSL <- MSL(densv) ## find the maximum of the log-likelihood surface
plot(currMSL); lines(seqx,liks/max(liks),col="red") # add true likelihood
currMSL <- refine(currMSL) ; lines(seqx,liks/max(liks),col="red")
if (inherits(currMSL$jointdens, "dMixmod")) plot(currMSL$jointdens, data=currMSL$logLs) 

if (FALSE) { # Instructive example of poor effect of (RF) projections
  bigparsp <- data.frame(theta=runif(3000,min=-4,max=4))
  bigsimuls <- add_reftable(Simulate="nips", parsTable=bigparsp)
  prth <- project("theta", stats="D",data=bigsimuls)
  prDobs <- project(Dobs,projectors=list(PRTH=prth))
  prsimuls <- project(bigsimuls,projectors=list(PRTH=prth))
  densv <- infer_SLik_joint(prsimuls,stat.obs=prDobs)
  slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
  plot(slik_j)
}
if (FALSE) { # More on example of poor effect of projections: standard workflow does not help much.
  parsp <- data.frame(theta=runif(npoints,min=-4,max=4))
  simuls <- add_reftable(Simulate="nips", parsTable=parsp)
  prth <- project("theta", stats="D",data=simuls)
  prDobs <- project(Dobs,projectors=list(PRTH=prth))
  prsimuls <- project(simuls,projectors=list(PRTH=prth))
  densv <- infer_SLik_joint(prsimuls,stat.obs=prDobs)
  slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
  plot(slik_j)
  slik_j <- refine(slik_j)
  slik_j <- refine(slik_j)
  slik_j <- refine(slik_j)
}
