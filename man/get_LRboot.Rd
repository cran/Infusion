\name{get_LRboot}
\alias{get_LRboot}
\alias{SLRT}
\title{
Summary likelihood ratio tests
}
\description{
\code{get_LRboot} provides a fast approximation to bootstrap distribution of likelihood ratio statistic.
The bootstrap distribution of the likelihood ratio (LR) statistic may be used to correct the tests based on its asymptotic chi-square distribution. However, the standard bootstrap involves resimulating the data-generating process, given the ML estimates on the original data. This function implements a fast approximation avoiding such simulation, instead drawing from the inferred distribution of (projected, if relevant) summary statistics, again given the maximum (summary-)likelihood estimates.   

\code{SLRT} computes likelihood ratio tests based on the summary-likelihood surface and optionally on \code{get_LRboot} results. Several correction of the basic likelihood ratio test may be reported, some more speculative than others.
}
\usage{
SLRT(object, h0, nsim=0L, BGP=NULL, ...)
get_LRboot(object, h0_pars = NULL, nsim = 100L, reset = TRUE, BGP=object$MSL$MSLE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik_j} object.
}
  \item{h0}{
Numeric named vector of tested parameter values.    
}
  \item{h0_pars}{
either \code{NULL} (the default), to approximate the distribution of the LR statistic for the full vector of estimated parameters; or a vector of names of a subset of this vector, to approximate the distribution of the profile LR statistic for this subset.    
}
  \item{nsim}{
Integer: number of bootstrap replicates. Values lower than the default are not recommended. Note that this will be ignored if the distribution has previously been simulated and \code{reset=FALSE}.
}
  \item{reset}{
Boolean: Whether to use any previously computed distribution (see Details) or not.
}
\item{BGP}{Named numeric vector of \dQuote{Bootstrap-Generating Parameters}. Ideally the distribution of the LR test statistic would be pivotal and thus the parameter values under which this distribution is simulated would not matter. In practice, simulating by default its distribution under the \dQuote{best} available information (the MSLE for \code{get_LRboot}, or the specifically tested hypothesis defined by the \code{h0} argument of \code{SLRT}) may be more accurate than under alternative parametric values. For \code{h0} being an incomplete parameter vector and \code{BGP} is NULL (the default), \code{SLRT} will simulate under a completed parameter vector using estimates of other parameters maximizing the likelihood profile for \code{h0}. 
}
  \item{\ldots}{
  For \code{SLRT}: further arguments passed to \code{get_LRboot}. For \code{get_LRboot}: further arguments controlling parallelization, including \code{nb_cores}. However, parallelization may be best ignored in most cases (see Details).
}
}
\details{
The result of calling \code{get_LRboot} (either directly or through \code{SLRT}) with given \code{h0_pars} is stored in the \code{object} (until the next \code{refine}), and this saved result is returned by a next call to \code{get_LRboot} with the same \code{h0_pars} if \code{reset=FALSE}. The default is however to recompute the distribution (\code{reset=TRUE}).

Parallelization is possible but maybe not useful because computations for each bootstrap replicate are fast relative to parallelization overhead. It will be called when the \ldots arguments include an \code{nb_cores}>1. The \ldots may include further arguments passed to \code{\link[spaMM]{dopar}}, but among the \code{dopar} arguments, \code{iseed} will be ignored, and \code{fit_env} should not be used.

  A raw bootstrap p-value can be computed from the simulated distribution as \code{(1+sum(t >= t0))/(N+1)} where \code{t0} is the original likelihood ratio, \code{t} the vector of bootstrap replicates and \code{N} its length. See Davison & Hinkley (1997, p. 141) for discussion of the adjustments in this formula. However, a sometimes more economical use of the bootstrap is to provide a Bartlett correction for the likelihood ratio test in small samples. According to this correction, the mean value \eqn{m} of the likelihood ratio statistic under the null hypothesis is computed (here estimated by simulation) and the original LR statistic is multiplied by \eqn{n/m} where \eqn{n} is the number of degrees of freedom of the test. Unfortunately, the underlying assumption that the corrected LR statistic follows the chi-square distribution does not always work well. 
}
\value{
\code{get_LRboot} returns a numeric vector representing the simulated distribution of the LR statistic, i.e. \bold{twice} the log-likelihood difference, as directly used in \code{pchisq()} to get the p-value.

\code{SLRT} returns a list with the following element(s), each being a one-row data frame:
\item{basicLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the p-value;}

and, if a bootstrap was performed: 
\item{BartBootLRT}{A data frame including values of the Bartlett-corrected likelihood ratio chi2 statistic, its degrees of freedom, and its p-value;}
\item{rawBootLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the raw bootstrap p-value;}
\item{safeBartBootLRT}{equal to \code{rawBootLRT} if the mean bootstrap value of the LR statistic is lower than the number of degrees of freedom, and to \code{BartBootLRT} otherwise.}
}
\references{
 Bartlett, M. S. (1937) Properties of sufficiency and statistical tests. Proceedings of the Royal Society (London) A 160: 268-282.

 Davison A.C., Hinkley D.V. (1997) Bootstrap methods and their applications. Cambridge Univ. Press, Cambridge, UK.
}
\examples{ 
## See help("example_reftable") for SLRT() examples;
## continuing from there, after refine() steps for good results:
# set.seed(123);mean(get_LRboot(slik_j, nsim=500, reset=TRUE)) # close to df=2 
# mean(get_LRboot(slik_j, h0_pars = "s2", nsim=500, reset=TRUE)) # close to df=1 

\dontrun{
### Simulation study of performance of the corrected LRTs: % bootstrap took 1314.42 s on pre v2.0

## Same toy example as in help("example_reftable"):
blurred <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    s <- exp(s/4)
    return(c(mean=mean(s),var=var(s)))
  }

## First build a largish reference table and projections to be used in all replicates
# Only the 600 first rows will be used as initial reference table for each "data"
#
set.seed(123)
#
parsp_j <- data.frame(mu=runif(6000L,min=2.8,max=5.2),
                      s2=runif(6000L,min=0.4,max=2.4),sample.size=40)
dsimuls <- add_reftable(,Simulate="blurred",par.grid=parsp_j,verbose=FALSE)
#
mufit <- project("mu",stats=c("mean","var"),data=dsimuls,verbose=TRUE)
s2fit <- project("s2",stats=c("mean","var"),data=dsimuls,verbose=TRUE)
dprojectors <- list(MEAN=mufit,VAR=s2fit)
dprojSimuls <- project(dsimuls,projectors=dprojectors,verbose=FALSE)

## Function for single-data analysis:
#
foo <- function(y, refine_maxit=0L, verbose=FALSE) {
  dSobs <- blurred(mu=4,s2=1,sample.size=40) 
  ## ----Inference workflow-----------------------------------------------
  dprojSobs <- project(dSobs,projectors=dprojectors)
  dslik <- infer_SLik_joint(dprojSimuls[1:600,],stat.obs=dprojSobs,verbose=FALSE)
  dslik <- MSL(dslik, verbose=verbose, eval_RMSEs=FALSE)
  if (refine_maxit) dslik <- refine(dslik, maxit=refine_maxit)
  ## ---- LRT-----------------------------------------------
  lrt <- SLRT(dslik, h0=c(s2=1), nsim=200)
  c(basic=lrt$basicLRT$p_value,raw=lrt$rawBootLRT$p_value,
    bart=lrt$BartBootLRT$p_value,safe=lrt$safeBartBootLRT$p_value)
}

## Simulations using convenient parallelization interface:
#
# library(doSNOW) # optional
#
bootreps <- spaMM::dopar(matrix(1,ncol=200,nrow=1),              # 200 replicates of foo()
  fn=foo, fit_env=list(blurred=blurred, dprojectors=dprojectors, dprojSimuls=dprojSimuls), 
  control=list(.errorhandling = "pass", .packages = "Infusion"),
  refine_maxit=5L,
  nb_cores=parallel::detectCores()-1L, iseed=123)
#
plot(ecdf(bootreps["basic",]))
abline(0,1)
plot(ecdf(bootreps["bart",]), add=TRUE, col="blue")
plot(ecdf(bootreps["safe",]), add=TRUE, col="red")
plot(ecdf(bootreps["raw",]), add=TRUE, col="green") 
#
# Note that refine() iterations are important for good performance.
# Without them, even a larger reftable of 60000 lines 
# may exhibit poor results for some of the CI types.
}
}
