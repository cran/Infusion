\name{confint.SLik}
\alias{confint.SLik}
\alias{confint.SLik_j}
\alias{confint.SLikp}
\alias{confint}
\alias{allCIs}
\title{
  Compute confidence intervals by (profile) summary likelihood 
}
\description{
  \code{confint} takes an \code{SLik} object (as produced by \code{\link{MSL}}) and deduces confidence bounds for each parameter, using a (profile, if relevant) likelihood ratio method, and optionally a bootstrap method. 
  
  \code{allCIs} calls \code{confint} for all fitted parameters and re-structure the results.
}
\usage{
allCIs(object, level=0.95, verbose=TRUE, ...)
\method{confint}{SLik_j}(object, parm, level=0.95, verbose=interactive(), fixed=NULL,
                         which=c(TRUE,TRUE), nsim=0L, reset=TRUE, 
                         cluster_args=NULL, nb_cores=NULL, type="perc",...)
\method{confint}{SLik}(object, parm, level=0.95, verbose=interactive(), fixed=NULL, 
                       which=c(TRUE,TRUE), ...)
%\method{confint}{SLikp}(object, parm,
%                       level=0.95, verbose=interactive(),
%                       fixed=NULL,which=c(TRUE,TRUE),...)
}
\arguments{
  \item{object}{
an \code{SLik} or \code{SLik_j} object
}
  \item{parm}{
  The parameter which confidence bounds are to be computed
}
  \item{level}{
  The desired coverage of the interval
 }
  \item{verbose}{
  Whether to print some information or not
 }
  \item{fixed}{When this is \code{NULL} the computed interval is a profile confidence interval over all parameters excluding \code{parm}.
     \code{fixed} allows one to set fixed values to some of these parameters.  
}
  \item{which}{
  A pair of booleans, controlling whether to compute respectively the lower and the upper CI bounds.
 }
  \item{nsim}{
Integer: number of bootstrap replicates. If >1, bootstrap interval(s) are computed as further controlled by the  \code{type} argument. Note that this will be ignored if the bootstrap has previously been run and \code{reset=FALSE}.
 }
   \item{reset}{
Boolean: Whether to use any previously computed distribution (see Details) or not.
}
   \item{cluster_args, nb_cores}{
Passed to parallelization wrappers such as \code{\link[spaMM]{dopar}}.
}
   \item{type}{
Character vector, ignored if \code{nsim} = 0: bootstrap CI type(s). Possible types are \code{"norm", "basic",  "perc"} (as handled by \code{\link[boot]{boot.ci}}), and \code{"Bartlett"} (where the interval bounds are defined by threshold values of the likelihood ratio statistics modified using a Bartlett correction).
}
  \item{\dots}{
further arguments passed to or from other methods. \code{allCIs} passes them to \code{confint}, so that, e.g.,  \code{nsim} can be passed through the \ldots. \code{confint} passes them to parallelization wrappers such as \code{\link[spaMM]{dopar}}.
}
}
\details{
\code{confint.SLik_j} results are stored in the \code{object} (until the next \code{refine}), including the result of the bootstrap simulations if it was performed. This distribution may then be reused by a next call to \code{confint} for the same \code{parm} if \code{reset=FALSE}. The default is however to recompute the distribution (\code{reset=TRUE}).

Bootstrap CIs computed using \code{boot.ci} are stored as distinct elements of the return list (see Value). However, for the \code{"Bartlett"} type of CI, the \code{interval} element of the return value is modified. 

}
\value{
Both functions modify the fit \code{object} as a side effect (see Details).

\code{confint} returns a list with sublists for each parameter, each sublist containing: the bounds of the one-dimensional confidence interval (element \code{interval}, a vector); the parameter point for the lower bound (element \code{lowerpar}, a vector including all parameters fitted in the \code{SLik} object), the full parameter point for the upper bound (element \code{upperpar}, formatted as \code{lowerpar}), and optionally if a bootstrap was run, the return value of a \code{boot::boot.ci} call (element \code{bootCI}) and the simulated distribution of parameter estimates (element \code{booreps}, 1-column matrix).

\code{allCIs} returns invisibly a list with elements including \code{CIs} (itself a list of \code{confint} results), \code{bounds} (a matrix made of bound points for all parameters), and some other elements.
}
\seealso{
\code{\link{SLRT}}
}
\examples{
if (Infusion.getOption("example_maxtime")>3) {
#### Provide fit for minimal toy example: 
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

# Infer surface:
densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
# Usual workflow using inferred surface:
slik_1 <- MSL(densv) ## find the maximum of the log-likelihood surface

####  Confidence interval calculations:
(ci1 <- confint(slik_1,"mu")) # basic likelihood ratio interval
(ci2 <- confint(slik_1,"mu", nsim=199L)) # Percentile interval added
(ci3 <- confint(slik_1,"mu", nsim=199L, type="Bartlett")) # 'interval' corrected

# Previous bootstrap computations are stored in the fit object,  
# and recycled if reset=FALSE *and* nsim > 0:
(ci4 <- confint(slik_1,"mu", nsim=199L, type= "Bartlett", reset=FALSE)) # = ci3
(ci5 <- confint(slik_1,"mu", nsim=199L, type= "perc", reset=FALSE)) # = ci2
}
}
