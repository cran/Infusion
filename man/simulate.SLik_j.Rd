\name{simulate.SLik_j}
\alias{simulate}
\alias{simulate.SLik_j}
\title{
Simulate method for an \code{SLik_j} object.
}
\description{
\code{\link[stats]{simulate}} method for \code{SLik_j} objects, by default simulating realizations of the vector of projected summary statistics, drawn from their inferred distribution, given the summary-ML estimates which are the default value of the \code{given} argument. 

For any non-default \code{given} argument, the sampling distribution is still deduced from the multivariate Gaussian mixture fit of the reference table, by conditioning it on \code{given} values. Any variable included in the mixture model may be included in \code{given}, allowing to simulate from other distributions than that of the vector of projected summary statistics.  

This usage should not be confused with simulating the sample-generating process, necessarily distinctly available to the user, and which does not rely on the mixture model stored in the fit object.
Simulations of the sample-generating process for \code{given} parameter values can be obtained by setting non-default option \code{SGP=TRUE}.  
}
\usage{
\method{simulate}{SLik_j}(object, nsim = 1, seed = NULL, given=object$MSL$MSLE, 
                          norm_or_t=.wrap_rmvnorm, SGP=FALSE, ...)
}
\arguments{
  \item{object}{
An object of class \code{SLik_j} as produced by the up-to-date workflow. 
}
  \item{nsim}{
number of response vectors of projected summary statistics to simulate. 
}
  \item{seed}{
Seed for the random number generator (RNG). Here this controls the \code{.Random.seed} in the global environment, as in \code{simulate.lm}. This means that if a non-NULL \code{seed} is specified, it controls the RNG during the \code{simulate} call, but the RNG is reset to its prior state on exit.  
}
  \item{given}{
The default is the summary-MLE, a full vector of fitted parameters; but Any variable included in the mixture fit of the referencetable may be included (see Description).  
}
  \item{norm_or_t}{
Controls the sampler in in cluster of the mixture. The default value is a trivial wrapper around the \code{\link[mvtnorm]{rmvnorm}} sampler (consistently with the fitted model), but this argument makes it possible to specify other samplers (e.g., \code{norm_or_t=Infusion:::.wrap_rmvt} to sample from \code{\link[mvtnorm]{rmvt}(., df=1)}; or used-defined samplers with the same interface).   
}
\item{SGP}{Boolean. Whether to sample from the sample-generating process.}
\item{\ldots}{
Additional arguments. Currently ignored, except when \code{SGP=TRUE}, in which case e.g. \code{control.Simulate} can be passed through the dots to control the sample simulator.
}
}
\value{
By default (\code{SGP=FALSE}), a matrix of size \code{nsim} times the number of \bold{projected} summary statistics; if \code{SGP=TRUE}, a data frame with columns for parameters, for \bold{raw} summary statistics, and optionally for latent variables if relevant.
}
\examples{
## Assuming an object 'slik_j' of class 'SLik_j':
# simulate(slik_j, nsim=3)
}