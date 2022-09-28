\name{summLik}
\alias{summLik}
\alias{summLik.SLik_j}
\alias{summLik.default}
\title{
Model density evaluation for given data and parameters
}
\description{
Evaluation of inferred probability density as function of parameters and of (projected) summary statistics is implemented as a generic function \code{summLik}.
Given the (projected) statistics for the data used to build the \code{SLik_j} object, and the fitted parameters, this returns the (log)likelihood as the generic
\code{logLik} extractor. However, parameters can be varied (providing the likelihood function), and the data too.

This documentation deals mostly with the method for objects of class \code{SLik_j} produced by the up-to-date version of the summary-likelihood workflow.
}
\usage{
summLik(object, parm, data, ...)

# S3 method for class 'SLik_j'
\method{summLik}{SLik_j}(object, parm, data=t(attr(object$logLs,"stat.obs")), 
                         log=TRUE, which="lik", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{SLik} or \code{SLik_j} object}
  \item{parm}{data frame or matrix, containing coordinates of parameter points for which (log) likelihoods will be computed}
  \item{data}{The (projected, if relevant) summary statistics for which the likelihood of given parameters is to be computed. By default, the (projected) statistics for the data used to build the \code{SLik_j} object}
  \item{log}{Boolean: whether to return log likelihood or raw likelihood. Better ignored.}
  \item{which}{character string: \code{"lik"} for (log) likelihood inferred directly from the gaussian mixture model for joint parameters and summary statistics.
  But \code{"safe"}, which deals with a possible problem of this direct computation (see Details), is used internally by Infusion in all maximizations of likelihood.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
An  object of class \code{SLik_j} contains a simulated joint distribution of parameters and (projected) summary statistics, and a fit of a multivariate gaussian mixture model to this simulated distribution, the \dQuote{jointdens}, from which a marginal density \dQuote{margpardens} of parameters can be deduced. The raw likelihood(P;D) is the probability of the data D given the parameters P, viewed as function the parameters and for fixed data. It is inferred as jointdens(D,P)/margpardens(P) (for different P, each of jointdens and margdens are probabilities from a single (multivariate) gaussian mixture model, but this is not so for their ratio). 

When margdens(P) is low, indicating that the region of parameter space around P has been poorly sampled in the simulation step, inference of likelihood is unreliable. Spuriously high likelihood may be inferred, which results notably in poor inference based on likelihood ratios. For this reason, it is often better to use the argument  \code{which="safe"} whereby the likelihood is penalized when margdens(P) is low. The penalization is of the form\cr 
\code{penalized= unpenalized * pmin(1,margpardens/object$thr_dpar)}, where \code{thr_dpar} is a marginal density threshold stored in the \code{SLik_j} object. The source code should be consulted for details, and is subject to changes without notice.  
}
\seealso{ 
\code{\link{logLik}}.
}
\value{
Numeric vector
}
\examples{
\dontrun{
## Using 'slik_j' object from the example in help("example_reftable") 
summLik(slik_j, parm=slik_j$MSL$MSLE+0.1)

# summLik() generalizes logLik():
summLik(slik_j, parm=slik_j$MSL$MSLE) == logLik(slik_j) # must be TRUE
}
}
