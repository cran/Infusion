\name{summLik}
\alias{summLik}
\alias{summLik.SLik_j}
\alias{summLik.default}
\title{
Model density evaluation for given data and parameters
}
\description{
Evaluation of inferred probability density as function of parameters and of (projected) summary statistics is implemented as a generic function \code{summLik}.
This documentation deals mostly with its method for objects of class \code{SLik_j} produced by the up-to-date version of the summary-likelihood workflow.

Given the (projected) statistics for the data used to build the \code{SLik_j} object, and the fitted parameters, this returns the (log)likelihood, as the generic
\code{logLik} extractor does. However, parameters can be varied (so that \code{summLik} provides the likelihood function rather than simply its maximum), the data can be varied too, and likelihood profiles (or even full new estimates) are computed when an incomplete parameter vector (or even NULL) is specified. 
}
\usage{
summLik(object, parm, data, ...)

\method{summLik}{SLik_j}(object, parm, data=t(get_from(object,"proj_data")), 
                         log=TRUE, which="safe", constr_tuning = Inf, 
                         newMSL=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{SLik} or \code{SLik_j} object}
  \item{parm}{Vector, data frame or matrix, containing coordinates of parameter points for which (log) likelihoods will be computed; or NULL. A profile will be computed if a single incomplete parameter vector is provided. A full new estimate will be computed in \code{parm} is NULL.}
  \item{data}{Matrix of (projected, if relevant) summary statistics for which the likelihood of given parameters is to be computed. By default, the (projected) statistics for the data used to build the \code{SLik_j} object}
  \item{log}{Boolean: whether to return log likelihood or raw likelihood. Better ignored.}
  \item{which}{character string: \code{"lik"} for (log) likelihood deduced from the multivariate gaussian mixture model for joint parameters and summary statistics, without further modifications. But the default, \code{"safe"}, may correct this result to deal with possible extrapolation artefacts (see Details of \code{\link{predict.SLik_j}}).
  }
  \item{constr_tuning}{Passed to \code{\link{predict.SLik_j}}.}
  \item{newMSL}{
  Boolean. If this is TRUE and a profile was computed, attributes are added to the result (see Value).
  }
  \item{\dots}{further arguments passed to or from other methods. Currently only passed to \code{predict.SLik_j} when no likelihood profile is computed.}
}
\seealso{ 
\code{\link{predict.SLik_j}}) for case without profiling; 
\code{\link{logLik}}, the standard extractor of likelihood for the model fitted to the original data.
}
\value{
Numeric vector, with optional attribute(s).

If no profile is computed, it may have attributes from the return value of \code{\link{predict.SLik_j}}. If a profile is computed, the returned value has attribute \code{"profpt"} giving the profile-maximizing parameter vector. Further, if \code{newMSL=TRUE}, the following attributes are added: \code{"newobs_MSL"}, a list with information about unconstrained summary-likelihood maximization (useful mainly when there are new \code{data}); and \code{"LRstat"}, the resulting log-likelihood ratio.  
}
\examples{
\dontrun{
## Using 'slik_j' object from the example in help("example_reftable") 
summLik(slik_j, parm=slik_j$MSL$MSLE+0.1)

# summLik() generalizes logLik():
summLik(slik_j, parm=slik_j$MSL$MSLE) == logLik(slik_j) # must be TRUE
}
}
