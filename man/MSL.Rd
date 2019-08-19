\name{MSL}
\alias{MSL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Maximum likelihood from an inferred likelihood surface
}
\description{
This computes the maximum of an object of class \code{SLik} representing an inferred (summary) likelihood surface}
\usage{
MSL(object, CIs = TRUE, level = 0.95, verbose = interactive(),
    eval_RMSEs = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{SLik} as produced by \code{\link{infer_surface.logLs}}
}
  \item{CIs}{
  If \code{TRUE}, construct one-dimensional confidence intervals for all parameters.
}
  \item{level}{
  Intended coverage probability of the confidence intervals.
}
  \item{verbose}{
  Whether to display some information about progress and results.
}
  \item{eval_RMSEs}{
  Logical: whether to evaluate prediction uncertainty for likelihoods/ likelihood ratios/ parameters. 
  %By default TRUE for \code{SLik} objects, and best kept so, as it is necessary for the automated iterative method. Default may be FALSE for other classes of objects, in particular \code{SLik_j}.  
}
  \item{\dots}{Further arguments passed from or to other methods. 
  }
}
\value{
  The \code{object} is returned invisibly, with added members 
  \describe{
    \item{\code{MSL}}{, itself with members \code{MSLE} and \code{maxlogL} that match the \code{par} and \code{value} returned by an \code{optim} call.}
    \item{\code{RMSEs}}{root mean square errors of the log-likelihood at its inferred maximum and of the log-likelihood ratios at the CI bounds.}
    \item{\code{par_RMSEs}}{root mean square errors of the CI bounds.}
  }
}
\details{
If Kriging has been used to construct the likelihood surface, \code{RMSEs} are computed using approximate formulas for prediction (co-)variances in linear mixed midels (see Details in \code{\link[spaMM]{predict}}). Otherwise, a more computer-intensive bootstrap method is used.   
\code{par_RMSEs} are computed from \code{RMSEs} and from the numerical gradient of profile log-likelihood at each CI bound. Only \code{RMSEs}, not \code{par_RMSEs}, are compared to \code{precision}.
}
\examples{
## see main documentation page for the package
}
