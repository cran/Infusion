\name{extractors}
\alias{extractors}
\alias{summary}
\alias{print}
\alias{logLik}
\alias{summary.logLs}
\alias{summary.SLik}
\alias{summary.SLik_j}
\alias{summary.SLikp}
\alias{print.SLik}
\alias{print.SLik_j}
\alias{print.logLs}
\alias{print.SLikp}
\alias{logLik.SLik}
\alias{logLik.SLik_j}
\title{
Summary, print and logLik methods for Infusion results.
}
\description{
\code{summary} prints information about the fit.
\code{print} is an alias for \code{summary}.
\code{logLik} extracts the log-likelihood (exact or approximated). 
}
\usage{
\method{summary}{SLik}(object, ...)
\method{print}{SLik}(x, ...)
\method{logLik}{SLik}(object, ...)
# and identical usage for 'SLik_j' objects
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{
An object of class \code{SLik} or  \code{SLik_j};
}
  \item{\dots}{  further arguments passed to or from other methods (currently without any specific effect). }
}
\value{
\code{logLik} returns the inferred likelihood maximum, with attribute \code{RMSE} giving its root means square error of estimation. 
\code{summary} and \code{summary} return the object invisibly. They print details of the fits in a convenient form.
}
\seealso{
See \code{\link{get_from}} for a more general interface for extracting elements from Infusion results, and \code{\link{summLik}} for using a fit object to evaluate the likelihood function for distinct parameter values and even distinct data.
}
\note{See workflow example in \code{\link{example_reftable}}.}
\examples{
# See Note
}
