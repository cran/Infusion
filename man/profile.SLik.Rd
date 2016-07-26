\name{profile.SLik}
\alias{profile.SLik}
\alias{profile.SLik_j}
\alias{profile}
\title{
  Compute profile summary likelihood 
}
\description{
  Predicts the profile likelihood for a given parameter value (or vector of such values) using predictions from an \code{SLik} object (as produced by \code{\link{MSL}}).
}
\usage{
\method{profile}{SLik}(fitted, value, fixed=NULL, return.optim=FALSE, ...)
\method{profile}{SLik_j}(fitted, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted}{
an \code{SLik} object.
}
  \item{value}{
  The parameter value (as a vector of named values) for which the profile is to be computed
}
  \item{fixed}{When this is \code{NULL} the computed interval is a profile confidence interval over all parameters excluding \code{value}.
     \code{fixed} allows one to set fixed values to some of these parameters.  
}
\item{return.optim}{If this is TRUE, and if maximization of likelihood given \code{value} and \code{fixed} is indeed required, then the full result of the \code{optim} call is returned.}
  \item{\dots}{
  For \code{SLik_j} method, arguments passed to \code{SLik} method.
  For \code{SLik_j} method, currently not used.
}
}
\value{
The predicted summary profile log-likelihood; or possibly the result of an \code{optim} call if \code{return.optim} is TRUE. 
}
\examples{
## see main documentation page for the package
}
