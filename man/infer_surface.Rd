\name{infer_surface}
\alias{infer_surface}
%\alias{infer_surface.default}
\alias{infer_surface.logLs}
\alias{infer_surface.tailp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Infer a (summary) likelihood or tail probability surface from inferred likelihoods
}
\description{
  The \code{logLs} method uses a standard smoothing method (prediction under linear mixed models, a.k.a. Kriging) to infer a likelihood surface, using as input likelihood values themselves inferred with some error for different parameter values. The \code{tailp} method use a similar approach for smoothing binomial response data, using the algorithms implemented in the spaMM package for fitting GLMMs with autocorrelated random effects.
}
\usage{
\method{infer_surface}{logLs}(object, method="REML",verbose=interactive(),allFix=NULL,...)
\method{infer_surface}{tailp}(object, method="PQL",verbose=interactive(),allFix,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{ 
  A data frame with attributes, containing independent prediction of logL or of LR tail probabilities for different parameter points, as produced by \code{\link{infer_logLs}} or \code{\link{infer_tailp}}.  
 }
 \item{method}{
  methods used to estimate the smoothing parameters. If \code{method="GCV"}, a generalized cross-validation procedure is used (for \code{logLs} method only). Other methods are as described in the \code{\link[spaMM]{HLfit}} documentation.   
 }
\item{verbose}{
  Whether to display some information about progress or not.
}
\item{allFix}{
  Fixed values in the estimation of smoothing parameters. For development purposes, not for routine use. For \code{infer_surface.logLs}, this should typically include values of all parameters fitted by \code{spaMM::corrHLfit} (\eqn{\rho,\nu,\phi,\lambda}, and \code{$etaFix=}\eqn{\beta}).
}
%\item{stat.obs}{
%  Observed value of summary statistics
%}
%\item{logLname}{
%  The name to be given to the log Likelihood in the return object, or the root of the latter name in case of conflict with other names in this object.
%}
  \item{\dots}{
further arguments passed to or from other methods (currently not used).
}
}
\value{
  An object of class \code{SLik} or \code{SLikp}, which is a list including an \code{HLfit} object as returned by \code{\link{corrHLfit}}, and additional members not documented here.
}
\examples{
## see main documentation page for the package
}
