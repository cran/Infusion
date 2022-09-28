\name{predict.SLik_j}
\alias{predict.SLik_j}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate log-likelihood for given parameters
}
\description{
As the Title says. Implemented as a method of the \code{predict} generic, for objects created by the up-to-date workflow using gaussian mixture modelling of a joint distribution of parameters and statistics (hence the \code{newdata} argument, shared by many \code{predict} methods; but these \code{newdata} should be parameter values, not data).
}
\usage{
\method{predict}{SLik_j}(
  object, newdata, log = TRUE, which = "lik", 
  tstat = t(attr(object$logLs, "stat.obs")), 
  solve_t_chol_sigma_lists = object$solve_t_chol_sigma_lists, 
  ...)
}
\arguments{
  \item{object}{
an object of class \code{SLik_j}, as produced by \code{\link{infer_SLik_joint}}.
}
  \item{newdata}{
A matrix, whose rows each contain a full vector of the fitted parameters; or a single vector. If parameter names are not provided (as column names in the matrix case), then the vector is assumed to be ordered as \code{object$colTypes$fittedPars}.
}
  \item{log}{
Boolean: whether to return log-likelihood or likelihood.
}
  \item{which}{
\code{"lik"} or \code{"safe"}. The latter protects against some artefacts of predictions beyond the regions of parameter space well sampled by the inference procedure. 
}
  \item{tstat}{
The data (as projected summary statistics). Defaults to the data input in the inference procedure (i.e., the projected statistics used as \code{stat.obs} argument of \code{infer_SLik_joint}). 
}
  \item{solve_t_chol_sigma_lists}{For programming purposes. Do not change this argument.}
  \item{\dots}{
For consistency with the generic. Currently ignored.
}
}
\value{
Numeric: a single value, or a vector of (log-)likelihoods for different rows of the input \code{newdata}.
}
\examples{
## see help("example_reftable")
}
