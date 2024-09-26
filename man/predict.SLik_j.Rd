\name{predict.SLik_j}
\alias{predict.SLik_j}
\title{
Evaluate log-likelihood for given parameters
}
\description{
As the Title says. As likelihood is obtained as a prediction from a statistical model for the likelihood surface, this has been implemented as a method of the \code{predict} generic, for objects created by the up-to-date workflow using gaussian mixture modelling of a joint distribution of parameters and statistics. Hence, it has a \code{newdata} argument, as shared by many \code{predict} methods; but these \code{newdata} should be parameter values, not data). You can use the alternative \code{\link{summLik}} extractor if you do not like this syntax.
}
\usage{
\method{predict}{SLik_j}(
  object, newdata, log = TRUE, which = "safe", 
  tstat = t(get_from(object,"stat.obs")), 
  solve_t_chol_sigma_lists = object$clu_params$solve_t_chol_sigma_lists, 
  constr_tuning= FALSE,
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
\code{"safe"} or \code{"lik"}. The default protects against some artefacts of extrapolation: see Details. 
}
  \item{tstat}{
The data (as projected summary statistics). Defaults to the data input in the inference procedure (i.e., the projected statistics used as \code{stat.obs} argument of \code{infer_SLik_joint}). 
}
  \item{solve_t_chol_sigma_lists}{For programming purposes. Do not change this argument.}
  \item{constr_tuning}{positive number, or FALSE: controls the effect of parameter constraints specified by the \code{constr_crits} argument of \code{infer_SLik_joint} on the evaluation of summary log-likelihood. When it is \code{FALSE} (or 0), no penalty is applied; when this is \code{Inf}, the log-likelihood of parameters violating constraints will be \code{-Inf}. Intermediate values allow an intermediate penalization (the source code should be consulted for further details), but their use is not recommended. 
}
  \item{\dots}{
For consistency with the generic. Currently ignored.
}
}
\value{
Numeric: a single value, or a vector of (log-)likelihoods for different rows of the input \code{newdata}. When \code{which="safe"}, a \code{"lowdens"} attribute will be added when at least on parameter points had a low \dQuote{pardens} (see Details).
}
\details{
An  object of class \code{SLik_j} contains a simulated joint distribution of parameters and (projected) summary statistics, and a fit of a multivariate gaussian mixture model to this simulated distribution, the \dQuote{jointdens}, from which a marginal density \dQuote{pardens} of parameters can be deduced. The raw likelihood(P;D) is the probability of the data D given the parameters P, viewed as function the parameters and for fixed data. It is inferred as jointdens(D,P)/pardens(P) (for different P, each of jointdens and pardens are probabilities from a single (multivariate) gaussian mixture model, but this is not so for their ratio). 

When pardens(P) is low, indicating that the region of parameter space around P has been poorly sampled in the reference table, inference of likelihood in such regions is unreliable. Spuriously high likelihood may be inferred, which results notably in poor inference based on likelihood ratios. For this reason, it is often better to use the argument  \code{which="safe"} whereby the likelihood may be penalized where pardens(P) is low. The penalization is applied to cases where the uncorrected likelihood is higher than the maximum one in the reference table, and pardens(P) is lower than a threshold also determined from the reference table. The source code should be consulted for further details.  
}
\examples{
## see help("example_reftable")
}
