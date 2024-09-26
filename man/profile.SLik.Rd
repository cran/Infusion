\name{profile.SLik}
\alias{profile.SLik}
\alias{profile.SLik_j}
\alias{profile}
\title{
  Compute profile summary likelihood 
}
\description{
  Predicts the profile likelihood for a given parameter value (or vector of such values) using predictions from an \code{SLik_j} (or older \code{SLik}) object (as produced by \code{\link{MSL}}).
}
\usage{
\method{profile}{SLik_j}(fitted, value, fixed=NULL, return.optim=FALSE, 
        init = "default", which="safe", 
        constr_crits=fitted$constr_crits, 
        eq_constr=NULL, ...)
\method{profile}{SLik}(fitted, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted}{an \code{SLik} object.}
\item{value}{
  The parameter value (as a vector of named values) for which the profile is to be computed. Setting in explictly to NULL will maximize likelihood subject only to optional constraints specified by other arguments.
}
\item{fixed}{This argument appears redundant with the \code{value} argument, so it will be deprecated and its use should be avoided. When it is non-\code{NULL}, the profile is computed for \code{value} updated to \code{c(value, fixed)}.
}
\item{return.optim}{If this is TRUE, and if maximization of likelihood given \code{value} and \code{fixed} is indeed required, then the full result of the optimization call is returned.}
\item{constr_crits}{
Inequality constraints, by default those provided in the first iteration of the workflow, if any. See \code{\link{constr_crits}} for details. 
}
\item{eq_constr}{
Optional equality constraints, provided in the same format as \code{\link{constr_crits}} (\bold{This feature is experimental}: in particular the procedure for finding initial values for maximization of likelihood does not use the \code{eq_constr} information). 
}
  \item{\dots}{
  For \code{SLik_j} method, arguments passed to \code{SLik} method.
  For \code{SLik_j} method, currently not used.
}
\item{init}{Better ignored. Either a named vector of parameter values (initial value for some optimizations) or a character string. The default is to call a procedure to find a good initial point from a set of candidates. The source code should be consulted for further details and is subject to change without notice.}
\item{which}{Better ignored (for development purpose).} % passed to predict.SLik_j, but predict functions are not part of API
}
\value{
If \code{return.optim} is FALSE (default): the predicted summary profile log-likelihood, with possible attribute \code{"solution"}, the optimization solution vector (named numeric vector, missing if no profiling was needed). if \code{return.optim} is TRUE, the result of an optimization call, a list including elements \code{solution} (solution vector) and \code{objective} (log-likelihood). 
}
\seealso{\code{\link{example_reftable}}}
\examples{
## see e.g. 'example_reftable' documentation
}
