\name{constr_crits}
\alias{constr_crits}
\alias{constraints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specificying arbitrary constraints on parameters
}
\description{
\code{constr_crits} is an argument of \code{\link{init_reftable}}, \code{\link{add_reftable}},  \code{\link{infer_SLik_joint}}, and \code{profile.SLik_j},
allowing to specify constraints on parameters, beyond the ones defined by the ranges over each parameter.  Depending on the function it controls, this argument will affect the generation of parameter points or the maximization of the summary-likelihood. The \code{\link{infer_SLik_joint}} return value contains the \code{constr_crits} information it was given, allowing subsequent \code{\link{refine}} calls to take constraints into account. 
 
The constraints are represented as a mathematical expression for a vector of quantities that should all be negative when the constraints are satisfied. For example, to incorporate the constraint \code{t1 < t3 && t2 < t3} between three time parameters \code{t1}, \code{t2} and \code{t3} present in the reference table, one can use\cr
  \code{constr_crits = quote({ c(t1-t3, t2-t3) })}.\cr
See Examples in \code{\link{infer_SLik_joint}} for a nice artificial toy example.

For computational efficiency, it may be better to avoid using this feature when the constraints can be represented as box constraints (i.e., independent ranges for each parameter) by some intuitive reparametrization.
}
