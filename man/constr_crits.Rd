\name{constr_crits}
\alias{constr_crits}
\alias{constraints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specificying arbitrary constraints on parameters
}
\description{
\code{constr_crits} is an argument of \code{\link{init_reftable}}, \code{\link{add_reftable}},  \code{\link{infer_SLik_joint}}, and \code{profile.SLik_j},
allowing to specify constraints on parameters, beyond the ones defined by the ranges over each parameter.  
Depending on the function it controls, this argument will affect the generation of parameter points or the maximization of the summary-likelihood. 

It is advised to use the \code{constr_crits} when calling \code{init_reftable} so that the initial reference table satisfies the constraints. 
But it is more important to use it when calling \code{\link{infer_SLik_joint}}, whose return value 
contains the \code{constr_crits} information it was given, allowing subsequent \code{\link{refine}} calls to take constraints into account. 
 
The constraints are represented as a mathematical expression for a vector of quantities that should all be negative when the constraints are satisfied. For example, to incorporate the constraint \code{t1 < t3 && t2 < t3} between three time parameters \code{t1}, \code{t2} and \code{t3} present in the reference table, one can use\cr
  \code{constr_crits = quote({ c(t1-t3, t2-t3) })}.\cr

For computational efficiency, it may be better to avoid using this feature when the constraints 
can be represented as box constraints (i.e., independent ranges for each parameter) by some intuitive reparametrization.
}
\examples{
## A toy example
if (Infusion.getOption("example_maxtime")>9) {
  myrnorm <- function(mu,s2,sample.size=40L) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    return(c(mean=mean(s),var=var(s)))
  } # simulate means and variances of normal samples of size 'sample.size'
  
  # The constraint (the final plot produced by refine() explains it):
  heart <- quote({ x <- 3*(mu-4.25);  y <- 3*(s2-0.75); x^2+(y-(x^2)^(1/3))^2-1})
  
  set.seed(123)
  Sobs <- myrnorm(mu=4,s2=1) 
  parsp_h <- init_reftable(lower=c(mu=2.8,s2=0.4), upper=c(mu=5.2,s2=2.4),
                           constr_crits=heart)
  simuls_h <- add_reftable(Simulate="myrnorm", parsTable=parsp_h)
  c_densv <- infer_SLik_joint(simuls_h, stat.obs=Sobs, constr_crits = heart)
  c_slik_j <- MSL(c_densv, CIs=FALSE, eval_RMSEs=FALSE) 
  c_slik_j <- refine(c_slik_j, target_LR=10, ntot=3000, eval_RMSEs=FALSE) 
}
}
