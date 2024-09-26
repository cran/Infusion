.constrOptim <- function(init, objfn, lower, upper, object, neg_ineq_constrfn, template, 
                         eq_constrfn, ...) { 
  if ("template" %in% names(formals(objfn))) {  # profile
    optr <- nloptr(x0=init, eval_f=objfn, eval_g_ineq = neg_ineq_constrfn,
                   eval_g_eq=eq_constrfn,
                   lb=lower,
                   ub=upper, 
                   template=template,
                   opts=list("algorithm"="NLOPT_LN_AUGLAG",
                             xtol_rel=1.0e-8,
                             local_opts=list(algorithm="NLOPT_LN_BOBYQA",
                                             xtol_rel=1.0e-8))
    )
  } else { # global max. The init must be a safe one that satisfy constraints
    optr <- nloptr(x0=init, eval_f=objfn, eval_g_ineq = neg_ineq_constrfn,
                   eval_g_eq=eq_constrfn,
                   lb=lower,
                   ub=upper, 
                   opts=list("algorithm"="NLOPT_LN_AUGLAG",
                             xtol_rel=1.0e-8,
                             local_opts=list(algorithm="NLOPT_LN_BOBYQA",
                                             xtol_rel=1.0e-8))
    )
  }
  names(optr$solution) <- names(init)
  optr
}


# A public interface for .safe_optim is profile.SLik()
.safe_optim <- function (init, objfn, lower, upper, template=NULL, neg_ineq_constrfn=NULL, 
                         eq_constrfn=NULL, ..., 
                         object) { 
  # special handling of template argument bc nloptr in passing all \dots it receives to the objective function,
  # so none may be 'ignored' (unless that function itself has \dots, but then its even worse: nloptr finds that 
  # "eval_f requires argument '...' but this has not been passed to the 'nloptr' function.")
  if (is.null(neg_ineq_constrfn) && is.null(eq_constrfn)) { 
    if (is.null(template)) {
      .safe_opt(init=init, objfn=objfn, lower=lower, upper=upper, ...)
    } else .safe_opt(init=init, objfn=objfn, lower=lower, upper=upper, template=template, ...)
  } else {
    if (is.character(init)) { # no init point was found satisfying constraints
      return(list(info=init, objective=Inf)) # ie logL=-Inf
    } else if (is.null(template)) {
      .constrOptim(init=init, objfn=objfn, lower=lower, upper=upper, 
                   object=object, neg_ineq_constrfn=neg_ineq_constrfn, 
                   eq_constrfn=eq_constrfn, ...)
    } else .constrOptim(init=init, objfn=objfn, lower=lower, upper=upper, 
                              template=template, object=object, neg_ineq_constrfn=neg_ineq_constrfn,
                        eq_constrfn=eq_constrfn, ...)
  } 
}


## Problem that this fn aims to solve: 
## an MSLE optimization result may be beyond the bounds by O(flop error).
## If we recycle this MSLE to init a later optimization, an error may result. Hence
.sanitize_optim_solution <- function(sol, lower, upper) {
  safelow <- lower[names(sol)]+1e-15
  risky <- (sol < safelow)
  sol[risky] <- safelow[risky]
  safeup <- upper[names(sol)]-1e-15
  risky <- (sol > safeup)
  sol[risky] <- safeup[risky]
  sol
}
