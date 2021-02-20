.safe_constrOptim <- function (theta, f, 
                               grad=NULL, # currently not used
                               ui, ci, mu = 1e-04, 
                               outer.iterations = 100, 
                               lower=NULL,upper=NULL,
                               outer.eps = 1e-05, ...) 
{
  fbar <- function(theta, gi.old, ...) {
    ui.theta <- ui %*% theta
    gi <- ui.theta - ci
    if (any(gi < 0)) return(NaN)
    #gi.old <- ui %*% theta.old - ci
    bar <- sum(gi.old * log(gi) - ui.theta)
    if (!is.finite(bar)) bar <- -Inf # then (if mu>0) return value is +Inf. optim handled this, nloptr seems to, 
                                     # but bobyqa first points are at boundsand it mis-performs in that case
    totCounts <<- totCounts+1L
    f(theta, ...) - mu * bar
  }
  # dR <- function(theta, gi.old, ...) {
  #   ui.theta <- ui %*% theta
  #   gi <- drop(ui.theta - ci)
  #   #gi.old <- drop(ui %*% theta.old - ci)
  #   dbar <- colSums(ui * gi.old/gi - ui)
  #   grad(theta, ...) - mu * dbar
  # }
  if (any(ui %*% theta - ci <= 0)) 
    stop("initial value is not in the interior of the feasible region")
  totCounts <- 1L
  obj <- f(theta, ...)
  r <- fbar(theta, gi.old=drop(ui %*% theta - ci), ...)
  fbar_obj <- function(theta) fbar(theta, gi.old, ...) # removed the ... from 'fun' args  bc not working with nloptr -> .checkfunargs()
  s.mu <- sign(mu)
  if (is.null(lower)) lower <- rep(1,length(theta))*-Inf
  if (is.null(upper)) upper <- rep(1,length(theta))*Inf
  for (i in seq_len(outer.iterations)) {
    obj.old <- obj
    r.old <- r
    #theta.old <- theta
    gi.old <- drop(ui %*% theta - ci)
    optr <- .safe_opt(init=theta, fbar_obj, lower=lower, upper=upper, verbose=FALSE, LowUp=list(), ...) 
    r <- optr$objective
    if (is.finite(r) && is.finite(r.old) && abs(r - r.old) < (0.001 + abs(r)) * outer.eps) 
      break
    theta <- optr$solution
    obj <- f(theta, ...)
    if (s.mu * obj > s.mu * obj.old) {
      if (obj > obj.old)  optr <- optr_old # .safe_opt() result worse than initial value: can occur by bobyqa() with fbar_obj() infinite at boundaries
      break
    }
    optr_old <- optr
  }
  resu <- list(optr=optr, outer.iterations=i, objective=f(optr$solution, ...), solution=optr$solution, 
               counts=c(totCounts+1L,0L) # consistent with optim() result.
               )
  resu$barrier.value <- optr$objective - resu$objective
  if (i == outer.iterations) {
    resu$convergence <- 7
    resu$message <- gettext("Barrier algorithm ran out of iterations and did not converge")
  }
  resu
}

if (FALSE) {
  # from help("constrOptim")
  fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }
  grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
      200 *      (x2 - x1 * x1))
  }
  
  #Box-constraint, optimum on the boundary
  constrOptim(c(-1.2,0.9), fr, grr, ui = rbind(c(-1,0), c(0,-1)), ci = c(-1,-1))
  Infusion:::.safe_constrOptim(c(-1.2,0.9), fr, ui = rbind(c(-1,0), c(0,-1)), ci = c(-1,-1), lower=c(-2,-2),upper=c(1,1))
}

