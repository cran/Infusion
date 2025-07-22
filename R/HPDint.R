#  not used
.safe_rparam <- function(object, nsim) {
  pardens_pts <- .sample_in_absol_constrs(object, nsim=nsim, object$colTypes$fittedPars, 
                                          density=object$pardens, 
                                          norm_or_t=.wrap_rmvnorm, unsort=TRUE)
  nr <- nrow(pardens_pts)
  # precompute all values for Metropolis algo
  pardensvals <- predict(object, newdata=pardens_pts, which="parvaldens",log=TRUE)
  thr_dpar <- .get_thr_info(object)$thr_dpar
  # loc_thr_dpar <- 2*thr_dpar -max(pardensvals) # ie max(pardensvals) - 2*(max(pardensvals)-thr_dpar)
  loc_thr_dpar <- 1.8*thr_dpar -0.8*max(pardensvals) # ie max(pardensvals) - 1.8*(max(pardensvals)-thr_dpar)
  
  ok <- (pardensvals>loc_thr_dpar)
  
  list(pts=pardens_pts[ok,,drop=FALSE],
       vals=pardensvals[ok])
}

.propose <- function(curr_pars,curr_clu,solve_t_chol_sigma_list, means,clu_probs,
                     OU_t) {
  .expm1.2 <- exp(-OU_t/2)
  .sqrt1mexpm1 <- sqrt(1-exp(-OU_t))
  nclu <- ncol(means)
  npar <- length(curr_pars)
  # sample y
  ## first sample a clu
  if (TRUE) {
    in_clu_probs <- numeric(nclu)
    for (it in seq_len(nclu)) {
      in_clu_probs[it] <- .fast_dmvnorm(curr_pars, mean=means[,it], 
                                        solve_t_chol_sigma=solve_t_chol_sigma_list[[it]], log = FALSE)
    }
    in_clu_probs <- in_clu_probs*clu_probs
    in_clu_probs <- in_clu_probs/sum(in_clu_probs)
    #
  } else in_clu_probs <-  rep(1/nclu,nclu)
  ## OU sampling in chosen cluster
  prop_clu <-  sample.int(nclu, size=1, prob=in_clu_probs)
  E_proposal <- means[,prop_clu] + .expm1.2 * (curr_pars-means[,prop_clu])
  prop_pars <- E_proposal + 
    .sqrt1mexpm1 * solve(solve_t_chol_sigma_list[[prop_clu]],rnorm(npar, mean=0,sd=1))
  prop_pars <- drop(prop_pars)
  #
  # all paths from x to y
  Txy_terms <- matrix(ncol=nclu,nrow=1L)
  for (it in seq_len(nclu)) {
    E_proposal_it <- means[,it] + .expm1.2 * (curr_pars-means[,it])
    Txy_terms[it] <- log(in_clu_probs[it]) + 
      .fast_dmvnorm(prop_pars, 
                    mean=E_proposal_it, 
                    solve_t_chol_sigma=solve_t_chol_sigma_list[[it]]/.sqrt1mexpm1, 
                    log=TRUE)
  }
  Txy <- .logsumexp(Txy_terms)
  
  # all paths from y to x
  Tyx_terms <- matrix(ncol=nclu,nrow=1L)
  if (TRUE) {
    in_clu_probs <- numeric(nclu)
    for (it in seq_len(nclu)) {
      in_clu_probs[it] <- .fast_dmvnorm(prop_pars, mean=means[,it], 
                                        solve_t_chol_sigma=solve_t_chol_sigma_list[[it]], log = FALSE)
    }
    in_clu_probs <- in_clu_probs*clu_probs
    in_clu_probs <- in_clu_probs/sum(in_clu_probs)
    #
  } else in_clu_probs <-  rep(1/nclu,nclu)
  for (it in seq_len(nclu)) {
    E_proposal_it <- means[,it] + .expm1.2 * (prop_pars-means[,it])
    Tyx_terms[it] <- log(in_clu_probs[it]) + 
      .fast_dmvnorm(curr_pars, 
                    mean=E_proposal_it, 
                    solve_t_chol_sigma=solve_t_chol_sigma_list[[it]]/.sqrt1mexpm1, 
                    log=TRUE)
  }
  Tyx <- .logsumexp(Tyx_terms)
  
  list(prop_clu=prop_clu,
       prop_pars=prop_pars,
       Txy=Txy,
       Tyx=Tyx)
}

.Hastings_sampler_unstuck <- function(object, 
                                   nsim, 
                                   max.stuck=nsim %/% 100L,
                                   dif_thr=log(max.stuck/2), 
                                   target_clu_size=1000L,
                                   OU_t=0.02) {
  clu_probs <- exp(object$clu_params$logproportions)
  nclu <- length(clu_probs)
  clu_means <- object$clu_params$pardens_means
  solve_t_chol_sigma_list <- object$clu_params$solve_t_chol_sigma_lists$pardens
  posterior_sample <- proposal_pts <- matrix(NA, nrow=nsim, ncol=nrow(clu_means))
  curr_clu <- which.max(clu_probs)
  posterior_sample[1,] <- curr_pars <- clu_means[,curr_clu]
  curr_logL <- predict(object, newdata=curr_pars, which="safe",log=TRUE)  # includes user prior
  curr_postdens <- curr_logL
  
  # set.seed(123)
  switches <- stuck <- stucks <- 0L 
  for (it in 2:nsim) {
    # if (it==2500L) browser()
    proposal <- .propose(curr_pars=curr_pars,
                         curr_clu=curr_clu,
                         solve_t_chol_sigma_list, 
                         means=clu_means,clu_probs=clu_probs,
                         OU_t=OU_t)
    proposal_pts[it,] <- prop_pars <- proposal$prop_pars
    prop_logL <- predict(object, newdata=prop_pars, which="safe",log=TRUE)
    prop_postdens <- prop_logL  # includes user prior
    alpha <- prop_postdens - proposal$Txy - curr_postdens + proposal$Tyx
    alpha <- exp(alpha)
    if (runif(1)<=alpha) { 
      # if (clu_curr != clu_prop) acc_clu   <- acc_clu+1L  
      # acc <- acc+1L
      if (is_really_stuck <- stuck > max.stuck) {
        stucks <- stucks + 1L
# if (stuck>1000L) browser()        
      }
      curr_pars <- prop_pars
      curr_clu  <- proposal$prop_clu
      curr_postdens <- prop_postdens
      stuck <- 0L
      switches <- switches+1L
    } else  stuck <- stuck+1L
    posterior_sample[it,] <- curr_pars
    
  }
  if (stucks) warning(paste("The chain got stuck",stucks,"times."), immediate. = TRUE)
  colnames(posterior_sample) <- colnames(proposal_pts) <- names(curr_pars)
  
  
  if (FALSE) {
    transition_facs <- postdens-pardensvals        
    
    # => high logLik, low parvaldens: high transition_fac, the MC gets stuck. We try to fix this:
    # ***This block is purely ad hoc for the present application, and should be ignored in toy examples***
    ord  <- order(transition_facs,decreasing = TRUE)
    dif <- - diff(transition_facs[ord])
    pardensordm1 <- (pardensvals[ord])[-nr]
    max_extreme  <- which.min(dif > dif_thr &
                                pardensordm1 < quantile(pardensordm1,0.05)
    ) -1L # which.min() returns position of first FALSE
    for (it in seq_len(max_extreme)) {
      transition_facs[ord[1:it]]  <- transition_facs[ord[1:it]] + dif_thr-dif[it] # amounts to reducing large logL for small pardens
    } # so that max dif for extreme values is reduced to dif_thr.
    # ? Benefits are more obvious for larger nsim. ?
    
  }
  
  resu <- list(posterior_sample=posterior_sample,
                     parm_pts_ord=proposal_pts, # only info for diagnostic plot; ordering consistent with that of the posterior_sample
                     switches=switches)
  resu
}

.HPDint_from_post_sample <- function(posterior_parm, # posterior_sample, 
                     object, level, verbose, plot.,
                     seed=Infusion.getOption("mixmodSeed"),
                     parm=colnames(posterior_parm)) {
  
  posterior_parm  <- as.data.frame(posterior_parm)
  locstrategy <- .Infusion.data$options$get_mixModstrategy(nc=1)
  nbcluster <- max(get_nbCluster_range(projdata=posterior_parm, verbose=FALSE)) # nc=1 implicitly
  locarglist <- list(data=posterior_parm, strategy=locstrategy, 
                    nbCluster=nbcluster, seed=seed)
  clu_post <- .mixmodCluster_w_re_retries(data=posterior_parm, locarglist = locarglist)
  clu_post@freq <- 0
  clu_post@varNames <- parm
  clu_post@statNames <- character(0)
  clu_post@Sobs_activeBoundaries <- NULL
  clu_post@simuls_activeBoundaries <- NULL
  
  # Finely discretized uniform exploration:
  unifparm <- seq(object$LOWER[[parm]], object$UPPER[[parm]], length.out=10000)
  unifdf  <- data.frame(X1=unifparm)
  colnames(unifdf)  <- parm
  
  clu_post_densv  <- predict(
    clu_post, newdata = unifdf, # unifdf !!!!!!!
    logproportions = log(clu_post@parameters["proportions"]),
    clu_means = t(clu_post@parameters["mean",]),
    solve_t_chol_sigma_list = lapply(clu_post@parameters["variance"], .solve_t_cholfn)
  )
  if (plot.) plot(unifparm, clu_post_densv)
  # range of HPD region:
  orderdens <- order(clu_post_densv, decreasing = TRUE)
  cumdens <- cumsum(clu_post_densv[orderdens])
  cumdens <- cumdens/tail(cumdens,1)
  idx <- which.min(cumdens < level)
  range(unifparm[orderdens[1:idx]])
}

# This is tested by D:/home/francois/travail/stats/Infusionplus/diy2inf_simuls/Harmonia/D_axy_8pars_logTh4/calc_HPD_intervals.R
.HPDint <- function(object, parms, nsim=10000L, verbose=TRUE, plot.=2L, level=0.95,
                    method="Hastings",
                    max.cor=0.05, lag=NULL, ..., debug=FALSE) {
  warn <- NULL
  if (method=="HMC") {
    if ( ! is.null(object$prior_logL) && is.null(object$grad_prior_logL)) {
      stop("THe gradient function for the prior_logL() is not found.")
    }
    glogPOSTERIOR <- .get_gradfn(object, init=object$MSL$MSLE, sign=1L)
    hmc <- .get_wrap("hmc", pack="hmclearn")
    hmcresu <- hmc(nsim,theta.init=object$MSL$MSLE,
        logPOSTERIOR = function(v) predict(object,newdata=v, log=TRUE),
        glogPOSTERIOR = glogPOSTERIOR,
        verbose=FALSE, ...)
    posterior_sample <- hmcresu$thetaCombined[[1]]
    parm_pts_ord <- do.call(rbind,hmcresu$theta.all[[1]])
    switches <- hmcresu$accept
  } else {
    mcmcresu <- .Hastings_sampler_unstuck(object=object, nsim=nsim, ...)
    posterior_sample <- mcmcresu$posterior_sample
    parm_pts_ord <- mcmcresu$parm_pts_ord
    switches <- mcmcresu$switches
  }
  nr <- nrow(parm_pts_ord)
  
  if (plot. > 1L) oldpar <- par(mfrow=c(2L,1L)) 
  
  resu <- vector("list", length(parms))
  names(resu) <- parms
  for (parm in parms) {
    posterior_parm <- posterior_sample[,parm,drop=FALSE]
    if (plot.) {
      plot(parm_pts_ord[,parm], col="grey")
      points(posterior_parm)
    }
    if (is.null(lag)) {
      if (max.cor < 1L) {
        mlog2cor <- - log(max.cor,2)
        if (mlog2cor <= 0) stop("'max.cor' must be 0 < . < 1 .")
        lag.max <- 2*mlog2cor *nr / switches
        while (lag.max < nr) {
          acfresu <- stats::acf(posterior_parm,plot = FALSE,lag.max = lag.max)
          lag  <- which.max(acfresu$acf < max.cor)
          if (lag>1L) break
          lag.max  <- lag.max*5
          # if (plot.) plot(acfresu)
        }
      } 
      if (problem <- is.null(lag)) { 
        warn <- "No lag found. Check chain."
        if (verbose ) warning(warn, immediate.=TRUE)
        lag <- nr/31L
      }
    } else problem <- FALSE
    
    sequ <- seq(to=nr, by=lag)[-1]
    if (length(sequ)<30L) {
      warn <- "Lag too large. Check chain."
      if (verbose) warning(warn, immediate.=TRUE)
      ci <- .HPDint(object=object, parms=parm, nsim=31L*lag, verbose=verbose, plot.=plot., 
                    level=level, method=method,
                    # max.cor=0.05, 
                    lag=lag, ..., debug=FALSE)[[parm]]$CI
      if (problem) ci <- structure(c(NA,NA), info=list(badCI=ci))
    } else {
      if (verbose) {
        cat(length(sequ),"points retained out of",nr,"points sampled (",switches,"switches)\n")
      }
      ci  <- .HPDint_from_post_sample(posterior_parm=posterior_parm[sequ,, drop=FALSE],
                      # posterior_sample, 
                      object=object, level=level, verbose=verbose, plot.=plot.)
    } 
    attr(ci,"warn") <- warn
    if (debug && anyNA(ci)) browser()
    resu[[parm]] <- list(CI=ci)
  }
  if (plot. > 1L) par(oldpar)
  resu
}

