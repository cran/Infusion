cat(cli::col_yellow("test inference with reparametrized reftables:\n"))

## Toy simulation function
hezsim <- function(logNe=parvec["logNe"],logmu=parvec["logmu"],parvec) {
  Ne <- 10^logNe
  mu <- 10^logmu
  Es <- Ne*mu
  Vs <- 1/log(1+Ne) 
  genom_s <- rgamma(5, shape=Es/Vs,scale=Vs) # 5 summary statistics
  names(genom_s) <- paste0("stat",seq(5))
  genom_s
} 


{ ## Analysis with 'canonical' parameters
  #
  ## simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- hezsim(logNe=4,logmu=-4) 
  #
  parsp <- init_reftable(lower=c(logNe=1,logmu=-5), 
                         upper=c(logNe=6,logmu=-2))
  init_reft_size <- nrow(parsp)
  simuls <- add_reftable(Simulate=hezsim, parsTable=parsp)
  
  { 
    plogNe <- project("logNe", data=simuls, stats=paste0("stat",seq(5)))
    plogmu <- project("logmu", data=simuls, stats=paste0("stat",seq(5)))

    dprojectors <- list(plogNe=plogNe,plogmu=plogmu)
    
    projSimuls <- project(simuls,projectors=dprojectors,verbose=FALSE)
    projSobs <- project(Sobs,projectors=dprojectors)
  }
  
  { ## Estimation: 
    ddensv <- infer_SLik_joint(projSimuls,stat.obs=projSobs)
    dslik_j <- MSL(ddensv, eval_RMSEs=FALSE) ## find the maximum of the log-likelihood surface
    refined_dslik_j <- refine(dslik_j, eval_RMSEs=FALSE, CIs=FALSE)
  }
}

{  #### Reparametrization to composite parameters
  { ## Using convenience function reparam_fit()
    
    locreparamfn <- function(object, ...) {
      logTh <- object[["logmu"]]+object[["logNe"]]
      if (inherits(object,"data.frame")) {
        data.frame(logTh=logTh,
                   logNe=object[["logNe"]])
      } else if (is.matrix(object)) {
        cbind(logTh=logTh,
              logNe=object[["logNe"]])
      } else c(logTh=logTh,
               logNe=object[["logNe"]])
    }
    
    rps <- reparam_fit(refined_dslik_j, to=c("logTh","logNe"),
                       reparamfn = locreparamfn)
    plot(rps)
    
    hezsim2 <- function(logNe=parvec["logNe"],logTh=parvec["logTh"],parvec) {
      hezsim(logNe=logNe,logmu=logTh-logNe)
    } 
    
    rps <- reparam_fit(refined_dslik_j, to=c("logTh","logNe"),
                       reparamfn = locreparamfn, raw=TRUE,
                       # statNames=names(get_from(refined_dslik_j,"raw_data")),
                       reftable_attrs=list(Simulate=hezsim2))
    plot(rps)
    refine(rps)
    
  }
  
  { ## Lower-level coding allowing refines()
    reparametrize <- function(v) {
      logTh <- v[["logmu"]]+v[["logNe"]]
      return(c(logTh=logTh,logNe=v[["logmu"]]))
    }
    
    { ##Composite analysis directly on FINAL reftable from canonical -space analysis: 
      final_raw_canon <- get_from(refined_dslik_j, "reftable_raw") 
      final_rpsimuls <- cbind(t(apply(final_raw_canon[,1:2], 1L, reparametrize)),
                              final_raw_canon[,3:7])
      apply(final_rpsimuls,2,range)
      attr(final_rpsimuls,"LOWER") <- c(logTh=-4,logNe=-5)
      attr(final_rpsimuls,"UPPER") <- c(logTh=4,logNe=-2)
      
      prlogTh <- project("logTh", data=final_rpsimuls, stats=paste0("stat",seq(5)))
      prlogmu <- project("logNe", data=final_rpsimuls, stats=paste0("stat",seq(5)))
      
      rpprojectors <- list(prlogTh=prlogTh, prlogmu=prlogmu)
      
      projrpSimuls <- project(final_rpsimuls,projectors=rpprojectors,verbose=FALSE)
      projrpSobs <- project(Sobs,projectors=rpprojectors)
      
      { ## Estimation: 
        final_rpdensv <- infer_SLik_joint(projrpSimuls,stat.obs=projrpSobs)
        rpslik_j <- MSL(final_rpdensv, eval_RMSEs=FALSE) ## find the maximum of the log-likelihood surface
        plot(rpslik_j)
      }
      
    }
    
    { ##Composite analysis on INITIAL reftable followed by refine(., newsimuls)
      rpsimuls <- simuls[]
      rpsimuls <- cbind(t(apply(rpsimuls[,1:2], 1L, reparametrize)),
                        rpsimuls[,3:7])
      attr(rpsimuls,"LOWER") <- c(logTh=-4,logNe=-5)
      attr(rpsimuls,"UPPER") <- c(logTh=4,logNe=-2)
      
      prlogTh <- project("logTh", data=rpsimuls, stats=paste0("stat",seq(5)))
      prlogmu <- project("logNe", data=rpsimuls, stats=paste0("stat",seq(5)))
      
      rpprojectors <- list(prlogTh=prlogTh, prlogmu=prlogmu)
      
      projrpSimuls <- project(rpsimuls,projectors=rpprojectors,verbose=FALSE)
      projrpSobs <- project(Sobs,projectors=rpprojectors)
      
      { ## initial Estimation: 
        ini_rpdensv <- infer_SLik_joint(projrpSimuls,stat.obs=projrpSobs)
        ini_rpslik_j <- MSL(ini_rpdensv, eval_RMSEs=FALSE) ## find the maximum of the log-likelihood surface
        plot(ini_rpslik_j)
      }
      
      { # refine
        raw_canon_refined <- get_from(refined_dslik_j, "reftable_raw") 
        # str(raw_canon_refined)
        rpsimuls <- cbind(t(apply(raw_canon_refined[,1:2], 1L, reparametrize)),
                          raw_canon_refined[,3:7])
        newsimuls <- rpsimuls[-(1:init_reft_size), ]
        # str(newsimuls)
        refined_rpslik_j <- refine(ini_rpslik_j, newsimuls=newsimuls, eval_RMSEs=FALSE)
      }
      
    }
  }
}


