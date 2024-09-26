cat(cli::col_yellow("test inference about latent variables:\n"))

{
  
  ##### A toy example motivated by some inference problem for genomic data #####
  # A model with three parameters logNe, Vs and Es is fitted.
  # Data from 100 loci are here summarized by three genome-wide summary statistics 
  # (slogNe, sVs and sEs), and one locus-specific statistic that will provide 
  # information about a locus-specific latent variable.
  
  ## Simulation function 
  genomloc <- function(logNe=parvec["logNe"],Es=parvec["Es"],Vs=parvec["Vs"], 
                       latent=TRUE, # returns the latent value by default
                       parvec) {
    slogNe <- rnorm(1,logNe, sd=0.3)
    genom_s <- rgamma(99, shape=Es/Vs,scale=Vs) # all loci except the focal one
    sEs <- mean(genom_s)
    sVs <- var(genom_s)
    latentv <- rgamma(1, shape=Es/Vs,scale=Vs) # locus-specific latent variable to predict
    sloc <- rnorm(1, mean=latentv-sEs,sd=latentv/3) # locus-specific statistic
    resu <- list(slogNe=slogNe,sEs=sEs,sVs=sVs, sloc=sloc)
    if (latent) resu$latentv <- latentv
    unlist(resu)
  } 
  #
  ## simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- genomloc(logNe=4,Es=0.05, Vs=0.1,latent=FALSE) ## no latent value in Sobs
  #
  workflow_design <- get_workflow_design(npar=3L, n_proj_stats=4L, n_latent=1L)
  parsp <- init_reftable(lower=c(logNe=2,Es=0.001,Vs=0.001), 
                         upper=c(logNe=6,Es=0.2,Vs=0.2),
                         nUnique=workflow_design$init_reft_size)
  simuls <- add_reftable(Simulate=genomloc, parsTable=parsp)
  
  simuls <- declare_latent(simuls,"latentv") 
  
  ## Projections are not necessary here since the number of statistics is minimal,
  # but will be discussed later.
  { ############ Without projections
    { ## Usual workflow for estimation: 
      densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
      slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
      slik_j <- refine(slik_j,maxit=2,update_projectors=TRUE)
      # plot1Dprof(slik_j) ## 1D profiles show parameter inference is OK
    }
    { ## inference about latent values:
      pplatent(slik_j)
      pplatent(slik_j, type="median")
      latint(slik_j, nsim=999, levels=c(0.025,0.5,0.975))
    }
    { ## Assessing prediction of latent variable:
      # Builds testing set:
      test_simuls <- t(replicate(1000, genomloc(logNe=4,Es=0.05, Vs=0.1)))
      test_data <- test_simuls[,-5]
      # Point prediction:
      pred <- pplatent(slik_j, sumstats = test_data)
      
      plot(test_simuls[,"latentv"], pred); abline(0,1) # prediction vs. true latent values
    }
  }
  
  { ########## Beyond standard use of projections for estimation of parameter values, 
    # projections can also be used when several individual-level statistics inform about 
    # the latent variable, to reduce them to a single summary statistic.
    # Projection will then be needed at the prediction step too.
    
    { # projection with latent variable as response:
      platent <- (project("latentv", data=simuls, stats=c("slogNe","sEs","sVs","sloc")))
      # (This example only serves to show the syntax since no dimention reduction occurs)
      
      dprojectors <- list(SLOC=platent,slogNe=NULL,sEs=NULL, sVs=NULL)
      
      # => As soon as one projection is used, The 'projectors' argument must include 
      # all projectors used for the inference, whether for parameters or for latent variables. 
      # NULL projectors should then be declared for raw statistics retained 
      # in the projected reference table.
      
      # Apply projections on simulated statistics and 'data':
      projSimuls <- project(simuls,projectors=dprojectors,verbose=FALSE)
      projSobs <- project(Sobs,projectors=dprojectors)
    }
    
    { ## Estimation: 
      ddensv <- infer_SLik_joint(projSimuls,stat.obs=projSobs)
      dslik_j <- MSL(ddensv) ## find the maximum of the log-likelihood surface
      dslik_j <- refine(dslik_j,maxit=2,update_projectors=TRUE)
      # plot1Dprof(dslik_j)
    }
    
    { ## Assessing prediction of latent variable: do not forget to project!
      
      test_simuls <- t(replicate(1000, genomloc(logNe=4,Es=0.05, Vs=0.1)))
      test_data <- test_simuls[,-5] # removing column of latent variable
      ptest_data <- project(test_data,projectors=dprojectors,verbose=FALSE) # Here!
      pred <- pplatent(dslik_j, sumstats = ptest_data)
      
      plot(test_simuls[,"latentv"], pred); abline(0,1)
    }
  }
  
  
}
