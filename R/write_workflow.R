write_workflow <- function(con=stdout(), # passed to writeLines
                           lower, upper, nUnique, Simulate,
                           simulator_args=NULL, 
                           ...) { # dots are used as refine() arguments...    note the "max" inside
  
  template <- "
  # Script generated on date_ by Infusion vversion_
  
  set.seed(123)
  parsp <- init_reftable(lower=lower_,
                         upper=upper_, 
                         nUnique = nbUnique_) 

  simuls <- add_reftable(, Simulate=\"Simulate_\" ,par.grid=parsp, verbose.=FALSE simulator_args__)
                            
  nThr <- parallel::detectCores() - 1L  # Projections may use parallelisation.                          
projections__                  
    
  # With appropriate use of the 'is_trainset' and 'use_oob' arguments:   
  projSimuls <- project(dsimuls, projectors=projectors, verbose=FALSE, is_trainset=TRUE)
  projSobs <- project(Sobs, projectors=projectors, use_oob=FALSE)
    
  slik_j <- infer_SLik_joint(projSimuls, stat.obs=projSobs, verbose=TRUE)
  slik_j <- MSL(slik_j)
  
  plot1Dprof(slik_j)
    
  slik_jn <- refine(slik_j dots__)

  "
  lower_ <- spaMM:::.DEPARSE(lower) 
  upper_ <- spaMM:::.DEPARSE(upper)
  nbUnique_ <- paste(nUnique) 
  
  simulator_args__ <- ""
  for (st in names(simulator_args)) {
    simulator_args__ <-  paste0(simulator_args__,",\n      ",st,"=", simulator_args[[st]])  
  }
  
  tmp1 <- gsub("date_", date(), template)
  tmp1 <- gsub("version_", packageVersion("Infusion"), tmp1)
  tmp1 <- gsub("lower_",lower_, tmp1)
  tmp1 <- gsub("upper_",upper_, tmp1)
  tmp1 <- gsub("nbUnique_",nbUnique_, tmp1)
  tmp1 <- gsub("simulator_args__",simulator_args__, tmp1)

  projections__ <- ""
  for (st in names(lower)) {
    projections__ <- paste0(projections__,"  ",st,"fit <- project(\"",st,"\", stats=names(Sobs), data=simuls, verbose=TRUE, 
                      methodArgs=list(num.threads=nThr))\n")
  }
  
  projnames <- paste0("p", names(lower))
  projfits <- paste0(names(lower),"fit")
  locstring <- paste0(paste0(projnames,"=",projfits, sep=", "),c("","","","\n                     "), collapse="")
  locstring <- substr(locstring,0,nchar(locstring)-2L)
  projections__ <- paste0(projections__,"  projectors <- list(",locstring,")")
  
  tmp1 <- gsub("projections__",projections__, tmp1)

  mc <- match.call()
  dotlist <- list(...)
  if (is.null(dotlist$nbCluster)) dotlist$nbCluster <-"max"
  if (is.null(dotlist$methodArgs)) dotlist$methodArgs <- list()
  dots__ <- ""
  count <- 0L
  for (st in names(dotlist)) {
    count <- count + 1L
    if ((count %% 4L)==0) {
      dots__ <-  paste0(dots__,",\n                    ",st,"=", spaMM:::.DEPARSE(mc[[st]]))  
    } else dots__ <-  paste0(dots__,", ",st,"=", spaMM:::.DEPARSE(mc[[st]]))  
  }

  tmp1 <- gsub("dots__",dots__, tmp1)
  
  writeLines(tmp1, con=con)
}

if (FALSE) {
  write_workflow(
    # arguments for init_reftable():
    lower=c(logTh1=-2,logTh2=-2,logTh3=-2,logTh4=-2, ar=0.01,logMu=-5, MEANP=0.01),
    upper=c(logTh1=1,logTh2=1,logTh3=1,logTh4=1, ar=0.99,logMu=-2,  MEANP=0.99),
    nUnique = 1000,
    # for add_reftable():
    Simulate="onesim_c",
    simulator_args= list(
      exe_path="'C:/home/francois/travail/stats/Infusionplus/caseStudies/Harmonia/HA_for_testing_infusion_RF'",
      simulator_call="'diyabc-RF-windows-v1.1.27.exe -p ./ -R \"ALL\" -r 1 -g 1'"
    ),
    # optional arguments for refine():    
    n=9000/7.2, CIs=TRUE, update_projectors=FALSE)
  
  write_workflow(
    ## arguments for init_reftable():
    lower=c(logTh1=-2,logTh2=-2,logTh3=-2,logTh4=-2, ar=0.01,logMu=-5, MEANP=0.01),
    upper=c(logTh1=1,logTh2=1,logTh3=1,logTh4=1, ar=0.99,logMu=-2,  MEANP=0.99),
    nUnique = 1000,
    #
    ## for add_reftable():
    Simulate="schtroumf",   # name of a user-defined R function
    simulator_args= list(   # Imagine that schtroumf() has arguments 'exe_path' and 'cmdline':
      exe_path="'path_to_smurf_executable'",
      cmdline="'smurf.exe -a -b -c -d'"
    ),
    #
    ## optional arguments for refine():    
    n=9000/7.2, CIs=TRUE, update_projectors=FALSE)
  
  
}

