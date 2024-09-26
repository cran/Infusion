init_reftable <- function (lower = c(par = 0), upper = c(par = 1), steps = NULL, 
                           nUnique = NULL, 
                           maxmin=(nUnique * length(lower)^2)<400000L, 
                           jitterFac = 0.5,
                           constr_crits=NULL, ...) {
  if (is.null(nUnique)) {
    np <- length(which(upper[names(lower)] != lower))
    nUnique <- get_workflow_design(np, ...)$init_reft_size
    if (missing(maxmin)) maxmin <- (nUnique * length(lower)^2)<400000L
  }
  parsTable <- init_grid(lower = lower, upper = upper, steps = steps, 
            nUnique = nUnique, maxmin=maxmin, jitterFac = jitterFac, nRepl = 0L)  
  if ( ! is.null(constr_crits)) { # selection of points that satisfy parameter constraints
    constrs <- apply(parsTable,1L, function(v) all(eval(constr_crits, envir = as.list(v))<0))
    parsTable <- parsTable[constrs,, drop=FALSE]
    while (nrow(parsTable)< nUnique) {
      message(paste0(sum( ! constrs)," points not satisfying constraints: sampling new points..."))
      morepts <- init_grid(lower = lower, upper = upper, steps = steps, 
                             nUnique = nUnique, maxmin=maxmin, jitterFac = jitterFac, nRepl = 0L)  
      constrs <- apply(morepts,1L, function(v) all(eval(constr_crits, envir = as.list(v))<0))
      morepts <- morepts[constrs,, drop=FALSE]
      parsTable <- rbind(parsTable, morepts)
    }
    if (nrow(parsTable)> nUnique) parsTable <- parsTable[1:nUnique,, drop=FALSE]
  }
  parsTable
}

`[.reftable` <- function (x, i, j, 
                          drop = TRUE, ## replicating default [.data.frame behaviour to allow selecting a single column as a vector, for plots etc 
                          rowdrop = FALSE, ## contrary to the data.frame default, which may produce a list-not-data-frame
                          dropcol=drop # resolves partial-match issues
) {
  class(x) <- "data.frame" ## avoids recursive call to `[.reftable
  if (missing(i)) {
    if (missing(j)) {
      resu <- x[ , , drop=dropcol] # by default replicating the effect of foo <- data.frame(x=1); foo[,,]
    } else resu <- x[, j, drop=dropcol]
  } else if (missing(j)) {
    resu <- x[i, , drop=rowdrop]
  } else {
    resu <- x[ ,j, drop=dropcol]
    if (is.null(dim(resu))) { # col dim has been dropped
      resu <- resu[i]
    } else resu <- resu[i, , drop=rowdrop] 
  }
  if (is.data.frame(resu)) {
    attrx <- attributes(x)
    attrx <- attrx[setdiff(names(attrx), c("names","row.names","class","na.action"))] # it make senses to remove "na.action" too, 
    # but `[.data.frame does not remove it (! this may defeat the purpose of the attribute: see Details of ?na.omit)
    for (st in names(attrx)) attr(resu,st) <- attrx[[st]]
    parNames <- names(attr(resu, "LOWER"))
    parKepts <- intersect(parNames,colnames(resu))
    attr(resu, "LOWER") <- attr(resu, "LOWER")[parKepts]
    attr(resu, "UPPER") <- attr(resu, "UPPER")[parKepts]
    class(resu) <- c("reftable", class(resu))
  } ## else vector
  return(resu)
} # Use unlist() to remove attributes from the return value

.warn_once_progressr <- function() {
    if ( ! .one_time_warnings$progressr_warned) {
      message("If the 'progressr' package were attached, a progress bar would be available.") 
      # attached, not simply installed => requires user action
      .one_time_warnings$progressr_warned <- TRUE
    } 
  }

.simulate_by_row <- function(Simulate, parsTable,
                             nRealizations,
                             verbose, nb_cores, packages,env,
                             control.Simulate,
                             cluster_args, cl_seed) {
  nsim <- nrow(parsTable)
  gridsimuls <- list()
  #
  cluster_args <- .set_cluster_type(cluster_args, nb_cores) # PSOCK vs FORK
  cores_info <- .init_cores(cluster_args=cluster_args)
  nb_cores <- cores_info$nb_cores
  #
  cl <- cores_info$cl
  if ( ! is.null(cl) && cluster_args$type!="FORK") {
    parallel::clusterExport(cl, "packages",envir=environment()) ## passes the list of packages to load
    # Infusion not leaded in child process !
    #parallel::clusterExport(cl, list("nRealizations"),envir=environment()) 
    #parallel::clusterCall(cl,fun=Infusion.options,nRealizations=nRealizations)
    if ( ! is.null(env)) parallel::clusterExport(cl=cl, ls(envir=env),envir=env)
  }
  as_one <- identical(names(nRealizations),"as_one")
  if (as_one) { 
    which_cl <- "param"
  } else {
    which_cl <- names(nb_cores) ## the name comes from explicit user input; this is doc'ed in ?add_simulation
    if (is.null(which_cl)) {
      if (nRealizations>1L) {
        if (cluster_args$spec>1L && cluster_args$type=="FORK") {
          warning("Using a FORK cluster with nRealizations>1 may be unreliable.")
        }
        which_cl <- "replic"
      } else which_cl <- "param"
    }
  }
  if (which_cl=="replic") {
    repl_cl <- cl ## the one used by simfn_per_par()
    if (length(repl_cl) > 1L && cluster_args$type!="FORK") {
      if (cores_info$has_doSNOW) {
        ## arguments such as Simulate are evaluated in the child envir so we pass definition under the "Simulate" name
        pb_char <- "N" # nested
      } else { ## arguments such as Simulate are evaluated in the parent envir to "myrnorm" so we pass def under "myrnorm" name
        pb_char <- "n" # nested
      }
      dotenv <- list2env(list(Simulate=Simulate) )
      parallel::clusterExport(cl=repl_cl, ls(envir=dotenv),envir=dotenv) ## exports eg "myrnorm": works for pbapply:: by not doSNOW     
    } else pb_char <- "s"
    param_cl <- NULL
  } else {
    param_cl <- cl
    if (length(param_cl) > 1L && cluster_args$type!="FORK") {
      if (cores_info$has_doSNOW) {
        pb_char <- "P"
      } else {
        pb_char <- "p"
      }
      dotenv <- list2env(list(Simulate=Simulate) )
      parallel::clusterExport(cl=param_cl, ls(envir=dotenv),envir=dotenv) ## exports eg "myrnorm": works for pbapply:: by not doSNOW     
    } else pb_char <- "s"
    repl_cl <- NULL 
  }
  #
  ############################################# simfn_per_par ######################################"
  simfn_per_par <- function(ii) {
    par <- parsTable[ii,,drop=FALSE] ## 1-row data frame: treated as list by do.call() or c(); but:
    parlist <- c(par,control.Simulate) ## now parlist is a list, not a data.frame
    if (length(repl_cl)>1L  && cores_info$has_doSNOW) { ## for the balancing only since we don't want the progress bar
      # R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
      # rdS_fn <- get("registerDoSNOW", asNamespace("doSNOW"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
      # do.call(rdS_fn,list(cl=repl_cl)) # this is what makes foreach see it and perform parallel computations
      # assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
      # if ( ! is.null(cl_seed) ) parallel::clusterSetRNGStream(cl = repl_cl, cl_seed) # that would mean several replicates for same parameter point => same simulation result
      foreach_blob <- foreach::foreach(i=1:length(repl_cl))
      abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
      foreach_args <- list(
        iii = seq_len(nRealizations), 
        .combine = "cbind",
        .packages= packages,
        .inorder = TRUE, .errorhandling = "pass"## "pass" to see error messages
      )
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      simuls <- foreach::`%dopar%`(foreach_blob, do.call(Simulate, parlist))
      if ( inherits(data.frame(simuls),"list") && ! is.null(mess <- simuls$message && is.list(mess) && is.character(mess[[1L]]))) {
        stop(paste0("The parallel call in add_simulation returned messages as\n",
                    mess[[1L]],
                    "If these indicate that some function could not be found,\n",
                    "try the 'env' argument, as in env=list2env(list(myfn=myfn))\n",
                    "to pass the 'myfn' function to all cores."))
      }
      # simuls <- foreach::`%dopar%`(foreach_blob, do.call("myrnorm",parlist)) #works
    } else if (as_one) { # case: no repl_cl parallelisation
      simuls <- do.call(Simulate,parlist)
    } else {
      # if ( ! is.null(cl_seed) ) parallel::clusterSetRNGStream(cl = repl_cl, cl_seed) # that would mean several replicates for same parameter point => same simulation result
      opb <- pboptions(type="none") ## silences the inner progress bar.S..
      on.exit(pboptions(opb)) ## ...but not the outer one
      simuls <- pbreplicate(n=nRealizations,do.call(Simulate,parlist),cl=repl_cl)
      if ((len<- length(dim(simuls)))>2L) {
        mess <- paste("The 'Simulate' function appears to return a ",len-1,"-dimensional table,\n",
                      " in a case where it should return a vector;\n",
                      " pbreplicate(n=nRealizations,do.call(Simulate,.),.) has dimensions",
                      paste("(",paste(dim(simuls),collapse=","),")"),"\n",
                      "  where nRealizations =",paste(deparse(nRealizations)),"\n",
                      "If you indeed want 'Simulate' to return all realizations in a table,\n",
                      "then use the named form nRealizations=c(as_one=.).\n",
                      "otherwise, check the return format of 'Simulate'.")
        stop(mess)
      }
    }
    #
    if (as_one) {
      if (is.null(ncol(simuls))) {
        stop(paste("The 'Simulate' value conflicts with 'nRealizations', because 'nRealizations' has name 'as_one',\n",
                   "  but 'Simulate' did not return a matrix."))
      } else if (ncol(simuls)!=nRealizations) stop(paste0("The 'Simulate' value conflicts with 'nRealizations', because nRealizations=c(as_one=",nRealizations,"),\n",
                                                          "  while 'Simulate' returned a matrix with ",ncol(simuls)," columns."))
    }
    if (is.null(dim(simuls))) { ## converts to 1-row matrix
      colName <- names(simuls[1])
      dim(simuls) <- c(length(simuls),1)
      colnames(simuls) <- colName
    } else if (nrow(simuls)>1L) simuls <- t(simuls)
    if(inherits(simuls,"numeric")) simuls <- matrix(simuls) ## if scalar summ stat.
    #colnames(simuls) <- stats
    if (is.null(colnames(simuls))) stop("The 'Simulate' function must provide names for the returned statistics.")
    attr(simuls,"par") <- par ## 1-row data frame (see comment in project.character())
    return(simuls)
  }
  ##################################################################################################"
  if (length(param_cl) > 1L) { ## ~ .run_cores, but object -> parsTable, list element -> ii index, and using param_cl
    if ( ! is.null(cl_seed) ) {
      ori <- RNGkind("L'Ecuyer-CMRG") # clusterSetRNGStream() re-calls it
      set.seed(cl_seed)
    }
    if (cluster_args$type=="FORK") {
      #if (is.null(mc.silent <- control$mc.silent)) mc.silent <- TRUE 
      mc.silent <- TRUE # no 'control' argument in contrast to spaMM::dopar
      #if (is.null(mc.preschedule <- control$mc.preschedule)) mc.preschedule <- TRUE 
      mc.preschedule <- TRUE # no 'control' argument in contrast to spaMM::dopar   (___F I X M E___?)
      has_progressr <- ("package:progressr" %in% search())
      seq_nr <- seq_len(nrow(parsTable))
      if (has_progressr) {
        # progressor is the only progress function that 'works' with mclapply
        # although not with load-balancing (mc.preschedule=FALSE)
        # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
        prog_fn <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
        with_fn <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
        with_fn({
          p <- prog_fn(steps=ceiling(nrow(parsTable)/nb_cores))
          p_simfn_per_par <- function(ii) {
            res <- simfn_per_par(ii)
            p() # p() call necessary for actual progress report 
            res
          }
          gridsimuls <- try(
            parallel::mclapply(seq_nr, FUN = p_simfn_per_par, mc.silent=mc.silent, mc.cores=nb_cores,
                               mc.preschedule = mc.preschedule)
          )
        })
      } else {
        .warn_once_progressr()
        gridsimuls <- try(
          parallel::mclapply(seq_nr, FUN = simfn_per_par, mc.silent=mc.silent,  mc.cores=nb_cores,
                             mc.preschedule = mc.preschedule)
        )
      }
      bootreps <- do.call(rbind,gridsimuls)
    } else { # PSOCK
      if (cores_info$has_doSNOW) {
        # loading (?) the namespace of 'snow' changes the *parent* RNG state (as well as sons' ones)! so we save and restore it 
        R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
        rdS_fn <- get("registerDoSNOW", asNamespace("doSNOW"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
        do.call(rdS_fn,list(cl=cl)) # this is what makes foreach see it and perform parallel computations
        assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
        if ( ! is.null(cl_seed) ) parallel::clusterSetRNGStream(cl = cl, cl_seed) 
        #
        show_pb <- (verbose && ! isTRUE(getOption('knitr.in.progress')))
        if (show_pb) {
          pb <- txtProgressBar(max = nrow(parsTable), style = 3, char=pb_char) ## doSNOW => 'long' progress bar
          progress <- function(n) setTxtProgressBar(pb, n)
          parallel::clusterExport(cl=param_cl, c("Simulate","progress") ,envir=environment()) ## slow! why?
          .options.snow <- list(progress = progress)
        } else .options.snow <- NULL
        ii <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'ii'
        foreach_args <- list(
          ii = seq_len(nrow(parsTable)), 
          .packages= packages,
          .options.snow = .options.snow,
          .inorder = TRUE, .errorhandling = #"remove"
            "pass"## "pass" to see error messages
        )
        foreach_blob <- do.call(foreach::foreach,foreach_args)
        gridsimuls <- foreach::`%dopar%`(foreach_blob, simfn_per_par(ii))
        if (show_pb) close(pb)
      } else {
        if ( ! is.null(cl_seed) ) parallel::clusterSetRNGStream(cl = cl, cl_seed) 
        if (verbose && ! isTRUE(getOption('knitr.in.progress'))) {
          pbopt <- pboptions(nout=min(100,2*nrow(parsTable)),type="timer", char=pb_char) ## pbapply:: 'short' progress bar
        } else pbopt <- pboptions(type="none")
        gridsimuls <- pbsapply(seq_len(nrow(parsTable)), simfn_per_par, simplify=FALSE, cl=param_cl) 
        pboptions(pbopt)
      }
    }
    if ( ! is.null(cl_seed) ) do.call("RNGkind", as.list(ori))
  } else { # serial
    if (verbose && ! isTRUE(getOption('knitr.in.progress'))) {
      pbopt <- pboptions(nout=min(100,2*nrow(parsTable)),type="timer", char=pb_char)
    } else pbopt <- pboptions(type="none")
    gridsimuls <- pbsapply(seq_len(nrow(parsTable)), simfn_per_par, simplify=FALSE, cl=param_cl) 
    pboptions(pbopt)
  } # a list of simulated distributions is returned in all cases.
  .close_cores(cores_info)
  gridsimuls
} 

add_simulation <- function(simulations=NULL, Simulate, parsTable=par.grid, par.grid=NULL,
                           nRealizations=Infusion.getOption("nRealizations"),
                           newsimuls=NULL, verbose=interactive(), nb_cores=NULL, packages=NULL,env=NULL,
                           control.Simulate=NULL,
                           cluster_args=list(), cl_seed=NULL, 
                           ...) { ## fn for ABC like simulation
  add_reftable(reftable=simulations, Simulate=Simulate, parsTable=parsTable,
               nRealizations=nRealizations,                          
               newsimuls=newsimuls, verbose=verbose, nb_cores=nb_cores, packages=packages,env=env,
               control.Simulate=control.Simulate,
               cluster_args=cluster_args, cl_seed=cl_seed, ...)
} 

# handling the various Simulate() formals
.wrap_Simulate <- function(Simulate, parsTable, reftable=NULL, 
                           nRealizations=1L, verbose=interactive(), nb_cores=NULL, packages=NULL,env=NULL,
                           control.Simulate=NULL,
                           cluster_args=list(), cl_seed=NULL) {
  if ("parsTable" %in% names(formals(Simulate))) {
    Simulate_input <- "parsTable"
    arglist <- control.Simulate
    arglist$parsTable <- parsTable
    gridsimuls <- do.call(Simulate, arglist)
    gridsimuls <- cbind(parsTable, gridsimuls)
    reftable <- rbind(reftable, gridsimuls)
  } else {
    Simulate_input <- "vector_as_args"
    gridsimuls <- .simulate_by_row(Simulate=Simulate, parsTable=parsTable,
                                   nRealizations=nRealizations, control.Simulate=control.Simulate, 
                                   verbose=verbose, # controls the progress bar... 
                                   nb_cores=nb_cores, packages=packages, env=env,
                                   cluster_args=cluster_args, cl_seed=cl_seed)
    if (nRealizations>1) {
      reftable <- c(reftable,gridsimuls) 
    } else { ## nRealizations=1 (~ABC reference table)
      gridsimuls <- do.call(rbind,gridsimuls) # matrix
      gridsimuls <- cbind(parsTable,gridsimuls) # data.frame bc called with parsTable data frame
      reftable <- rbind(reftable,gridsimuls)
    }
  }
  reftable # data frame (in up to date workflow)
}


add_reftable <- function(reftable=NULL, Simulate, parsTable=par.grid, par.grid=NULL,
                         nRealizations=1L,
                         newsimuls=NULL, verbose=interactive(), nb_cores=NULL, packages=NULL,env=NULL,
                         control.Simulate=NULL,
                         cluster_args=list(), cl_seed=NULL, 
                         constr_crits=NULL,
                         ...) {
  if (is.null(control.Simulate)) control.Simulate <- list(...)
  if ( ! is.null(reftable)) {
    lowersList <- list(old=attr(reftable,"LOWER"))  ## that gives the parameter names...
    uppersList <- list(old=attr(reftable,"UPPER"))  
    if (is.null(parsTable) && is.null(lowersList$old)) {
      stop("No information available to determine the parameter names; add a 'LOWER' attribute to the 'reftable' or 'simulations'")
    }
  } else lowersList <- uppersList <- list()
  if ( ! is.null(parsTable)) {
    if ( ! inherits(parsTable,"data.frame")) stop("'parsTable' argument is not a data.frame")
    #if ( ! inherits(Simulate,"character")) stop("'Simulate' must be a character string")
    
    if ( ! is.null(constr_crits)) { # Possible selection of points that satisfy parameter constraints. HOWEVER
      # The doc explains that it is better to use constr_crits when calling init_reftable() than 
      # add_reftable() [similarly, refine() use them at the .rparam step]
      constrs <- apply(parsTable,1L, function(v) all(eval(constr_crits, envir = as.list(v))<0))
      parsTable <- parsTable[constrs,, drop=FALSE]
    }

    if ( inherits(Simulate,"character")) Simulate <- eval(parse(text=Simulate))   # __F I X M E__ or get() ?
    reftable <- .wrap_Simulate(Simulate, control.Simulate=control.Simulate, parsTable=parsTable, 
                               reftable=reftable, nRealizations=nRealizations, 
                               verbose=verbose, # controls the progress bar... 
                               nb_cores=nb_cores, packages=packages, env=env,
                               cluster_args=cluster_args, cl_seed=cl_seed)
    chk_LOWER <- sapply(parsTable,min)
    lowersList$pargrid <- attr(parsTable,"LOWER") ## not necessarily present
    if (is.null(lowersList$pargrid)) {
      lowersList$pargrid <- chk_LOWER
    } else if (any(chk_LOWER<lowersList$pargrid)) warning(
      "The reftable appears not to satisfy its 'LOWER' attribute.\n Inefficiencies and confusion may result.",
      immediate. = TRUE)
    
    chk_UPPER <- sapply(parsTable,max)
    uppersList$pargrid <- attr(parsTable,"UPPER")
    if (is.null(uppersList$pargrid)) {
      uppersList$pargrid <- chk_UPPER
    } else if (any(chk_UPPER>uppersList$pargrid)) warning(
      "The reftable appears not to satisfy its 'UPPER' attribute.\n Inefficiencies and confusion may result.",
      immediate. = TRUE)
  }
    
  if ( ! is.null(newsimuls)) { ## user input 'newsimuls'
    reftable <- rbind(reftable,newsimuls)
    if (inherits(newsimuls,"list")) {
      newsimuls_pars <- do.call(rbind,lapply(newsimuls,attr,which="par")) ## allows the following check 
      if (is.null(newsimuls_pars)) {
        stop("'par' attribute appears to be missing from each of the new simulations.")
      }
      lowersList$new <- attr(newsimuls,"LOWER") 
      if (is.null(lowersList$new)) lowersList$new <- sapply(newsimuls_pars,min)
      uppersList$new <- attr(newsimuls,"UPPER")
      if (is.null(uppersList$new)) uppersList$new <- sapply(newsimuls_pars,max)
    } else { ## also for SLik_j: otherwise the parNames info cannot be certain...
      lowersList$newLOWER <- LOWER <- attr(newsimuls,"LOWER") 
      if (is.null(LOWER)) stop('Please provide attr(newsimuls,"LOWER")')
      if (is.null(names(LOWER))) stop('attr(newsimuls,"LOWER") should be a named vector (and UPPER too).')
      lowersList$new <- sapply(newsimuls[ ,names(LOWER)],min) # values in LOWER are ingored; names are used
      uppersList$newUPPER <- attr(newsimuls,"UPPER")
      #if (is.null(uppersList$newUPPER)) stop('Please provide attr(newsimuls,"UPPER")') #not necessary
      lowersList$new <- sapply(newsimuls[ ,names(LOWER)],max)
    }
  }
  attr(reftable,"LOWER") <- do.call(pmin,lowersList) # pmin over successive parameters across lowersList members ($pargrid, $new...)
  attr(reftable,"UPPER") <- do.call(pmax,uppersList)
  if ( ! missing(Simulate)) attr(reftable,"Simulate") <- Simulate
  attr(reftable,"control.Simulate") <- control.Simulate
  # attr(reftable,"Simulate_input") <- Simulate_input
  attr(reftable,"packages") <- packages
  attr(reftable,"env") <- env ## an environment!
  attr(reftable,"workflow_env") <- list2env(list(cl_seed=cl_seed)) # cl_seed created once, *not* modified afterwards
  if (nRealizations==1L) { 
    if (inherits(chk <- try(data.matrix(reftable), silent=TRUE),"try-error")) {
      condmess <- attr(chk,"condition")$message
      warning(paste(
        cli::style_underline("Format of reference table is suspect and likely to lead to a later error;\n"),
        cli::style_underline("   trying to convert the table to a matrix led to the error\n   "),
        cli::style_inverse(condmess),
        cli::style_underline("\n    It may be worth checking the simulation program and simulation arguments.")
      ),
      immediate. = TRUE)
    }
    class(reftable) <- c("reftable",class(reftable))
  } else if ( ! inherits(reftable,"EDFlist")) class(reftable) <- c("EDFlist",class(reftable))
  return(reftable)
}  





