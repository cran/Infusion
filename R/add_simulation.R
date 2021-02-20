
add_reftable <- function(...) { ## fn for ABC like simulation
  resu <- add_simulation(...,nRealizations=1L)
  class(resu) <- c("reftable",class(resu))
  resu
} 

`[.reftable` <- function (x, i, j, 
                             drop = TRUE ## replicating default [.data.frame behaviour to allow selecting a single column as a vector, for plots etc 
) {
  class(x) <- "data.frame" ## avoids recursive call to `[.reftable
  resu <- x[i,j,drop=drop]
  if (is.data.frame(resu)) {
    attrx <- attributes(x)
    attrx <- attrx[setdiff(names(attrx), c("names","row.names","class"))]
    for (st in names(attrx)) attr(resu,st) <- attrx[[st]]
    parNames <- names(attr(resu, "LOWER"))
    parKepts <- intersect(parNames,colnames(resu))
    attr(resu, "LOWER") <- attr(resu, "LOWER")[parKepts]
    attr(resu, "UPPER") <- attr(resu, "UPPER")[parKepts]
    class(resu) <- c("reftable", class(resu))
  } ## else vector
  return(resu)
} # Use unlist() to remove attributes from the return value


add_simulation <- function(simulations=NULL, Simulate, par.grid=NULL,
                           nRealizations=NULL,
                           newsimuls=NULL, verbose=interactive(), nb_cores=NULL, packages=NULL,env=NULL,
                           control.Simulate=NULL,
                           cluster_args=list(),
                           ...) {
  if (is.null(control.Simulate)) control.Simulate <- list(...)
  old_nRealizations <- Infusion.getOption("nRealizations")
  if (is.null(nRealizations)) {
    nRealizations <- old_nRealizations
  } else {
    Infusion.options(nRealizations=nRealizations)
  }
  if ( ! is.null(simulations)) {
    lowersList <- list(old=attr(simulations,"LOWER"))  ## that gives the parameter names...
    uppersList <- list(old=attr(simulations,"UPPER"))  
  } else lowersList <- uppersList <- list()
  if ( ! is.null(par.grid)) {
    if ( ! inherits(par.grid,"data.frame")) stop("'par.grid' argument is not a data.frame")
    #if ( ! inherits(Simulate,"character")) stop("'Simulate' must be a character string")
    if ( inherits(Simulate,"character")) Simulate <- eval(parse(text=Simulate))   
    prevmsglength <- 0L
    nsim <- nrow(par.grid)
    gridsimuls <- list()
    if (is.null(cluster_args$spec)) cluster_args$spec <- nb_cores # which means that cluster_args$spec overrides nb_cores
    cores_info <- .init_cores(cluster_args=cluster_args)
    cl <- cores_info$cl
    if ( ! is.null(cl)) {
      parallel::clusterExport(cl, "packages",envir=environment()) ## passes the list of packages to load
      # Infusion not leaded in child process !
      #parallel::clusterExport(cl, list("nRealizations"),envir=environment()) 
      #parallel::clusterCall(cl,fun=Infusion.options,nRealizations=nRealizations)
      if ( ! is.null(env)) parallel::clusterExport(cl=cl, ls(env),envir=env)
    }
    nb_cores <- cores_info$nb_cores
    as_one <- identical(names(nRealizations),"as_one")
    if (as_one) { 
      which_cl <- "param"
    } else {
      which_cl <- names(nb_cores) ## the name comes from explicit user input; this is doc'ed in ?add_simulation
      if (is.null(which_cl)) {
        if (nRealizations>1L) {
          which_cl <- "replic"
        } else which_cl <- "param"
      }
    }
    if (which_cl=="replic") {
      repl_cl <- cl ## the one used by simfn_per_par()
      if (length(repl_cl) > 1L) {
        if (cores_info$has_doSNOW) {
          ## arguments such as Simulate are evaluated in the child envir so we pass definition under the "Simulate" name
          pb_char <- "N" # nested
        } else { ## arguments such as Simulate are evaluated in the parent envir to "myrnorm" so we pass def under "myrnorm" name
          pb_char <- "n" # nested
        }
        dotenv <- list2env(list(Simulate=Simulate) )
        parallel::clusterExport(cl=repl_cl, ls(dotenv),envir=dotenv) ## exports eg "myrnorm": works for pbapply:: by not doSNOW     
      } else pb_char <- "s"
      param_cl <- NULL
    } else {
      param_cl <- cl
      if (length(param_cl) > 1L) {
        if (cores_info$has_doSNOW) {
          pb_char <- "P"
        } else {
          pb_char <- "p"
        }
        dotenv <- list2env(list(Simulate=Simulate) )
        parallel::clusterExport(cl=param_cl, ls(dotenv),envir=dotenv) ## exports eg "myrnorm": works for pbapply:: by not doSNOW     
      } else pb_char <- "s"
      repl_cl <- NULL 
    }
    #
    ############################################# simfn_per_par ######################################"
    simfn_per_par <- function(ii) {
      par <- par.grid[ii,,drop=FALSE] ## 1-row data frame: treated as list by do.call() or c(); but:
      parlist <- c(par,control.Simulate) ## now parlist is a list, not a data.frame
      #
      if (length(repl_cl)>1L  && cores_info$has_doSNOW) { ## for the balancing only since we don't want the progress bar
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
    if (length(param_cl) > 1L) { ## ~ .run_cores, but object -> par.grid, list element -> ii index, and using param_cl
      if (cores_info$has_doSNOW) {
        show_pb <- (verbose && ! isTRUE(getOption('knitr.in.progress')))
        if (show_pb) {
          pb <- txtProgressBar(max = nrow(par.grid), style = 3, char=pb_char) ## doSNOW => 'long' progress bar
          progress <- function(n) setTxtProgressBar(pb, n)
          parallel::clusterExport(cl=param_cl, c("Simulate","progress") ,envir=environment()) ## slow! why?
          .options.snow <- list(progress = progress)
        } else .options.snow <- NULL
        ii <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'ii'
        foreach_args <- list(
          ii = seq_len(nrow(par.grid)), 
          .packages= packages,
          .options.snow = .options.snow,
          .inorder = TRUE, .errorhandling = #"remove"
                                          "pass"## "pass" to see error messages
        )
        foreach_blob <- do.call(foreach::foreach,foreach_args)
        gridsimuls <- foreach::`%dopar%`(foreach_blob, simfn_per_par(ii))
        if (show_pb) close(pb)
      } else {
        if (verbose && ! isTRUE(getOption('knitr.in.progress'))) {
          pbopt <- pboptions(nout=min(100,2*nrow(par.grid)),type="timer", char=pb_char) ## pbapply:: 'short' progress bar
        } else pbopt <- pboptions(type="none")
        gridsimuls <- pbsapply(seq_len(nrow(par.grid)), simfn_per_par, simplify=FALSE, cl=param_cl) 
        pboptions(pbopt)
      }
    } else { 
      if (verbose && ! isTRUE(getOption('knitr.in.progress'))) {
        pbopt <- pboptions(nout=min(100,2*nrow(par.grid)),type="timer", char=pb_char)
      } else pbopt <- pboptions(type="none")
      gridsimuls <- pbsapply(seq_len(nrow(par.grid)), simfn_per_par, simplify=FALSE, cl=param_cl) 
      pboptions(pbopt)
    } # a list of simulated distributions is returned in all cases.
    .close_cores(cores_info)
    if (nRealizations>1) {
      simulations <- c(simulations,gridsimuls) 
    } else { ## nRealizations=1 (~ABC reference table)
      gridsimuls <- do.call(rbind,gridsimuls)
      gridsimuls <- cbind(par.grid,gridsimuls)
      simulations <- rbind(simulations,gridsimuls)
    }
    if (verbose && prevmsglength>0L) cat("\n")
    lowersList$pargrid <- attr(par.grid,"LOWER") ## not necessarily present
    if (is.null(lowersList$pargrid)) lowersList$pargrid <- apply(par.grid,2,min)
    uppersList$pargrid <- attr(par.grid,"UPPER")
    if (is.null(uppersList$pargrid)) uppersList$pargrid <- apply(par.grid,2,max)
  } 
  if ( ! is.null(newsimuls)) { ## user input 'newsimuls'
    simulations <- rbind(simulations,newsimuls)
    if (inherits(newsimuls,"list")) {
      newsimuls_pars <- do.call(rbind,lapply(newsimuls,attr,which="par")) ## allows the following check 
      if (is.null(newsimuls_pars)) {
        stop("'par' attribute appears to be missing from each of the new simulations.")
      }
      lowersList$new <- attr(newsimuls,"LOWER") 
      if (is.null(lowersList$new)) lowersList$new <- apply(newsimuls_pars,2,min)
      uppersList$new <- attr(newsimuls,"UPPER")
      if (is.null(uppersList$new)) uppersList$new <- apply(newsimuls_pars,2,max)
    } else { ## also for SLik_j: otherwise the parNames info cannot be certain...
      lowersList$newLOWER <- LOWER <- attr(newsimuls,"LOWER") 
      if (is.null(LOWER)) stop('Please provide attr(newsimuls,"LOWER")')
      lowersList$new <- apply(newsimuls[,names(LOWER)],2L,min) # values in LOWER are ingored; names are used
      uppersList$newUPPER <- attr(newsimuls,"UPPER")
      #if (is.null(uppersList$newUPPER)) stop('Please provide attr(newsimuls,"UPPER")') #not necessary
      lowersList$new <- apply(newsimuls[,names(LOWER)],2L,max)
    }
  }
  attr(simulations,"LOWER") <- do.call(pmin,lowersList) # pmin over sucessive parameters across lowersList members ($pargrid, $new...)
  attr(simulations,"UPPER") <- do.call(pmax,uppersList)
  if ( ! missing(Simulate)) attr(simulations,"Simulate") <- Simulate
  attr(simulations,"control.Simulate") <- control.Simulate
  attr(simulations,"packages") <- packages
  attr(simulations,"env") <- env ## an environment!
  if ( nRealizations>1 && ! inherits(simulations,"EDFlist")) class(simulations) <- c("EDFlist",class(simulations))
  Infusion.options(nRealizations=old_nRealizations)
  return(simulations)
}  





