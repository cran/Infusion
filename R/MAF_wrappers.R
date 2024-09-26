.py_MAF_handle <- list2env(list(is_set=FALSE),parent = emptyenv())

.get_py_MAF_handle <- function(reset=FALSE, torch_device=.Infusion.data$options$torch_device,
                               GPU_mem=.Infusion.data$options$GPU_mem # , silent_ignore=FALSE
                               ) {
  if ( ! missing(torch_device)) .Infusion.data$options$torch_device <- torch_device
  # if (.py_MAF_handle$is_set && ! reset && ! silent_ignore) {
  #   warning("Python environment already exists. Use reset=TRUE to overwrite it\n (possibly losing previous results).")
  # }
  .py_MAF_handle <- mafR::get_py_MAF_handle(envir=.py_MAF_handle, 
                                            reset=reset, torch_device=torch_device,
                                            GPU_mem=GPU_mem)
  if ( (! .one_time_warnings$GPU_warned) && torch_device=="cpu" && ! reset) {
    message(paste("GPU will be ignored. Use config_mafR() to force its use."))
    .one_time_warnings$GPU_warned <- TRUE
  }
  .py_MAF_handle
}


config_mafR <- function(torch_device, 
                        #reset=FALSE, 
                        ...) {
  if (missing(torch_device)) torch_device <- .Infusion.data$options$torch_device
  if ( ( ! .one_time_warnings$GPU_warned ) && torch_device=="cpu" ) {
    warning("GPU will be ignored. See config_mafR(torch_device=.) to force its use.")
    .one_time_warnings$GPU_warned <- TRUE 
  }
  py_MAF_handle <- .get_py_MAF_handle(reset=TRUE, # reset, 
                                      torch_device=torch_device, ...)
  message("torch device set to '", py_MAF_handle$device$type,"'.")
  invisible(NULL)
}

# I should avoid using this fn ?
# .r_to_torch <- function(x, py_handle=.get_py_MAF_handle(), 
#                         device=py_handle$device) {
#   mafR::r_to_torch(x=x, py_handle=py_handle, device=device)
# }       

# A wild speculation. Allows much larger batch size than usually considered. 
.MAF_batchsize_mem_hack <- function(py_MAF_handle, trainsize) {
  if ( is.null(gpu_memory_B <- py_MAF_handle$gpu_memory)) {
    ram_kB <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))
    memory_B <- 1000 * ram_kB
  } else memory_B <- gpu_memory_B[2]
  batch_size <- 32L*min( max(1L, 
                             as.integer(memory_B/(trainsize*trainsize)),
                             as.integer(trainsize/ 320L) # => batch size ~trainsize/10 at most
                             ),
                         64L)
  batch_size # between 32L and min(2048L, O(trainsize/10))
} # MUST Be integer

# .MAF_validasize <- function(nr, ...) {as.integer(sqrt(nr))}

.design_MAF <- function(nr,nc, py_MAF_handle, transforms, 
                        design_fac, 
                        ... # to override elements of default return value
                        ) {
  hidden_units <- .Infusion.data$options$design_hidden_layers(nr=nr,nc=nc, transforms=transforms,
                                           design_fac=design_fac)
  patience <- .Infusion.data$options$MAF_patience
  valid_size <- .Infusion.data$options$MAF_validasize(nr)
  trainsize <- nr-valid_size
  batch_size <- .Infusion.data$options$MAF_batchsize(py_MAF_handle, trainsize)
  val_indices <- sample(nr, valid_size)
  train_indices <- seq_len(nr)
  train_indices[val_indices] <- NA_integer_
  train_indices <- na.omit(train_indices)
  design <- list(hidden_units=hidden_units, batch_size=batch_size,
       train_indices=train_indices, val_indices=val_indices,
       patience=patience)
  dotlist <- list(...)
  design[names(dotlist)] <- dotlist
  
  if (.is_devel_session()) {
    is_cpu <- py_MAF_handle$device$type=="cpu"
    if (is_cpu) {xpu <- "CPU"} else {xpu <- "GPU"}
    locmess <- with(design,
                    paste0("On ",xpu,": train size: ",length(train_indices),
                           "; valida.set size: ", length(val_indices),
                      "; batch size: ", batch_size,
                      "; hidden units: ",  paste0(hidden_units, collapse=","),
                      "; patience: ", patience, "\n")
    )
    cat(cli::col_green(locmess))
  }
  
  design
}

# the dots are passed to .design_MAF where they override elements in the returned list.
# so e.g passing patience through the dots overrides the default $options$MAF_patience 
.wrap_MAF <- function(data, inferredvars, statNames, verbose=TRUE,
                      py_MAF_handle=.get_py_MAF_handle(), 
                      transforms=.Infusion.data$options$MAF_auto_layers, 
                      design_fac=.Infusion.data$options$MAF_design_fac, 
                      Adam_learning_rate=.Infusion.data$options$Adam_learning_rate, 
                      ...) {
  allVars <- c(inferredvars, statNames)
  data <- as.matrix(data[,allVars, drop=FALSE])
  design <- .design_MAF(nr=nrow(data),nc=ncol(data), py_MAF_handle, transforms,
                        design_fac=design_fac,
                        ...)
  
  train <- data[design$train_indices, , drop=FALSE]
  train <- reticulate::r_to_py(train)
  valida_set <- data[design$val_indices, allVars, drop=FALSE] 
  valida_set <- reticulate::r_to_py(valida_set)
  
  captured <- switch(paste(as.integer(verbose)), "2"=c(), c("stdout"))
  time1 <- Sys.time()
  py_output <- reticulate::py_capture_output(
    MAF_object <- py_MAF_handle$MAF_density_estimation(y_train=train, 
                                                      y_test=valida_set, 
                                                      features=ncol(data), 
                                                      transforms=transforms,      
                                                      hidden_features=design$hidden_units, 
                                                      randperm=FALSE, 
                                                      max_epochs=500L, 
                                                      batch_size=design$batch_size,
                                                      device= .py_MAF_handle$device,
                                                      patience=design$patience,
                                                      learning_rate=Adam_learning_rate)
    , type = captured)
  attr(MAF_object, "train_time") <- round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1) ## spaMM:::.timerraw(time1)
  if (verbose==1L) {
    py_output <- strsplit(py_output, split="\n", fixed=TRUE)[[1]]
    cat(paste0("Output from MAF_density_estimation():\n",
                py_output[1],"\n  (...)\n", 
                paste(tail(py_output,3)[1:2],collapse="\n"),"\n") )
  }
  class(MAF_object) <- c("MAF", class(MAF_object))
  attr(MAF_object, "colTypes") <- list(allVars=allVars)
  MAF_object
}

# the dots are passed to .design_MAF where they override elements in the returned list.
# so e.g passing patience through the dots overridesthe default $options$MAF_patience 
.wrap_MAF_cond <- function(data, outVars, givenVars, verbose=TRUE,
                           py_MAF_handle=.get_py_MAF_handle(), 
                           transforms=.Infusion.data$options$MAF_auto_layers,
                           design_fac=.Infusion.data$options$MAF_design_fac, 
                           Adam_learning_rate=.Infusion.data$options$Adam_learning_rate, 
                           ...) {
  data <- as.matrix(data)
  design <- .design_MAF(nr=nrow(data),
                        nc=length(c(outVars,givenVars)), 
                        py_MAF_handle, transforms=transforms, 
                        design_fac=design_fac, ...)
  train_indices <- design$train_indices
  val_indices <- design$val_indices
  
  train_x <- as.matrix(data[train_indices, givenVars, drop=FALSE]) # conditions
  train_x <- reticulate::r_to_py(train_x)
  train_y <- as.matrix(data[train_indices, outVars, drop=FALSE])
  train_y <- reticulate::r_to_py(train_y)
  valida_set_x <- as.matrix(data[val_indices,givenVars, drop=FALSE])
  valida_set_x <- reticulate::r_to_py(valida_set_x)
  valida_set_y <- as.matrix(data[val_indices,outVars, drop=FALSE])
  valida_set_y <- reticulate::r_to_py(valida_set_y)
  
  captured <- switch(paste(as.integer(verbose)), "2"=c(), c("stdout"))
  time1 <- Sys.time()
  py_output <- reticulate::py_capture_output(
    MAF_object <- py_MAF_handle$MAF_conditional_density_estimation(train_y, train_x, valida_set_y, valida_set_x, 
                                                    features=length(outVars), 
                                                    context=length(givenVars), 
                                                    transforms=transforms, 
                                                    hidden_features=design$hidden_units, 
                                                    randperm=FALSE, 
                                                    max_epochs=500L, batch_size=design$batch_size,
                                                    device = py_MAF_handle$device,
                                                    patience=design$patience,
                                                    learning_rate=Adam_learning_rate)
    , type = captured)
  attr(MAF_object, "train_time") <- round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1) ## spaMM:::.timerraw(time1)
  
  if (verbose==1L) {
    py_output <- strsplit(py_output, split="\n", fixed=TRUE)[[1]]
    cat(paste0("Output from MAF_conditional_density_estimation():\n",
               py_output[1],"\n  (...)\n", 
               paste(tail(py_output,3)[1:2],collapse="\n"),"\n") )
  }
  
  class(MAF_object) <- c("MAF", class(MAF_object))
  attr(MAF_object, "colTypes") <- list(outVars=outVars, givenVars=givenVars)
  MAF_object
}

str.MAF <- function(object, ...) {
  # class(object) <- setdiff(class(object), "MAF") # to call other str() methods
  # str(object, ...)
  attrs <- attributes(object)
  extraAttrs <- setdiff(names(attrs),c("class","py_object"))
  cat("'MAF' object with attributes: ")
  str(attrs[c("py_object",extraAttrs,"class")])
}

.calc_all_MAFs <- function(data, statNames, latentVars, fittedPars, inferredVars, 
                           verbose, using="MAFmix") {
  locreftable <- na.omit(data)   
  if (nNAlines <- length(attr(locreftable,"na.action"))) {
    message(paste0(nNAlines," lines of reference table still contained NA's, and were removed. see ?handling_NAs for details."))
  }
  if (verbose$most) {
    if (length(grep("c.mafR",using))) {
      cat(paste0("Inferring densities [(S;P); (S,P); P instrumental; P|S posterior]: ",
                 nrow(data)," points\n"))
    } else cat(paste0("Inferring densities [joint(S,P); P instrumental; P|S posterior]: ",
                      nrow(data)," points\n"))
  }
  
  py_MAF_handle <- .get_py_MAF_handle()
  # reticulate::py_run_string("del complconddens") # no clear effect
  py_main <- py_MAF_handle$py_main  # Module(__main__)
  py_main$complconddens <- py_main$completedens <- py_main$pardens <- 
    py_main$postdens <- py_main$jointdens <- py_main$conddens <- NULL
  # py_MAF_handle$gc$collect()
  py_MAF_handle$torch$cuda$empty_cache() # _F I X M E_ ?
  
  if (length(grep("c.mafR",using))) {
    complconddens <- # pdf(S,L|P)
      .wrap_MAF_cond(locreftable, outVars=c(statNames,latentVars), 
                     givenVars=fittedPars, verbose=verbose$MAF,
                     py_MAF_handle=.get_py_MAF_handle()) 
  } else complconddens <- NULL
  py_main$complconddens <- complconddens
  
  completedens <- 
    py_main$completedens <- 
    .wrap_MAF(locreftable, inferredvars=inferredVars, statNames=statNames, verbose=verbose$MAF,
              py_MAF_handle=.get_py_MAF_handle())
  
  pardens <- 
    py_main$pardens <- 
    .wrap_MAF(locreftable, inferredvars=fittedPars, statNames=c(), verbose=verbose$MAF,
              py_MAF_handle=.get_py_MAF_handle())
  
  # postdens used by .inits_from_postdens() which is needed in maximizations... 
  # I really need a an object to simulate from, not some way to predict from the intrumental posterior,
  # so I cannot simply use a jointdens and a statdens.
  if ( ! length(grep("MAFmix",using))) {
    if (length(grep("u.mafR",using))) { # ___F I X M E___
      design_fac <- .Infusion.data$options$MAF_design_fac/2 # trying to minimize the conditional computation
      patience <- .Infusion.data$options$MAF_patience %/% 2L
    } else { # "c.mafR"
      design_fac <- .Infusion.data$options$MAF_design_fac # trying to minimize the conditional computation
      patience <- .Infusion.data$options$MAF_patience
    } 
    postdens <- 
      py_main$postdens <- # instrumental posterior  P|S
      .wrap_MAF_cond(locreftable, outVars=fittedPars, givenVars=statNames, verbose=verbose$MAF,
                     py_MAF_handle=.get_py_MAF_handle(),design_fac=design_fac,
                     # dots:
                     patience=patience)
    attr(postdens,"which") <- "I_postdens" 
  } else postdens <- NULL
  attr(pardens,"which") <- "pardens"
  if (length(latentVars)) { # then refit without the latent vars
    py_main$jointdens <- 
      .wrap_MAF(locreftable, inferredvars=fittedPars, statNames=statNames, verbose=verbose$MAF,
                py_MAF_handle=.get_py_MAF_handle())
    if (length(grep("c.mafR",using))) {
      conddens <- # pdf(S|P)
        .wrap_MAF_cond(locreftable, outVars=statNames, givenVars=fittedPars, verbose=verbose$MAF,
                       py_MAF_handle=.get_py_MAF_handle()) # !!! not tested with latentVars
      attr(complconddens,"which") <- "S,L|P"
      attr(conddens,"which") <- "S|P"
      attr(jointdens,"which") <- "jointdens"
      attr(completedens,"which") <- "completedens"
    } else conddens <- NULL
    py_main$conddens <- conddens
    resu <- list(jointdens=jointdens, conddens=conddens, completedens=completedens,
                 complconddens=complconddens, pardens=pardens, postdens=postdens)    
  } else {
    if (length(grep("c.mafR",using))) attr(complconddens,"which") <- "S|P"
    attr(completedens,"which") <- "jointdens"
    resu <- list(jointdens=completedens, conddens=complconddens,
                      pardens=pardens, postdens=postdens)  
  }
  py_MAF_handle$py_main <- py_main 
  resu
}


.predict_MAF <- function(X, # parameters only 
                         object, 
                         tstat.obs, # 1-row matrix as otherwise more cases should be considered for cbind'ing
                         # log, 
                         which, # "safe" version ignores, by correcting, spuriously high logL in area of low parameter density.
                         thr_info=.get_thr_info(object),
                         ...) {
  if (is.null(dim(X))) {
    # dim(X) <- c(1L, length(X)) drops dimnames (sigh)
    X <- t(X)
    if (which=="jointvaldens" ||
        length(grep("u.mafR|MAFmix",object$using))) newjointX <- cbind(X,tstat.obs) # cbind two 1-row matrices
  } else {
    if (length(intersect(colnames(X),colnames(tstat.obs)))) stop("'X' should contain only parameters, not summary statistics")
    if (which=="jointvaldens" ||
        length(grep("u.mafR|MAFmix",object$using))) newjointX <- cbind(X,tstat.obs[rep(1,nrow(X)),,drop=FALSE]) # cbind two nrow(X)-row matrices
  }
  nr <- nrow(X)
  
  # if (.py_MAF_handle$device$type=="cuda") {
  #   if (.py_MAF_handle$torch$cuda$memory_allocated() >
  #       100000000) {
  #     .py_MAF_handle$torch$cuda$empty_cache()
  #     .py_MAF_handle$gc$collect()
  #     }
  # }
  
  if (which %in% c("parvaldens","safe") ||
      length(grep("u.mafR|MAFmix",object$using))
      ) {
    X <- reticulate::r_to_py(as.matrix(X)) 
    parvaldens <- 
      .py_MAF_handle$MAF_predict_nocond(object$pardens, X, .py_MAF_handle$device)
    if (which=="parvaldens") return(parvaldens)
  }
  if (which=="jointvaldens" ||
      length(grep("u.mafR|MAFmix",object$using))) {
    newjointX <- reticulate::r_to_py(as.matrix(newjointX)) 
    jointvaldens <- 
      .py_MAF_handle$MAF_predict_nocond(object$jointdens, newjointX, .py_MAF_handle$device)
    return(jointvaldens)
  }
  # else "lik" or "safe"
  sobs <- reticulate::r_to_py(tstat.obs[rep(1,nr),, drop=FALSE]) 
  if (which=="lik") X <- reticulate::r_to_py(as.matrix(X)) 
  if (length(grep("u.mafR|MAFmix",object$using))) {
    condvaldens <- jointvaldens-parvaldens
  } else {
    condvaldens <- 
      .py_MAF_handle$MAF_predict_cond(object$conddens, 
                                      sobs, # 'Y', but I cannot name it
                                      X, # 'cond'... 
                                      .py_MAF_handle$device)
  }
  if (which=="safe") {
    ddens <- parvaldens-thr_info$thr_dpar
    negddens <- (ddens<0)
    if (any(negddens)) {
      dlogl <- condvaldens-thr_info$reft_maxlogl # perhaps by setting a sightly higher thr
      # high would allow a bit more extrapol and thus exploration? bc currently it may induce self-reinforcing 
      # high pardens in suboptimal regions.... but the effect seems opposite on exploration...
      posdlogl <- dlogl>0
      #lowextrapol <- negddens & ! posdlogl
      #condvaldens[lowextrapol] <- condvaldens[lowextrapol]+ ... # not clear what to do
      highextrapol <- negddens & posdlogl # ...................... potentially conservative intervals.............
      # make it continuous wrt to dlogl, but strongly compensating in most cases
      condvaldens[highextrapol] <- condvaldens[highextrapol]+ sqrt(dlogl[highextrapol])*ddens[highextrapol]
      # extrapol <-  posdlogl & ! negddens
      # condvaldens[extrapol] <- condvaldens[extrapol]- sqrt(dlogl[highextrapol])
      attr(condvaldens,"lowdens") <- negddens
    }
  }
  return(condvaldens)
  # if ( ! log) density <- exp(density)
  # return(density)
}

.simulate.MAF <- function(density, nsim, 
                          given, # =t(get_from(object,"stat.obs"))
                          object, # I need 'given' or 'object'
                          batchsize= 4000L
) {
  colTypes <- attr(density,"colTypes")
  which. <- attr(density,"which")
  # .py_MAF_handle$torch$cuda$empty_cache()
  if (which. %in% c("jointdens","pardens","completedens")) {
    trypoints <- 
      density()$sample(reticulate::tuple(as.integer(nsim)))$cpu()$numpy()
    colnames(trypoints) <- colTypes$allVars
  } else { # all conditional simulations
    givenVars <- colTypes$givenVars
    coln_given <- colnames(given) # expect 1-row *matrix*
    if (which.=="I_postdens") { # INSTRUMENTAL postdens
      if ( ! setequal(givenVars, coln_given)) {
        if (length(setdiff(coln_given,colTypes$outVars))) { 
          stop("Wrong 'given' value in .simulate.MAF()")
        } else {
          # patch for sampling instrumental postdens for given parameter values
          given <- t(get_from(object,"stat.obs"))
        }
      }
    } else if (which. %in% c("S|P","S,L|P")) { # then 'given' should be params, not stat.obs...
      if ( ! setequal(givenVars, coln_given)) {
        stop("Wrong 'given' value in .simulate.MAF()")
      }
    } else stop("'which.' value not handled in .simulate.MAF()" ) # might be feasible with MAFcond but...
    if (nrow(given)>1L) {
      stop(".simulate.MAF() handles only a single 'given' vector")
      # The Python batching code does not yet handle other cases. See further comments there.
    }
    given <- reticulate::r_to_py(given[,givenVars,drop=FALSE]) # controls order
    trypoints <- 
      .py_MAF_handle$MAF_simulate_cond(density, 
                                       reticulate::tuple(as.integer(nsim)),
                                       given, # cond
                                       .py_MAF_handle$device,
                                       batchsize
      )
    colnames(trypoints) <- colTypes$outVars
  }
  trypoints
}
