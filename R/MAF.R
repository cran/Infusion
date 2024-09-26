
save_MAFs <- function(object, ext="_MAF.pkl", prefix="") {
  if (inherits(object,"SLik_j")) {
    for (st in grep("dens",names(object),value=TRUE)) {
      if (inherits(object[[st]],"MAF")) {
        reticulate::py_save_object(object[[st]], 
                                   file=paste0(prefix,st,ext))
      }
    }
    
    if ( ! is.environment(load_MAFs_info <- object$load_MAFs_info)) {
      object$load_MAFs_info <- 
        list2env(list(ext=ext,prefix=prefix,pwd=getwd()), parent = emptyenv())
      warning("Old-style fit object: the return value of load_MAFs() must be kept.",
              immediate. = TRUE)
    } else {
      load_MAFs_info$ext <- ext
      load_MAFs_info$prefix <- prefix
      load_MAFs_info$pwd <- getwd()
    }
  }
  object
}

load_MAFs <- function(object, ext="_MAF.pkl", prefix="", set_path_only=FALSE) {
  if (inherits(object,"SLik_j")) {
    for (st in grep("dens",names(object), value=TRUE)) {
      if (inherits(object[[st]],"MAF") && 
          reticulate::py_is_null_xptr(attr(object[[st]],"py_object"))
          ) {
        MAF_file <- paste0(prefix,st,ext)
        if ( ! file.exists(MAF_file)) stop(paste(MAF_file, "needed but not available from pwd =\n",getwd()))
        if ( ! set_path_only ) attr(object[[st]],"py_object") <- reticulate::py_load_object(file=MAF_file)
      }
    }
    if (set_path_only || is.null(object$load_MAFs_info$pwd)) {
      if ( ! is.environment(load_MAFs_info <- object$load_MAFs_info)) {
        object$load_MAFs_info <- 
          list2env(list(ext=ext,prefix=prefix,pwd=getwd()), parent = emptyenv())
      } else {
        load_MAFs_info$ext <- ext
        load_MAFs_info$prefix <- prefix
        load_MAFs_info$pwd <- getwd()
      }
    } 
  }
  object # the MAF objects will be accessible only through the return object as they are not in an envir
}

.design_hidden_layers_MGM_like <- function(
    projdata, nr= nrow(projdata), nc=ncol(projdata),
    transforms, # autoregressive layers: the K ~ 5 of PPM17; cf MAF_auto_layers for default value
    n_hid_layers=2L, # The L ~ 2 of PPM17 ans PSM19
    design_fac, 
    ...) { 
  # too high value => very poor fit => not good candidates when sampling pardens
  # occurred with nr= ~ 6000, nc=4, hidden_units= rep(3000L,3)
  n1 <- (Infusion.getOption("maxnbCluster"))(nr=nr,nc=nc)
  n2 <-  nr^0.31
  nbClu <- max(1L, as.integer(min(n1,n2))) # effectively the # of clusters of MGM fit
  npMGM <-  ((nc*(nc+3L))%/%2L+1L)*nbClu-1L # exact # of param of MGM fit
  n_units_per_layer <- (- 3*nc + sqrt(9*nc*nc+8*npMGM/transforms) )/2 # O(npMGM)
  n_units_per_layer <- n_units_per_layer*design_fac # allows ad-hoc adjustment
  n_units_per_layer <- 2^as.integer(3+ # so that ~ npMAF is power of two,  O(8 to 16 times npMGM)
                                      log(n_units_per_layer,base=2))
  n_units_per_layer <- as.integer(n_units_per_layer +1/(100*npMGM)) # safe rounding
  resu <- rep(n_units_per_layer,n_hid_layers) # zuko.flows.MAF()'s hidden_features argument:
  ## => the number of elements gives the number of hidden layers.
  
  ## Info: approximate number of parameters according to Table 1 SM PapamakariokPM17:
  npMAF <- (3L*transforms*nc*n_units_per_layer + 
              transforms*(n_hid_layers-1)*n_units_per_layer*n_units_per_layer
  )/2L
  attr(resu,"info") <- c(npMGM=npMGM, npMAF= npMAF)
  
  resu
}

.design_hidden_layers_PSM19 <- function(
    projdata, nc=ncol(projdata),
    transforms, # autoregressive layers: the K ~ 5 of PPM17; cf MAF_auto_layers for default value
    n_hid_layers=2L, # The L ~ 2 of PPM17 ans PSM19
    n_units_per_layer=50L,
    ...
) { 
  resu <- rep(n_units_per_layer, n_hid_layers) 
  ## Info: approximate number of parameters according to Table 1 SM PapamakariokPM17:
  npMAF <- (3L*transforms*nc*n_units_per_layer + 
              transforms*(n_hid_layers-1)*n_units_per_layer*n_units_per_layer
  )/2L
  attr(resu,"info") <- c(npMGM=NA_integer_, npMAF= npMAF)
  resu
}

MAF.options <- function(template="zuko-like", ...) {
  if (template=="zuko-like") {
    optns <- list(
      design_hidden_layers = .design_hidden_layers_MGM_like,
      MAF_patience=30, # (2017,p.7) # PSM 2019 used 20
      MAF_auto_layers=3L,
      Adam_learning_rate=1e-3 # but PSM used 1e-4 
    )
  } else if (template=="PSM19") {
    optns <- list( 
      design_hidden_layers= .design_hidden_layers_PSM19,
      MAF_patience=20,
      MAF_auto_layers=5L,
      Adam_learning_rate=1e-4 
    )
  } else if (is.list(template)) {
    optns <- template
  } else if (is.null(template)) {
    optns <- list()
  } else stop("Unknown 'template'")
  dotlist <- list(...)
  for (st in names(dotlist)) optns[names(dotlist)] <- dotlist
  Infusion.options(optns)
}
