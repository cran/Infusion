# default value of option 'mixturing_errorfn' ; also used in .rparam...
# Use Infusion.options(is_devel_session=TRUE) to force it.
.is_devel_session <- local({  
  return_value <- NULL
  function(opt=.Infusion.data$options$is_devel_session) { # 'opt' allows global control, or local non-API control
    if (is.null(opt)) {
      if (is.null(return_value)) {  
        knitr_i_p <- isTRUE(getOption('knitr.in.progress')) # assuming that the editor starts a new R session when calling knitr. 
        onMyMachine <- (Sys.getenv("_LOCAL_TESTS_")=="TRUE")
        if (onMyMachine && ! knitr_i_p) {
          user <- Sys.info()["user"]
          onMyMachine <- (user=="francois.rousset")
          if (onMyMachine) {
            if (inInfusionProject <- .inRstudio(bool=TRUE)) {
              activeproject <- .inRstudio(bool=FALSE) 
              inInfusionProject <- ( ! is.null(activeproject) &&
                                       tail(strsplit(activeproject,"/")[[1]],1)=="Infusion")
              if ( ! inInfusionProject) {
                mess <- paste("\nSys.info()['user']='francois.rousset', but not in an Infusion devel session.\n", 
                               "Set   Infusion.options(is_devel_session=TRUE)\n", 
                               "to activate private messages, and to call the browser on .densityMixmod() errors.\n")
                cat(cli::col_green(mess))
              }
            } # else we may be running the checks... '___F I X M E___' it would be nice to handle that case
              # An ex-tempo may be to add Infusion.options(is_devel_session=TRUE) in test-all.R
          } else if (length(grep("francois|rousset", user))) {
            ## 'onMyMachine' is currently FALSE and 'inInfusionProject' is yet undefined
            if (inInfusionProject <- .inRstudio(bool=TRUE)) {
              activeproject <- .inRstudio(bool=FALSE) 
              inInfusionProject <- ( ! is.null(activeproject) &&
                                       tail(strsplit(activeproject,"/")[[1]],1)=="Infusion")
              if ( inInfusionProject) {
                mess <- paste("\nIt looks like we're in an Infusion devel session,\n", 
                              "and user name matches to grep('francois|rousset', .) but is not 'francois.rousset',\n", 
                              "Use   Infusion.options(is_devel_session=TRUE)\n", 
                              "to activate private messages, and to call the browser on .densityMixmod() errors.\n")
                cat(cli::col_green(mess))
              }
            } 
          }
        } 
        return_value <<- onMyMachine && ! knitr_i_p && inInfusionProject
      } 
      return_value
    } else return(opt)
  }
})

.one_time_warnings <- list2env(list(GPU_warned=FALSE,
                                    caret_warned=FALSE,
                                    update_warned=FALSE,
                                    progressr_warned=FALSE,
                                    nb_cores_warned=FALSE,
                                    cores_avail_warned=FALSE,
                                    doSNOW_warned=FALSE,
                                    RF_warned=FALSE,
                                    nThr_warned=FALSE),
                               parent = emptyenv())

.Infusion.data <- new.env(parent = emptyenv())
.Infusion.data$Constants <- list(Version = NA)
.Infusion.data$options <- list(
  ## should be documented: 
  LRthreshold= - qchisq(0.999,df=1)/2, 
  nRealizations=1000,
  seq_nbCluster= function(projdata, nr=nrow(projdata)) {seq(ceiling(nr^0.31))}, # M O D I F Y doc if this is modified ! 
  #### 
  ## Cf get_workflow_design for the older version(s) of maxnbCluster()
  ## Full MGM model has P = (nc*(nc+3)/2 +1)G-1) params 
  ## so G ~ (nr+1)/(nc*(nc+3)/2 +1) for a saturated model P ~ nr
  ## We compare P to 4 nr to set the maximum G such that P ~ nr/4:
  ## so G ~ (nr/4+1)/(nc*(nc+3)/2 +1) = (nr+4)/(2*nc*(nc+3) +4)
  ## test: MVNcovmat_identif/summaries_list.1_200.reclu15.v2.1.186.2 (in v185- P ~ nr/6)
  ## argument 'projdata' to emphasize that this is not the row reftable
  maxnbCluster= function(projdata, nr= nrow(projdata), nc=ncol(projdata)) { 
    (nr+4L)%/%(nc*(nc+3L)*2L+4L) 
  },
  ## Corresponds to zuko-like:
  design_hidden_layers = .design_hidden_layers_MGM_like,
  MAF_batchsize = function(...) 100L, # PSM19 used 100, (2017) used 50
  # MAF_batchsize = .MAF_batchsize_mem_hack, 
  MAF_patience=30, # (2017,p.7) # PSM 2019 used 20
  MAF_auto_layers=3L,
  Adam_learning_rate=1e-3, # but PSM used 1e-4 
  ##
  MAF_validasize = function(nr, ...) {nr %/% 20L}, 
  MAF_design_fac=1L,
  #
  gof_nstats_fn=function(nr,nstats) floor(nr^(1/3)),  # heuristically balancing gdim and max nbCluster
  mixturing="Rmixmod", # alternatives: "xLLiM", "mclust"
  #infer_logL_method="infer_logL_by_Rmixmod", # "infer_logL_by_mclust", ## string for clusterExport !
  mixmodGaussianModel="Gaussian_pk_Lk_Ck", # all free vs # "Gaussian_pk_Lk_Dk_A_Dk", # shape is constant but volume and orientation are free
  mclustModel="VVV", ## "VEV", ## equivalent to default mixmodGaussianModel
  # Rmixmod controls
  criterion="AIC", # criterio for selection of number of clusters (not BIC by default!)
  #strategy=quote(.do_call_wrap("mixmodStrategy",arglist=list(),pack="Rmixmod")), ## not doc'ed...
  #strategy=quote(Rmixmod::mixmodStrategy()), ## given Rmixmod in Suggests
  global_strategy_args=list(nbIterationInAlgo=1000L), # list of args for Rmixmmod::mixmodStrategy()
  get_mixModstrategy=.get_mixModstrategy, # default fn is internal fn... must return a Strategy object
  mixmodSeed=123,
  #
  precision=0.25, # threshold for stopping iterations
  #
  train_cP_size=quote(.train_cP_size_fn(method_string,stats)),
  trainingsize=quote(.trainingsize_fn(method_string,stats)),
  projKnotNbr=1000,
  use_oob=TRUE, # for .predictWrap()
  target_size_factor=1000L, # scaling factor for selection of points for projection
  upd_proj_subrows_thr=40000L,
  subrows_overlap_thr=0.9, # control of re-projection. 0.9 is value used in MS (____F I X M E____ rethink?)
  proj_methodArgs=list(num.trees=1000L), # session-level defaults 
  #      1000 trees is Infusion::project() default if missing here; but set here for ref by abcrf...
  #
  repeat_stat_thr=5L, ## For data autopsy
  #
  example_maxtime=2.5,
  ## partly doc'ed and does not need more
  nb_cores=NULL,
  ## undocumented
  constrOptim=FALSE,
  binningExponent=0.5,
  zeromargin=0.1,
  logLname="logL",
  tailNames=c("qsup","qinf"),
  EDFestLevelName="EDFestLevel", # tailp code
  HLmethod="PQL/L", ## fast; REPQL is fast also; ML et REML notably slower
  useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), ## does not seem to be used
  CIweight=1, ## used in rparam to enhance (or not) sampling near CI bounds
  # Allows choice of function:
  rparamFn=.rparam_SLik_j_in_out, ## for SLik_j
  rparamFn_SLik="rparam", ## for primitive SLik
  # .rparam controls
  freq_runif_weights_thr=25, # control of adj_target_size: reduce it to allow smaller subblocks
  default_ceil_size_fac=10, # control of target_ceil_size
  freq_full_rejection_bnd=0.01, # control of ceil_size 
  expansion=1,
  target_LR_min_value=quote(qchisq(level,df=2)/2), # eval() in contexts where 'level' and (speculatively) 'npar' are available.
  safe_prof4rparam=TRUE, # FALSE allows a specific use of nlminb (but not very useful)
  penal_facs=c(pardens=2, # ___F I X M E__ topic for eternal rethink but instr_dens 0<.<1 seems best.
               instr_dens=0.5), # 0.25 not clearly better
  exploration_fac=1/3, # controls w_l/w_u correction in .calc_filltop_weights()
  norm_or_t=.wrap_rmvt, # debug(mvtnorm::rmvt) to see that it is used.
  #
  thr_info_ctrl=0L, # quote(log(length(logls)))
  #
  maxeval=quote(10^(3+(log(length(initvec))-log(5))/log(4))), # nloptr; *modified for bobyqa (which recommends > 10 * npar^2)
  SLRTopts=c(), # c("gamma") to activate gamma fit of bootstrap distribution of LR
  RMSE_nsim=10L, # for the RMSE bootstrap
  max_base=Inf, # control a number of points in .safe_init()
  chol_error=identity, # default is thus to do nothing... was previously warning, but not helpful so far
  #
  torch_device="cpu", 
  is_devel_session=NULL, # controls .is_devel_session()
  mixturing_errorfn=.is_devel_session # here a function. Return value controls what to do
  # in case of problem (here controls whether to dump frames)
)

Infusion.options <- function(...) {
  if (nargs() == 0) return(.Infusion.data$options)
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.Infusion.data$options[arg]),  ## return here for eg ... = "NUMAX"
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(.Infusion.data$options)
  argnames <- names(temp)
  if (is.null(argnames)) stop("options must be given by name")
  old <- .Infusion.data$options[argnames]
  names(old) <- argnames ## bc names are not valid for previously absent elements
  .Infusion.data$options[argnames] <- temp
  invisible(old)
}

Infusion.getOption <- function (x) {Infusion.options(x)[[1]]}


".onAttach" <- function (lib, pkg) {
  version <- utils::packageVersion("Infusion")
  packageStartupMessage("Infusion (version ", version, 
                        ") is loaded.", 
                        "\nType 'help(\"Infusion\")' for a short introduction,",
                        "\nand see https://gitlab.mbb.univ-montp2.fr/francois/Infusion",
                        "\nfor more documentation.")
}


".onLoad" <- function (lib, pkg) {
  .Infusion.data$Constants$Version <- utils::packageVersion("Infusion")
}  

.safe_tempdir <- function(temp_dir, sep = NULL, sub = NULL) {
  if (is.null(sep)) {
    if (.Platform$OS.type == "windows")
      sep <- "\\"
    else
      sep <- "/"
  }
  if (is.null(temp_dir))
    temp_dir <- tempdir()
  temp_dir <- normalizePath(temp_dir, sep)
  #temp_dir <- gsub("\\","/", temp_dir, fixed=TRUE)
  if (substring(temp_dir, nchar(temp_dir)) != sep) {
    temp_dir <- paste0(temp_dir, sep)
  }
  if (!is.null(sub))
    temp_dir <- paste0(temp_dir, sub, sep)
  if (!dir.exists(temp_dir))
    dir.create(temp_dir)
  temp_dir
} 

