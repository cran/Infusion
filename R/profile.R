.Infusion.data <- new.env(parent = emptyenv())
.Infusion.data$Constants <- list(Version = NA)
.Infusion.data$options <- list(## should be documented: 
                               LRthreshold= - qchisq(0.999,df=1)/2, 
                               nRealizations=1000,
                               seq_nbCluster= function(projdata, nr=nrow(projdata)) {seq(ceiling(nr^0.3))}, # M O D I F Y doc if this is modified ! 
                               maxnbCluster= function(projdata, nr= nrow(projdata), nc=ncol(projdata)) { # 'projdata' to emphasize that these are not the user-level data (but rather, projected ones)
                                 nr_corr <- (nr*2L) %/% 3L #4 clusters for 200 rows, 6 cols => 171 parameters sometimes fails. 
                                 #                                      So we correct nrow adhocly so that it appear to be less than 171
                                 floor((nr_corr+1L)/(1L+nc*(nc+1L))) # compar nrow to param of free model ((gd+1)*gd+1)*G-1 (cov mat*G + means*G + G-1 proportions)
                               },
                               gof_nstats_fn=function(nr,nstats) floor(nr^(1/3)),  # heuristically balancing gdim and max nbCluster
                               mixturing="Rmixmod", # alternatives: "xLLiM", "mclust"
                               #infer_logL_method="infer_logL_by_Rmixmod", # "infer_logL_by_mclust", ## string for clusterExport !
                               mixmodGaussianModel="Gaussian_pk_Lk_Ck", # all free vs # "Gaussian_pk_Lk_Dk_A_Dk", # shape is constant but volume and orientation are free
                               mclustModel="VVV", ## "VEV", ## equivalent to default mixmodGaussianModel
                               precision=0.1,
                               #
                               train_cP_size=quote(.train_cP_size_fn(method_string,stats)),
                               trainingsize=quote(.trainingsize_fn(method_string,stats)),
                               projKnotNbr=1000,
                               use_oob=TRUE, # for .predictWrap()
                               nodesize=5L, ## that the default for regression by ranger and randomForest # F I X M E see comments where this is used
                               trainsize=quote(nraw/log(nraw)), ## for UNUSED code in update_projectors block 
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
                               EDFestLevelName="EDFestLevel",
                               HLmethod="PQL/L", ## fast; REPQL is fast also; ML et REML notably slower
                               useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), ## does not seem to be used
                               CIweight=1, ## used in rparam to enhance (or not) sampling near CI bounds
                               rparamfn="rparam", ## allows choice of function
                               # Rmixmod controls
                               criterion="AIC", # criterio for selection of number of clusters (not BIC by default!)
                               #strategy=quote(.do_call_wrap("mixmodStrategy",arglist=list(),pack="Rmixmod")), ## not doc'ed...
                               strategy=quote(Rmixmod::mixmodStrategy()), ## given Rmixmod in Suggests
                               mixmodSeed=123,
                               clu_optimizer="",
                               # trypoints controls
                               samplingType=c(default=1,posterior=0),
                               expansion=1,
                               maxeval=quote(10^(3+(log(length(initvec))-log(5))/log(4))), # nloptr; *modified for bobyqa (which recommends > 10 * npar^2)
                               SLRTopts=c() # c("gamma") to activate gamma correction
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
                        "\nType 'help(Infusion)' for a short introduction,",
                        "\nand see https://gitlab.mbb.univ-montp2.fr/francois/Infusion",
                        "\nfor more documentation.")
}


".onLoad" <- function (lib, pkg) {
  .Infusion.data$Constants$Version <- utils::packageVersion("Infusion")
  #abyss <- suppressMessages(delaunayn(matrix(1,nrow=2,ncol=1))) # *sigh*
}  
