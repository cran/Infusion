.Infusion.data <- new.env(parent = emptyenv())
.Infusion.data$Constants <- list(Version = NA)
.Infusion.data$options <- list(## should be documented: 
                               LRthreshold= - qchisq(0.999,df=1)/2, 
                               nRealizations=1000,
                               nbCluster= quote(seq(ceiling(nrow(data)^0.3))), # M O D I F Y doc if this is modified ! # note that functiosn such as goftest handle only this format, not a list
                               mixturing="Rmixmod",
                               #using="Rmixmod", ## or ""mclust
                               #infer_logL_method="infer_logL_by_Rmixmod", # "infer_logL_by_mclust", ## string for clusterExport !
                               mixmodGaussianModel="Gaussian_pk_Lk_Ck", # all free vs # "Gaussian_pk_Lk_Dk_A_Dk", # shape is constant but volume and orientation are free
                               mclustModel="VVV", ## "VEV", ## equivalent to default mixmodGaussianModel
                               precision=0.1,
                               #
                               projTrainingSize=quote(.trainsize_fn(method,stats)),
                               knotnbr=quote(.knotnbr_fn(method,stats)),
                               projKnotNbr=1000,
                               oob=TRUE, # for .predictWrap()
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
                               fitmeCondition=quote(TRUE),
                               # Rmixmod controls
                               criterion="AIC", # criterio for selection of number of clusters (not BIC by default!)
                               #strategy=quote(.do_call_wrap("mixmodStrategy",arglist=list(),pack="Rmixmod")), ## not doc'ed...
                               strategy=quote(Rmixmod::mixmodStrategy()), ## given Rmixmod in Suggests
                               mixmodSeed=123,
                               clu_optimizer="",
                               # trypoints controls
                               samplingType=c(default=1,posterior=0),
                               expansion=1,
                               maxeval=quote(10^(3+(log(length(initvec))-log(5))/log(4))) # nloptr; *modified for bobyqa (which recommends > 10 * npar^2)
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
                        "\nand see https://kimura.univ-montp2.fr/~rousset/Infusion.htm",
                        "\nfor more documentation.")
}


".onLoad" <- function (lib, pkg) {
  .Infusion.data$Constants$Version <- utils::packageVersion("Infusion")
  #abyss <- suppressMessages(delaunayn(matrix(1,nrow=2,ncol=1))) # *sigh*
}  
