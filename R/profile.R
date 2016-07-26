.Infusion.data <- new.env(parent = emptyenv())
.Infusion.data$Constants <- list(Version = NA)
.Infusion.data$options <- list(## should be documented: 
                               LRthreshold= - qchisq(0.999,df=1)/2, 
                               nRealizations=1000,
                               nbCluster= quote(seq(ceiling(nrow(data)^0.3))), # M O D I F Y doc if this is modified !
                               #nbCluster= quote(ceiling(nrow(data)^0.3)), 
                               mixmodGaussianModel="Gaussian_pk_Lk_Dk_A_Dk",
                               precision=0.1,
                               projTrainingSize=200,
                               projKnotNbr=300,
                               example_maxtime=2.5,
                               ## undocumented
                               binningExponent=0.5,
                               zeromargin=0.1,
                               logLname="logL",
                               tailNames=c("qsup","qinf"),
                               EDFestLevelName="EDFestLevel",
                               HLmethod="PQL/L", ## fast; REPQL is fast also; ML et REML notably slower
                               useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), ## does not seem to be used
                               CIweight=1, ## used in rparam to enhance (or not) sampling near CI bounds
                               rparamfn="rparam", ## allows choice of function
                               fitmeCondition=quote(TRUE)
                               )


Infusion.options <- function(...) {
  if (nargs() == 0) return(.Infusion.data$options)
  current <- .Infusion.data$options
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.Infusion.data$options[arg]),  ## return here for eg ... = "NUMAX"
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  current[n] <- temp
  .Infusion.data$options <- current
  invisible(current)
}

Infusion.getOption <- function (x) {Infusion.options(x)[[1]]}


".onAttach" <- function (lib, pkg) {
  version <- utils::packageVersion("Infusion")
  packageStartupMessage("Infusion (version ", version, 
                        ") is loaded.", 
                        "\nType 'help(Infusion)' for a short introduction,",
                        "\nand see http://kimura.univ-montp2.fr/~rousset/Infusion.htm",
                        "\nfor more documentation.")
}


".onLoad" <- function (lib, pkg) {
  .Infusion.data$Constants$Version <- utils::packageVersion("Infusion")
  abyss <- suppressMessages(delaunayn(matrix(1,nrow=2,ncol=1))) # *sigh*
}  
