\name{project.character}
\alias{project}
\alias{project.character}
\alias{project.default}
\alias{get_projector}
\alias{get_projection}
\alias{project.default}
\alias{neuralNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Learn a projection method for statistics and applies it 
}
\description{
\code{project} is a generic function with two methods. If the first argument is a parameter name, 
\code{project.character} (alias: \code{get_projector}) defines a projection function from several statistics to an output statistic predicting  
this parameter. \code{project.default}  (alias: \code{get_projection}) produces a vector of projected statistics using such a projection. \code{project} is particulary useful to reduce a large number of summary statistics to a vector of projected summary statistics, with as many elements as parameters to infer. This dimension reduction can substantially speed up subsequent computations. 
The concept implemented in \code{project} is to fit a parameter to the various statistics available, using machine-learning or mixed-model prediction methods. All such methods can be seen as nonlinear projection to a one-dimensional space. \code{project.character} is an interface that allows different projection methods to be used, provided they return an object of a class that has a defined \code{predict} method with a \code{newdata} argument (as expected, see \code{\link[stats]{predict}}).     
}
\usage{
project(x,...)

## S3 method for building the projection 
\method{project}{character}(x, stats, data, 
             trainingsize= eval(Infusion.getOption("projTrainingSize")),
             knotnbr=  eval(Infusion.getOption("knotnbr")),
             method,methodArgs=list(),verbose=TRUE,...)
get_projector(...) # alias for project.character

## S3 method for applying the projection
\method{project}{default}(x, projectors,...)
get_projection(...) # alias for project.default
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The name of the parameter to be predicted, or a vector/matrix/list of matrices of summary statistics.  
}
  \item{stats}{
  Statistics from which the predictor is to be predicted 
}
  \item{data}{
  A list of simulated empirical distributions, as produced by \code{\link{add_simulation}}, or a data frame with all required variables.
}
  \item{trainingsize}{
  For REML only: size of random sample of realizations from the \code{data} from which the smoothing parameters are estimated.  
}
  \item{knotnbr}{
  Size of random sample of realizations from the \code{data} from which the predictor is built given the smoothing parameters.
}
 \item{method}{
  character string: \code{"REML"}, \code{"GCV"}, or the name of a suitable projection function. The latter may be defined in another package, e.g. \code{"ranger"} or \code{"randomForest"}, or predefined by \verb{Infusion} (function \code{"nnetwrap"}), or defined by the user. See Details for predefined functions and for defining new ones. The default method is \code{"randomForest"} if this package is installed, and \code{"REML"} otherwise. Defaults may change in later versions, so it is advised to provide an explicit \code{method} to improve reproducibility. 
 }
 \item{methodArgs}{
   A list of arguments for the projection method. One may not need to provide arguments in the following cases, where 
   \code{project} kindly (tries to) assign values to the required arguments if they are absent from \code{methodArgs}: 
   
   If \code{"REML"} or \code{"GCV"} methods are used (in which case \code{methodArgs} is completely ignored); or
   
   if the projection method uses \code{formula} and \code{data} arguments (in particular if the formula is of the form 
   \code{response ~ var1 + var2 + ...}; otherwise the formula should be provided through \code{methodArgs}). This works 
   for example for methods based on \code{nnet}; or
   
   if the projection method uses \code{x} and \code{y} arguments. This works for example for \code{randomForest} (though not with the generic function \code{method="randomForest"}, but only with the internal function 
   \code{method="randomForest:::randomForest.default"}).
 }
  \item{projectors}{
  A list with elements of the form \code{<name>=<project result>}, where the \code{<name>} must differ from any name of \code{x}. 
  \code{<project result>} may indeed be the return object of a \code{project} call. 
}
  \item{verbose}{
Whether to print some information or not. In particular, \code{TRUE}, true-vs.-predicted diagnostic plots will be drawn for projection methods \dQuote{known} by Infusion (notably \code{"ranger"}, \code{"keras::keras.engine.training.Model"}, \code{"randomForest"}, \code{"GCV"}, \code{caret::train}).
}
  \item{\dots}{
Further arguments passed to or from other functions.
}
}
\details{
  Prediction can be based on a linear mixed model (LMM) with autocorrelated random effects,
  internally calling the \code{\link[spaMM]{corrHLfit}} function with formula
  \code{<parameter> ~ 1+ Matern(1|<stat1>+...+<statn>)}. This approach allows in principle to produce arbitrarily 
  complex predictors (given sufficient input) and avoids overfitting in the same way as restricted likelihood methods avoids overfitting in LMM. REML methods are then used by default to estimate the smoothing parameters. However, faster methods may be required, and random forests by \code{ranger} is the default method, if that package is installed (\code{randomForest} can also be used). Some other methods are possible (see below).    
  
  Machine learning methods such as random forests overfit, \emph{except if} out-of-bag predictions are used. When they are not, the bias is manifest in the fact that using the same simulation table for learning the projectors and for other steps of the analyses tend to lead to too narrow confidence regions. This bias disappears over iterations of \code{\link{refine}} when the projectors are kept constant. \code{Infusion} avoid this bias by using out-of-bag predictions when relevant, when \code{ranger} and \code{randomForest} are used. But it provides no code handing that problem for other machine-learning methods. Then, users should cope with that problems, and at a minimum should not update projectors in every iteration (the \dQuote{\href{https://kimura.univ-montp2.fr/~rousset/Infusion/InfusionIntro.pdf}{Gentle Introduction to Infusion} may contain further information about this problem}).       
  
  To keep REML computation reasonably fast, the \code{trainingsize} and \code{knotnbr} arguments determine respectively the size of the subset used to estimate the smoothing parameters and the size of the subset defining the predictor given the smoothing parameters.  REML fitting is already slow for data sets of this size (particularly as the number of predictor variables increase).
  
  If \code{method="GCV"}, a generalized cross-validation procedure (Golub et al. 1979) is used to estimate the smoothing parameters. This is faster but still slow, so a random subset of size \code{knotnbr} is still used to estimate the smoothing parameters and generate the predictor.
  
  Alternatively, various machine-learning methods can be used (see e.g. Hastie et al., 2009,  for an introduction). A random subset of size \code{knotnbr} is again used, with a larger default value bearing the assumption that these methods are faster. Predefined methods include \code{"ranger"}, \code{"randomForest"}, and \code{method="neuralNet"} which interfaces a neural network method, using  the \code{train} function from the \code{caret} package.
  
  In principle, any object suitable for prediction could be used as one of the \code{projectors}, and \code{Infusion} implements their usage so that in principle unforeseen projectors could be used. That is, if predictions of a parameter can be performed using an object \code{MyProjector} of class \code{MyProjectorClass}, 
  \code{MyProjector} could be used in place of a \code{project} result 
  if \code{predict.MyProjectorClass(object,newdata,...)} is defined. However, there is no guarantee that this will work on unforeseen projetion methods, as each method tends to have some syntactic idiosyncrasies. For example, if the learning method that generated the projector used
  a formula-data syntax, then its \code{predict} method is likely to request names for its \code{newdata}, that need to be provided through \code{attr(MyProjector,"stats")} (these names cannot be assumed to be in the \code{newdata} when \code{predict} is called through \code{optim}). 
  
}
\value{
  \code{project.character} returns an object of class returned by the \code{method} (methods \code{"REML"} and \code{"GCV"} will call \code{\link[spaMM]{corrHLfit}} which returns an object of class \code{spaMM})
  \code{project.default} returns an object of the same class and structure as the input \code{x}, containing the projected statistics inferred from the input summary statistics.
}
\references{
  Golub, G. H., Heath, M. and Wahba, G. (1979) Generalized Cross-Validation as a method for choosing a good ridge parameter. 
  Technometrics 21: 215-223.

  T. Hastie, R. Tibshirani, and J. Friedman. The Elements of Statistical Learning: Data Mining, Inference,
and Prediction. Springer, New York, 2nd edition, 2009.
}
\examples{
##########
if (Infusion.getOption("example_maxtime")>250) {
## Transform normal random deviates rnorm(,mu,sd)
## so that the mean of transformed sample is not sufficient for mu,
## and that variance of transformed sample is not sufficient for sd,
blurred <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  s <- exp(s/4)
  return(c(mean=mean(s),var=var(s)))
}

set.seed(123)
dSobs <- blurred(mu=4,s2=1,sample.size=20) ## stands for the actual data to be analyzed

## Sampling design as in canonical example 
parsp <- init_grid(lower=c(mu=2.8,s2=0.4,sample.size=20),
                      upper=c(mu=5.2,s2=2.4,sample.size=20))
# simulate distributions
dsimuls <- add_simulation(,Simulate="blurred", par.grid=parsp) 

## Use projection to construct better summary statistics for each each parameter 
mufit <- project("mu",stats=c("mean","var"),data=dsimuls)
s2fit <- project("s2",stats=c("mean","var"),data=dsimuls)

## additional plots for some projection method
if (inherits(mufit,"HLfit")) mapMM(mufit,map.asp=1,
  plot.title=title(main="prediction of normal mean",xlab="exp mean",ylab="exp var"))
if (inherits(s2fit,"HLfit")) mapMM(s2fit,map.asp=1,
  plot.title=title(main="prediction of normal var",xlab="exp mean",ylab="exp var"))

## apply projections on simulated statistics
corrSobs <- project(dSobs,projectors=list("MEAN"=mufit,"VAR"=s2fit))
corrSimuls <- project(dsimuls,projectors=list("MEAN"=mufit,"VAR"=s2fit))

## Analyze 'projected' data as any data (cf canonical example)
densb <- infer_logLs(corrSimuls,stat.obs=corrSobs) 
} else data(densb)
#########
if (Infusion.getOption("example_maxtime")>10) {
slik <- infer_surface(densb) ## infer a log-likelihood surface
slik <- MSL(slik) ## find the maximum of the log-likelihood surface
}
if (Infusion.getOption("example_maxtime")>500) {
slik <- refine(slik,10, update_projectors=TRUE) ## refine iteratively
}
}
