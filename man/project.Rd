\name{project.character}
\alias{project}
\alias{project.character}
\alias{project.default}
\alias{neuralNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Learn a projection method for statistics and applies it 
}
\description{
\code{project} is a generic function with two methods. If the first argument is a parameter name, 
\code{project.character} defines a projection function from several statistics to an output statistic predicting  
this parameter. \code{project.default} produces a vector of projected statistics using such a projection. \code{project} is particulary useful to reduce a large number of summary statistics to a vector of projected summary statistics, with as many elements as parameters to infer. This dimension reduction can substantially speed up subsequent computations. 
The concept implemented in \code{project} is to fit a parameter to the various statistics available, using machine-learning or mixed-model prediction methods. All such methods can be seen as nonlinear projection to a one-dimensional space. \code{project.character} is an interface that allows different projection methods to be used, provided they return an object of a class that has a defined \code{predict} method with a \code{newdata} argument (as expected, see \code{\link[stats]{predict}}).     
}
\usage{
project(x,...)

## S3 method for building the projection 
\method{project}{character}(x, stats, data, 
             trainingsize= if (method=="REML") 
               {Infusion.getOption("projTrainingSize")} else {NULL},
             knotnbr= if (method \%in\% c("REML","GCV"))  {
                        Infusion.getOption("projKnotNbr")
                      } else {floor(1000*log2(length(stats)+1))},
             method="REML",methodArgs=list(),verbose=TRUE,...)

## S3 method for applying the projection
\method{project}{default}(x, projectors,...)
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
  character string: \code{"REML"}, \code{"GCV"}, or the name of a suitable projection function. The latter may be defined in another package, e.g. \code{"randomForest"}, or predefined by \verb{Infusion} (function \code{"nnetwrap"}), or defined by the user. See Details for predefined functions and for defining new ones.
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
Whether to print some information or not. In particular, \code{TRUE}, true-vs.-predicted diagnostic plots will be drawn if any of the following methods have been used: \code{"REML"}, \code{"GCV"}, or a call to \code{caret::train}.
}
  \item{\dots}{
further arguments passed to or from other methods (currently not used).
}
}
\details{
  Prediction can be based on a linear mixed model (LMM) with autocorrelated random effects,
  internally calling the \code{\link[spaMM]{corrHLfit}} function with formula
  \code{<parameter> ~ 1+ Matern(1|<stat1>+...+<statn>)}. This approach allows in principle to produce arbitrarily 
  complex predictors (given sufficient input) and avoids overfitting in the same way as restricted likelihood methods avoids overfitting in LMM. REML methods are then used by default to estimate the smoothing parameters. However, faster methods may be required, and method \code{"neuralNet"} interfaces a neural network approach.    
  
  The data may involve hundreds of thousands of realizations of the summary statistic, and REML fitting is already slow for much smaller data sets, which is why faster alternative methods may be worth considering, and why random subset(s) of the data may be considered at various steps. The default size of these subsets aim to ensure that the computations can be performed in reasonable time. 
  
  For REML, the \code{trainingsize} and \code{knotnbr} arguments determine respectively the size of the subset used to estimate the smoothing parameters and the size of the subset defining the predictor given the smooothing parameters.
  
  If \code{method="GCV"}, a generalized cross-validation procedure (Golub et al. 1979) is used to estimate the smoothing parameters. This is faster but still slow, so a random subset of size \code{knotnbr} is still used to estimate the smoothing parameters and generate the predictor.
  
  Alternatively, various machine-learning methods can be used (see e.g. Hastie et al., 2009,  for an introduction). A random subset of size \code{knotnbr} is again used, with a larger default value bearing the assumption that these methods are faster. \code{method="neuralNet"} interfaces a neural network method. It calls the \code{\link[caret]{train}} function from the \code{caret} package.
  
  In principle, any object suitable for prediction could be used as one of the \code{projectors}. 
  That is, if predictions of a parameter can be performed using an object \code{MyProjector} of class \code{MyProjectorClass}, 
  \code{MyProjector} could be used in place of a \code{project} result 
  if \code{predict.MyProjectorClass(object,newdata,...)} is defined. However, if the learning method that generated the projector used
  a formula-data syntax, then its \code{predict} method is likely to request names for its \code{newdata}, that need to be provided through
  \code{attr(MyProjector,"stats")} (these names cannot be assumed to be in the \code{newdata} when \code{predict} is called through \code{optim}).
  
}
\value{
  \code{project.character} returns an object of class returned by the \code{method} (methods \code{"REML"} and \code{"GCV"} will call \code{\link[spaMM]{corrHLfit}} which return an object of class \code{spaMM})
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

## plots
mapMM(mufit,map.asp=1,
  plot.title=title(main="prediction of normal mean",xlab="exp mean",ylab="exp var"))
mapMM(s2fit,map.asp=1,
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
slik <- refine(slik,10) ## refine iteratively
}
}
