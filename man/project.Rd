\name{project.character}
\alias{project}
\alias{project.character}
\alias{project.default}
\alias{get_projector}
\alias{get_projection}
\alias{project.default}
\alias{deforest_projectors}
\alias{neuralNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Learn a projection method for statistics and apply it 
}
\description{
\code{project} is a generic function with two methods. If the first argument is a parameter name, 
\code{project.character} (alias: \code{get_projector}) defines a projection function from several statistics to an output statistic predicting  
this parameter. \code{project.default}  (alias: \code{get_projection}) produces a vector of projected statistics using such a projection. \code{project} is particularly useful to reduce a large number of summary statistics to a vector of projected summary statistics, with as many elements as parameters to infer. This dimension reduction can substantially speed up subsequent computations. 
The concept implemented in \code{project} is to fit a parameter to the various statistics available, using machine-learning or mixed-model prediction methods. All such methods can be seen as nonlinear projection to a one-dimensional space. \code{project.character} is an interface that allows different projection methods to be used, provided they return an object of a class that has a defined \code{predict} method with a \code{newdata} argument.

\code{deforest_projectors} is an utility to reduce the saved size of objects containing \pkg{ranger} objects (\code{\link{reproject}} can be used to reverse this). 
}
\usage{
project(x,...)

## S3 method for building the projection 
\method{project}{character}(x, stats, data, 
             trainingsize=  eval(Infusion.getOption("trainingsize")),
             train_cP_size= eval(Infusion.getOption("train_cP_size")), method, 
             methodArgs=eval(Infusion.getOption("proj_methodArgs")), 
             verbose=TRUE, keep_data=TRUE, ...)
get_projector(...) # alias for project.character

## S3 method for applying the projection
\method{project}{default}(x, projectors, use_oob=Infusion.getOption("use_oob"), 
                          is_trainset=FALSE, methodArgs=list(), ext_projdata, ...)
get_projection(...) # alias for project.default

##
deforest_projectors(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The name of the parameter to be predicted, or a vector/matrix/list of matrices of summary statistics.  
}
  \item{stats}{
  Statistics from which the predictor is to be predicted 
}
  \item{use_oob}{
  Boolean: whether to use out-of-bag predictions for data used in the training set, when such oob predictions are available (i.e. for random forest methods). Default as controlled by the same-named package option, is TRUE. This by default involves a costly check on each row of the input \code{x}, whether it belongs to the training set, so it is better to set it to FALSE if you are sure \code{x} does not belong to the training set (for true data in particular). Alternatively the check can be bypassed if you are sure that \code{x} was used as the training set.
}
\item{is_trainset}{
  Boolean. In a \code{project} call, set it to TRUE if \code{x} was used as the training set, to bypass a costly check (see \code{use_oob} argument). The same logic may apply in a \code{plot_proj} call, except that it is not immediately obvious for users whether the full reference table in an \code{object} was used as the training set, so trying to save time by setting \code{is_trainset=TRUE} requires more insight.
}
  \item{data}{
  A list of simulated empirical distributions, as produced by \code{\link{add_simulation}}, or a data frame with all required variables.
}
  \item{trainingsize, train_cP_size}{ Integers;
  for most projection methods (excluding \code{"REML"} but including \code{"ranger"}) only  \code{trainingsize} is taken into account: it gives the maximum size of the training set (and is infinite by default for \code{"ranger"} method). If the \code{data} have more rows the training set is randomly sampled from it. For the \code{"REML"} method, \code{train_cP_size} is the maximum size of the data used for estimation of smoothing parameters, and \code{trainingsize} is the maximum size of the data from which the predictor is built given the smoothing parameters.
}
 \item{method}{
  character string: \code{"REML"}, \code{"GCV"}, or the name of a suitable projection function. The latter may be defined in another package, e.g. \code{"ranger"} or \code{"randomForest"}, or predefined by \verb{Infusion}, or defined by the user. See Details for predefined functions and for defining new ones. The default method is \code{"ranger"} if this package is installed, and \code{"REML"} otherwise. Defaults may change in later versions, so it is advised to provide an explicit \code{method} to improve reproducibility. 
 }
 \item{methodArgs}{
   A list of arguments for the projection method. For \code{project.character}, the \code{ranger} method is run with some default argument if no \code{methodArgs} are specified. Beware that a NULL \code{methodArgs$splitrule} is interpreted as the \code{"extratrees"} \code{splitrule}, whereas in a direct call to \code{ranger}, this would be interpreted as the \code{"variance"} \code{splitrule}. For \code{project.default}, the only \code{methodArgs} element handled is \code{num.threads} passed to \code{predict.ranger} (which can also be controlled globally by \code{Infusion.options(nb_cores=.)}).   
   
   For other methods, \code{project} kindly (tries to) assign values to the required arguments if they are absent from \code{methodArgs}, according to the following rules:
   
   If \code{"REML"} or \code{"GCV"} methods are used (in which case \code{methodArgs} is completely ignored); or
   
   if the projection method uses \code{formula} and \code{data} arguments (in particular if the formula is of the form 
   \code{response ~ var1 + var2 + ...}; otherwise the formula should be provided through \code{methodArgs}). This works 
   for example for methods based on \code{nnet}; or
   
   if the projection method uses \code{x} and \code{y} arguments. This works for example for the (somewhat obsolete) method \code{randomForest} (though not with the generic function \code{method="randomForest"}, but only with the internal function 
   \code{method="randomForest:::randomForest.default"}).
 }
  \item{projectors}{
  A list with elements of two possible forms: (1) \code{<name>=<project result>}, where the \code{<name>} must differ from any name of \code{x} and \code{<project result>} is the return object of a \code{project} call; or (2) \code{<name>=NULL} where \code{<name>} is the name of a variable (raw summary statistic) in \code{x} (such explicit NULLs are \bold{needed} for any raw statistic to be retained in the projected data; see Value). 
}
  \item{verbose}{
Whether to print some information or not. In particular, \code{TRUE}, true-vs.-predicted diagnostic plots will be drawn for projection methods \dQuote{known} by Infusion (notably \code{"ranger"}, \code{"fastai.tabular.learner.TabularLearner"}, \code{"keras::keras.engine.training.Model"}, \code{"randomForest"}, \code{"GCV"}, \code{caret::train}).
}
  \item{keep_data, ext_projdata}{
(experimental, and only when \pkg{ranger} is used). Setting \code{keep_data=FALSE} allows the input data to be removed from the return object of \code{project.character} (where they are otherwise part of its \code{call} element). This may be useful to save memory when multiple projections are based on the same data. However, as this data information is sometimes used, it must then be manually added as element \code{projdata} to the return value of \code{infer_SLik_joint}, and provided to \code{project.default} calls through the \code{ext_projdata} argument.
}
  \item{object}{
An object of class \code{SLik_j}.
}
  \item{\dots}{
Further arguments passed to or from other functions. Currently, they are passed to \code{plot}.
}
}
\details{
  The preferred \code{project} method is non-parametric regression by (variants of) the random forest method as implemented in \pkg{ranger}. It is the default method, if that package is installed. Alternative methods have been interfaced as detailed below, but the functionality of most interfaces is infrequently tested. 
  
  By default, the ranger call through \code{project} will use the split rule \code{"extratrees"}, with some other controls also differing from the \pkg{ranger} package defaults. If the split rule \code{"variance"} is used, the default value of \code{mtry} used in the call is also distinct from the \pkg{ranger} default, but consistent with Breiman 2001 for regression tasks. 

  Machine learning methods such as random forests overfit, \emph{except if} out-of-bag predictions are used. When they are not, the bias is manifest in the fact that using the same simulation table for learning the projectors and for other steps of the analyses tend to lead to too narrow confidence regions. This bias disappears over iterations of \code{\link{refine}} when the projectors are kept constant. \code{Infusion} avoid this bias by using out-of-bag predictions when relevant, when \code{ranger} and \code{randomForest} are used. But it provides no code handling that problem for other machine-learning methods. Then, users should cope with that problems, and at a minimum should not update projectors in every iteration (the \dQuote{\href{https://gitlab.mbb.univ-montp2.fr/francois/Infusion/-/blob/master/documents/InfusionIntro.pdf}{Gentle Introduction to Infusion} may contain further information about this problem}).       

  Prediction can be based on a linear mixed model (LMM) with autocorrelated random effects,
  internally calling the \code{\link[spaMM]{fitme}} function with formula
  \code{<parameter> ~ 1+ Matern(1|<stat1>+...+<statn>)}. This approach allows in principle to produce arbitrarily 
  complex predictors (given sufficient input) and avoids overfitting in the same way as restricted likelihood methods avoids overfitting in LMM. REML methods are then used by default to estimate the smoothing parameters. However, faster methods are generally required.
  
  To keep REML computation reasonably fast, the \code{train_cP_size} and \code{trainingsize} arguments determine respectively the size of the subset used to estimate the smoothing parameters and the size of the subset defining the predictor given the smoothing parameters.  REML fitting is already slow for data sets of this size (particularly as the number of predictor variables increase).
  
  If \code{method="GCV"}, a generalized cross-validation procedure (Golub et al. 1979) is used to estimate the smoothing parameters. This is faster but still slow, so a random subset of size \code{trainingsize} is still used to estimate the smoothing parameters and generate the predictor.
  
  Alternatively, various machine-learning methods can be used (see e.g. Hastie et al., 2009,  for an introduction). A random subset of size \code{trainingsize} is again used, with a larger default value bearing the assumption that these methods are faster. Predefined methods include
  \itemize{
\item \code{"ranger"}, the default, a computationally efficient implementation of random forest;
\item \code{"randomForest"}, the older default, probably obsolete now;
\item \code{"neuralNet"}, a neural network method, using  the \code{train} function from the \code{caret} package (probably obsolete too);
\item \code{"fastai"} deep learning using the \code{fastai} package;
\item \code{"keras"} deep learning using the \code{keras} package.
}
The last two interfaces may yet offer limited or undocumented control: using deep learning seems attractive but the benefits over \code{"ranger"} are not clear (notably, the latter provide out-of-bag predictions that avoid overfitting).   

  In principle, any object suitable for prediction could be used as one of the \code{projectors}, and \code{Infusion} implements their usage so that in principle unforeseen projectors could be used. That is, if predictions of a parameter can be performed using an object \code{MyProjector} of class \code{MyProjectorClass}, 
  \code{MyProjector} could be used in place of a \code{project} result 
  if \code{predict.MyProjectorClass(object,newdata,...)} is defined. However, there is no guarantee that this will work on unforeseen projection methods, as each method tends to have some syntactic idiosyncrasies. For example, if the learning method that generated the projector used
  a formula-data syntax, then its \code{predict} method is likely to request names for its \code{newdata}, that need to be provided through \code{attr(MyProjector,"stats")} (these names cannot be assumed to be in the \code{newdata} when \code{predict} is called through \code{optim}). 
  
}
\value{
  \code{project.character} returns an object of the class returned by the called \code{method} (by default, a \code{ranger} object for the up-to-date workflow).
  
  \code{project.default} returns an object of the same class and structure as the input \code{x}, containing the variables named in the \code{projectors} argument, each variable being a projected statistics inferred from the input summary statistics, or a summary statistic copied from the input \code{x} (if an \bold{explicit} NULL projector was included for this statistic in the \code{projectors} argument).
  
  \code{deforest_projectors} is used for its side effect  (the contents of an environment within the input object been modifed), and returns a character string emphasizing this feature.
}
\references{
  Breiman, L. (2001). Random forests. Mach Learn, 45:5-32. <doi:10.1023/A:1010933404324>

  Golub, G. H., Heath, M. and Wahba, G. (1979) Generalized Cross-Validation as a method for choosing a good ridge parameter. 
  Technometrics 21: 215-223.

  T. Hastie, R. Tibshirani, and J. Friedman. The Elements of Statistical Learning: Data Mining, Inference,
and Prediction. Springer, New York, 2nd edition, 2009.
}
\seealso{\code{\link{plot_proj}} and \code{\link{plot_importance}} for diagnostic plots of the projections.}
\note{See workflow examples in \code{\link{example_reftable}} and \code{\link{example_raw_proj}}.}
\examples{
  ## see Note for links to examples.
}
