\name{NEWS}
\title{Infusion News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
\section{CHANGES UP TO VERSION 2.1.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item add_reftable() allows a more extended interface for the \code{Simulate} function, which may now handle a data frame of parameter values.
\item New default refine(., nbCluster) value, based on a more definite guess of a good value and aiming at faster analyses. 
\item New convenient wrapper init_reftable() around the init_grid() function imported from the 'blackbox' package. Infusion now depends on versions >= 1.1.41 of that package, which include a faster init_grid() function. 
\item plot2Dprof() can use parallelisation. 
\item plot2Dprof(., pars) syntax extended for more flexible specification of sets of profiles. 
\item refine.default() gets new argument 'CIs' for more transparent control of CI computations. 
\item refine.default() can now request parallel execution of ranger(), according to refine()'s 'nb_cores' argument by default, but with refine()'s 'cluster_args' argument further allowing independent control of ranger()'s 'num.threads' independently from 'nb_cores'.  
\item project.default() gets new arguments 'use_oob' and 'is_trainset' to bypass a costly step in specific (but commonplace) cases; and also argument 'methodArgs'. 
\item nbCluster="max" syntax allowed in some contexts for fast gaussian mixture modelling by fixing the tried number of clusters to a single value automatically generated. 
\item The 'xLLiM' package (new 'Suggested' dependency) can be used as an alternative to 'Rmixmod' for joint density modelling, by calling infer_SLik_joint(., using="xLLiM").
\item New convenience function write_workflow().
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item The structure of 'SLik_j' objects has been modified, meaning that the result of running infer_SLik_joint() has to be re-generated to be compatible with functions of which this object is an argument. 
\item add_reftable() and add_simulation() gain argument 'parsTable' with is an alias (and now the preferred name) for pre-existing argument 'par.grid'. 
}
}
\subsection{UP-TO-DATE WORKFLOW}{
\itemize{
\item The documentations for add_reftable() and add_simulation() have been revised, to better emphasize the inference workflow based on add_reftable() over the primitive workflow based on add_simulation(), and split in different files for clarity. The definitions of the functions have also been revised, mostly in a backward-compatible way, although subtle changes may result in the primitive workflow.    
}
}
\subsection{DEPENDENCY CHANGES}{
\itemize{
\item 'matrixStats' is an imported package.
}
}
}
\section{CHANGES UP TO VERSION 2.0.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item New function SLRT() for "summary-likelihood ratio tests", including bootstrap-corrected ones. For the latter correction, another function get_LRboot() provides a fast approximation to bootstrap distribution of likelihood ratio statistic. 
\item Better choice of initial value of optimization operations in likelihood profile computations, giving more correct 1D profile plots and likelihood ratio tests (and possibly improving all subsequent operations in iterative workflows). 
\item More thorough implementation of parallelization in add_simulation() (and add_reftable()), allowing forking, and control of random number generator.
\item Finer control of parallelisation in refine(), through new features of the 'cluster_args' argument. 
\item project() accepts (barely documented) "fastai" and "keras" methods, interfacing the packages of the same name (themselves interfacing python libraries).
\item New default values of various controls in calls to ranger::ranger.  
\item infer_SLik_joint() and goftest() will check for linear dependencies among the summary statistics (which would cause bugs in goftest()), thanks to the new convenience function check_raw_stats().
\item MSL() now computes the hessian of summary likelihood at its maximum, to check parameter identifiability (with limited success). This hessian is included in the return value of MSL() and then in the 'MSL' environment stored in e.g.) objects of class 'SLik_j'. 
\item New function focal_refine() to refine the likelihood surface for given parameter values.  
\item The previously private function predict.SLik_j() that evaluates the likelihood for specific parameter points is now part of the API.
\item New argument 'plot.slices' for plot.SLik_j() and plot.SLik().  
\item New argument 'decorations' for plot1Dprof() and plot2Dprof().  
\item Improved control of number of clusters in Gaussian mixture modelling to avoid fitting more parameters than data. 
\item projections using 'ranger' now calls it with argument importance="permutation" which may be used to select raw statistics for the goodness of fit test.
\item Check by add_reftable() of its result, to detect and more clearly warn about otherwise obscure problems that may occur at later steps.
\item profile.SLik() gets an 'init' argument, though it's more for programming purposes than intended as a user-level feature. 
\item New 'summLik' extractor. summLik() uses a fit object to evaluate the summary-likelihood function for distinct parameter values and even for new data.
\item plot_proj(), a hastily written convenience function to plot a diagnostic plot for a projection from an SLik_j object.
\item Several code fixes for more efficient handling of large reference tables (notably, resolving some memory issues when using "ranger" projection results).
\item New functions get_nbCluster_range() (and seq_nbCluster(), less directly useful) to help in controlling the number of clusters used in Gaussian mixture modelling.
\item get_from() gets new argument 'force' to force computation of elements that may be missing from the fit object.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Formals of add_reftable() modified (backward incompatibility if the former first argument 'simulations' was named in the function call).
\item Infusion.getOption("nb_cores") now controls ranger::predict and ranger::ranger, with the expected effect according to the general effect of 'nb_cores' in Infusion (a single CPU is used if 'nb_cores' is NULL, while ranger functions use multithreading by default, i.e., when their 'num.threads' argument is NULL).  
\item Backward-incompatible changes in project() arguments controlling training sizes.
\item Backward-incompatible changes in global options controlling the number of clusters used in Gaussian mixture modelling.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Incorrect reference to randomForest in documentation for project().
}
}
\subsection{DEPENDENCY CHANGES}{
\itemize{
\item 'caret' is back in suggested packages, but for its findLinearCombos() function rather than for any modelling method.
\item 'ranger' is now an imported package.
}
}
}
\section{CHANGES UP TO VERSION 1.5.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item New goftest() function, which may return a test of goodness of fit.
\item Improved sampling of parameter space in refine().
\item Several functions have a new argument 'cluster_args', passed to parallel::makeCluster(), for better control of the parallel computations. These functions' 'nb_cores' argument now acts as a shortcut for cluster_args$spec, in a backward-compatible way.
\item add_reftable() results can be more easily subsetted thanks to a new '[' method.
\item For better control of 1D profiles, plot1Dprof() gets new options for argument 'type', and new argument 'control'. 
\item Refined use of function optimizers, which generally performs better than the previously used base optimizer (optim()). 
\item plot.SLik_j(), plot1Dprof() and plot2Dprof() now re-run MSL(., eval_RMSEs=FALSE, CIs=FALSE) when they detect a new likelihood maximum, and discard pre-existing CI/RMSE information as these depend on the maximum. plot1Dprof() and plot2Dprof() now have a return value indicating whether MSL() was re-run.
\item New extractor get_from() to extract elements from summary-likelihood objects in a backward-compatible way.
\item Better control of verbosity in refine().
\item Tenfold larger default 'knotNbr' for ranger() projection method than previously.
\item New 'logLik' extractor.
}
}
\subsection{DEPENDENCY CHANGES}{
\itemize{
\item 'Rmixmod', which had disappeared from the dependencies in recent versions, is back in Suggested packages.
\item 'nloptr' and 'minqa' now used for optimization.
\item 'crayon' added in Suggests, with no visible effects yet in standard use.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
%\item Both default value of global option 'nbCluster' and new precise meaning of this option. Results of all methods using Rmixmod (the default) are modified. 
\item Default model for Gaussian mixture modelling changed.
\item add_simulation() argument 'Simulate' can now be a function (rather than a name) and the return value's attr(.,"Simulate") is always the function rather than its name. 
\item New Infusion.options "mixturing" (which replaces some previously undocumented options). 
\item Slightly modified parameter sampling in refine(), so that all ensuing results will be slightly modified. 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item add_simulation() could fail to pass the definition of 'Simulate' to child processes in parallel computations.
\item project.default() could fail on one-dimensional statistics (a toy case).
}
}
}
\section{CHANGES UP TO VERSION 1.4.1}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Out-of-bag projections used whenever relevant in projections using random-forests methods (which definitely improves the results over "in-bag" projections).
\item New 'eval_RMSEs' argument of function refine.default().
\item New argument 'update_projectors' in refine.default(), effective only for the reference table method.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
%\item Both default value of global option 'nbCluster' and new precise meaning of this option. Results of all methods using Rmixmod (the default) are modified. 
\item 'ranger' package included in Suggests; random forest by 'ranger' is the new default projection 'method' in project.character().
\item infer_SLik_joint() can now (automatically) provide correct results when there is *one* probability mass in the distribution of the summary statistics.
\item infer_SLik_joint() can now handle 'using="mclust"' (but this will still provide incorrect results when there is a probability mass in the distribution of the summary statistics).
\item Objects of class 'SLik_j' using projection now retain a cumulative table of unprojected simulations, from which new projections can be generated.
\item Revised documentation including examples for the inference method based on reference tables.
\item Internal methods in infer_SLik_joint() are modified by default so that numerical results are slightly modified. See its new argument 'marginalize' to control this.
\item Defaults Infusion options 'projTrainingSize' and 'projKnotNbr' for projection by REML have been increased, because the called procedures are more efficiently implemented in versions of spaMM >= 2.7.0.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item RMSE computation on SLik_j objects could fail following Rmixmod::mixmodCluster()'s failure to fit a bootstrap replicate. 
\item Potential problems with figure rendering when Rstudio's graphic device was *not* used. 
}
}
}
\section{CHANGES UP TO VERSION 1.3.0}{
\subsection{NEW FEATURES}{
\itemize{
\item New 'nb_cores', 'packages' and 'env' arguments for more widely allowing parallel simulation in add_simulation() and refine().
\item New '...' (or 'control.Simulate') argument in add_simulation() for controlling the 'Simulate' arguments.
\item Dependencies changed for parallel computations; ''doSNOW' usage has been reinstated optionally (though only 'foreach' appears in the package dependencies). See help("Infusion.options") for details.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item plot.SLik() failed for more than two parameters.
\item project() did not work on reference tables. 
}
}
}
\section{CHANGES UP TO VERSION 1.2.0}{
\subsection{MAJOR DEPENDENCY CHANGES}{
\itemize{
\item DESCRIPTION now Suggests 'Rmixmod' rather than Imports it. This means that, even if 'Rmixmod' were to be archived from CRAN (which has nearly occurred in March 2018), 'Infusion' could still be installed and run without 'Rmixmod', using 'mclust' as an alternative for model-based clustering. 'Rmixmod' is still the preferred method. 
\item 'mclust', 'caret' and 'ks' have been removed from DESCRIPTION (but can still be used without changes in user's code). 
\item Imports package 'pbapply' to draw progress bars, including in parallel computations; no longer Suggests the package 'foreach' (and optionally 'doSNOW') for the same effect.
}
}
\subsection{NEW FEATURES}{
\itemize{
\item infer_SLik_joint() can now use procedures from the 'mclust' package (controlled by new argument 'using'). 
\item The calls to 'mclust' procedures in other Infusion functions have been revised to be comparable to calls to 'Rmixmod' procedures. This means in particular that the 'G' and 'modelNames' arguments of 'mclust' procedures are set to values equivalent to the corresponding arguments of calls to 'Rmixmod' procedures, and that Infusion uses AIC to select the number of gaussian clusters in both cases (although selection by AIC is implemented in neither of the clustering packages). 
}
}
}
\section{CHANGES UP TO VERSION 1.1.8}{
\subsection{NEW FEATURES}{
\itemize{
\item infer_logLs() can now analyze the input distributions in parallel. This is controlled by its new 'nb_cores' argument, and may make use of another new argument, 'packages'. 
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item Infusion.options(<new values>) did not return old values.
}
}
}
\section{CHANGES UP TO VERSION 1.1.0}{
\subsection{NEW FEATURES}{
\itemize{
\item New function infer_SLik_joint() to infer likelihood surfaces from a simulation table where each simulated data set is drawn for a distinct (vector-valued) parameter, as is usual for reference tables in Approximate Bayesian Computation.
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item 'Infusion' now depends on package 'blackbox' for several functions moved there from package 'spaMM'. 
}
}
}
