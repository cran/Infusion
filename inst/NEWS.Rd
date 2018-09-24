\name{NEWS}
\title{Infusion News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
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
