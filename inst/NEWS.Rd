\name{NEWS}
\title{Infusion News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
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
