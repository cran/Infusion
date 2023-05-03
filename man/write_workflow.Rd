\name{write_workflow}
\alias{write_workflow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Workflow template
}
\description{
code{write_workflow} writes a workflow script for inference. Beyond possibly saving some typing, this suggests what may be a reasonably good starting workflow. One should not expect to control all options of the workflow through the \code{write_workflow} arguments.
}
\usage{
write_workflow(con = stdout(), lower, upper, nUnique, Simulate, simulator_args=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
  A connection object or a character string. Passed to \code{\link{writeLines}}.
}
  \item{lower}{
A named numeric vector of lower bounds for parameter space. 
}
  \item{upper}{
A named numeric vector of upper bounds for parameter space. 
}
  \item{nUnique}{
Number of simulations of the process (i.e. of rows of the reference table) in the first iteration. 
}
  \item{Simulate}{
  Function, or function name as a string. Sets the same-named \code{add_reftable} argument.
}
  \item{simulator_args}{
  list of arguments for the simulator. Sets the \dots in the \code{add_reftable} call.
}
  \item{\dots}{
  Sets additional arguments in the \code{refine} call.
}
}
\value{
No return value. Used for the side-effect text, written to the connection.
}
\examples{
  write_workflow(
    ## arguments for init_reftable():
    lower=c(logTh1=-2,logTh2=-2,logTh3=-2,logTh4=-2,ar=0.01,logMu=-5,MEANP=0.01),
    upper=c(logTh1=1, logTh2=1, logTh3=1, logTh4=1, ar=0.99,logMu=-2,MEANP=0.99),
    nUnique = 1000,
    #
    ## for add_reftable():
    Simulate="schtroumf",   # name of a user-defined R function
    simulator_args= list(   # Imagine that schtroumf() require arguments 'exe_path' and 'cmdline':
      exe_path="'path_to_smurf_executable'",      
      cmdline="'smurf.exe -a -b -c -d'"
    ),    # Do check the quotation marks in the output...
    #
    ## optional arguments for refine():    
    n=8000/3.2, CIs=TRUE, update_projectors=FALSE)
}
