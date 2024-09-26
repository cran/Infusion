\name{handling_NAs}
\alias{handling_NAs}
\alias{NA_handling}
\alias{boundaries-attribute}
\title{
Discrete probability masses and NA/NaN/Inf in distributions of summary statistics.
}
\description{
This section was written for the primitive workflow and may be largely irrelevant for the up-to-date one. 
It explains the use of the \code{boundaries} attribute of observed statistics to handle (1) values of the summary statistics that can occur with some probability mass; (2) special values (NA/NaN/Inf) in distributions of summary statistics. This further explains why \code{Infusion} handles special values by removing affected distributions unless the \code{boundaries} attribute is used.   
}
\details{
  Special values may be encountered in an analysis. For example, trying to estimate a regression coefficient when the predictor variable is constant may return a NaN. Since functions such as \code{refine} automatically add simulated distributions, this problem must be automatically handled by the user's simulation function or by the package functions, rather than by user's tinkering with the Infusion procedures.   
  
  The user must consider what s-he would do if actual data also included NA/NaN/Inf values. If such data would not be subject to a statistical analysis, then the simulation procedure must reflect that, otherwise the analysis will be biased. The processing of reference tables by Infusion functions applies \code{na.omit()} on the tables so any line containing NA's will be removed. The drawbacks are that the number of informative simulations is reduced and that inference will be difficult if the data-generating parameters were indeed prone to induce data that would not be subject to statistical analysis. Thus, it may be necessary to simulate alternative data until no special values are obtained and the target size of the simulated distribution is reached. One solution is for the user to write a simulation function that calls itself recursively until a valid summary statistic is produced. Care is then needed to avoid infinite recursion (which might well indicate unlikely parameter values). 
  
  Alternatively, if one considers that special values are informative about parameters (in the above example of a regression coefficient, if a constant predictor variable says something about the parameters), then NA/NaN/Inf must be replaced by a (fixed) dummy numerical value which is flagged to be distinctly handled, using the \code{boundaries} attribute of the observed summary statistics. The simulation function should return statistic \code{foo=-1} (say) instead of \code{foo=NaN}, and one should then set \code{attr(<observed>,"boundaries") <- c(foo=-1)}. 
  
    The boundary attribute is also useful to handle all values of the summary statistics that can occur with some probability mass. For example if the estimate \code{est_p} of a probability takes values 0 or 1 with positive probability, one should set \code{attr(<observed>,"boundaries") <- c(p_est=0,p_est=1)}.

}
