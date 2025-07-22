cat(cli::col_yellow("test 'constr_crits':\n"))

if (Infusion.getOption("example_maxtime")>9) { 
  myrnorm <- function(mu,s2,sample.size=40L) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    return(c(mean=mean(s),var=var(s)))
  } # simulate means and variances of normal samples of size 'sample.size'
  set.seed(123)
  # simulated data with stands for the actual data to be analyzed:  
  Sobs <- myrnorm(mu=4,s2=1) 
  # With constraints (the final plot explains the constraint):
  heart <- quote({ x <- 3*(mu-4.25);  y <- 3*(s2-0.75); x^2+(y-(x^2)^(1/3))^2-1})
  parsp_h <- init_reftable(lower=c(mu=2.8,s2=0.4), upper=c(mu=5.2,s2=2.4)
                           #, constr_crit=heart
                           )
  simuls_h <- add_reftable(Simulate="myrnorm", parsTable=parsp_h)
  c_densv <- infer_SLik_joint(simuls_h, stat.obs=Sobs, constr_crits = heart)
  c_slik_j <- MSL(c_densv, CIs=FALSE, eval_RMSEs=FALSE) 
  c_slik_j <- refine(c_slik_j, target_LR=10, ntot=3000, eval_RMSEs=FALSE) 
}

