cat(cli::col_yellow("test Infusion.Rd long examples (primitive workflow):\n"))

data(densv)
data(myrnorm) # 'Test are run in an environment that inherits from the package's namespace environment'
set.seed(456)
summliksurf <- infer_surface(densv) ## infers a log-likelihood surface
currMSL <- MSL(summliksurf) ## find the maximum of the log-likelihood surface
currMSL <- refine(currMSL, map.asp=1) ## refine iteratively

#expect_equal(currMSL$MSL$MSLE,c("mu"=4.1161537,"s2"=0.8982683),tolerance=1e-4)
