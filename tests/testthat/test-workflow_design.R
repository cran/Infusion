cat(cli::col_yellow("test get_workflow_design():\n"))

# with .... , version="2.3.0" the first gave non-integer sizes and the second gave a  -Inf
x1 <- get_workflow_design(npar=2L, test_fac=1.666); x1$fac <- NULL
x2 <- get_workflow_design(npar=2L, test_fac=1/(1/1.666)); x2$fac <- NULL
crit <- identical(x1,x2)
testthat::test_that(
  'get_workflow_design() robust to floating point precision.', 
  testthat::expect_true(crit)
) 

# Future version:
(foo <- sapply(seq(20)*2L, function(npar) get_workflow_design(npar)$subblock_nbr))
testthat::test_that(
  'get_workflow_design()$subblock_nbr is back compatible.', 
  testthat::expect_true({identical(
    foo[1:12],
    c(8L, 8L, 10L, 14L, 18L, 22L, 26L, 30L, 34L, 38L, 40L, 40L))
  })
)

(foo <- sapply(seq(20)*2L, function(npar) get_workflow_design(npar)$subblock_size) )
testthat::test_that(
  'get_workflow_design()$subblock_size is back compatible.', 
  testthat::expect_true({identical(
    foo[1:12],
    c(125L, 375L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 525L, 575L))
  })
)

