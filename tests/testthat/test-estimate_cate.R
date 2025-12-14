test_that("estimate_cate works with OLS learner", {
  set.seed(123)
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- rbinom(n, 1, 0.5)
  Y <- 1 + 2 * Z + X[,1] + rnorm(n)
  
  cate <- estimate_cate(
    Y, Z, X,
    learner = learner_ols
  )
  
  expect_type(cate, "double")
  expect_length(cate, n)
  expect_true(abs(mean(cate) - 2) < 0.5)
})

test_that("estimate_cate rejects bad learner output", {
  bad_learner <- function(...) 1:3
  
  expect_error(
    estimate_cate(
      Y = 1:10,
      Z = rep(0:1, 5),
      X = matrix(rnorm(20), 10, 2),
      learner = bad_learner
    )
  )
})




