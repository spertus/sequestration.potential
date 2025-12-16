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

test_that("estimate_cate returns CIs with OLS", {
  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- rbinom(n, 1, 0.5)
  Y <- 2 * Z + rnorm(n)
  
  res <- estimate_cate(
    Y, Z, X,
    learner = learner_ols,
    num.trees = 500
  )
  
  expect_true(all(c("cate", "se", "lower", "upper") %in% names(res)))
  expect_equal(nrow(res), n)
  expect_true(all(res$upper >= res$lower))
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


test_that("estimate_cate returns CIs with RF learner", {
  skip_if_not_installed("grf")
  
  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- rbinom(n, 1, 0.5)
  Y <- 2 * Z + rnorm(n)
  
  res <- estimate_cate(
    Y, Z, X,
    learner = learner_rf,
    num.trees = 500
  )
  
  expect_true(all(c("cate", "se", "lower", "upper") %in% names(res)))
  expect_equal(nrow(res), n)
  expect_true(all(res$upper >= res$lower))
})




