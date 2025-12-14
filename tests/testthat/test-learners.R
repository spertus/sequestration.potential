# OLS learner
test_that("learner_ols returns correct length", {
  n <- 50
  X <- matrix(rnorm(n * 3), n, 3)
  Z <- rbinom(n, 1, 0.5)
  Y <- Z + rnorm(n)
  
  cate <- learner_ols(Y, Z, X, X)
  
  expect_length(cate, n)
  expect_true(is.numeric(cate))
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



# Random forest learner
test_that("learner_rf returns numeric vector of correct length", {
  skip_if_not_installed("grf")
  
  set.seed(123)
  n <- 100
  X <- matrix(rnorm(n * 3), n, 3)
  Z <- rbinom(n, 1, 0.5)
  Y <- 1 + 2 * Z + X[,1] + rnorm(n)
  
  cate <- learner_rf(
    Y = Y,
    Z = Z,
    X = X,
    x_pop = X,
    num.trees = 500
  )
  
  expect_type(cate, "double")
  expect_length(cate, n)
})


test_that("learner_rf recovers average treatment effect in simple DGP", {
  skip_if_not_installed("grf")
  
  set.seed(42)
  n <- 300
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- rbinom(n, 1, 0.5)
  
  tau <- 3
  Y <- X[,1] + tau * Z + rnorm(n)
  
  cate <- learner_rf(
    Y = Y,
    Z = Z,
    X = X,
    x_pop = X,
    num.trees = 1000
  )
  
  expect_true(abs(mean(cate) - tau) < 0.75)
})


test_that("learner_rf is reproducible with fixed seed", {
  skip_if_not_installed("grf")
  
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- rbinom(n, 1, 0.5)
  Y <- Z + rnorm(n)
  
  set.seed(1)
  cate1 <- learner_rf(Y, Z, X, X, num.trees = 500)
  
  set.seed(1)
  cate2 <- learner_rf(Y, Z, X, X, num.trees = 500)
  
  expect_equal(cate1, cate2)
})


test_that("learner_rf errors if grf is not installed", {
  withr::local_envvar(R_TESTS_SKIP_GRF = "true")
  
  if (requireNamespace("grf", quietly = TRUE)) {
    skip("grf is installed; skipping missing-package test")
  }
  
  expect_error(
    learner_rf(1:10, rep(0:1, 5), matrix(rnorm(20), 10, 2), matrix(rnorm(20), 10, 2))
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



