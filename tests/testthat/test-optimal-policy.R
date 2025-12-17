# test optimal policy estimator

test_that("estimate_optimal_policy returns policy of correct length", {
  set.seed(1)
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- rbinom(n, 1, 0.5)
  tau <- 0 # true treatment effect
  Y <- tau * Z + X[,1] + rnorm(n)
  
  res <- estimate_optimal_policy(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner_ols
  )
  
  expect_type(res, "list")
  expect_true("policy" %in% names(res))
  expect_length(res$policy, n)
  expect_true(all(res$policy %in% c(0,1)))
})


test_that("policy assigns treatment when true CATE is positive", {
  set.seed(2)
  n <- 200
  X <- matrix(rnorm(n), n, 1)
  Z <- rbinom(n, 1, 0.5)
  
  tau <- 3
  Y <- tau * Z + rnorm(n)
  
  res <- estimate_optimal_policy(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner_ols
  )
  
  # Most units should be treated
  expect_true(mean(res$policy) > 0.8)
})


test_that("oracle value dominates estimated policy value", {
  set.seed(3)
  n <- 150
  X <- matrix(rnorm(n), n, 1)
  Z <- rbinom(n, 1, 0.5)
  
  tau <- rnorm(n) # ATE is 0
  y_0 <- rnorm(n)
  y_1 <- y_0 + tau
  Y  <- ifelse(Z == 1, y_1, y_0)
  
  res <- estimate_optimal_policy(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner_ols,
    y_pop = cbind(y_1, y_0)
  )
  
  expect_true(res$oracle_value >= res$policy_value)
  expect_true(res$regret >= 0)
})





